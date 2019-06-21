%% Description
% @file 	saturation.m 
%           Dependencies: Scenario.m, Antenna.m, MDC.m, VM.m, euclidian_sm_LTE.m
% @author	Jonathan Mendes de Almeida (MSc Student) & Marcelo A. Marotta, PhD (MSc Co-Adivisor)
% @email	jonathanalmd at gmail dot com / jonathan at aluno dot unb dot br / marcelo dot marotta at unb dot br 
% @page     jonathanalmd.github.io
% @date     06/10/2019 
% @info     MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
% @brief	MatLab code for the saturation problem


%% Scenario Parameters
scenario = Scenario;
scenario = scenario.start();
display(scenario);

%% Basic problem inputs

% Number of areas (hexagons) -- Scenario length
n_sites = 7;

% Time slots
T = 24; % 

%% Antennas
% Number of MacroCells antennas per covered area
mc_antennas_per_area = 1;

% Number of SmallCells clusters per covered area
sc_clusters_per_area = 1;
% Number of SmallCells antennas per cluster
sc_antennas_per_cluster = 4;

% Number of Antennas 
M_macrocell = n_sites * mc_antennas_per_area;
M_smallcell = n_sites * sc_clusters_per_area * sc_antennas_per_cluster;
M = M_macrocell + M_smallcell;

%% MDCs    
% https://www.ec2instances.info/?selected=a1.medium,c4.8xlarge
% Number of machine classes
I = 3;

% Number of MDC's
S_macrocell = n_sites * mc_antennas_per_area % one MDC per MC antenna
S_smallcell = n_sites * sc_clusters_per_area % one MDC per SC cluster
S = S_macrocell + S_smallcell

% https://aws.amazon.com/ec2/instance-types/c5/
% https://www.microway.com/knowledge-center-articles/detailed-specifications-of-the-skylake-sp-intel-xeon-processor-scalable-family-cpus/
% Processor cycles (for each machine class i - from S U S')
P_is = [2 3;
        3 4;
        4 5
        ]; 

% Number of cores
N_is = [4   8;
        8   16;
        16  32
        ];

% Processor efficiency (for each machine class i - from S U S')
Ef_is = [2 4;
         4 8;
         8 16
         ];

% Machine pricing
A_is = [30 20;
        50 30;
        90 40
        ];


%% Workload
% Decoder: linear complexity
% Number of decoder recursions
decoder_recursions = 7;
% Number of decoder instructions
decoder_instructions = 200;
% Number of operations for each bit
W = decoder_recursions * decoder_instructions;

% Transmited data (from normal distribution) - Gamma_st (antenna saturation)
close all;  
transmited_data_mt = zeros([M T]);
for i = 1:M
    %Macrocells: Residential area 
    if i <= 3
        time_var = i;
        probabilities = normpdf(-1.96:1.98/(T/2):1.96, 0, 1)*2;
        time=[6+time_var:24 1:5+time_var];
        transmited_data_mt(i,time) = norminv(probabilities,  460, 350);
    %Macrocells: Urban area 
    elseif i <= 7
        time_var = 7 - i;
        probabilities = normpdf(-1.96:1.98/(T/2):1.96, 0, 1)*2;
        time=[19+time_var:24 1:18+time_var];
        transmited_data_mt(i,time) = norminv(probabilities,  180, 130);
    %Smallcells: Residential area 
    elseif i <= 19
        time_var = 10 - i;
        probabilities = normpdf(-1.96:1.98/(T/2):1.96, 0, 1)*2;
        time=[6:24 1:5];
        transmited_data_mt(i,time) = norminv(probabilities,  460, 350);
    %Smallcells: Urban area 
    else
        time_var = 28 - i;
        probabilities = normpdf(-1.96:1.98/(T/2):1.96, 0, 1)*2;
        time=[6:24 1:5];
        transmited_data_mt(i,time) = norminv(probabilities,  180, 130);
    end
end

bar(1:24, transmited_data_mt(1,:));


% Workload: W * Gamma

%% Vertical allocation constraint variables
% Block length (worst case - LTE)
block_len = 6114; % bits
% Processor cycles (for each machine class i)
% Processor efficiency (for each machine class i)
% Number of operations 
n_operations = W * block_len;


%% Round-trip Delay (RTD) variables
% Speed of light: 299.792 km/s or 299792458 m/s
c = 299.792;
% Distance between an antenna m and an MDC s (km) -> FIX: use euclidian distance
d_sm = 1;
% 1) Propagation Delay
prop_delay = (3 * d_sm) / c;

% Block length (block_len)
% Fiber-optic flow rate (Gbit/s)
fiber_flow = 10;
% 2) Transmission Delay
trans_delay = block_len / fiber_flow;

% Hops distance (km)
d_hops = 50;
% 3) Hops Delay
hop_delay = floor(d_sm / (d_hops/10));

% Round-trip Delay - i, s, m (machine, mdc, antenna)
RTD = prop_delay + trans_delay + hop_delay;

% Time Constraint -- b_ismt (seconds)
% RTD < sigma
sigma = 0.003;




%% Saturation Problem
function [vec, fval, answer, resume, output_a, output_b] = opt_assignment( )
    %% Description
    % Function Output: vec, fval, answer, resume, output_a, output_b
    % Function Parameters: none
    
    %% Navigation 
    % Anonynous function to navigate vectors
    % Reshape cannot work with common reasoning since it works with rows first
    % followed by collumns, use this definition
    % Lines first 
    navA = @(i,s,t) sub2ind([I,M,T],i,m,t);
    navB = @(i,s,m,t) sub2ind([I,S,M,T],i,s,m,t) + I*M*T;
    navC = @(i,s,m,t) sub2ind([I,S,M,T],i,s,m,t) + I*M*T + I*S*M*T;
    
%     nav2d = @(m, n) (s-1)*N+m;

    %% A and b constraints matrixes
    % One line for each constraint
    n_constr_lines = I*M*T + I*S*M*T + I*S*M*T; % forall t in T, forall i in I, forall m in M U M'
    % I*S*M*T columns 
    n_constr_cols = I*S*M*T;
    
    % A definition
    A = zeros([n_constr_lines, n_constr_cols]);
    % Right side of the equation 
    b = zeros([n_constr_lines, 1]);
    
    %% Map the constraints
    % Head declaration for matrix navigation
    ihead = 1;
    
    % First constraint: horizontal alocation
    for m = 1:M
        for t = 1:T
            for i = 1:I
                for s = 1:S
                    A(ihead, navB(i,s,m,t)) = transmited_data_mt(m,t) * W;
                end
                A(ihead, navA(i,s,t)) = -(P_is(i) * N_is(i,s) * Ef_is(i,s));
                b(ihead) = 0;
                ihead = ihead + 1;
            end
        end
    end
    
    % Second constraint: vertical alocation
    for m = 1:M
        for t = 1:T
            for s = 1:S
                for i = 1:I
                    t_proc = (W * block_len) / (P_im(i,m) * Ef_im(i,m)) 
                    %(ihead, navB(i,s,m,t) = t_proc + (RTD(i,s,m) * 2);
                    b(ihead) = sigma;
                    ihead = ihead + 1;
                end
            end
        end
    end
    
    % Third constraint: migration cost
    
    
    %% Map the objective function
    % (?)
    
    %% Set upper and lower bounds
    
    
    
    
    %% Get optimal solution
    % [vec, fval, answer, resume] = intlinprog(f,1:N*N, [], [], Aeq,beq,l_bound,u_bound);

    
    
    
    
    
    %% Transform solution into matrix (NxN)
    % output_a = reshape(vec, [N,N]);
    
    %% Get final solution: solution matrix (binary matrix) x weight matrix 
    % output_b = output_a .* a;
    
end