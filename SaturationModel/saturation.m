%% Description
% @file 	saturation.m 
%           Dependencies: Scenario.m, Antenna.m, MDC.m, VM.m, euclidian_sm_LTE.m
% @author	Jonathan Mendes de Almeida (MSc Student) & Marcelo A. Marotta, PhD (MSc Co-Adivisor)
% @email	jonathanalmd at gmail dot com / jonathan at aluno dot unb dot br / marcelo dot marotta at unb dot br 
% @page     jonathanalmd.github.io
% @date     06/10/2019 
% @info     MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
% @brief	MatLab code for the saturation problem


%% Create Scenario
% % Number of sites (hexagons) -- Scenario length
% n_sites = 7;
% % Time slots
% T = 24; % 
% % Antennas
% % Number of MacroCells antennas per covered area
% mc_antennas_per_site = 1;
% % Number of SmallCells clusters per covered area
% mc_clusters_per_site = 1;
% % Number of SmallCells antennas per cluster
% sc_antennas_per_cluster = 4;
% 
% % Workload
% % Decoder: linear complexity
% % Number of decoder recursions
% decoder_recursions = 7;
% % Number of decoder instructions
% decoder_instructions = 200;
% 
% % Vertical allocation constraint variables
% % Block length (worst case - LTE)
% block_len = 6114; % bits
% 
% % Round-trip Delay (RTD) variables
% % Speed of light: 299.792 km/s or 299792458 m/s
% c = 299.792;
% fiber_flow = 10;
% % Hops distance (km)
% d_hops = 50;
% % RTD < sigma
% sigma = 0.003;
% 
% dmacromacro = 1000; % 500
% dmacrocluster = 105;
% dsmallsmall = 20;
% dropradius_mc = 250;
% dropradius_sc = 500;
% dropradius_sc_cluster = 50; 
% scenario = Scenario(n_sites, T, mc_antennas_per_site, mc_clusters_per_site, ...
%                     sc_antennas_per_cluster, decoder_recursions, decoder_instructions, ...
%                     block_len, c, fiber_flow, d_hops, sigma, dmacromacro, dmacrocluster, ...
%                     dsmallsmall, dropradius_mc, dropradius_sc, dropradius_sc_cluster ...
%                     );
scenario = Scenario;
scenario = scenario.start();
display(scenario);


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
    navA = @(i,s,t) sub2ind([I,S,T],i,s,t);
    navB = @(i,s,m,t) sub2ind([I,S,M,T],i,s,m,t) + I*S*T;
    navC = @(s,m,t) sub2ind([S,M,T],s,m,t) + I*S*T + I*S*M*T;
    navD = @(s,m) sub2ind([s,m],s,m) + I*S*T + I*S*M*T + S*M*T;
%     nav2d = @(m, n) (s-1)*N+m;

    %% A and b constraints matrixes
    % One line for each constraint
    n_constr_lines = I*S*T + I*S*M*T + S*M*T + S*M; % forall t in T, forall i in I, forall m in M U M'
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
    for s = 1:s
        for t = 1:T
            for i = 1:I
                for m = 1:M
                    A(ihead, navB(i,s,m,t)) = transmited_data_mt(m,t) * W;
                end
                % A(ihead, navA(i,s,t)) = -(P_is(i,s) * N_is(i,s) * Ef_is(i,s));
                % scenario.mdcs(s).vms(i)
                P_is = scenario.mdcs(s).vms(i).cycles;
                N_is = scenario.mdcs(s).vms(i).n_cores;
                Ef_is = scenario.mdcs(s).vms(i).efficiency;
                A(ihead, navA(i,s,t)) = -(P_is * N_is * Ef_is);
                
                b(ihead) = 0;
                
                ihead = ihead + 1;
            end
        end
    end
    
    % Second constraint: vertical alocation
    for s = 1:S
        for t = 1:T
            for m = 1:M
                for i = 1:I
                    P_is = scenario.mdcs(s).vms(i).cycles;
                    Ef_is = scenario.mdcs(s).vms(i).efficiency;
                    t_proc = (W * block_len) / (P_is * Ef_is); 
                    
                    % 1) Propagation Delay
                    prop_delay = (3 * scenario.d_sm(s,m)) / obj.c;
                    % 2) Transmission Delay
                    trans_delay = scenario.block_len / scenario.fiber_flow;
                    % 3) Hops Delay
                    hop_delay = floor(scenario.d_sm(s,m) / (scenario.d_hops));
                    % Delays
                    delays = prop_delay + trans_delay + hop_delay;
            
                    A(ihead, navB(i,s,m,t) = t_proc + (delays * 2);
                    
                    b(ihead) = sigma;
                    
                    ihead = ihead + 1;
                end
            end
        end
    end
    
    % Third constraint: migration cost
    for s = 1:s
        for t = 1:T
            for m = 1:M
                for s1 = 1:S
                    if s1 ~= s
                        A(ihead, navB(s,m,t)) = 0 % 0 or 1 (?) -> b_sm (t-1)
                    end
                end
                
                b(ihead) = 1;
                
                ihead = ihead + 1;
            end
        end
    end
    
    
    
    %% Map the objective function
    % (?)
    
    %% Set upper and lower bounds
    % Min 0 Max 1
    u_bound(1:n_constr_lines, 1:1) = 99999999;
    l_bound = zeros([n_constr_lines, 1]);

    
    
    
    %% Get optimal solution
    % [vec, fval, answer, resume] = intlinprog(f,1:N*N, A, b, [], [],l_bound,u_bound);

    
    
    
    
    
    %% Transform solution into matrix (NxN)
    % output_a = reshape(vec, [N,N]);
    
    %% Get final solution: solution matrix (binary matrix) x weight matrix 
    % output_b = output_a .* a;
    
end