%% Description
% @file saturation.m 
%           Dependencies: Scenario.m, Antenna.m, MDC.m, VM.m, euclidian_sm_LTE.m
% @author	Jonathan Mendes de Almeida (MSc Student) & Marcelo A. Marotta, PhD (MSc Co-Adivisor)
% @email	jonathanalmd at gmail dot com / jonathan at aluno dot unb dot br / marcelo dot marotta at unb dot br 
% @page     jonathanalmd.github.io
% @date     06/10/2019 
% @info     MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
% @brief	MatLab code for the saturation problem

%% Create Scenario
%scenario = Scenario;
%scenario = scenario.start();
%display(scenario);

%% Saturation Problem
function [vec, fval, answer, resume, output_a, output_b] = opt_assignment(scenario)
    %% Description
    % Function Output: vec, fval, answer, resume, output_a, output_b
    % Function Parameters: scenario obj
    
    %% Navigation 
    % Anonynous function to navigate vectors
    % Reshape cannot work with common reasoning since it works with rows first
    % followed by collumns, use this definition
    % Lines first 
    I = scenario.I;
    S = scenario.S;
    M = scenario.M;
    T = scenario.T;
    navA = @(i,s,t) sub2ind([I,S,T],i,s,t);
    navB = @(i,s,m,t) sub2ind([I,S,M,T],i,s,m,t) + I*S*T;
    navC = @(i,s,m,t) sub2ind([I,S,M,T],i,s,m,t) + I*S*T + I*S*M*T;

    %% A and b constraints matrixes
    % One line for each constraint
    n_constr_lines = I*S*T + I*S*M*T + I*S*M*T; % forall t in T, forall i in I, forall m in M U M'
    % I*S*M*T columns 
    n_constr_cols = I*S*T + I*S*M*T + I*S*M*T;
    
    % A definition
    A = zeros([n_constr_lines, n_constr_cols]);
    % Right side of the equation 
    b = zeros([n_constr_lines, 1]);
    
    %% Map the constraints
    % Head declaration for matrix navigation
    ihead = 1;
    
    % First constraint: horizontal alocation
    for s = 1:S
        for t = 1:T
            for i = 1:I
                for m = 1:M
                    A(ihead, navB(i,s,m,t)) = scenario.transmited_data_mt(m,t) * scenario.W;
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
                    t_proc = (scenario.W * scenario.block_len) / (P_is * Ef_is); 
                    
                    % 1) Propagation Delay
                    prop_delay = (3 * scenario.d_sm(s,m)) / scenario.c;
                    % 2) Transmission Delay
                    trans_delay = scenario.block_len / scenario.fiber_flow;
                    % 3) Hops Delay
                    hop_delay = floor(scenario.d_sm(s,m) / (scenario.d_hops));
                    % Delays
                    delays = prop_delay + trans_delay + hop_delay;
            
                    A(ihead, navB(i,s,m,t)) = t_proc + (delays * 2); % RTD_ism
                    
                    b(ihead) = scenario.sigma;
                    
                    ihead = ihead + 1;
                end
            end
        end
    end
    
    % Third constraint: migration cost
    for s = 1:S
        for t = 2:T
            for m = 1:M
                for i = 1:I
                    for s1 = 1:S
                        if s1 ~= s
                            A(ihead, navB(i,s,m,t-1)) = -1; 
                        else
                            A(ihead, navB(i,s,m,t-1)) = 1; 
                        end
                        A(ihead, navC(i,s,m,t)) = 1; 
                    end
                    b(ihead) = 1;

                    ihead = ihead + 1;
                end
            end
        end
    end
    
    
    %% Map the objective function
    f = zeros([I*S*T + I*S*M*T + I*S*M*T, 1]);
    for i = 1:I
        for s = 1:s
            for t = 1:T
                f(navA(i,s,t)) = scenario.mdc(s).vm(i).price;
                for m = 1:M
                    K = scenario.mdc(s).vm(i).price * 0.5;
                    f(navC(i,s,m,t)) = K;
                end
            end
        end
    end
                    
    %% Set upper and lower bounds
    % Min 0 Max 1
    u_bound(1:I*S*T) = 99999999;
    u_bound(I*S*T+1 : I*S*T + I*S*M*T + I*S*M*T) = 1;
    l_bound = zeros([I*S*T + I*S*M*T + I*S*M*T, 1]);
    
    %% Get optimal solution
    [vec, fval, answer, resume] = intlinprog(f,1 : I*S*T + I*S*M*T + I*S*M*T, A, b, [], [], l_bound, u_bound);
        
    %% Transform solution into matrix (NxN)
    output_a = reshape(vec(1 : I*S*T), [I,S,T]);
    
    %% Get final solution: solution matrix (binary matrix) x weight matrix 
    output_b = reshape(vec(I*S*T+1 : I*S*T + I*S*M*T), [I,S,M,T]);
    
end