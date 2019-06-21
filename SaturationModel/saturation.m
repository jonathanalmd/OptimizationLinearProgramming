%% Description
% @file 	saturation.m 
%           Dependencies: Scenario.m, Antenna.m, MDC.m, VM.m, euclidian_sm_LTE.m
% @author	Jonathan Mendes de Almeida (MSc Student) & Marcelo A. Marotta, PhD (MSc Co-Adivisor)
% @email	jonathanalmd at gmail dot com / jonathan at aluno dot unb dot br / marcelo dot marotta at unb dot br 
% @page     jonathanalmd.github.io
% @date     06/10/2019 
% @info     MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
% @brief	MatLab code for the saturation problem


%% Creating Scenario
scenario = Scenario;
scenario = scenario.start();
display(scenario);


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
        time=[19:24 1:18];
        transmited_data_mt(i,time) = norminv(probabilities,  180, 130);
    end
end

bar(1:24, transmited_data_mt(1,:));



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
    navC = @(i,s,m,t) sub2ind([I,S,M,T],i,s,m,t) + I*S*T + I*S*M*T;
    
%     nav2d = @(m, n) (s-1)*N+m;

    %% A and b constraints matrixes
    % One line for each constraint
    n_constr_lines = I*S*T + I*S*M*T + I*S*M*T; % forall t in T, forall i in I, forall m in M U M'
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
    % [vec, fval, answer, resume] = intlinprog(f,1:N*N, A, b, [], [],l_bound,u_bound);

    
    
    
    
    
    %% Transform solution into matrix (NxN)
    % output_a = reshape(vec, [N,N]);
    
    %% Get final solution: solution matrix (binary matrix) x weight matrix 
    % output_b = output_a .* a;
    
end