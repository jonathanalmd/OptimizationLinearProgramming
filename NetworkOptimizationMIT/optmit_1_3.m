%% Description
% @file 	opmit_1_3.m
% @author	Marcelo A. Marotta, PhD (MSc Adivisor) & Jonathan Mendes de Almeida (MSc Student)
% @email	marcelo dot marotta at unb dot br / jonathanalmd at gmail dot com / jonathan at aluno dot unb dot br
% @page     jonathanalmd.github.io
% @date     05/20/2019 
% @info     MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
% @brief	MatLab code for the problem formalization (Chapter 1, Example 1.3): The Max-Flow Problem
%           Network Optimization: Continuous and Discrete Models, Dimitri P. Bertsekas, Massachusetts Institute of Technology (MIT)

%% Max-Flow - Linear Programming (LP)
function [vec, fval, answer, resume, output_a] = max_flow( )
    %% Description
    % Function Output: vec, fval, answer, resume, output_a
    % Function Parameters: none
    %% Problem formulation
    %% Basic problem inputs
    N = 8; 
    % Capacity matrix (N x N) 
    a = [
         0 3 2 2 0 0 0 0;
         0 0 0 0 5 1 0 0;
         0 0 0 0 1 3 1 0;
         0 0 0 0 0 1 0 0;
         0 0 0 0 0 0 0 4;
         0 0 0 0 0 0 0 2;
         0 0 0 0 0 0 0 4;
         0 0 0 0 0 0 0 0
         ];
    source = 1;
    dest = 8;
    %     N = 6;
    %     a = [
    %          0 16 13  0  0  0;
    %          0  0 10 12  0  0;
    %          0  4  0  0 14  0;
    %          0  0  9  0  0 20;
    %          0  0  0  7  0  4;
    %          0  0  0  0  0  0;
    %         ];
    %     source = 1;
    %     dest = 6;
    
    % Head declaration for matrix navigation
    ihead = 1;
    
    %% Navigation 
    % Anonynous function to navigate vectors
    % Reshape cannot work with common reasoning since it works with rows first
    % followed by collumns, use this definition
    % Lines first 
    nav2d = @(m, n) (n-1)*N+m; % Gamma function (annonymous)

    %% Equality matrix (Aeq / beq)
    n_lin = N
    n_col = N * N
    Aeq = zeros([n_lin, n_col]);
    % Right side of the equation 
    beq = zeros([n_lin, 1]);
    
    %% Map the constraints
    % First constraint: flow conservation
    % N constraints (for each i)       
    for i=1:N
        for j=1:N
 
            if i ~= source && i ~= dest
                % Each sum
                % -> in
                Aeq(ihead, nav2d(i,j)) = 1; 
                % <- out
                Aeq(ihead, nav2d(j,i)) = -1;
            end
        end
        ihead=ihead+1;
    end
    %     Aeq = [
    %            0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 1 0 0 0 0 0 0  0 0 0 0 -1 -1 0 0  0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 1 0 0 0 0 0  0 0 0 0 0 0 0 0    0 0 0 0 -1 -1 -1 0  0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 1 0 0 0 0  0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0     0 0 0 0 0 -1 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0  0 0 0 0 1 0 0 0    0 0 0 0 1 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0  0 0 0 0 0 1 0 0    0 0 0 0 0 1 0 0     0 0 0 0 0 1 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0    0 0 0 0 0 0 1 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0
    %           ];
    %     Aeq = [
    %            0 -1 -1 -1 0 0 0 0  0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  1 0 0 0 0 0 0 0;
    %            0 1 0 0 0 0 0 0     0 0 0 0 -1 -1 0 0  0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 1 0 0 0 0 0     0 0 0 0 0 0 0 0    0 0 0 0 -1 -1 -1 0  0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 1 0 0 0 0     0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0     0 0 0 0 0 -1 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0     0 0 0 0 1 0 0 0    0 0 0 0 1 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0     0 0 0 0 0 1 0 0    0 0 0 0 0 1 0 0     0 0 0 0 0 1 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0    0 0 0 0 0 0 1 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 -1 0 0 0 0 0 0 0 0;
    %            0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0    0 0 0 0 0 0 0 0     0 0 0 0 0 0 0 0   0 0 0 0 0 0 0 1  0 0 0 0 0 0 0 1  0 0 0 0 0 0 0 1  -1 0 0 0 0 0 0 0
    %           ];
    
    %% Map the objective function
    f = zeros([1, N*N]);
    % f(nav2d(1,8)) = a(8,1);
    for j = 1:N
        % Maximization problem to minimization problem (MatLab)
        if a(1,j) > 0
            %f(j)=-a(1,j);
            f(j) = -a(1,j);
        end
    end

    %% Set upper and lower bounds
    l_bound = zeros([N*N, 1]);
    
    % Second constraint: capacity constraint
    A = a;
    reshape(A,1,[]);
    transpose(A);
    A.';
    A(:);
    u_bound = reshape(A.',[],1);
    
    %% Get optimal solution
    [vec, fval, answer, resume] = intlinprog(f,1:N*N, [], [], Aeq,beq,l_bound,u_bound);

    %% Transform solution into matrix (NxN)
    %output_a = reshape(vec, [N,N]);
    output_a = vec2mat(vec,N)

    %% Display the solution
    display(output_a);
    display(sum(vec));
   
end