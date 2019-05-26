%% Description
% @file 	opmit_1_1.m
% @author	Marcelo A. Marotta, PhD (MSc Adivisor) & Jonathan Mendes de Almeida (MSc Student)
% @email	marcelo dot marotta at unb dot br / jonathanalmd at gmail dot com / jonathan at aluno dot unb dot br
% @page     jonathanalmd.github.io
% @date     05/20/2019 
% @info     MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
% @brief	MatLab code for the problem formalization (Chapter 1, Example 1.1): Shortest Path Problem
%           Network Optimization: Continuous and Discrete Models, Dimitri P. Bertsekas, Massachusetts Institute of Technology (MIT)

%% Optimal Shortest Path - Linear Programming (LP)
function [vec, fval, answer, resume, output_a, output_b] = optimal_algorithm1( )
    %% Prototype
    %% Problem formulation
    N = 4;
    dest = N;
    source = 1;
    % Navigation 
    % Anonynous function to navigate vectors
    % Reshape cannot work with common reasoning since it works with rows first
    % followed by collumns, use this definition
    % Lines first 
    nav2d = @(m, n) (n-1)*N+m; % Gamma function (annonymous)
    nav2d = @(m, n) (n-1)*N+m;

    a_ij= [
          99 10  3 99
          99 99  1 20
          99 10 99 15
          99 99 99 99
    ];
    %     a_ij= [
    %          99  5  3 99 99 99 99 99
    %          99 99 99  1  2 99 99 99
    %          99  1 99  1 99 99  5 99
    %          99 99 99 99  2  3  4 99
    %          99  1 99 99 99  3 99  6
    %          99 99 99 99  1 99  1  2
    %          99 99 99 99 99  2 99  3
    %          99 99 99 99 99 99 99 99
    %     ];
    % N constraints (1 constraint for each graph node)
    % NxN colunns 
    Aeq = zeros(N, N*N);
    beq = zeros(N, 1);
    ihead = 1;
    for i=1:N
        for j=1:N

            nav2d(i,j)
            % Each sum
            % ->
            Aeq(ihead, nav2d(i,j)) = 1; 
            % <- 
            Aeq(ihead, nav2d(j,i)) = -1;

        end
        % For each edge i-j:
        % If Source (i = s) -> 1
        if i == source
            beq(ihead) = 1;
        else
            % If Sink (i = t) -> -1
            if i == dest
                beq(ihead) = -1;
            % Otherwise -> 0
            else
                beq(ihead) = 0;
            end
        end
        ihead=ihead+1;
    end
    
    % f = a_ij (4x4 matrix to vector[16])
    f = reshape(a_ij, [N*N,1]);
    
    % Min 0 Max 1
    u_bound = ones([N*N, 1]);
    l_bound = zeros([N*N, 1]);

    [vec, fval, answer, resume] = intlinprog(f,1:N*N, [], [], Aeq, beq,l_bound,u_bound);

    output_a = reshape(vec, [N,N]);

    output_b = output_a .* a_ij;

    display(output_a);
    display(output_b);
    display(Aeq);
    display(beq);
end