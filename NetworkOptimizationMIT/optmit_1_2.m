% @file 	opmit_1_2.m
% @author	Marcelo A. Marotta, PhD & Jonathan Mendes de Almeida
% @email	jonathanalmd@gmail.com / jonathan@aluno.unb.br
% @page     jonyddev.github.io
% @date     05/20/2019 
% @info     MSc Research at Computer Networks Lab (COMNET) -- University of Brasília (UnB)
% @brief	MatLab code for the problem formalization (Chapter 1, Example 1.4): The Assignment Problem 
%           Network Optimization: Continuous and Discrete Models, Dimitri P. Bertsekas, Massachusetts Institute of Technology (MIT)

function [vec, fval, answer, resume, output_a, output_b] = optimal_association( )
    %% Description
    % Function Output: vec, fval, answer, resume, output_a, output_b
    % Function Parameters: none
    %% Problem formulation
  
    %% Parameters
    % N persons and N gifts 
    N = 4; 
    
    % Weight matrix (m x n)
    a = [
         8     8     2     3;
         9     4     7     2;
         10    3    10     5;
         1     4     5     9
         ];
    % Head declaration for matrix navigation
    head = 1;

    %% Navigation 
    % Anonynous function to navigate vectors
    % Reshape cannot work with common reasoning since it works with rows first
    % followed by collumns, use this definition
    % Lines first 
    nav2d = @(m, n) (n-1)*N+m; % Gamma function (annonymous)

    %% Equality matrix (Aeq / beq)
    % N lines for each constraint (total = 2*N)
    % m*n (N*N) columns 
    Aeq = zeros([2*N, N*N]);
    % Right side of the equation 
    beq = ones([2*N, 1]);

    % First constraint
    % For all n in {1,...,N}
    for n = 1:N
        % For all m in {1,...,M}
        for m = 1:N
            Aeq(head, nav2d(m,n)) = 1;
        end
        % Increment Head according to n
        head = head+1;
    end
    % Second constraint
    for m = 1:N
        for n = 1:N
            Aeq(head, nav2d(m,n)) = 1;
        end
        head = head+1;
    end

    % Function
    % Transform a_mn matrix to vector 
    f = ones([1, N*N]);
    for m = 1:N
        for n = 1:N
            % Maximization problem to minimization problem (MatLab)
            f(nav2d(m,n))=-a(m,n);
        end
    end

    % Set upper and lower bounds
    % Binary
    u_bound = ones([N*N, 1]);
    l_bound = zeros([N*N, 1]);

    % Get optimal solution
    [vec, fval, answer, resume] = intlinprog(f,1:N*N, [], [], Aeq,beq,l_bound,u_bound);

    % Transform solution into matrix (NxN)
    output_a = reshape(vec, [N,N]);
    % Get final solution: solution matrix (binary matrix) x weight matrix 
    output_b = output_a .* a;
end