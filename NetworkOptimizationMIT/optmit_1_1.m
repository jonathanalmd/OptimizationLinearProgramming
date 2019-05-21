function [vec, fval, answer, resume, output_a, output_b] = optimal_algorithm1( )
    %% Prototype
    %% Problem formulation
    %Anonynous function to navigate the vectors 4d and 3d
    % Reshape cannot work with the logical below since it works with rows first
    % followed by collumns, use the second definition
    % nav4d = @(vec, s, m, n, t) vec((s-1)*M*N*T+(m-1)*N*T+(n-1)*T+t );
    %navCloud = @(s, t) (t-1)*S+s + S*M*N*T;
    %nav4d = @(s, m, n, t) (t-1)*M*N*S+(n-1)*M*S+(m-1)*S+s;
    % nav3d = @(m, n, t) (t-1)*M*N+(n-1)*M+m;
    N = 4;
    dest = N;
    source = 1;
    nav2d = @(m, n) (n-1)*N+m;
    % Number of variable
    %A[s,m,n,t] + B[s,t]
    % Lower bound
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
    A = zeros(N, N*N);
    b = zeros(N, 1);
    ihead = 1;
    for i=1:N
        for j=1:N

            nav2d(i,j)
            % Each sum
            % ->
            A(ihead, nav2d(i,j)) = 1; 
            % <- 
            A(ihead, nav2d(j,i)) = -1;

        end
        % For each edge i-j:
        % If Source (i = s) -> 1
        if i == source
            b(ihead) = 1;
        else
            % If Sink (i = t) -> -1
            if i == dest
                b(ihead) = -1;
            % Otherwise -> 0
            else
                b(ihead) = 0;
            end
        end
        ihead=ihead+1;
    end
    
    % f = a_ij (4x4 matrix to vector[16])
    f = reshape(a_ij, [N*N,1]);
    
    % Min 0 Max 1
    u_bound = ones([N*N, 1]);
    l_bound = zeros([N*N, 1]);

    [vec, fval, answer, resume] = intlinprog(f,1:N*N, A, b, [],[],l_bound,u_bound);

    output_a = reshape(vec, [N,N]);

    output_b = output_a .* a_ij;

    display(output_a);
    display(output_b);
    display(A);
    display(b);
end
