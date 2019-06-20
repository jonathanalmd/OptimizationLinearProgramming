function d = euclidian_sm_LTE(RRH, UE)
    N = size(UE, 2); % Size of the second dimension of UE (1 x number of UEs)
    M = size(RRH, 2); % Size of the second dimension of RRH (1 x number of RRHs)
    d = zeros(N, M);

    for i = 1:N
        for j = 1:M
            d(i,j) = ((UE(i).x - RRH(j).x)^2 + (UE(i).y - RRH(j).y)^2)^0.5;
        end
    end
end