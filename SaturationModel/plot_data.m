scenario = Scenario;
scenario = scenario.start();

[vec, fval, answer, resume, n_ismt, output_a] = opt_assignment(scenario);

total_cost = zeros([1,24]);
mc_cost = zeros([1,24]);
sc_cost = zeros([1,24]);
allocated_cp = zeros([1,24]);

I = scenario.I;
S = scenario.S;
M = scenario.M;
T = scenario.T;
ihead = 1;

for t = 1:T
    for m = 1:M
        for s = 1:S
            for i = 1:I
                if output_a(i,s,m,t) == 1 
                    P_is = scenario.mdcs(s).vms(i).cycles;
                    Ef_is = scenario.mdcs(s).vms(i).efficiency;
                    N_is = scenario.mdcs(s).vms(i).n_cores;
                    A_is = scenario.mdcs(s).vms(i).price;
                    cur_cost = (A_is * n_ismt(i,s,m,t));
                    total_cost(ihead) = total_cost(ihead) + cur_cost;
                    if s <= 7
                        mc_cost(ihead) = mc_cost(ihead) + cur_cost;
                    else
                        sc_cost(ihead) = sc_cost(ihead) + cur_cost;
                    end 

                    proc = P_is * Ef_is * n_ismt(i,s,m,t);
                    allocated_cp(ihead) = allocated_cp(ihead) + proc;
                end
            end
        end
    end
    ihead = ihead + 1; 
end


plot(1:24, total_cost(1,:));
plot(1:24, mc_cost(1,:));
hold on
plot(1:24, sc_cost(1,:));
hold off
% 
gamma = sum(scenario.transmited_data_mt) * scenario.W;
plot(1:24, allocated_cp(1,:));
hold on
plot(1:24, gamma(1,:));
hold off


