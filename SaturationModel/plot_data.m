scenario = Scenario;
scenario = scenario.start();

[vec, fval, answer, resume, n_ismt, output_a] = opt_assignment(scenario);


% P_is = scenario.mdcs(s).vms(i).cycles;
% Ef_is = scenario.mdcs(s).vms(i).efficiency;
% N_is = scenario.mdcs(s).vms(i).n_cores;
% A_is = scenario.mdcs(s).vms(i).price;

total_cost = zeros([1,24]);
mc_cost = zeros([1,24]);
sc_cost = zeros([1,24]);

I = scenario.I;
S = scenario.S;
M = scenario.M;
T = scenario.T;
ihead = 1;

for t = 1:T
    for m = 1:M
        for s = 1:S
            for i = 1:I
                A_is = scenario.mdcs(s).vms(i).price;
                cur_cost = (A_is * output_a(i,s,m,t));
                total_cost(ihead) = cost(ihead) + cur_cost;
                if M == 1
                    mc_cost(ihead) = mc_cost(ihead) + cur_cost;
                else
                    sc_cost(ihead) = sc_cost(ihead) + cur_cost;
                end 
            end
        end
    end
    ihead = ihead + 1; 
end



