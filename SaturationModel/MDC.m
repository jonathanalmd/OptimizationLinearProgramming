classdef MDC
    properties
        index = 1;
        % Position
        x = 0;
        y = 0;
        % machine class i
        type = 1;
        % Processor cycles (P_is)
        cycles = 1;
        % Number of cores (N_is)
        n_cores = 1;
        % Processor eficiency (Ef_is)
        eficiency = 1; 
        % Pricing (A_is)
        price = 1;
        

    end
    
    methods
        function obj = MDC(type, position, index)
            obj.type = type;
            obj.x = position(1);
            obj.y = position(2);
            obj.index = index;
            if type == 1
                obj.bandwidth = 10;       %MHz
                obj.frequency = 2.0;      %GHz Frequency: 2.0 GHz
                obj.maxpower = 30;        %dBm
                %Outdoor 0, otherwise: 20+0.5*rand()*min(25, distance)
                obj.height = 10;          %meters
                obj.gain = 5;             %dBi Antenna gain + connector loss
                obj.tti = 0.01;            %TTI of 10 ms
            end
        end
        
    end
    
end
