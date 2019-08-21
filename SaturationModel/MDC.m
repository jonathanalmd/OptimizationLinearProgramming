classdef MDC
    properties
        index = 1;
        % Position
        x = 0;
        y = 0;
        % Antenna Type (MC = 1, SC = 2)
        antenna_type = 1;
        % VMs
        vms = [];
    end
    
    methods
        function obj = MDC(antenna_type, position, index)
            obj.antenna_type = antenna_type;
            obj.x = position(1);
            obj.y = position(2);
            obj.index = index;
            % If MacroCell
            if antenna_type == 1
                % machineclass, cycles, n_cores, efficiency, price, index
                obj.vms = [obj.vms VM(1, 3.9*10^9, 16, 16, 10.680, index)];
                obj.vms = [obj.vms VM(2, 3.9*10^9, 36, 16, 1.530, index)];
                obj.vms = [obj.vms VM(3, 3.9*10^9, 48, 16, 2000.040, index)];
            % Else SmallCell
            else
                obj.vms = [obj.vms VM(1, 3.6*10^9, 2, 8, 10.085, index)];
                obj.vms = [obj.vms VM(2, 3.6*10^9, 4, 8, 0.170, index)];
                obj.vms = [obj.vms VM(3, 3.6*10^9, 8, 8, 1000.340, index)];
            end
        end        
    end
    
end
