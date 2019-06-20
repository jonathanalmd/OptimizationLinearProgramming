classdef MDC
    properties
        index = 1;
        % Position
        x = 0;
        y = 0;
        % Antenna Type (MC = 1, SC = 2)
        antenna_type = 1;
        % VMs
        mdc_vms = [];
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
                obj.mdc_vms = [obj.mdc_vms VM(1, 3*10^9, 8, 4, 20)];
                obj.mdc_vms = [obj.mdc_vms VM(2, 4*10^9, 16, 8, 30)];
                obj.mdc_vms = [obj.mdc_vms VM(3, 5*10^9, 32, 16, 40)];
            % Else SmallCell
            else
                obj.mdc_vms = [obj.mdc_vms VM(1, 2*10^9, 4, 2, 30)];
                obj.mdc_vms = [obj.mdc_vms VM(2, 3*10^9, 8, 4, 50)];
                obj.mdc_vms = [obj.mdc_vms VM(3, 4*10^9, 16, 8, 90)];
            end
        end
        
    end
    
end
