classdef UE
    %UE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        height = 1.5;                %meters
        tti = 0.01;                  %TTI of 10 ms
        type = 0;                    %0 outdoor, 1 otherwise
        x = 0;
        y = 0;
        sector = 0;
        antenna = 0;
        shadowing = 3;
        noise_figure = 9;
        gain = 0;
        nantennas = 1;
        powerClass = 23 - rand();
        P0 = -53;                   % P0 = -106 dB 3GPP 36.814 Table: A.2.2-1
        buffer = 1024;
    end
    
    methods
        
        function obj = UE(position)
            obj.x = position(1);
            obj.y = position(2);
        end
        
        %% Calculating power from 3GPP TS 36.213 p13
        function p = signalPower(obj, antennas, alpha, prbs)
            d = euclidian(obj, antennas)/1000;
            %PCmax_i = MIN {PEMAX ? DTC, PPowerClass ? MAX(MPR + A-MPR, P-MPR) ? DTC}
            %Pemax to be used in the target cell. If absent the UE applies the
            %maximum power according to the UE capability 3GPP TS 36.331 version 10.7.0 Release 10
            %?TC - Allowed operating band edge transmission power relaxation. =0 for band 33 1.9GHz ~ 1.92GHz
            %MPR <= 1 for PRB <=50 and QPSK modulation
            %A-MPR used for mimo transmissions
            %P-MPR important for networks not in accordance to 3GPP
            % In the end PCmax_i = PowerClass - rand(); 
            
            % alpha = 1  3GPP 36.814 Table: A.2.2-1
            
            pClass = obj.powerClass; 
            Pmax = pClass - rand();
            % P = min{Pmax, P0+10?log10 MPUSCH+??L+?mcs+f(?i)}
            % P = min{Pmax, 10log10(MPUSCH,c)+PO_PUSCH,c(j)+?.PLc+?TF,c(i) + fc(i)}
            p = zeros(length(antennas));
%             p(ant) = max(-43, min(Pmax, P0 + 10 * log10(prbs) + alpha * obj.pathloss(antennas(ant), d(ant)) + obj.penetration(d(ant))));
            p = min(Pmax, obj.P0 + 10 * log10(prbs) + alpha * obj.pathloss(antennas, d));
            %     ch = 3.2*(log10(11.75*1.5))^2-4.97;
            %     s = gain_t + gain_r + 69.55+26.16*log10(frequency) - 13.82*log10(32) - ch + (44.9-6.55*log10(32))*log(d);
        end
        
        %% return SNR in dB
        function p = sinr(obj, antenna, alpha, prbs)
            d = euclidian(obj, antenna);
            thermal = thermal_noise(300, prbs * 180000);
            IoT =  (thermal + obj.noise_figure)/thermal;
            %p(ant) = max(-43, min(Pmax, P0 + 10 * log10(prbs) + alpha * obj.pathloss(antennas(ant), d(ant)) + obj.penetration(d(ant))));
            p = obj.P0 + 10 * log10(prbs) + (alpha - 1) * obj.pathloss(antenna, d) - obj.penetration(d) - IoT - thermal - awgn(9,0.5);
            %ch = 3.2*(log10(11.75*1.5))^2-4.97;
            %p = gain_t + gain_r + 69.55+26.16*log10(frequency) - 13.82*log10(32) - ch + (44.9-6.55*log10(32))*log(d);
        end
        
        %% return the pathloss in DB
        function pl = pathloss(obj, antenna, d)
            %Calculating the path loss for macrocell NLOS - 3GPP TR 36.814 
            % UMa NLOS - Table B.1.2.1-1
            if antenna.type == 0
                pl = 161.04 - 7.1*log10(20) + 7.5 * log10(20) ...
                    -(24.37 - 3.7*(20/antenna.height)^2) * log10(antenna.height) ...
                    + (43.42-3.1*log10(antenna.height)) * (log10(d)-3) ...
                    + 20*log10(antenna.frequency) - (3.2*(log10(11.75*obj.height))^2-4.97);
            else
            %Calculating the path loss for microcell NLOS - 3GPP TR 36.814 
            % UMi NLOS - Table B.1.2.1-1
                pl = 36.7*log10(d) + 22.7 + 26*log10(antenna.frequency);                
            end
        end
        function pe = penetration(obj, d)
            %Calculating the penetration outdoor/indoor - 3GPP TR 36.872
            if obj.type == 0    %outdoor
                pe=0;
            else                %indoor
                pe = 20 + 0.5 * rand() * min(d, 25); 
            end
        end
    end
    
end

