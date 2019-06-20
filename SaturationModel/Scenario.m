classdef Scenario
    %% System parameters
    % Scenario 1 - 36.872 Annex A
    % Hexagonal grid, 3 sectors per site, case 1
    % Both 19 Macro sites and 7 Macro sites can be used. Companies should 
    % indicate whether 19 or 7 sites are used when presenting the results.
    % Clusters uniformly random within macro geographical area; 
    properties
        dmacromacro = 1000;
%         dmacromacro = 500;
        dmacroue = 35;
        dmacrocluster = 105;
        dsmallue = 5;
        dsmallsmall = 20;
        
        dropradius_mc = 250;
        dropradius_sc = 500;
        dropradius_sc_cluster = 50;
        dropradius_ue_cluster = 70;
        
        
        n_sites = 7;
        n_clusters = 1;         %1, 2, optional 4
        n_antennas = 10;             %4, 10
        n_ues = 60;            %Per macrocell area
        
        uesindoor = 0.2;       %Type ratio
        uessmallcell = 2/3;    %Dropping ratio
        
        
        index = 1;
        macrocells;
        smallcells;
        clusters;
        cluster_centers;
        center;
        ues;
        cluster_ues;
        dmacroradius;
        R;
    end
    
    methods

        %% Main loop for starting the scenario
        function obj = start(obj)
            obj.dmacroradius = obj.dmacromacro * 0.425;
            % Creating macrocell antennas positioning
            % Centralized antenna
            obj.center = Antenna(0, [obj.dmacroradius*obj.n_sites/2, obj.dmacroradius*obj.n_sites/2], obj.index);
%             obj.center = Antenna(0, [750, 750], obj.index);
            obj.index = obj.index+1;
            [obj.macrocells, obj.index] = obj.hexaCluster(obj.dmacromacro, obj.center, obj.index);
            
            % Spawning rrhs    
            obj.smallcells = [];
            for i = obj.macrocells
                for j = 1:obj.n_clusters
                    [scs_aux, cluster_center, obj.index] = obj.spawnAntennas(i, 1, obj.index);
                    obj.cluster_centers = [obj.cluster_centers; cluster_center];
                    obj.smallcells = [obj.smallcells scs_aux];
                end
            end
            
            obj.clusters = reshape(obj.smallcells, obj.n_sites, obj.n_clusters, obj.n_antennas);
            obj.cluster_centers = reshape(obj.cluster_centers, obj.n_sites, obj.n_clusters);
            
            %Spawning UEs 
            for i = 1:length(obj.macrocells)
                ues_aux = obj.spawnUEs(obj.macrocells(i), obj.clusters(i,:), obj.cluster_centers(i,:));
%                 plot([ues_aux.x], [ues_aux.y], '.', 'MarkerSize',4);
                obj.ues = [obj.ues ues_aux];
            end
            obj.cluster_ues = reshape(obj.ues, obj.n_sites, obj.n_ues);
%             hold off;
        end
        
        %% Creates a hexagon map structure for an antenna
        function [cluster, index] = hexaCluster(obj, radius, center, index)
             cluster = center;
%             for i = [-1, 1, 2, 3, 5, 6] 45degree
             for i = 2*(0:(obj.n_sites-2))+1 
                position = [center.x + radius * cos(i*pi/6), center.y + radius * sin(i*pi/6)];
                ant = Antenna(center.type, position, index);
                index = index+1;
                cluster = [cluster ant];
            end
        end
        
        %% Random uniform distribution of antennas in a cluster 
        function [antennas, cluster_center, index] = spawnAntennas(obj, center, type, index)
            antennas = [];
            cluster_center = obj.dropObject(center, obj.dmacromacro*0.425, obj.dmacrocluster);
            reset = 0;
            count = 1;
            positions = [];
            while count <= obj.n_antennas
                if reset > 1000
                    count=1;
                    reset=0;
                    cluster_center = obj.dropObject(center, obj.dmacromacro*0.425, obj.dmacrocluster);
                    positions = [];
                end
                %Positioning
                position = obj.dropObject(cluster_center, obj.dropradius_sc_cluster, 0);
                d = [euclidian_LTE(position, positions) euclidian_LTE(position, [obj.smallcells])];
                if  isempty(d) || isempty((d(d < obj.dsmallsmall)))  
                    count=count+1;
                    positions = [positions; position];
                else
                    reset = reset + 1;
                end
            end
            
            for position = 1:length(positions)
                antenna = Antenna(type, [positions(position).x,positions(position).y], index);
                antenna.parent = center.index;
                antennas = [antennas antenna];
                index = index+1;
            end
        end
        %% Randon uniform distribution of ues in clusters 
        function ues = spawnUEs(obj, macrocell, clusters, clustercenter)
            ues = [];
            reset = 0;
            count = 1;
            positions = [];
            while count <= obj.n_ues
                if reset > 1000
                    count=1;
                    reset=0;
                    positions = [];
                end

                if rand() < 0.6666
                    %Positioning
                    cluster = clustercenter(randi([1 obj.n_clusters]));
                    position = obj.dropObject(cluster, obj.dropradius_ue_cluster, 0);
                else
                    position = obj.dropObject(macrocell, obj.dmacroradius, obj.dmacroue);
                end
                d = euclidian_LTE(position, clusters(:));
                
                if ~isempty(d(d < obj.dsmallue)) 
                    reset = reset + 1;
                else
                    count=count+1;
                    positions = [positions; position];
                end
            end
            
            for position = 1:length(positions)
                ue = UE([positions(position).x,positions(position).y]);
                ue.type = rand() < 0.80;
                ues = [ues ue];
            end
        end
        function d = dropObject(~, center, radius, min_distance)
            not_done = true;
            while not_done
                d.x = radius * (1 - 2 * rand()) + center.x;
                d.y = radius * (1 - 2 * rand()) + center.y;
                not_done = euclidian_LTE(d, center) < min_distance;
            end
        end
        function p = signalMatrix(prbs)
            %% Calculating power from 3GPP TS 36.213 p13 
            %PCmax_i = MIN {PEMAX ? DTC, PPowerClass ? MAX(MPR + A-MPR, P-MPR) ? DTC}
            %Pemax to be used in the target cell. If absent the UE applies the
            %maximum power according to the UE capability 3GPP TS 36.331 version 10.7.0 Release 10
            %?TC - Allowed operating band edge transmission power relaxation. =0 for band 33 1.9GHz ~ 1.92GHz
            %MPR <= 1 for PRB <=50 and QPSK modulation
            %A-MPR used for mimo transmissions
            %P-MPR important for networks not in accordance to 3GPP
            % In the end PCmax_i = PowerClass - rand(); 
            % P0 = -106 dB 3GPP 36.814 Table: A.2.2-1
            % alpha = 1  3GPP 36.814 Table: A.2.2-1
            P0 = -106;
            Pmax = PowerClass - rand(size(prbs));
            % P = min{Pmax, P0+10?log10 MPUSCH+??L+?mcs+f(?i)}
            % P = min{Pmax, 10log10(MPUSCH,c)+PO_PUSCH,c(j)+?.PLc+?TF,c(i) + fc(i)}
            p = max(Pmin, min(Pmax, P0 + 10 * log10(prbs) + Pmax + alpha * L ));
        %     ch = 3.2*(log10(11.75*1.5))^2-4.97;
        %     s = gain_t + gain_r + 69.55+26.16*log10(frequency) - 13.82*log10(32) - ch + (44.9-6.55*log10(32))*log(d);
        end
        
        %%  Return a matrix of distances - Antennas x UEs
        function d_mn = distances(obj)
            d_mn =  euclidian_LTE([obj.macrocells obj.smallcells], obj.ues);
        end
        
        %%  Return a matrix of distances - Antennas x MDCs (smallcells = cluster center / macrocells = same place)
        function d_sm = distances_sm(obj)
            d_sm =  euclidian_sm_LTE([obj.macrocells obj.smallcells], obj.ues);
        end
        
        %% Plot the scenario considering the antennas and UEs placement
        function plotScenario(obj, sectors)
            if exist('sectors','var') == 0 
                sectors = 0;
            end
            
            figure;
            plot([obj.macrocells.x], [obj.macrocells.y],'k^','MarkerSize',8);
            hold on;
            if sectors == 1
                plot([obj.macrocells(1).x, obj.macrocells(1).x + 0.575*obj.dmacromacro * cos(pi/3*(1))], [obj.macrocells(1).y, obj.macrocells(1).y + 0.575*obj.dmacromacro*sin(pi/3*(1))],'-.k');
                plot([obj.macrocells(1).x, obj.macrocells(1).x + 0.575*obj.dmacromacro * cos(pi/3*(3))], [obj.macrocells(1).y, obj.macrocells(1).y + 0.575*obj.dmacromacro*sin(pi/3*(3))],'-.k');
                plot([obj.macrocells(1).x, obj.macrocells(1).x + 0.575*obj.dmacromacro * cos(pi/3*(5))], [obj.macrocells(1).y, obj.macrocells(1).y + 0.575*obj.dmacromacro*sin(pi/3*(5))],'-.k');
            end
            for i = 2:obj.n_sites
                plot(obj.macrocells(i).x + 0.575*obj.dmacromacro * cos(pi/3*(0:6)), obj.macrocells(i).y + 0.575*obj.dmacromacro*sin(pi/3 * (0:6)),'--');
                if sectors == 1
                    plot([obj.macrocells(i).x, obj.macrocells(i).x + 0.575*obj.dmacromacro * cos(pi/3*(1))], [obj.macrocells(i).y, obj.macrocells(i).y + 0.575*obj.dmacromacro*sin(pi/3*(1))],'-.k');
                    plot([obj.macrocells(i).x, obj.macrocells(i).x + 0.575*obj.dmacromacro * cos(pi/3*(3))], [obj.macrocells(i).y, obj.macrocells(i).y + 0.575*obj.dmacromacro*sin(pi/3*(3))],'-.k');
                    plot([obj.macrocells(i).x, obj.macrocells(i).x + 0.575*obj.dmacromacro * cos(pi/3*(5))], [obj.macrocells(i).y, obj.macrocells(i).y + 0.575*obj.dmacromacro*sin(pi/3*(5))],'-.k');
                end
            end
            plot([obj.smallcells.x], [obj.smallcells.y], 'o', 'MarkerSize',4);
            plot([obj.cluster_centers.x], [obj.cluster_centers.y],'r*');
            plot([obj.ues.x], [obj.ues.y], '.', 'MarkerSize',8);
            hold off;  
            if obj.n_sites > 1
                ylim([0 obj.dmacromacro*floor(obj.n_sites/2) ])
                xlim([0 obj.dmacromacro*floor(obj.n_sites/2) ])
            else
                ylim([0 obj.dmacromacro/2])
                xlim([0 obj.dmacromacro/2])
            end
            display('plotted');
        end
    end
    
end