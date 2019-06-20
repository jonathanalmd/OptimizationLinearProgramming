classdef Scenario
    %% System parameters
    % Scenario 1 - 36.872 Annex A
    % Hexagonal grid, 3 sectors per site, case 1
    % Both 19 Macro sites and 7 Macro sites can be used. Companies should 
    % indicate whether 19 or 7 sites are used when presenting the results.
    % Clusters uniformly random within macro geographical area; 
    properties
        %% Basic problem inputs
        % Number of sites (hexagons) -- Scenario length
        n_sites = 7;
        % Time slots
        T = 24; % 

        %% Antennas
        % Number of MacroCells antennas per covered area
        mc_antennas_per_site = 1;
        % Number of SmallCells clusters per covered area
        mc_clusters_per_site = 1;
        % Number of SmallCells antennas per cluster
        sc_antennas_per_cluster = 4;
        % Number of Antennas 
        % M_macrocell = n_sites * mc_antennas_per_site;
        M_macrocell;
        % M_smallcell = n_sites * mc_clusters_per_site * sc_antennas_per_cluster;
        M_smallcell;
        % M = M_macrocell + M_smallcell;
        M;

        %% MDCs    
        % https://www.ec2instances.info/?selected=a1.medium,c4.8xlarge
        % Number of machine classes
        I = 3;
        % Number of MDC's
        % S_macrocell = n_sites * mc_antennas_per_site % one MDC per MC antenna
        S_macrocell;
        % S_smallcell = n_sites * mc_clusters_per_site % one MDC per SC cluster
        S_smallcell;
        % S = S_macrocell + S_smallcell
        S;
        % MC MDCs
        mc_mdcs;
        % SC MDCs
        sc_mdcs;
        % List of all MDCs
        mdcs;
%         % MDCs positions
%         S_xy;
%         
%         % https://aws.amazon.com/ec2/instance-types/c5/
%         % https://www.microway.com/knowledge-center-articles/detailed-specifications-of-the-skylake-sp-intel-xeon-processor-scalable-family-cpus/
%         % Processor cycles (for each machine class i - from S U S')
%         P_is = [2 3;
%                 3 4;
%                 4 5
%                 ]; 
%         % Number of cores
%         N_is = [4   8;
%                 8   16;
%                 16  32
%                 ];
%         % Processor efficiency (for each machine class i - from S U S')
%         Ef_is = [2 4;
%                  4 8;
%                  8 16
%                  ];
%         % Machine pricing
%         A_is = [30 20;
%                 50 30;
%                 90 40
%                 ];

        %% Workload
        % Decoder: linear complexity
        % Number of decoder recursions
        decoder_recursions = 7;
        % Number of decoder instructions
        decoder_instructions = 200;
        % Number of operations for each bit
        % W = decoder_recursions * decoder_instructions;
        W;
        %transmited_data_mt = zeros([M T]);
        transmited_data_mt;
        % Workload = transmited_data_mt * W
        
        %% Vertical allocation constraint variables
        % Block length (worst case - LTE)
        block_len = 6114; % bits
        % Processor cycles (for each machine class i)
        % Processor efficiency (for each machine class i)
        % Number of operations 
        % n_operations = W * block_len;
        n_operations;

        %% Round-trip Delay (RTD) variables
        % Speed of light: 299.792 km/s or 299792458 m/s
        c = 299.792;
        % Distance between an antenna m and an MDC s (km) 
        % d_sm = zeros([S M]);
        d_sm;
        % 1) Propagation Delay
        % prop_delay = (3 * d_sm) / c;
        prop_delay
        % Block length (block_len)
        % Fiber-optic flow rate (Gbit/s)
        fiber_flow = 10;
        % 2) Transmission Delay
        % trans_delay = block_len / fiber_flow;
        trans_delay;
        % Hops distance (km)
        d_hops = 50;
        % 3) Hops Delay
        % hop_delay = floor(d_sm / (d_hops/10));
        hop_delay;
        % Round-trip Delay - i, s, m (machine, mdc, antenna)
        %RTD_ism = prop_delay + trans_delay + hop_delay;
        RTD_ism;

        % Time Constraint -- b_ismt (seconds)
        % RTD < sigma
        sigma = 0.003;

        %% ORIGINAL
        dmacromacro = 1000; % 500
        dmacroue = 35;
        dmacrocluster = 105;
        dsmallue = 5;
        dsmallsmall = 20;
        
        dropradius_mc = 250;
        dropradius_sc = 500;
        dropradius_sc_cluster = 50;
        % dropradius_ue_cluster = 70;
             
        %n_sites = 7;
        %n_clusters = 1;         %1, 2, optional 4
        % n_antennass = 10;             %4, 10
        % n_ues = 60;            %Per macrocell area 
        % uesindoor = 0.2;       %Type ratio
        % uessmallcell = 2/3;    %Dropping ratio
           
        index = 1;
        macrocells;
        smallcells;
        clusters;
        cluster_centers;
        center;
        % ues;
        % cluster_ues;
        dmacroradius;
        % R;
    end
    %% Methods
    methods
        %% Main loop for starting the scenario
        function obj = start(obj)
            % Number of Antennas 
            obj.M_macrocell = obj.n_sites * obj.mc_antennas_per_site;
            obj.M_smallcell = obj.n_sites * obj.mc_clusters_per_site * obj.sc_antennas_per_cluster;
            obj.M = obj.M_macrocell + obj.M_smallcell;

            % Number of MDC's
            obj.S_macrocell = obj.n_sites * obj.mc_antennas_per_site % one MDC per MC antenna
            obj.S_smallcell = obj.n_sites * obj.mc_clusters_per_site % one MDC per SC cluster
            obj.S = obj.S_macrocell + obj.S_smallcell

            % Number of operations for each bit
            obj.W = obj.decoder_recursions * obj.decoder_instructions;

            obj.transmited_data_mt = obj.antennaSaturationNorm();
            % Workload = transmited_data_mt * W

            % Number of operations 
            obj.n_operations = obj.W * obj.block_len;

            % Distance between an antenna m and an MDC s (km) -> FIX: use euclidian distance
            obj.d_sm = zeros([obj.S obj.M]);
            % 1) Propagation Delay
            obj.prop_delay = (3 * obj.d_sm) / obj.c;

            % 2) Transmission Delay
            obj.trans_delay = obj.block_len / obj.fiber_flow;

            % 3) Hops Delay
            obj.hop_delay = floor(obj.d_sm / (obj.d_hops/10));
            % Round-trip Delay - i, s, m (machine, mdc, antenna)
            obj.RTD_ism = obj.prop_delay + obj.trans_delay + obj.hop_delay;
 
            
            % ORIGINAL
            obj.dmacroradius = obj.dmacromacro * 0.425;
            % Creating macrocell antennas positioning
            % Centralized antenna
            position = [obj.dmacroradius*obj.n_sites/2, obj.dmacroradius*obj.n_sites/2]
            obj.center = Antenna(0, position, obj.index);
            % MC MDC (center)
            obj.mc_mdcs = MDC(1, [position(1)+50 position(2)+50], obj.index);
            obj.index = obj.index+1;
            % Create mc antennas and MDCs
            [obj.macrocells, obj.mc_mdcs, obj.index] = obj.hexaCluster(obj.dmacromacro, obj.center, obj.mc_mdcs, obj.index);
            
            % Spawning rrhs    
            obj.smallcells = [];
            for i = obj.macrocells
                for j = 1:obj.mc_clusters_per_site
                    [scs_aux, cluster_center, obj.index] = obj.spawnAntennas(i, 1, obj.index);
                    obj.cluster_centers = [obj.cluster_centers; cluster_center];
                    obj.smallcells = [obj.smallcells scs_aux];
                    % SC MDCs
                    obj.sc_mdcs = [obj.sc_mdcs MDC(2, [cluster_center.x cluster_center.y], obj.index)];
                end
            end
            
            % All MDCs
            obj.mdcs = [obj.mc_mdcs obj.sc_mdcs];
            
            
            
            obj.clusters = reshape(obj.smallcells, obj.n_sites, obj.mc_clusters_per_site, obj.sc_antennas_per_cluster);
            obj.cluster_centers = reshape(obj.cluster_centers, obj.n_sites, obj.mc_clusters_per_site);
            
            %             % Spawning UEs 
            %             for i = 1:length(obj.macrocells)
            %                 ues_aux = obj.spawnUEs(obj.macrocells(i), obj.clusters(i,:), obj.cluster_centers(i,:));
            %                 % plot([ues_aux.x], [ues_aux.y], '.', 'MarkerSize',4);
            %                 obj.ues = [obj.ues ues_aux];
            %             end
            %             obj.cluster_ues = reshape(obj.ues, obj.n_sites, obj.n_ues);
            %             % hold off;
        end
        
        %% Creates a hexagon map structure for an antenna
        function [cluster, mdcs, index] = hexaCluster(obj, radius, center, centermdc, index)
             cluster = center;
             mdcs = centermdc;
           % for i = [-1, 1, 2, 3, 5, 6] 45degree
             for i = 2*(0:(obj.n_sites-2))+1 
                position = [center.x + radius * cos(i*pi/6), center.y + radius * sin(i*pi/6)];
                % mc MDC: antenna_type, position, index
                mdc = MDC(1, [position(1)+50 position(2)+50], index);
                ant = Antenna(center.type, position, index);
                index = index+1;
                cluster = [cluster ant];
                mdcs = [mdcs mdc];
            end
        end
        
        %% Random uniform distribution of antennas in a cluster 
        function [antennas, cluster_center, index] = spawnAntennas(obj, center, type, index)
            antennas = [];
            cluster_center = obj.dropObject(center, obj.dmacromacro*0.425, obj.dmacrocluster);
            reset = 0;
            count = 1;
            positions = [];
            while count <= obj.sc_antennas_per_cluster
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
        
        %% Drop Object
        function d = dropObject(~, center, radius, min_distance)
            not_done = true;
            while not_done
                d.x = radius * (1 - 2 * rand()) + center.x;
                d.y = radius * (1 - 2 * rand()) + center.y;
                not_done = euclidian_LTE(d, center) < min_distance;
            end
        end
        
        %% Creating Transmited data matrix (from normal distribution) - Gamma_mt (antenna saturation)
        function transmited_data_mt = antennaSaturationNorm(obj)
            transmited_data_mt = zeros([obj.M obj.T]);
            for i = 1:obj.M
                %Macrocells: Residential area 
                if i <= 3
                    time_var = i;
                    probabilities = normpdf(-1.96:1.98/(obj.T/2):1.96, 0, 1)*2;
                    time=[6+time_var:24 1:5+time_var];
                    transmited_data_mt(i,time) = norminv(probabilities,  460, 350);
                %Macrocells: Urban area 
                elseif i <= 7
                    time_var = 7 - i;
                    probabilities = normpdf(-1.96:1.98/(obj.T/2):1.96, 0, 1)*2;
                    time=[19+time_var:24 1:18+time_var];
                    transmited_data_mt(i,time) = norminv(probabilities,  180, 130);
                %Smallcells: Residential area 
                elseif i <= 19
                    time_var = 10 - i;
                    probabilities = normpdf(-1.96:1.98/(obj.T/2):1.96, 0, 1)*2;
                    time=[6:24 1:5];
                    transmited_data_mt(i,time) = norminv(probabilities,  460, 350);
                %Smallcells: Urban area 
                else
                    time_var = 28 - i;
                    probabilities = normpdf(-1.96:1.98/(obj.T/2):1.96, 0, 1)*2;
                    time=[6:24 1:5];
                    transmited_data_mt(i,time) = norminv(probabilities,  180, 130);
                end
            end

            bar(1:24, transmited_data_mt(1,:));
        end
        
        
        %%  Return a matrix of distances - Antennas x MDCs (smallcells = cluster center / macrocells = same place)
        function d_sm = distances_sm(obj)
            d_sm =  euclidian_sm_LTE(obj.ues, [obj.macrocells obj.smallcells]);
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
            % plot([obj.ues.x], [obj.ues.y], '.', 'MarkerSize',8);
            plot([obj.mdcs.x], [obj.mdcs.y], 'square', 'MarkerSize', 8);
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