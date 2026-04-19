clc;clearvars;close all; format shortG;
%% ------------------------------------------------------------------------
%% This codes analyses the metric specific data and creates the graphs
%% Input: Data from GRAINS/CHRONO, of the dynamical systems at the time steps
%% Output: The results desired, and the data for the edge elimination process
%% ------------------------------------------------------------------------

tic 
format_ast = {"Kleopatra.obj"};   % Set the name of the folder where the text files containing the data from GRAINS are found. - Kleopatra, Geographos, Bennnu
densitiesall  = {"rho 2000\","rho 3000\","rho 2400\","rho 2000\","rho 2000\"};   % Set the densities you want to consider for the network analysis
rotationsall = {"Rotation Period 1.980\", "Rotation Period 0.720\", "Rotation Period 1.980\","Rotation Period 3.060\", "Rotation Period 2.520\"};  % Set the rotation you would like to consider for the network analysis

linestyle = {"-", ":", "--", "-.", "-*"};

%% Initialize DATA

for format_i = 1:length(format_ast)
    for r_i = 1:length(rotationsall)
            for d_i = 1:length(densitiesall)


                r_i = d_i;   % Imposed only if the cases are set like this in the two structures

                density_folder = densitiesall{d_i};
                rotation_folder = rotationsall{r_i};

                photo_folder = "SavePhotos";
                mkdir(photo_folder);

                % Load the data saved from GRAINS in txt file
                timeAll.positions = readmatrix("Positions.txt");
                timeAll.velocities = readmatrix("Velocities.txt");
                timeAll.mass = readmatrix("Mass.txt");
                timeAll.radius = readmatrix("Radius.txt");
                timeAll.angmomentum = readmatrix("AngularMomentum.txt");   % This is the angular momentum
                timeAll.inertia = readmatrix("Inertia.txt");               % Moment of inertia of each body around itself in the form [x,y,z] - > diag([x,y,z])
                timeAll.accforces = readmatrix("AccumulatedForce.txt");    
                timeAll.contforces = readmatrix("ContactForces.txt");      % The sum of all contact forces on each body 
                timeAll.contactpairs = readlines("ContactPairs.txt");
                timeAll.eachforce = readlines("EachForce.txt");            % The contact forces for each contact pair
                timeAll.eachforcesave = timeAll.eachforce;
                timeAll.contactpairssave = timeAll.contactpairs;
               
                % Load the data from the initial aggregate simulation -
                % when the initial body has been created
                
                sim_inputs_path = append("C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\",core{core_i},"results\", format_ast{doc_i},density_folder,rotation_folder,"simInputs.txt");
                field = "Universal gravity constant G = ";
                if (format_ast{doc_i} == "Geographos.obj\")
                    sim_inputs_path  = append("C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\spherical_large_no_core_fixed\results\Geographos Radar-based, mid-res.obj\",density_folder,rotation_folder,"simInputs.txt");
                end

                Grav = readdata(sim_inputs_path, field);
                % Create colors for different plots
                hex_colors = { "#FF0000";  "#00FF00";   "#0000FF";   "#00FFFF";  "#FF00FF";   "#FFFF00";  "#000000";   "#FFFFFF"; "#0072BD";   "#D95319";  "#EDB120";   "#7E2F8E"; "#77AC30";   "#4DBEEE";  "#A2142F"; "#1A2B3C"; "#4F5E6D"; "#A1B2C3"; "#D4E5F6";"#FF5733"; "#33FF57";"#5733FF"; "#C0C0C0";"#800000"; "#008000";"#000080"; "#FFA500";"#4B0082"; "#EE82EE";"#4682B4"; "#20B2AA";"#DC143C"; "#8B0000";"#556B2F2"; "#2F4F4F"};

                %% Load the data into Matlab format
                
                % Check to see other distance multipliers
                threshold.distance = 10 * (timeAll.radius(1) + timeAll.radius(2));
                sizes.time = length(timeAll.mass(:,1));
                sizes.bodies = length(timeAll.mass(1,:));
                A.distance = zeros(sizes.bodies,sizes.bodies);
                %Initialize data
                centralities_distance.closeness = zeros(sizes.time, sizes.bodies);
                centralities_distance.closesum = zeros(sizes.time, sizes.bodies);
                
                % Check the number of instances that can be created, there are times when
                % errors come from GRAINS and the forces are not in equilibrium, making it
                % not suitable for networks symmilarity
                % Initialize vectors for the simulation to extract the
                % data at each time frame
                pairs.time_i = zeros(sizes.bodies,sizes.bodies, 96);
                pairs.forces_time_i = zeros(sizes.bodies,sizes.bodies);
                pairs.forces_time_All = zeros(sizes.bodies, sizes.bodies,96);
                pairs.time_All = zeros(sizes.bodies,sizes.bodies, 96);
                sumtimeAll = zeros(sizes.bodies,sizes.time);
                time_vec = round(linspace(1,sizes.time,96));
                time_ii = 1;
                
                % Create a system tolerance
                tol = 1e-6;                
                
                % Transform Adjacency List into Adjacency Matrix
                for time_i = 1:size(timeAll.positions,1)
                
                    % Save the lines that have more values in common in pairs
                    savelines = [];
                
                    if(any(time_vec == time_i))
                        length_max_local = 0;
                        for pairs_i = 1:sizes.bodies
                            length_max_local =  max(length_max_local, length(str2num(timeAll.contactpairs(pairs_i))));
                        end
                    
                        length_max = length_max_local;
                    
                        if (length_max > 2)
                            for i = 1:sizes.bodies
                                line = str2num(timeAll.contactpairs(i)) + 1;
                                timell(time_ii) = line(1);
                                line(1:2) = [];
                                line2 = str2num(timeAll.eachforce(i));
                                id_line_init = line2(2) + 1;
                                line2(1:2) = [];        
                                if(length(line) >= 1)
                                    pairs.time_All(i,line(line > i),time_ii) = 1;
                                    pairs.time_All2(i,line,time_ii) = 1;
                                    % Put in also the contact forces
                                    for line_i = 1:length(line)
                                        line3 = str2num(timeAll.eachforce(line(line_i)));
                                        line4 = str2num(timeAll.contactpairs(line(line_i))) + 1;
                                        id_line_comp = line4(2);
                                        line3(1:2) = [];
                                        line4(1:2) = [];
                    
                                        count = 0;
                                        for line_j = 1:length(line2)
                                            if (any ( abs(line2(line_j) - line3) < tol))     
                                                count = count + 1;
                                                save_j = line_j;
                                            end
                                        end
                                        if (count > 1)
                                            savelines = [savelines; [id_line_init id_line_comp time_ii]];
                                            continue;
                                        else
                                            pairs.forces_time_All(i,line(line_i),time_ii) = line2(save_j);
                                        end
                    
                                    end
                                end
                            end
                        end
                
                   % Check the pairs between the combinations of pairs that have not been
                   % imposed yet
                         for k = 1:size(savelines,1)
                    
                            i = savelines(k,1);
                            j = savelines(k,2);
                            t = savelines(k,3);
                    
                            % Skip if already assigned
                            if pairs.forces_time_All(i,j,t) ~= 0
                                continue;
                            end
                    
                            % Companion already knows → inherit
                            if pairs.forces_time_All(j,i,t) ~= 0
                                pairs.forces_time_All(i,j,t) = pairs.forces_time_All(j,i,t);
                                continue;
                            end
            
                            % Extract force lists
                            f_i = str2num(timeAll.eachforce(i)); f_i(1:2) = [];
                            f_j = str2num(timeAll.eachforce(j)); f_j(1:2) = [];
                    
                            used_i = nonzeros(pairs.forces_time_All(i,:,t));
                            used_j = nonzeros(pairs.forces_time_All(j,:,t));
                    
                            rem_i = setdiff(f_i, used_i);
                            rem_j = setdiff(f_j, used_j);
                    
                            % Find common force (tolerance-based)
                            common = rem_i(ismembertol(rem_i, rem_j, tol));
                    
                             if isempty(common)
                                 continue;   % NOTHING written → safe
                             end
                    
                          % If multiple identical forces, sum or take one (your choice)
                                force_value = sum(common);
                    
                               pairs.forces_time_All(i,j,t) = force_value;
                                pairs.forces_time_All(j,i,t) = force_value;
                          end
                    
                    
                             sumtimeAll(1:sizes.bodies,time_ii) = sum(pairs.time_All(:,:,time_ii),2);
                             time_ii = time_ii + 1;
                             pairs.time_i = timeAll.contactpairs(1:sizes.bodies);
                             pairs.forces_time_i = timeAll;
                    end
                        timeAll.contactpairs(1:sizes.bodies) = [];
                        timeAll.eachforce(1:sizes.bodies) = [];
                end
                % Check symmetry of created matrices
                %%
                for i = 1:length(time_vec)
                    pairs.time_All(:,:,i) = pairs.time_All(:,:,i) + transpose(pairs.time_All(:,:,i));
                    pairs.time_All(:,:,i) = abs(pairs.time_All(:,:,i));
                    if (issymmetric(pairs.forces_time_All(:,:,i)) && issymmetric(pairs.time_All(:,:,i)))
                        disp(true);
                    else
                        disp(i);
                        disp(false);
                    end
                end
                toc

                %% Create an initialization
                E_time = zeros(sizes.bodies,numel(time_vec));
                Ltotal = zeros(1,numel(time_vec));
                L_total_vec = zeros(sizes.bodies,3);
                L_total = zeros(numel(time_vec),3);
                Adj = cell(1,numel(time_vec));
                
                for time_i = 1:length(time_vec)
                
                    %% --- INITIALIZATION ---
                    clear E.c E.pe E.r E.t E.escape
                    tic
                
                    time_ii = time_vec(time_i);
                    disp(time_ii)
                
                    % Extract physical quantities at current time
                    positions  = transpose(reshape(timeAll.positions(time_ii,:), [3, sizes.bodies]));
                    velocities = transpose(reshape(timeAll.velocities(time_ii,:), [3, sizes.bodies]));
                    radius     = transpose(timeAll.radius(time_ii,:));
                    Mass.all   = transpose(timeAll.mass(time_ii,:));
                
                    acc_forces      = transpose(reshape(timeAll.accforces(time_ii,:), [3, sizes.bodies]));
                    contact_forces  = transpose(reshape(timeAll.contforces(time_ii,:), [3, sizes.bodies]));
                
                    acc_forces_norm     = vecnorm(acc_forces,2,2);
                    contact_forces_norm = vecnorm(contact_forces,2,2);
                
                    omega   = transpose(reshape(timeAll.angmomentum(time_ii,:), [3, sizes.bodies]));
                    inertia = transpose(reshape(timeAll.inertia(time_ii,:), [3, sizes.bodies]));
                
                
                    %% --- DISTANCE & RELATIVE VELOCITY MATRICES ---
                    % Compute pairwise distances and velocity differences
                
                    R = zeros(sizes.bodies);
                    V = zeros(sizes.bodies);
                
                    for i = 1:sizes.bodies
                        for j = i+1:sizes.bodies
                            R(i,j) = norm(positions(i,:) - positions(j,:));
                            R(j,i) = R(i,j);
                
                            V(i,j) = norm(velocities(i,:) - velocities(j,:));
                            V(j,i) = V(i,j);
                        end
                    end
                
                
                    %% --- GRAPH CONSTRUCTION ---
                    % Build adjacency matrices and corresponding graphs
                
                    % Contact graph
                    A.pairs = pairs.time_All(:,:,time_i);
                    G.pairs = graph(A.pairs, string(1:sizes.bodies));
                
                    density.pairs(time_i) = 2 * numedges(G.pairs) / ...
                        (numnodes(G.pairs) * (numnodes(G.pairs) - 1));
                
                
                    %% --- DISTANCE GRAPH ---
                    A.distance = zeros(sizes.bodies);
                
                    for i = 1:sizes.bodies
                        for j = i+1:sizes.bodies
                            dist_ij = norm(positions(i,:) - positions(j,:));
                
                            A.distance(i,j) = dist_ij;
                            A.distance(j,i) = dist_ij;
                        end
                    end
                
                    G.distances = graph(A.distance, string(1:sizes.bodies));
                
                    density.distances(time_i) = 2 * numedges(G.distances) / ...
                        (numnodes(G.distances) * (numnodes(G.distances) - 1));
                
                
                    %% --- FORCE GRAPH ---
                    A.forces = abs(pairs.forces_time_All(:,:,time_i));
                    G.forces = graph(A.forces, string(1:sizes.bodies));
                
                    density.forces(time_i) = 2 * numedges(G.forces) / ...
                        (numnodes(G.forces) * (numnodes(G.forces) - 1));
                
                
                    %% --- GRAVITY GRAPH ---
                    A.gravity = zeros(sizes.bodies);
                
                    for i = 1:sizes.bodies
                        for j = i+1:sizes.bodies
                            dist_ij = norm(positions(i,:) - positions(j,:));
                
                            % Contact condition (percolation-like threshold)
                            if dist_ij <= radius(i) + radius(j)
                                Contact(i,j) = 1;
                                Contact(j,i) = 1;
                            end
                
                            % Gravitational interaction
                            A.gravity(i,j) = Grav * Mass.all(i) * Mass.all(j) / dist_ij^2;
                            A.gravity(j,i) = A.gravity(i,j);
                        end
                    end
                
                    G.gravity = graph(A.gravity, string(1:sizes.bodies));
                
                    density.gravity(time_i) = 2 * numedges(G.gravity) / ...
                        (numnodes(G.gravity) * (numnodes(G.gravity) - 1));
                
                    toc
                
                
                    %% --- ANGULAR VELOCITY GRAPH ---
                    A.omega = zeros(sizes.bodies);
                
                    for i = 1:sizes.bodies
                        for j = i+1:sizes.bodies
                            A.omega(i,j) = norm(omega(i,:) - omega(j,:));
                            A.omega(j,i) = A.omega(i,j);
                        end
                    end
                
                    G.omega = graph(A.omega, string(1:sizes.bodies));
                
                
                    %% --- ENERGY CALCULATIONS ---
                    Mass.total = sum(Mass.all);
                
                    r_com        = sum(Mass.all .* positions) / Mass.total;
                    velocity_com = sum(Mass.all .* velocities) / Mass.total;
                
                    L = zeros(sizes.bodies,3);
                
                    % Angular momentum per particle
                    for i = 1:sizes.bodies
                        L(i,:) = cross((positions(i,:) - r_com), ...
                                       Mass.all(i)*(velocities(i,:) - velocity_com));
                    end
                
                    % Initialize energies
                    E.c = zeros(sizes.bodies,1); % kinetic
                    E.pe = zeros(sizes.bodies,1); % potential
                    E.r = zeros(sizes.bodies,1); % rotational
                
                    for i = 1:sizes.bodies
                        % Kinetic energy
                        E.c(i) = 0.5 * Mass.all(i) * norm(velocities(i,:))^2;
                
                        % Gravitational potential
                        for j = 1:sizes.bodies
                            if i ~= j
                                E.pe(i) = E.pe(i) - Grav * Mass.all(i)*Mass.all(j) / ...
                                          norm(positions(i,:) - positions(j,:));
                            end
                        end
                
                        % Self potential
                        U_ii = -3/5 * Grav * Mass.all(i) / radius(i);
                
                        % Rotational energy
                        E.r(i) = 0.5 * omega(i,:) * diag(inertia(i,:)) * omega(i,:)';
                
                        % Total energy
                        E.t(i) = E.c(i) + E.r(i) + U_ii + E.pe(i);
                
                        % Escape energy
                        E.escape(i) = Grav * Mass.total / norm(positions(i,:));
                    end
                
                    % Energy statistics
                    E.total(time_i) = sum(E.t);
                    E_time(:,time_i) = E.t';
                
                
                    %% --- ENERGY-BASED GRAPHS ---
                    % Energy similarity (used for percolation-like clustering)
                
                    sigma = median(pdist(E.t'));
                    A.energies = exp(-(E.t' - E.t).^2 / (2*sigma^2));
                    A.energies(1:sizes.bodies+1:end) = 0;
                
                    G.energies = graph(abs(A.energies), string(1:sizes.bodies));
                
                    % Potential difference graph
                    A.potentials = abs(E.pe' - E.pe);
                    G.potentials = graph(A.potentials, string(1:sizes.bodies));
                
                    % Kinetic difference graph
                    A.kinetics = abs(E.c' - E.c);
                    G.kinetics = graph(A.kinetics, string(1:sizes.bodies));
                
                
                    %% --- ENTROPY (NETWORK COMPLEXITY) ---
                    H.pairs(time_i,r_i,d_i)    = calculate_entropy(A.pairs,"normal");
                    H.force(time_i,r_i,d_i)    = calculate_entropy(A.forces,"normal");
                    H.energies(time_i,r_i,d_i) = calculate_entropy(A.energies,"normal");
                
                    H.pairs_lambda(time_i,r_i,d_i)    = calculate_entropy(A.pairs,"random walker");
                    H.energies_lambda(time_i,r_i,d_i) = calculate_entropy(A.energies,"random walker");
                
                    H.energies_kolmogorov(time_i,r_i,d_i) = ...
                        calculate_entropy(A.energies,"kolmogorov");
                
                
                    %% --- PERCOLATION ANALYSIS ---
                    % Degree distribution and critical thresholds
                
                    k = sum(A.pairs,2);
                
                    [counts, edges] = histcounts(k, 'Normalization','probability');
                    k_vals = edges(1:end-1);
                
                    mean_degree   = mean(k);
                    second_moment = sum(k_vals.^2 .* counts);
                
                    % Clustering coefficient
                    num = trace(A.pairs^3);
                    den = sum(0.5 * k .* (k-1));
                
                    C_global(time_i,d_i,r_i) = num / den;
                
                    % Critical probabilities
                    pc(time_i,d_i,r_i)   = (mean_degree) / (second_moment - mean_degree) ...
                                          / (1 - C_global(time_i,d_i,r_i));
                
                    pc_c(time_i,d_i,r_i) = (mean_degree) / (second_moment - mean_degree);
                
                    pc_pois(time_i,d_i,r_i) = 1 / mean_degree;
                
                    % Generating function approach
                    p_c(time_i,d_i,r_i) = generatingfunction(G.pairs);
                
                    if C_global(time_i,d_i,r_i) <= 0.1
                        pc_active(time_i,d_i,r_i) = ...
                            (1/(1-C_global(time_i,d_i,r_i))) * generatingfunction(G.pairs);
                    else
                        pc_active(time_i,d_i,r_i) = generatingfunction(G.pairs);
                    end
                
                    toc
                
                    %% --- STORE DATA ---
                    A_all_times{time_i} = A;
                    E.pe_time(time_i,d_i) = sum(E.pe);
                
                end
            end
    end
end


%% Plot
markerstyle = {"o", "x", "square", "v", "pentagram"};

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(pc_active,2)
        pc_save(:,i) = pc_active(:,i,i);
        pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Percolation Threshold $p_c$[-]','Interpreter','latex','FontSize',18)
ylim([0.15 1.65])
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});


 %% Plot clusterring coefficient global

 figure() 
hold on;
if (numel(densitiesall) > 1)
for i = 1:size(pc_active,2)
    C_global_save(:,i)  = C_global(:,i,i);
    pp(i) = plot(time_vec,C_global_save(:,i),'LineWidth',3,"LineStyle", linestyle{i});
end
else
for i = 1:size(pc,3)
    C_global_save(:,i)  = C_global(:,i,i);
    pp(i) = plot(time_vec,C_global_save(:,i),'LineWidth',3);
end
end

hold off;


xlabel('Time [h]','Interpreter','latex','FontSize',16)
ylabel('Global Clustering Coefficient $c$[-]','Interpreter','latex','FontSize',16)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});
if (numel(densitiesall) > 1)
 legend([pp(1), pp(2), pp(3), pp(4), pp(5)], ...
        {"Case A", "Case B", "Case C", "Case D", "Case E"}, ...
        'Interpreter', 'latex','FontSize',18);
 %       {'$\rho = 1200 \ \mathrm{kg/m^3}$', '$\rho = 1600 \ \mathrm{kg/m^3}$', '$\rho = 2000 \ \mathrm{kg/m^3}$', '$\rho = 2400 \ \mathrm{kg/m^3}$', '$\rho = 2800 \ \mathrm{kg/m^3}$', '$\rho = 3000 \ \mathrm{kg/m^3}$'}, ...
        
else
 legend([pp(2), pp(3), pp(4), pp(1), pp(5), pp(6)],...a
    {"$P = 21.8$ hrs", "$P = 5$ hrs","$P = 3.2$ hrs", "$P = 2.5$ hrs", "$P = 2.1$ hrs", "$P = 1.9$ hrs" }, 'Interpreter','latex','FontSize',18);
end
 

%% Plot Pcs density cases
for i = 1:5
    pcc(:,i) = pc(:,i,i);
end
figure()
plot(time_vec, pcc,'LineWidth',3)
xlabel('Time[H]','Interpreter','latex','FontSize',16)
ylabel('Percolation Threshold $p_c$[-]','Interpreter','latex','FontSize',16)
xticks([0.25,24,48,72,96])
xticklabels({'0.25','2','4','6','8'});
legend({'$\rho = 2000 \ \mathrm{kg/m^3}$','$\rho = 1200 \ \mathrm{kg/m^3}$','$\rho = 1600 \ \mathrm{kg/m^3}$','$\rho = 2400 \ \mathrm{kg/m^3}$','$\rho = 2800 \ \mathrm{kg/m^3}$'},'Interpreter','latex');
%% Plot the averages
figure()
hold on;
if (numel(densitiesall) > 1)
for i = 1:size(avg_degree.pairs,2)
    pp(i) = plot(time_vec, avg_degree.pairs(:,i,i),'LineWidth',3,"LineStyle", linestyle{i});
end
else
for i = 1:size(avg_degree.pairs,3)
    pp(i) = plot(time_vec, avg_degree.pairs(:,i,i),'LineWidth',3,"LineStyle", linestyle{i});
end
end

hold off;
xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Average Degree of the system [-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});

%% Plot Graph Entryopy 
figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(H.pairs,2)
        H.pairs2(:,i) = H.pairs_lambda(:,i,i);
        pp(i) = plot(time_vec, H.pairs2(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec, pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Pairs $h_P$[-]','Interpreter','latex','FontSize',18)
xticks([3,24,48,72,96])
xticklabels({'0.25','2','4','6','8'});

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(H.energies,2)
        H.energies2(:,i) = H.energies_lambda(:,i,i);
        pp(i) = plot(time_vec,H.energies2(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;


xlabel('Time [h]','Interpreter','latex','FontSize',24)
ylabel('Entropy Energy $h_E$[-]','Interpreter','latex','FontSize',24)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 2:size(H.pairs,2)
        H.pairs5(:,i) = H.pairs_von(:,i,i);
        pp(i) = plot(time_vec,H.pairs5(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 2:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Pairs $h_P$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});

%%
figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 2:size(H.force,2)
        H.force5(:,i) = H.force_von(:,i,i);
        pp(i) = plot(time_vec,H.force5(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 2:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Force $h_F$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(H.energies,2)
        H.energies5(:,i) = H.energies_von(:,i,i);
        pp(i) = plot(time_vec,H.energies5(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;

xlabel('Time [h]','Interpreter','latex','FontSize',24)
ylabel('Entropy Energy $h_E$[-]','Interpreter','latex','FontSize',24)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(H.energies,2)
        H.distance(:,i) = H.distance(:,i,i);
        pp(i) = plot(time_vec,H.distance(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;


xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Distance $h_D$[-]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});

%% Kolmogorov

figure()
hold on;
if (numel(densitiesall) > 1)
    for i = 1:size(H.energies_kolmogorov,2)
        H.energies5(:,i) = H.energies_kolmogorov(:,i,i);
        pp(i) = plot(time_vec, H.energies5(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
else
    for i = 1:size(pc_active,3)
     pc_save(:,i) = pc_active(:,:,i);
     pp(i) = plot(time_vec,pc_save(:,i),'LineWidth',3, "LineStyle", linestyle{i});
    end
end
hold off;

xlabel('Time [h]','Interpreter','latex','FontSize',18)
ylabel('Entropy Energies $h_E$[nats]','Interpreter','latex','FontSize',18)
xticks([0,24,48,72,96])
xticklabels({'0','2','4','6','8'});


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value = readdata(path, field)
% READDATA Reads the wanted parameter from GRAIN file.
%
% INPUTS:
%   path  - String or char array representing the path to the text file.
%   field - String or char array representing the field name to search for.
%
% OUTPUT:
%   value - Numeric value corresponding to the specified field. The function
%           scans each line of the file, and when a line starts with the given
%           field name, it extracts and converts the remaining part into a number.
%
% DESCRIPTION:
%   This function reads a file line by line, searches for a line that begins
%   with a given field name, and extracts the numeric value that follows it.

    data = readlines(path);              % Read all lines from the file into a string array
    lengthoffield = length(char(field)); % Determine the length of the field name
    
    for i = 1:length(data)               % Loop through each line in the file
        data_char = char(data(i));       % Convert the current line to a character array
        
        if (length(data_char) > lengthoffield) % Ensure the line is long enough to contain the field
            % Check if the beginning of the line matches the field name
            if (data_char(1:lengthoffield) == char(field))
                % Extract the substring after the field and convert it to a number
                value = str2num(data_char(lengthoffield+1:end));
                break;                  % Stop once the field is found
            end
        end
    end
end


function [entropy] = calculate_entropy(A,entropy_case)
% CALCULATE_ENTROPY Computes different types of entropy for a graph.
%
% INPUTS:
%   A            - Adjacency matrix representing the graph.
%   entropy_case - String specifying the type of entropy to compute:
%                  "normal"         : Shannon entropy based on node strengths
%                  "Von Newmann"    : Von Neumann entropy using Laplacian spectrum
%                  "random walker"  : Entropy based on spectral radius
%                  "kolmogorov"     : Kolmogorov-Sinai entropy (MERW-based)
%
% OUTPUT:
%   entropy - The computed entropy value.
%
% DESCRIPTION:
%   This function calculates different entropy measures depending on the
%   selected case. Each method captures different structural or dynamical
%   properties of the graph represented by adjacency matrix A.

    entropy = 0;                        % Initialize entropy value
    
    switch entropy_case                % Choose entropy type
        
        case "normal"
        % Compute probability distribution based on node strengths
        p_i = sum(A,2)/sum(A,'all');   % Normalize row sums
        
        % Compute Shannon entropy while handling p_i = 0 safely
        for i = 1:numel(p_i)
            if (p_i(i) == 0)           % Skip zero probabilities (limit case)
                continue;
            else
                entropy = entropy - p_i(i) * log(p_i(i)); % Shannon entropy formula
            end
        end
        
        case "Von Newmann"
        % Compute graph Laplacian
        L = laplacian(graph(A));
        
        if (trace(L) == 0)             % Handle empty graph case
            entropy = 0;
            return;
        end
        
        rho = L/trace(L);              % Normalize Laplacian (density matrix)
        lambda = eig(rho);             % Compute eigenvalues
        
        % Compute Von Neumann entropy
        for i = 1:numel(lambda)
            if(lambda(i) <= 1e-3)      % Ignore very small eigenvalues
                continue;
            else
                entropy = entropy - lambda(i) * log(lambda(i));   
            end
        end
        
        case "random walker"
            % Entropy based on the largest eigenvalue (spectral radius)
            
            if isempty(find(A, 1))     % Check if matrix has no edges
                entropy = 0;           % No edges implies zero entropy
            else
                lambda_max = eigs(sparse(A), 1); % Largest eigenvalue
                
                % Entropy is logarithm of spectral radius
                entropy = log(lambda_max);
            end
            
        case "kolmogorov"
            % Compute Kolmogorov-Sinai entropy using MERW
            
            % Step 1: Eigen decomposition of adjacency matrix
            [V, D] = eig(A);
            [lambda_max, idx] = max(diag(D)); % Largest eigenvalue
            v = V(:, idx);                   % Corresponding eigenvector
            
            v = abs(v);                      % Ensure positivity (Perron-Frobenius)
        
            % Step 2: Construct transition matrix P
            n = size(A,1);                  % Number of nodes
            P = zeros(n,n);                 % Initialize transition matrix
            
            for i = 1:n
                for j = 1:n
                    if A(i,j) ~= 0          % Only consider existing edges
                        % MERW transition probability
                        P(i,j) = (A(i,j) * v(j)) / (lambda_max * v(i));
                    end
                end
            end
        
            % Step 3: Compute stationary distribution
            [W, D_P] = eig(P');            % Eigen decomposition of transpose
            [~, idx_p] = min(abs(diag(D_P) - 1)); % Eigenvalue closest to 1
            pi_vec = W(:, idx_p);          % Corresponding eigenvector
            pi_vec = pi_vec / sum(pi_vec); % Normalize distribution
            
            % Step 4: Compute entropy
            H_KS = 0;
            for i = 1:n
                if pi_vec(i) > 1e-12       % Ignore negligible probabilities
                    row_entropy = 0;
                    for j = 1:n
                        if P(i,j) > 0
                            % Shannon entropy contribution per transition
                            row_entropy = row_entropy + P(i,j) * log(P(i,j));
                        end
                    end
                    % Weight by stationary probability
                    H_KS = H_KS + pi_vec(i) * row_entropy;
                end
            end
            
            entropy = -H_KS;               % Apply negative sign
        
        % Note: For symmetric A, entropy equals log(lambda_max)
    end
end



function p_c = generatingfunction(G)
% GENERATINGFUNCTION Computes the critical probability using generating functions.
%
% INPUT:
%   G - Graph object.
%
% OUTPUT:
%   p_c - Critical probability derived from the generating function of the
%         degree distribution.
%
% DESCRIPTION:
%   This function computes the generating function of the degree distribution
%   of a graph and uses its derivative to estimate the critical threshold p_c,
%   commonly used in percolation theory.

    d = degree(G);                      % Get degree of each node
    k_values = min(d):max(d);           % Possible degree values

    Pk = histcounts(d, 'Normalization', 'probability'); % Degree distribution

    syms x;                             % Define symbolic variable
    g0 = sum(Pk .* x.^k_values);        % Generating function G0(x)
    
    g0_p = diff(g0,x);                  % First derivative of G0
    g1 = (g0_p/subs(g0_p,x,1));         % Normalize to obtain G1
    
    if (g0_p == 0)                      % Handle degenerate case
        g1 = 0;
    end
    
    g1_p = matlabFunction(diff(g1));    % Convert derivative to MATLAB function
    p_c = 1/g1_p(1);                    % Critical probability = 1 / G1'(1)

end

