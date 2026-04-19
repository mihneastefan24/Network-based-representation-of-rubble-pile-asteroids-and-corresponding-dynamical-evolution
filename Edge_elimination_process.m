%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Martin M.S.
%% 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code simulates the asteroid break-up through percolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

%% Load WorkSpace

format_ast = { "Kleopatra.obj\"};
core = {"spherical_large_no_core_fixed\"}; %,"spherical_large_no_core_fixed\"};
density_folder = "rho 2000\";
rotation_folder = "Rotation Period 2.520\";
density_folder = char(density_folder);
rotation_folder = char(rotation_folder);
formatast = char(format_ast{1});

% load(append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/",core{1},formatast(1:end-5),"_",density_folder(5:8),"_",rotation_folder(end-5:end-1),".mat"));
load("C:\Users\mihne\Documents\GitHub\Chrono_Projects\files\spherical_large_core_fixed\Bennu_data.mat")

% Set the timestep at which you want to consider the matrix tstep = 1;
G.pairs = graph(pairs.time_All(:,:,1));

%% Load the data
N = sizes.bodies;
placement = transpose(reshape(timeAll.positions(time_vec(1),:),[3,sizes.bodies]));
%%

positions_rand = rand(N,3); 
steps = 1000; % Create a monte carlo of 1000
D = 0.01; % What is D?

trajectory = zeros(N,3,steps);
trajectory(:,:,1) = positions_rand;
for t = 2:steps
    trajectory(:,:,t) = trajectory(:,:,t-1) + D * randn(N,3);
end
% What is the radius?
radius = 0.5;

edge_counts = zeros(N);
for t =1:steps
    pos = trajectory(:,:,t);
    for i = 1:N
        for j = i+1 :N
            if(norm(pos(i,:) - pos(j,:)) < radius)
                edge_counts(i,j,t) = edge_counts(i,j) + 1;
                edge_counts(j,i,t) = edge_counts(j,i) + 1;
            end
        end
    end
end

edge_probs = rand([sizes.bodies, sizes.bodies]);

% Extract the edge list as numeric indeices
edgeList = G.pairs.Edges.EndNodes;

% Assuming edge probs is a matix of probabilities for each edge (i,j)
% Get probabilities for the corresponding edges
edge_probs_vals = edge_probs(sub2ind(size(edge_probs), edgeList(:,1), edgeList(:,2)));

% Generate random values and keep edges where rand  < prob
keep = rand(size(edgeList,1),1) < edge_probs_vals;

% Remove edges that are not kept
G_perc = rmedge(G.pairs, edgeList(~keep,1), edgeList(~keep,2));

% Find conncted components in the resulting graph
bins = conncomp(G_perc);

num_fragments = max(bins);
fragment_sizes = histcounts(bins,1:num_fragments + 1);

figure()
colormap("cool");
plot(G_perc, "XData", positions(:,1), "YData",positions(:,2), "ZData", positions(:,3), ...
    'NodeCData',bins, 'MarkerSize',6);
title(['Fragments: ', num2str(num_fragments)]);
colorbar;

%% Chanage threshold based on your needs
threshold = 2;
[i_add, j_add] = find(triu(scores,1) > threshold);
H = addedge(G.pairs, i_add, j_add, ones(length(i_add),1));
A.pairs_update = adjacency(H);

figure()
hold on;
scatter3(positions(:,1),positions(:,2),positions(:,3),'filled')
edges = double(string(G.pairs.Edges.EndNodes));
for i = 1:size(edges,1)
    x = [positions(edges(i,1),1) positions(edges(i,2),1)];
    y = [positions(edges(i,1),2) positions(edges(i,2),2)];
    z = [positions(edges(i,1),3) positions(edges(i,2),3)];
    plot3(x,y,z,'Color',[0 0 0 + 0.1])
end

for  i = 1:size(i_add,1)
    x = [positions(i_add(i),1) positions(j_add(i),1)];
    y = [positions(i_add(i),2) positions(j_add(i),2)];
    z = [positions(i_add(i),3) positions(j_add(i),3)];
    plot3(x,y,z,'Color',[1 0 0])
end
hold off;

%% Check the centre of the graph - the points that corresponds to the most centred region of the graph
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 


% Calcualte the closest node to the actual centre of mass of the system
Mass.total = sum(Mass.all);
r_com = sum(Mass.all .* placement) / Mass.total;
[~,idx_centre] = min(vecnorm(placement - r_com,2,2));

% Calculate the points with highest centrality(type is chosen by you), this
% are the points most centered in the graph
closeCent = centrality(G.distances, 'closeness', 'Cost',G.distances.Edges.Weight);
[~, centerNode] = max(closeCent);
% Move on from the centre of mass to the other points not the so called
% centre of nodes

numbernodes = neighbors(G.forces, idx_centre);

% Create a matrix that represents the stress imposed to the edges

A.stress = Adj{1}.pairs .* Adj{1}.acc_forces - Adj{1}.forces;
G.stress = graph(A.stress);
% Continue the percolation process
% Continue through by doing also MCMC - Monte Carlo Markov Chain
% division contains certain degrees of simplification based on maximum
% length reduction in the system - line 141
% division = [25, 20, 15, 10];
division = 15;
edgeCounts_structure = cell(1,4);
time_parallel = zeros(1,4);
for i = 1:numel(division)
    clear A.stress G.stress;
    A.distance_min = A_all_times{1}.distance < max(max(A_all_times{1}.distance))/division(i);
    A.stress = (0.5 .* A_all_times{1}.gravity / max(A_all_times{1}.gravity,[],'all') + +0.5 .* A_all_times{1}.forces/max(A_all_times{1}.forces,[],'all')) .* A.distance_min;
    G.stress = graph(A.stress);
    tic
    [edgeCounts_structure{i},vec1]= simulate_edge_breaks1(G.stress, 500,A.stress);
    time_parallel(i) = toc;
end
%%

G.contact = graph(A_all_times{1}.forces .* 0.5 + 0.5 .* A_all_times{1}.gravity);

Mass.total = sum(Mass.all);
r_com = sum(Mass.all .* placement) / Mass.total;
[~,idx_centre] = min(vecnorm(placement - r_com,2,2));
[edgeCounts22,vec]= simulate_edge_breaks(G.contact, 1000, idx_centre);
limit = round(0.10 * length(edgeCounts22));
kmax_val = maxk(edgeCounts22,limit);

idx_EdgeCounts = find(edgeCounts22 < kmax_val(end));

G_new = rmedge(G.contact,idx_EdgeCounts);

positions = transpose(reshape(timeAll.positions(time_vec(1),:),[3,sizes.bodies]));
% Take the clusters of the final state of the system for comparison
clustervalue1 = clusterpairs(A.pairs, zeros(0, numedges(G_new)));

% Calculate the clusters for the matrix obtained from the percolation
A_new = full(adjacency(G_new));
clustervalue = clusterpairs(A_new, zeros(0, numedges(G_new)));

%%
% Plot the new clusters in the first time step with the new clusters
figure()
    hold on;
    grid on;
    num_clusters = length(clustervalue);
    colors = hsv(num_clusters);
    for i = 1:length(clustervalue)
         if length(clustervalue{i}) > 1
            scatter3(positions(clustervalue{i},1),positions(clustervalue{i},2),positions(clustervalue{i},3),'filled','MarkerFaceColor',colors(i,:))
            % Plot the edges
            for j = 1:length(clustervalue{i})
                pairs_j = find(A_new(clustervalue{i}(j),:));
                for pairs_i = 1:length(pairs_j)
                    x = [positions(clustervalue{i}(j),1), positions(pairs_j(pairs_i),1)];
                    y = [positions(clustervalue{i}(j),2), positions(pairs_j(pairs_i),2)];
                    z = [positions(clustervalue{i}(j),3), positions(pairs_j(pairs_i),3)];
    
                    plot3(x, y, z, 'Color',[0 0 0]+0.1);  % 'k' = black line
                end
            end
        else
            body = clustervalue{i};
            scatter3(positions(body,1),positions(body,2),positions(body,3),'filled','MarkerFaceColor',hex_colors{end})
        end
    end
    hold off;


limit1 = round(0.50 * length(edgeCounts1));
kmax_val1 = maxk(edgeCounts1,limit1);

idx_EdgeCounts1 = find(edgeCounts1 > kmax_val1(end));

G_stress_new = rmedge(G.stress,idx_EdgeCounts1);
positions = transpose(reshape(timeAll.positions(time_vec(1),:),[3,sizes.bodies]));

% Calculate the clusters for the matrix obtained from the percolation
A_stress_new = full(adjacency(G_stress_new));
clustervalue_stress = clusterpairs(A_stress_new, zeros(0, numedges(G_stress_new)));


% Plot the new clusters in the first time step with the new clusters
figure()
    hold on;
    grid on;
    num_clusters = length(clustervalue_stress);
    colors = hsv(num_clusters);
    for i = 1:length(clustervalue_stress)
         if length(clustervalue_stress{i}) > 1
            scatter3(positions(clustervalue_stress{i},1),positions(clustervalue_stress{i},2),positions(clustervalue_stress{i},3),'filled','MarkerFaceColor',colors(i,:))
            % Plot the edges
            for j = 1:length(clustervalue_stress{i})
                pairs_j = find(A_stress_new(clustervalue_stress{i}(j),:));
                for pairs_i = 1:length(pairs_j)
                    x = [positions(clustervalue_stress{i}(j),1), positions(pairs_j(pairs_i),1)];
                    y = [positions(clustervalue_stress{i}(j),2), positions(pairs_j(pairs_i),2)];
                    z = [positions(clustervalue_stress{i}(j),3), positions(pairs_j(pairs_i),3)];
    
                    plot3(x, y, z, 'Color',[0 0 0]+0.1);  % 'k' = black line
                end
            end
        else
            body = clustervalue_stress{i};
            scatter3(positions(body,1),positions(body,2),positions(body,3),'filled','MarkerFaceColor',hex_colors{end})
        end
    end
    hold off;

%% Make the comparison
positions = transpose(reshape(timeAll.positions(time_vec(4),:),[3,sizes.bodies])); 

figure()
    hold on;
    grid on;
    num_clusters = length(clustervalue);
    colors = hsv(num_clusters);
    for i = 1:length(clustervalue)
         if length(clustervalue{i}) > 1
            scatter3(positions(clustervalue{i},1),positions(clustervalue{i},2),positions(clustervalue{i},3),'filled','MarkerFaceColor',colors(i,:))
            % Plot the edges
            for j = 1:length(clustervalue{i})
                pairs_j = find(A_new(clustervalue{i}(j),:));
                for pairs_i = 1:length(pairs_j)
                    x = [positions(clustervalue{i}(j),1), positions(pairs_j(pairs_i),1)];
                    y = [positions(clustervalue{i}(j),2), positions(pairs_j(pairs_i),2)];
                    z = [positions(clustervalue{i}(j),3), positions(pairs_j(pairs_i),3)];
    
                    plot3(x, y, z, 'Color',[0 0 0]+0.1);  % 'k' = black line
                end
            end
        else
            body = clustervalue{i};
            scatter3(positions(body,1),positions(body,2),positions(body,3),'filled','MarkerFaceColor',hex_colors{end})
        end
    end
    title("Pairs Cluster mine")


%%                          Functions                              %%
%%%---------------------------------------------------------------%%%
%%%---------------------------------------------------------------%%%
%%%---------------------------------------------------------------%%%
%%%---------------------------------------------------------------%%%

function clusters = clusterpairs(A, delta)
% CLUSTERPAIRS Identifies connected clusters (components) in a graph.
%
% INPUTS:
%   A     - Adjacency matrix (NxN), where A(i,j) = 1 indicates an edge
%           between node i and node j.
%   delta - (Unused parameter in current implementation; reserved for
%           thresholding or percolation criteria).
%
% OUTPUT:
%   clusters - Cell array where each cell contains a vector of node indices
%              belonging to the same connected cluster.
%
% DESCRIPTION:
%   This function detects clusters (connected components) in a network
%   represented by adjacency matrix A. It iteratively groups nodes that
%   are directly or indirectly connected, mimicking cluster formation
%   in percolation processes.
%
% STEPS:
%   1. Iterate through each node and identify its direct neighbors.
%   2. Merge nodes into clusters if they share connections.
%   3. Merge overlapping clusters (transitive closure).
%   4. Remove empty clusters.

    clusters = {};  % Initialize empty list of clusters
    
    % --- Step 1: Build initial clusters from adjacency ---
    for i = 1:size(A, 1)
        contacts = find(A(i, :) == 1);   % Find neighbors of node i
        contacts = [i, contacts];        % Include node itself
        
        found = false; % Flag to check if cluster already exists
        
        % --- Step 2: Merge with existing clusters if overlap exists ---
        for j = 1:length(clusters)
            if any(ismember(contacts, clusters{j}))
                % Merge current contacts into existing cluster
                clusters{j} = unique([clusters{j}, contacts]);
                found = true;
                break;
            end
        end
        
        % --- Step 3: Create new cluster if no overlap ---
        if ~found
            clusters{end+1} = contacts;
        end
    end

    % --- Step 4: Merge overlapping clusters (transitive connections) ---
    for clu_i = length(clusters):-1:2
        for clu_j = clu_i-1 :-1:1
            if any(ismember(clusters{clu_i}, clusters{clu_j}))
                clusters{clu_j} = unique([clusters{clu_j}, clusters{clu_i}]);
                clusters{clu_i} = []; % Mark for deletion
            end
        end
    end
    
    % --- Step 5: Remove empty clusters ---
    for i = length(clusters):-1:1
        if (numel(clusters{i}) == 0)
            clusters(i) = [];
        end
    end
end
%% 
function [edgeCounts, vecAll] = simulate_edge_breaks(G, numSim, nodes)
% SIMULATE_EDGE_BREAKS Simulates edge removal via biased random walks.
%
% INPUTS:
%   G       - Graph object with weighted edges.
%   numSim  - Number of independent simulations (Monte Carlo runs).
%   nodes   - (Currently unused) optional parameter for node selection.
%
% OUTPUTS:
%   edgeCounts - Vector (size = number of edges) counting how many times
%                each edge was removed across all simulations.
%   vecAll     - Cell array storing visited nodes for each simulation.
%
% DESCRIPTION:
%   This function models a percolation-like process where edges are removed
%   sequentially via a biased random walk. At each step:
%     - The walker moves probabilistically based on edge weights.
%     - The traversed edge is removed.
%   This mimics progressive network fragmentation.
%
% STEPS:
%   1. For each simulation, initialize graph and starting node.
%   2. Perform a weighted random walk.
%   3. Remove edges as they are traversed.
%   4. Stop when no valid moves remain or all nodes are visited.
%   5. Aggregate edge removal statistics.

    vecAll = cell(numSim,1);           % Store visited nodes per simulation
    removeEdgesAll = cell(numSim,1);  
    edgeCounts = zeros(numedges(G),1); % Global edge removal counts

    parfor sim = 1:numSim
        node = 351;                    % Fixed starting node
        disp(sim);
        
        Gtemp = G;                     % Copy graph
        visited = false(numnodes(G), 1);
        SaveRemoveEdges = [];
        vec = [];                      % Track visited nodes
        
        edgeCount = zeros(numedges(G),1); % Local edge removal counter

        % --- Random walk with edge removal ---
        while true
            neighborsList = neighbors(Gtemp, node);
            
            % Only allow moves to nodes with degree > 1 (avoid dead ends)
            eligible = neighborsList(degree(Gtemp, neighborsList) > 1);

            if(isempty(eligible))
                break; % Stop if no valid moves
            else
                % --- Compute transition probabilities ---
                edgeIdxs = findedge(Gtemp, node * ones(size(eligible)), eligible);
                weights = Gtemp.Edges.Weight(edgeIdxs);
                probs = weights / sum(weights); % Normalize
                
                % --- Choose next node ---
                if(isscalar(eligible))
                    nextNode = eligible;
                else
                    nextNode = randsample(eligible, 1, true, probs);
                end
            end
           
            % --- Remove traversed edge ---
            edgeToRemove = findedge(Gtemp, node, nextNode);
            edgeCount(edgeToRemove) = edgeCount(edgeToRemove) + 1;
            
            Gtemp = rmedge(Gtemp, node, nextNode); % Edge deletion
            SaveRemoveEdges = [SaveRemoveEdges; edgeToRemove];

            % --- Move walker ---
            nodeOld = node;
            node = nextNode;
            
            vec = unique([vec,nextNode]); % Track visited nodes
            
            % Stop if all nodes have been visited
            if(numel(vec) == numnodes(G))
                break;
            end
        end
        
        % --- Aggregate results ---
        edgeCounts = edgeCounts + edgeCount;
    end
end

function [edgeCounts, vecAll] = simulate_edge_breaks1(G, numSim, A)
% SIMULATE_EDGE_BREAKS1 Simulates biased edge removal using structural + stochastic rules.
%
% INPUTS:
%   G       - Graph object (network structure).
%   numSim  - Number of simulations.
%   A       - Adjacency matrix used to bias transitions.
%
% OUTPUTS:
%   edgeCounts - Vector counting how often each edge is removed.
%   vecAll     - Cell array storing visited nodes per simulation.
%
% DESCRIPTION:
%   This function extends the previous random walk model by introducing
%   a *biased transition mechanism* combining:
%       - Structural weights (adjacency matrix A)
%       - Exposure (node degree influence)
%       - Random noise
%
%   This models *biased percolation*, where edge failure is influenced
%   by both topology and stochastic perturbations.
%
% STEPS:
%   1. Initialize random starting node.
%   2. Compute transition probabilities using weighted combination.
%   3. Perform random walk and remove edges.
%   4. Track edge removals and visited nodes.
%   5. Aggregate results over simulations.

    edgeCounts = zeros(numedges(G),1);
    vecAll = cell(numSim,1);
    removeEdgesAll = cell(numSim,1);  
    
    parfor sim = 1:numSim
        disp(sim);

        Gtemp = G;
        visited = false(numnodes(G), 1);
        SaveRemoveEdges = [];
        vec = [];
        
        node = randi(numnodes(G)); % Random starting node
        edgeCount = zeros(numedges(G),1);

        while true
            neighborsList = neighbors(Gtemp, node);
            eligible = neighborsList(degree(Gtemp, neighborsList) > 1);

            if(isempty(eligible))
                break;
            else
                % --- Base weights ---
                edgeIdxs = findedge(Gtemp, node * ones(size(eligible)), eligible);
                weights = Gtemp.Edges.Weight(edgeIdxs);

                % --- Bias parameters ---
                alpha = 0.9; beta = 0.0; gamma = 0.1;
                
                % Structural influence (from adjacency matrix)
                structural = abs(A(node, eligible))' / abs(sum(A(node,eligible')));
                
                % Exposure term (penalizes high-degree nodes)
                exposure = 1 ./ (degree(G,node) + degree(G, eligible));
                exposure = exposure / sum(exposure);
                
                % Random noise
                noise = rand(length(eligible),1);
                noise = noise / sum(noise);

                % --- Combined transition probability ---
                probs = alpha * structural + beta * exposure + gamma * noise;

                % --- Choose next node ---
                if(isscalar(eligible))
                    nextNode = eligible;
                else
                    nextNode = randsample(eligible, 1, true, probs);
                end
            end
           
            % --- Remove edge ---
            edgeToRemove = findedge(Gtemp, node, nextNode);
            edgeCount(edgeToRemove) = edgeCount(edgeToRemove) + 1;

            Gtemp = rmedge(Gtemp, node, nextNode);
            SaveRemoveEdges = [SaveRemoveEdges; edgeToRemove];

            % --- Move walker ---
            nodeOld = node;
            node = nextNode;
            
            vec = unique([vec,nextNode]);

            if(numel(vec) == numnodes(G))
                break;
            end
        end

        % --- Aggregate results ---
        edgeCounts = edgeCounts + edgeCount;
    end
end
