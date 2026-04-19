%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Martin M.S.
%% 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This code is divided in two parts:
%%%  1) Shape the rubble-pile (GRAINS) into the desire geometry.
%%%  2) Adapt and size the parameters of the rubble pile to fit the initial input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clearvars; clc; close all;
format("shortG")

% This is a Matlab add-on available to extract the points that are present
% inside the obj shape.
addpath(fullfile(getenv('APPDATA'), ...
    'MathWorks\MATLAB Add-Ons\Functions\inpolyhedron - are points inside a triangulated volume_'))

%%%
%%% This part is handling the cutting of the asteroid into shape
%%% 

%% Initialize data for the cutting

% Input the directories and files of the data

particles_format = "spherical_large_no_core_fixed";
text_names = "_Fixed_tstep_small_NSC";
directory_load = append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/",particles_format,"/");
% Location of the .txt files on my laptop, change accordingly
cd(directory_load)
position_txt = append("Position",text_names,".txt"); velocity_txt = append("Velocity",text_names,".txt"); mass_txt = append("Mass",text_names,".txt");
% Extract the positions, velocity and masses
positions = readmatrix(position_txt);
velocity = readmatrix(velocity_txt);
mass = readmatrix(mass_txt);

% Calculate the barycentre
distance_each = (positions(:,4).^2 + positions(:,2).^2 + positions(:,3).^2).^0.5;
distance_mass = (mass(:,2) .* positions(:,2:4));
CoM = sum(distance_mass,1)/sum(mass(:,2));

% Check the particle radius to match the one from C++
rhocpp = 3000;               % density used in the cpp simulation of the first aggregate
volume =[mass(:,1), mass(:,2) ./ rhocpp];
radius_each_large_agg = [volume(:,1), ((3 * volume(1:end,2)) / (4 * pi)).^(1/3)];
radius_all = sum(radius_each_large_agg);

%% Introduce the fomat script to check the overlapping of the particles in chrono, and adapt their size->volume->sf

radius = 0.999 * round(radius_each_large_agg(:,2));

% Compute pairwise distances and check for overlaps
overlapping_pairs = []; % Store indices of overlapping spheres
numSpheres = length(radius);
k = 1;
tic
for i = 1:numSpheres
    for j = i+1:numSpheres
        % Compute distance between sphere centers
        dist_ij = norm(positions(i,2:4) - positions(j,2:4));

        % Check if spheres overlap
        if dist_ij < (radius(i) + radius(j))
                overlapping_pairs(k,:) = [i, j, dist_ij, (radius(i) + radius(j))];  
                k = k + 1;
        end
    end
end


max_overlap = min(overlapping_pairs(:,3)./ overlapping_pairs(:,4));
% Create a radius copy which doesn't change during the for loop
radius_copy = radius; %max_overlap;
% Check the mass change and update mass
mass = rhocpp * 4/3 *pi *radius_each_large_agg.^3;
t2 = toc;

% Asteroid list and properties
fig_name = ["Bennu", "Kleopatra","Geographos","Apophis","Arrokoth"];
filename = ["Bennu_v20_200k.obj","Kleopatra.obj", "Geographos Radar-based, mid-res.obj", "Apophis Model 1.obj","Arrokoth Stern 2019.obj"];  % Replace with your OBJ file path

% Bulk densities cases
bulk_densities = [1200, 1600, 2000, 2400, 2800, 3000];

% Create a vector to match the id and the mass which will not be changed
% inside the for loop
mass_copy = [[1:1e4]',  mass];
% Make a copy of the initial positions which will not be changed
% inside the for loop
positions_copy = positions;

for file_i = 1:length(filename)

    % Initialize data for each file
    mass = mass_copy; 
    positions = positions_copy;
    radius_each = radius_copy;

    %% Obj Cutout Shape
    % Select the kind of the document in order to extract the data change
    % accordingly

    kind = [1 2 2 2 1];
    fid = fopen(filename(file_i), 'r');

    if fid == -1
        error('Cannot open the file.');
    end
            
    C = textscan(fid, '%c %f %f %f');

    % Initialize vector to save the vertices and faces
    vertices = [];
    faces = [];
    % Implement method for fid 2
    if filnr == 2 
        while ~feof(fid)
            line = fgetl(fid);
            if startsWith(line,'v')
                vData = textscan(line, '%c %f %f %f');
                vertices = [vertices; [vData{2} vData{3} vData{4}]];
            elseif startsWith(line,'f')
                fData = textscan(line, '%c %f %f %f');
                faces = [faces; [fData{2} fData{3} fData{4}]];
            end
        end
    % Implement method for fid 1
    elseif filnr == 1 
        data = [C{2} C{3} C{4}];
        faces = data(C{1} == 'f',:);
        vertices = data(C{1} == 'v',:);
    end
    fclose(fid);


    % Check to see if the extract shape is the desired one
    figure(); % figure1
    fig = fig + 1;
    ast = patch('Faces',faces,'Vertices',vertices);
    set(ast,'Edgealpha',0.5);
    set(ast,'EdgeColor','white');
    lighting gouraud;
    axis equal;


        %% 
    
    % Save the scaled mesh
    id = fopen(filename(file_i));
    C=textscan(id,'%c %f %f %f');
    fclose(id);
    data = [C{2} C{3} C{4}];
    face=data(C{1}=='f',:);
    vertex = data(C{1}=='v',:);
    fv.vertices = vertices;
    fv.faces = faces;
    %plot the mesh
    figure(fig); %figure2
    fig = fig + 1;
    ast=patch('Vertices',vertex,'Faces',face);
    set(ast,'Edgealpha',0.5);
    set(ast,'EdgeColor','white');
    lighting gouraud;
    axis equal;
    
    figure(fig); %figure3
    fig = fig + 1;
    scatter3(100 * vertices(:,1), 100 * vertices(:,2), 100 * vertices(:,3), 'filled');
    axis equal;
    grid on;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('3D Points from OBJ');
    fig_title = append("figure_",fig_name(file_i),".png");
    saveas(gcf,fig_title);


   %%  Impose the limits from obj to the positions from chrono but first adapt the size of the obj to that of chrono
      
   % Calculate the distance from the GRAINS reference system for each
   % particle
    for dis_i = 1:size(positions,1)
        distance(dis_i) = norm(positions(dis_i,2:4));
    end
            
    % Eliminate the particles that were not able to aggregate based on
    % percentile
    data = distance; 
    
    % Compute Q1 (25th percentile) and Q3 (75th percentile)
    Q1 = prctile(data, 25);
    Q3 = prctile(data, 75);
            
    % Compute interquartile range
    IQR = Q3 - Q1;
            
    % Define a threshold (typically 1.5 * IQR above Q3 is an outlier)
    upper_bound = Q3 + 1.5 * IQR;
            
    % Remove outliers
    filtered_data = data(data <= upper_bound);
    filtered_data_sorted_up = sort(filtered_data,"ascend");
    filtered_data_sorted_down = sort(filtered_data,"descend");
    distance_data_sorted_up = sort(distance,"ascend");
    distance_data_sorted_down = sort(distance,"descend");

    % Update the data vectors: position, velocity and mass
    [Lia,Locb] = ismember(distance,filtered_data);
    mass = mass(Lia,:);
    positions = positions(Lia,:);
    radius_each = radius_each(Lia,:);
            
    % Update the CoM
    distance_mass = (mass(:,2) .* positions(:,2:4));
    CoM = sum(distance_mass,1)/sum(mass(:,2));
    mass_tot = sum(mass(:,2));

    % Update the positions: Create the positions with respect to the new
    % CoM

    positions_CoM = positions(:,2:4) - CoM;
    positions_vec_1 = positions_CoM(:,1);
    positions_vec_2 = positions_CoM(:,2);
    positions_vec_3 = positions_CoM(:,3);


    % The aggregate is not perfectly spherical, thus a change in the axis
    % position is made to maximize the number of partiles kept during the
    % cutting process
    positions_CoM = [positions_vec_3, positions_vec_2, positions_vec_1];
    %% Update the shape of OBJ to match that of the asteroid
            
    figure(); %figure4
    fig = fig + 1;
    scatter3(positions(:,2),positions(:,3),positions(:,4),'filled');
            
    % Calculate the volume of each particle and the total volume without
    % the porosity
    volume =[mass(:,1), mass(:,2) ./ rhocpp]; % include particle id
    radiuse = [volume(:,1),(3 * volume(:,2)/(4 * pi)).^(1/3)];
    volume_all = sum(volume(:,2));
    M_tot = sum(mass(:,2));
            
    NB.Rb = 10;
    NB.Ras_min = 10;
    NB.Ras_max = 100000000000000;
    Mdensity = 3000; % 3000
    NB.Nb =  size(positions,1);
    NB.Vmat = M_tot/Mdensity;
    NB.Vmat = volume_all;
    
    % Calculate the maximum possible stretch in each direction for the
    % asteroid

    xL0 = abs(max(positions_CoM(:,1))) + abs(min(positions_CoM(:,1)));
    yL0 = abs(max(positions_CoM(:,2))) + abs(min(positions_CoM(:,2)));
    zL0 = abs(max(positions_CoM(:,3))) + abs(min(positions_CoM(:,3)));
    
    % Calculate the maximum possible stretch in each direction for the
    % obj
    obj_xL = abs(max(vertices(:,1))) + abs(min(vertices(:,1)));
    obj_yL = abs(max(vertices(:,2))) + abs(min(vertices(:,2)));
    obj_zL = abs(max(vertices(:,3))) + abs(min(vertices(:,3)));
            
    % Calculate the ratio in each direction and take the minimum value, in
    % this way in one direction the obj and asteroid will have the same
    % length; 
    % A safety factor can be implemented to eliminate the possible loose
    % partiles that are not strongly aggregated
    dim_sf = min(([xL0 yL0 zL0]./[obj_xL obj_yL obj_zL])) .* [1 1 1];      % 0.98 * min([xL0 yL0 zL0]./[obj_xL obj_yL obj_zL]) .* [1 1 1]
            
    %% Check also visually if the two images overlap accordingly
    figure(fig); %figure6
    fig = fig + 1;
    scatter3(positions_CoM(:,1),positions_CoM(:,2),positions_CoM(:,3),'filled')
    hold on;
    scatter3(dim_sf(1) * vertices(:,1),dim_sf(2) * vertices(:,2),dim_sf(3) * vertices(:,3));
            
    % Implement the expansion ration to all vertices
    vertices = vertices .* dim_sf; 
    distance_v = zeros(size(vertices,1),1);

    %% Verify and eliminate the partiles outside the .obj

    % It uses the matlab add-on inpolyhedron, a second function was created
    % to test the performance of inpolyhedron and it is adequate
    fv.vertices = vertices;
    inside = inpolyhedron(fv, positions_CoM);
    figure(fig), hold on, view(3)        % Display the result
    fig = fig + 1;
    patch(fv,'FaceColor','g','FaceAlpha',0.2)
    plot3(positions_CoM(inside,1),positions_CoM(inside,2),positions_CoM(inside,3),'bo','MarkerFaceColor','b')
    hold off;

    clear pos0_i vert_i;
    positions_inside = positions_CoM(inside,:);
    pos_in = positions(inside,2:4);

    % Check to verify the shape
    figure(); %figure7
    scatter3(positions_inside(:,1),positions_inside(:,2),positions_inside(:,3));
    hold on;
    scatter3(vertices(:,1), vertices(:,2), vertices(:,3));
    xlim([-6 6]);
    ylim([-6 6]);
    zlim([-6 6]);
    hold off;
    
    figure() %figure8
    scatter3(positions_inside(:,1),positions_inside(:,2),positions_inside(:,3),'filled');
    fig_title = append("figure_",fig_name(file_i),".png");
    saveas(gcf,fig_title);
    
   
    %% Save the data into txt files
    
    
    save_directory = append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/spherical_large_no_core_fixed/results/",filename(file_i), "/");   
    mkdir(save_directory);

    % Check the id of the particles that are inside the shape
    [~, idpos] = ismember(positions_inside, positions_CoM,'rows');
    masses_inside = [idpos,masses(idpos,:)];
    radius_inside = radius_each(idpos,:); 
    
    fileID_pos = fopen(append(save_directory, "Positions.txt"),'w');
    fprintf(fileID_pos, "Positions x y z \n");
    for pos_i = 1:length(positions_inside)
        fprintf(fileID_pos,"%.6e %.6e %.6e \n",positions_inside(pos_i,:));
    end
    fclose(fileID_pos);
    fileID_mass = fopen(append(save_directory, "Mass.txt"),'w');
    fprintf(fileID_mass, "Masses \n");
    for mass_i = 1:length(masses_inside)
        fprintf(fileID_mass,"%.10e \n",masses_inside(mass_i,2));
    end
    fclose(fileID_mass);
    fileID_radius = fopen(append(save_directory, "Radius.txt"),'w');
    fprintf(fileID_radius, "Radii \n");
    for rad_i = 1:length(radius_inside)
        fprintf(fileID_radius,"%.6e \n",radius_inside(rad_i,1));
    end
    fclose(fileID_radius);

    % Check the number of particles are in each shape
    disp(size(radius_inside));

%% Second part size the physical parameters to be in accordance with the data created so far
%%% 
%%% This part is created by I. Fodde
%%%

    rmatSpec = '%.f';
    
    % Initial setup of config file for spinup(IMPORTANT, SOME VALUES WILL BE CHANGED BY THE PROGRAM)
    spinConfig.G_univ = 6.67430e-11 ;   % [L^3/kg*s^2]   universal gravitational constant (scaled to L)
    spinConfig.velcase = 2;             % [-]           1-no velocity, 2-angular momentum, 3-load
    spinConfig.initialSpin = 7.5e-04; % [rad/s]      initial angular spin of the aggregate (only if velcase=2)
    spinConfig.load_directory = '';
    spinConfig.density = 3000;  % [kg/L^3]  set it to zero to load density info
    spinConfig.densityIC = 3000;      % [kg/L^3] density of inner core (set it to zero to not have a heterogeneous density distribution)
    spinConfig.rIC = 1;        % [L] radius of the inner core (only if densityIC!=0)
    % surface properties %%%%
    spinConfig.restitution = 0;         % [-]   normal restitution coefficient. Default =0.4
    spinConfig.friction = 0.6 ;         % [-]   Set both static friction and kinetic friction at once, with same value.
    spinConfig.adhesion = 0;            % [N]   constant cohesion force. Dafault = 0
    spinConfig.adhesionMult = 0;        % [?]   adhesion multiplier. Default = 0
                            %               DMT adhesion model:
                            %                       adhesion = adhesionMultDMT * sqrt(R_eff).
                            %               Given the surface energy, w,
                            %                       adhesionMultDMT = 2 * CH_C_PI * w * sqrt(R_eff).
                            %               Given the equilibrium penetration distance, y_eq,
                            %               adhesionMultDMT = 4.0 / 3.0 * E_eff * powf(y_eq, 1.5)
    spinConfig.Kn = 200000;             % [?]   normal stiffness coefficient. Default=200000
    spinConfig.Kt = 200000;             % [?]   tangential stiffness coefficient. Default=200000
    spinConfig.Gn = 40;                 % [?]   normal damping coefficient. Default=40
    spinConfig.Gt = 20;                 % [?]   tangential damping coefficient. Default=20
    % numerical integration
    spinConfig.maxT = 8;             % [h]              maximum simulation time
    spinConfig.IntTimeStep = 1;      % [s]            integration time step -> to be tuned according to collision characteristic time (depends on relative velocities), can be tuned by looking at totH
    spinConfig.deltaSpin = 0;         % [rad/s]       amount of spin increase every 'spinupT' time steps
    spinConfig.spinupT = 0;             % [-]           perform spinup every 'spinupT' time steps (set spinupT=0 to not spinup the aggregate)
    spinConfig.binX = 1;                % [-]   number of bins for parallel simulation in x direction (such to have 2-4 objects per bin per direction)
    spinConfig.binY = 1;               % [-]   number of bins for parallel simulation in y direction (such to have 2-4 objects per bin per direction)
    spinConfig.binZ = 1;                % [-]   number of bins for parallel simulation in z direction (such to have 2-4 objects per bin per direction)
    spinConfig.minThreads = 20;         % [-]   number of parallel threads to run the simulations
    % force computation
    spinConfig.gravT = 1;               % [-]           compute gravity every 'gravT' time steps (set gravT=0 to not compute gravity)
    spinConfig.doparallel = false;      % [bool]        do computation of forces using BH-GPU (true) or N2 method (false)
    % save results
    spinConfig.SaveResults = 3;         % [-]           save results: 0-no, 1-subdirectories, 2-one file per body, 3-single file
    spinConfig.SaveT = 1;             % [-] save results every 'SaveT' time steps for all bodies (set SaveT=0 to save no results in time)
    spinConfig.SaveVideo = 1;         % [-] save more data to reproduce video: one system file each 'SaveVideo' time steps (set SaveVideo=0 to save no video data)
    spinConfig.save_directory = '';
    
    angular_velocities = [2e-4, 3.5e-4, 5.5e-4 , 7e-4, 8.5e-4, 9e-4];
    time_vel = angular_velocities * 3600;
    % Get the volumes of each asteroid or at least the mass
    volumes = 608974.358974; % Corresponds to dataset volume
    
    real_volume = 608974.358974;

    masses = masses_inside(:,2:3); mass = masses_inside(:,2:3);
    positions = positions_inside; radius_each = radius_inside(:,1);
    for rho_i = 1:length(bulk_densities)
        
        real_bulk_density = bulk_densities(rho_i);
        
        for w_i = 1:length(angular_velocities)
            angular_velocity = angular_velocities(w_i);
            real_period = (2*pi/angular_velocity)/3600;
                    
            NB.prorositycom = 1;

            % CHANGE: fill in here the directory from which it loads the shape and where it stores the results
            spinConfigShort = spinConfig;
            % The shape is defined below
            spinConfigShort.loadf_directory = append("../" , filename(file_i),"/");
            spinConfigShort.save_directory = append("../", filename(file_i),"/results_singleStep/");

            % Set final time of simulation to one timestep
            spinConfigShort.maxT = spinConfigShort.IntTimeStep / (3600.0);

            % Fill the config file using the struct
            fileID = fopen('../config_preprocess.txt', 'w');
            fields = fieldnames(spinConfigShort);
            for k = 1:numel(fields)
                % Print field name and its corresponding value to the file
                fprintf(fileID, '%s = %s\n', fields{k}, num2str(spinConfigShort.(fields{k})));
            end
            fclose(fileID);
             Nbodies = length(mass);
            M_tot = sum(mass(:,2));
            rhocpp = 3000;
            volume =[mass(:,1), mass(:,2) ./ rhocpp];
            volume_all = sum(volume(:,2));
           

            NB.Rb = 10;
            NB.Ras_min = 10;
            NB.Ras_max = 100000000000000;
            Mdensity = rhocpp;
            NB.Nb =  size(positions,1);
            NB.Vmat = M_tot/Mdensity;
            NB.Vmat = volume_all;
            RW.rhoB_real = real_bulk_density;
            RW.T_real = real_period;
            RW.Vt_real = real_volume;
            SM.useVmin = 0;
            
            if(NB.prorositycom == 1)
                Pend = zeros(NB.Nb,3);
                for j = 1:NB.Nb
                    Pend(j,:) = positions(j,:);
                end
            end

            xL=abs(max(Pend(:,1)))+abs(min(Pend(:,1)));
            yL=abs(max(Pend(:,2)))+abs(min(Pend(:,2)));
            zL=abs(max(Pend(:,3)))+abs(min(Pend(:,3)));
            xL = xL + 2 * NB.Rb;yL = yL + 2*NB.Rb; zL = zL + 2*NB.Rb;

            %%          Check the volume of the shape with the alphavol function
            [Vol_min,~] = alphavol(Pend,NB.Ras_min);
            [Vol_max,~] = alphavol(Pend,NB.Ras_max);
            
            % Check the porosity of the rubble pile
            porosity_min=1-NB.Vmat/Vol_min;
            porosity_max=1-NB.Vmat/Vol_max;

        
        %% RP model with real world data
            % "_real"  = real world d             
            G = 6.67408e-11;
            % shape model in simulation units
            % Here there is a large difference between the sf_obj and the sf
            %Vt_sh = volume_all;
            if SM.useVmin
                Vt_sh=Vol_min;    % [L^3] total volume (min envelope)
            else
                Vt_sh=Vol_max;    % [L^3] total volume (max envelope)
            end

            Vm_sh = NB.Vmat;
            stff = RW.Vt_real/Vt_sh;
            sf = (stff)^(1/3);

            % simulation parameters
            rhoM_sim = RW.rhoB_real*RW.Vt_real/Vm_sh;
            AngVel_sim = 2*pi/(RW.T_real*3600);
            G_sim = G/stff;

            % Change the directory from which it loads the shape and the results one
            spinConfigFull = spinConfig;
            spinConfigFull.load_directory = append("../Shapes/Asteroids/", filename(file_i),"/");
            file_directory = append(filename(file_i),"/rho ", num2str(real_bulk_density ,formatSpec), "/Rotation Period ", num2str(angular_velocity * 3600,'%.3f'));
            spinConfigFull.save_directory = append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/spherical_large_no_core_fixed/results/",file_directory, "/");
            % Change the scaled factors
            spinConfigFull.G_univ = G_sim;
            spinConfigFull.initialSpin = AngVel_sim;
            spinConfigFull.density = rhoM_sim;

            % One timestep
            spinConfigFull.maxT = spinConfigFull.IntTimeStep / (3600.0);

            % Fill the config file
            mkdir(append('C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/spherical_large_no_core_fixed/config/config_shapes/', file_directory));
            conf_file = append('C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/spherical_large_no_core_fixed//config/config_shapes/', file_directory, '/config_shape_', filename(file_i), '.txt');
            fileID = fopen(conf_file, 'w');
            fields = fieldnames(spinConfigFull);
            for k = 1:numel(fields)
                % Print field name and its corresponding value to the file
                fprintf(fileID, '%s = %s\n', fields{k}, num2str(spinConfigFull.(fields{k})));
            end
            fclose(fileID);

            % Extra file that fills in the results of this script, check it for consistency
            mkdir(spinConfigFull.save_directory);
            fileID = fopen(append(spinConfigFull.save_directory, "simInputs.txt"),'w');
            fprintf(fileID, "Characteristic lengths (x,y,z) = (%.0f,%.0f,%.0f) \n",xL,yL,zL);
            fprintf(fileID, "Geometrical slenderness ratio = %.3f\n",min([xL,yL,zL])/max([xL,yL,zL]));
            fprintf(fileID, "Volume filled with material = %.4e\n",NB.Vmat);
            fprintf(fileID, "Min-max volume of the aggregate = %.3e-%.3e\n",volume_all);
            fprintf(fileID, "Min-max porosity of the aggregate = %.3f-%.3f\n\n",porosity_min,porosity_max);
            fprintf(fileID, "Total volume = %.4e \n", RW.Vt_real);
            fprintf(fileID, "Bulk density = %.2f \n", RW.rhoB_real);
            fprintf(fileID, "Spin period = %.2f \n\n", RW.T_real);
            fprintf(fileID, "Scaling factor^3 = %f \n", stff);
            fprintf(fileID, "Scaling factor = %f \n\n", sf);
            fprintf(fileID, "Material density = %.2f \n", rhoM_sim);
            fprintf(fileID, "Spin rate = %.4e \n", AngVel_sim);
            fprintf(fileID, "Universal gravity constant G = %.4e \n\n", G_sim);
            M_tot= rhoM_sim * NB.Vmat;    % [kg]      total mass of the aggregate
            fprintf(fileID, "Total mass = %.3e \n",M_tot);
            fprintf(fileID, "Material density = %.3f \n",rhoM_sim/stff);
            fprintf(fileID, "Characteristic lengths (x,y,z) = (%.0f,%.0f,%.0f) \n",xL*sf,yL*sf,zL*sf);
            fprintf(fileID, "Geometrical slenderness ratio = %.3f\n",min([xL,yL,zL])/max([xL,yL,zL]));
            fprintf(fileID, "Total volume = %.3e \n",Vt_sh*stff);
            porosity = 1-NB.Vmat/Vt_sh;        % [-] porosity
            fprintf(fileID, "Porosity of the aggregate = %.3f\n",porosity);
            Bdensity = M_tot/Vt_sh/stff;    % [kg/m^3]  bulk density
            fprintf(fileID, "Bulk density of the aggregate = %.3f\n",Bdensity);
            AngVel_n = AngVel_sim/sqrt(4*pi*G*Bdensity/3); % [-] normalized spin rate
            fprintf(fileID, "Normalized angular velocity = %.3f\n\n",AngVel_n);
            fclose(fileID);


            % Check the maximum time step possible and take the
            % minimum value obtained; The solution was obtained by
            % taking as reference velocities the maximum velocity
            % encoutered in similar simulations
            max_time_Step_1(rho_i, w_i, file_i) = 1/(2 * sqrt(G_sim * rhoM_sim/stff));
            delta = 1e-3;
            R = 2;

            filename1 = ["Bennu_v20_200k.obj","Kleopatra.obj", "Geographos.obj"];  % Replace with your OBJ file path

            dens = ["rho 1200/", "rho 1600/", "rho 2000/", "rho 2400/", "rho 3000/"];
            pers = ["Rotation Period 0.288/", "Rotation Period 1.260/", "Rotation Period 1.980/", "Rotation Period 2.520/", "Rotation Period 3.060/", "Rotation Period 3.240/"];
            x_dot = readmatrix(append("C:/Users/mihne/Documents/GitHub/Chrono_Projects/files/spherical_large_no_core_fixed/Chrono_Reintroduction/CutOut_simulation_build/Release/ResultsSave/", filename1(file_i), "\", pers(w_i), dens(rho_i),"Velocities.txt"));
            max_x_dot =  max(vecnorm(x_dot(:,2:4),2,2));

            max_time_step_2(rho_i, w_i, file_i) = 2 * sqrt(4 * R * delta - delta^2 ) / max_x_dot;
        end
    end
end


