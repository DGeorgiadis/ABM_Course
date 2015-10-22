% Code for Axelrod's (1997) benchmark model of cultural dissemination.

%% Clean.
clc;close all; clear 

%% Step 1: Define matrix of cultural traits.

% 1.1 Define variables of benchmark model.

C1 = 2;
n = 10; % Nr. regions per row.
tn = n^2; % Total nr. regions.
f = 5; % Nr. features for each region = f.
t = 10; % Nr. different traits for each region = t.
a = 0; % Minimum value of trait.
b = 9; % Maximum value of trait.

% Save a video of the simulation.
makeVideo = 1;

% Every x ticks the display is updated.
DISPLAY_UPDATE_INTERVAL = 100;

% 1.2 Generate matrix of traits.

% Randomly generate a matrix of (n*f)*(n*f) = (10*5)*(10*5) of uniformly
% distributed numbers with values in the range [a,b] = [0,9]
% where x = a + (a+b)*rand(n*f,n*f).
traits = a + (a+b)*rand(n,n*f);
% Round up to get rid of decimals.
traits = floor(traits);

% Original matrix.
original_matrix = traits;


% 1.3 Special regions (either corner or edge).
vector_positions = 1 : tn;
matrix_positions = reshape(vector_positions,n,n);
matrix_positions = transpose(matrix_positions);
all_special_positions = [matrix_positions(1,:) ...
                         matrix_positions(n,:) ...
                         matrix_positions(2:n-1,1)' ...
                         matrix_positions(2:n-1,n)'];
all_special_positions = sort(all_special_positions);

corner_positions = [1, n, (n-1)*n+1, n*n];
edge_positions = setdiff(all_special_positions,corner_positions);

% Corner positions
corner_NW = 1;
corner_NE = n;
corner_SW = (n-1)*n+1;
corner_SE = n*n;
% Edge positions
edge_West = matrix_positions(2:n-1,1)';
edge_North = matrix_positions(1,2:n-1);
edge_East = matrix_positions(2:n-1,n)';
edge_South = matrix_positions(n,2:n-1);


%% Step 2: Model the interaction between regions.

% This consists in: 
% 2.1 Select randomly one active region (from 1 to 100)
% 2.2 Select randomly a neighbor of the active region
%     (from 1 to 4, where 1=West, 2=North, 3=East, 4=South)
% 2.3 Compare the active region with its neighbor, compute similarity
%     percentage (sim = nr active trains in common/nr traits)
% 2.4 If sim>0, select randomly one different trait and replace its value
%     with the value of the neighbor's trait

% Probability of interaction with neighbor 
% belongs to [0, 0.2, 0.4, 0.6, 0.8, 1].
matrix_probability_interaction = zeros(f+1,f);
matrix_probability_interaction(2,1) = 1;
matrix_probability_interaction(3,1:2) = 1;
matrix_probability_interaction(4,1:3) = 1;
matrix_probability_interaction(5,1:4) = 1;
matrix_probability_interaction(6,1:5) = 1;

% Vector position traits.
vector_position_traits = 1 : f;

% Repeat NSim times.
NSim = 150000;
neighbor_vec = zeros(NSim,1);
vector_start_sum = 1 : 5 : tn;
vector_stop_sum = 5 : 5 : tn;
tsum = zeros(n,n);
color_limits = [0 40];
%M = zeros(NSim,1);

if (makeVideo)
    % Get the handle of the figure
    h = figure();
    % Prepare the new file.
    vidObj = VideoWriter('video.avi');
    open(vidObj);
end

for sim = 1 : NSim
    
    % Print sim step.
    sim
    
    % Select radomly a region (from 1 to 100).
    region_active = floor(1 + tn*rand(1));
    % Where is this region placed in the traits? 
    rest_position_divided_by_0 = rem(region_active,n);
    if rest_position_divided_by_0 == 0
        % Column.
        column_active_region = n;
        % Row.
        row_active_region = floor(region_active/n);
    else
        % Column.
        column_active_region = rest_position_divided_by_0;
        % Row.
        row_active_region = floor(region_active/n) + 1;
    end
    

    % What are the traits of the active region? (A vector of length 5).
    traits_active = traits(row_active_region,(column_active_region-1)*f+1:column_active_region*f);
    
    %% Select randomly a neighbor of the active region
    % Is the active region in a special position?
    indie_special = isempty(intersect(all_special_positions,region_active));
    if indie_special ==  1
        neighbor_region = floor(1 + 4*rand(1));
        % Is the neighbor in the West? (Same row, column to the left).
        if neighbor_region == 1 
            traits_neighbor = traits(row_active_region,(column_active_region-2)*f+1:(column_active_region-1)*f);
        % Is the neighbor in the North? (Row above, same column).
        elseif neighbor_region == 2 
            traits_neighbor = traits(row_active_region-1,(column_active_region-1)*f+1:column_active_region*f);
        % Is the neighbor in the East? (Same row, column to the right).
        elseif neighbor_region == 3
            traits_neighbor = traits(row_active_region,column_active_region*f+1:(column_active_region+1)*f);
        % Is the neighbor in the South? (Row below, same column).
        else
            traits_neighbor = traits(row_active_region+1,(column_active_region-1)*f+1:column_active_region*f);
        end
        
        
        
        % If the active region is in a special position, account
        % differently for the possible neighbors.
        
    %%%%%%%% Is the active region on the CORNER?
    elseif isempty(intersect(corner_positions,region_active)) == 0
        
        neighbor_region = floor(1 + 2*rand(1));
        
        if region_active == corner_NW
            % Is the neighbor in the South? (Row below, same column)
            if neighbor_region == 1 
                traits_neighbor = traits(row_active_region+1,(column_active_region-1)*f+1:column_active_region*f);
            % Is the neighbor in the East? (Same row, next column)
            else
                traits_neighbor = traits(row_active_region,column_active_region*f+1:(column_active_region+1)*f);
            end
        elseif region_active == corner_NE
            % Is the neighbor in the South? (Row below, same column)
            if neighbor_region == 1 
                traits_neighbor = traits(row_active_region+1,(column_active_region-1)*f+1:column_active_region*f);
            % Is the neighbor in the West? (Same row, previous column)
            else
                traits_neighbor = traits(row_active_region,(column_active_region-2)*f+1:(column_active_region-1)*f);
            end
        elseif region_active == corner_SE
            % Is the neighbor in the North? (Row above, same column)          
            if neighbor_region == 1
                traits_neighbor = traits(row_active_region-1,(column_active_region-1)*f+1:column_active_region*f);
            % Is the neighbor to the West? (Same row, previous column)
            else 
                traits_neighbor = traits(row_active_region,(column_active_region-2)*f+1:(column_active_region-1)*f);
            end
        elseif region_active == corner_SW
            % Is the neighbor in the North? (Row above, same column)          
            if neighbor_region == 1
                traits_neighbor = traits(row_active_region-1,(column_active_region-1)*f+1:column_active_region*f); 
            % Is the neighbor in the East? (Same row, next column)
            else
                traits_neighbor = traits(row_active_region,column_active_region*f+1:(column_active_region+1)*f);
            end 
        end
            

    %%%%%%%%%%% Otherwise the active region is on the EDGE     
    else
        indie_edge_West = isempty(intersect(edge_West,region_active));
        indie_edge_North = isempty(intersect(edge_North,region_active));
        indie_edge_East = isempty(intersect(edge_East,region_active));
        indie_edge_South = isempty(intersect(edge_South,region_active));
        
        
        neighbor_region = floor(1 + 3*rand(1));
        % Is the active region on the West edge?
        if indie_edge_West == 0
            % Is the neighbor in the North? (Row above, same column)
            if neighbor_region == 1
                traits_neighbor = traits(row_active_region-1,(column_active_region-1)*f+1:column_active_region*f);
            % Is the neighbor in the East? (Same row, column to the right)
            elseif neighbor_region == 2
                traits_neighbor = traits(row_active_region,column_active_region*f+1:(column_active_region+1)*f);
            % Is the neighbor in the South? (Row below, same column)
            else
               traits_neighbor = traits(row_active_region+1,(column_active_region-1)*f+1:column_active_region*f);
            end
            
        elseif indie_edge_North == 0
            % Is the neighbor in the West? (Same row, column to the left)
            if neighbor_region == 1
                traits_neighbor = traits(row_active_region,(column_active_region-2)*f+1:(column_active_region-1)*f);
            % Is the neighbor in the East? (Same row, column to the right)
            elseif neighbor_region == 2
                traits_neighbor = traits(row_active_region,column_active_region*f+1:(column_active_region+1)*f);
            % Is the neighbor to the South? (Row below, same column)
            else
                traits_neighbor = traits(row_active_region+1,(column_active_region-1)*f+1:column_active_region*f);
            end
            
        elseif indie_edge_East == 0
            % Is the neighbor in the West? (Same row, column to the left)
            if neighbor_region == 1
                traits_neighbor = traits(row_active_region,(column_active_region-2)*f+1:(column_active_region-1)*f);
            % Is the neighbor in the North? (Row above, same column)
            elseif neighbor_region == 2
                traits_neighbor = traits(row_active_region-1,(column_active_region-1)*f+1:column_active_region*f);
            % Is the neighbor to the South? (Row below, same column)
            else
                traits_neighbor = traits(row_active_region+1,(column_active_region-1)*f+1:column_active_region*f);
            end
            
        else 
            % Is the neighbor in the West? (Same row, column to the left)
            if neighbor_region == 1
                traits_neighbor = traits(row_active_region,(column_active_region-2)*f+1:(column_active_region-1)*f);
            % Is the neighbor in the North? (Row above, same column)
            elseif neighbor_region == 2
                traits_neighbor = traits(row_active_region-1,(column_active_region-1)*f+1:column_active_region*f);   
            % Is the neighbor in the East? (Same row, column to the right)
            else 
                traits_neighbor = traits(row_active_region,column_active_region*f+1:(column_active_region+1)*f);    
            end
        end
    end
    
    
    
    %% What is the degree of similarity between active and neighbor?
    nr_common_traits = sum(traits_active == traits_neighbor);
    nr_different_traits = f - nr_common_traits;
    degree_sim = nr_common_traits/f;
        
    % Based on the degree of similarity, determine if they will
    % actually interact or not. For this, select randomly one element
    % from the vector probability interaction.
    vector_probability_interaction = matrix_probability_interaction(nr_common_traits+1,:);
    indicator_interaction = vector_probability_interaction(floor(1 + 4*rand(1)));
        
    %% Model the interaction and the change in active's traits
    if (indicator_interaction > 0) && (degree_sim < 1)
        % On which position (from 1 to f) are the different traits?
        indie_difference = traits_active ~= traits_neighbor;  % the mask of different traits
        different_positions = vector_position_traits(indie_difference); 
        random_position_to_change = floor(1 + nr_different_traits*rand(1));
        position_to_change = different_positions(random_position_to_change);
        traits_active(position_to_change) = traits_neighbor(position_to_change)+C1*rand(1);
        if length(different_positions)==0
            a_c = 1; b_c = 8;
            position_to_change = (b_c-a_c).*rand(1) + a_c;
            traits_active(position_to_change) = C1*rand(1);
        end
        traits_active(traits_active<1) = 1;
        traits_active(traits_active>8) = 8;
    end
          
    %% Save changes of active region into the traits matrix 
    traits(row_active_region,(column_active_region-1)*f+1:column_active_region*f) = traits_active;
            
    
    %% Visualize graphically every 20 steps.
    
    if mod(sim, DISPLAY_UPDATE_INTERVAL) == 0        
        % For each region, compute the sum of its total traits
        % (to obtain a n*n matrix).
        for ii = 1 : n
            indie_start_sum = vector_start_sum(ii);
            indie_stop_sum = vector_stop_sum(ii); 
            tsum(:,ii) = sum(traits(:,indie_start_sum:indie_stop_sum),2);       
        end
        imagesc(tsum, color_limits)
        % colorbar
        
        if (makeVideo)
            currFrame = getframe(h);
            writeVideo(vidObj,currFrame);
        end
        pause(0.01)
   end
   
end

if (makeVideo)
    % SAVE VIDEO FRAMES TO FILE
    close(vidObj);
end
