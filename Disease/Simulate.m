% Modeling and Simulating Social Systems with MATLAB
% http://www.soms.ethz.ch/matlab
% Author: Stefano Balietti and Karsten Donnay, 2012
clear; clc
% Simulate disease spreading on a 2D grid

% Set parameter values
N=200;              % Grid size (NxN)
beta=0.2;           % Infection rate
gamma=0.1;          % Immunity rate
delta= 0.01;        % Re-infection rate
epsilon = 0.01;      %  death rate
Rinit = 5;

% define grid
x = zeros(N, N);    % The grid x, is coded as:  0=Susceptible, 
                    % 1=Infected, 2=Removed, 3=Dead

% Set the initial grid, x with a circle of infected individuals in the
% center of the grid, and with a radius of Rinit cells.
for i=1:N
    for j=1:N
        dx = i-N/2;
        dy = j-N/2;
        R = sqrt(dx*dx+dy*dy);
        if R <= Rinit
            x(i,j) = 1;
        end
    end
end

figure(1)
% Animate
clf                                 % Clear figure
imagesc(x, [0 2])                   % Display grid
colormap([0 0 1; 1 0 0; 0 1 0]);    % Define colors: Red, Green, Blue
pause(0.01)                         % Pause for 0.01 s



% Define the Moore neighborhood, i.e. the 8 nearest neighbors
neigh = [-1 -1; 0 -1; 1 -1; 1 0; 1 1; 0 1; -1 1; -1 0];

% Go back to the docked figure
figure(1)

% main loop, iterating the time variable, t
for t=1:100000
 
    % iterate over all cells in grid x, for index i=1..N and j=1..N
    for i=1:N
        for j=1:N
            
            if x(i,j) == 3
                % do nothing
            elseif x(i,j) == 0 % check if the cell is Susceptible, if not move along
                % Iterate over the neighbors and spread the disease
                for k=1:8
                    i2 = i+neigh(k, 1);
                    j2 = j+neigh(k, 2);
                    % Check that the cell is within the grid boundaries
                    if ( i2>=1 && j2>=1 && i2<=N && j2<=N )
                        % if cell is in state Susceptible and neighboring cell
                        % Infected => Spread infection with probability beta
                        if x(i2,j2) == 1 && rand<beta
                            x(i,j) = 1;
                        end
                        
                    end
                end
                
            
            % If infected => Recover from disease with probability gamma
            elseif ( x(i,j)==1 )
                if rand<gamma
                    x(i,j) = 2;
                elseif rand<epsilon
                    x(i,j) = 3;
                end
            elseif ( x(i,j)==2 && rand<delta )
                x(i,j) = 0;
            end            
            
        end
    end
    
    % If no more infected => Stop the simulation
    if ( sum(x==1)==0 )
        x(x==2) = 0;
        flag=1;
    end

    % Animate
    clf                                 % Clear figure
    imagesc(x, [0 3])                   % Display grid
    colormap([1 1 1; 1 0 0; 0 1 0; 0 0 0]);    % Define colors: Red, Green, Blue
    pause(0.00001)                         % Pause for 0.01 s
    
    if flag==1; break; end
end
