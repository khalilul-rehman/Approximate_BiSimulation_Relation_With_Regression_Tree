clc; clear; close all;

%% Parameters
n = 3;  % Grid size (rows)
m = 3;  % Grid size (columns)
dt = 0.1;  % Time step
T = 10;  % Simulation duration
A = [-1.2 0.1; 
    0.1 -1.2];  % System matrix
B = -A; % Input matrix (identity for direct control of velocity)
C = eye(2); % Output matrix
D = zeros(2); % No direct feedthrough

%% Convert Continuous System to Discrete using c2d
sys_c = ss(A, B, C, D);  % Create continuous state-space system
sys_d = c2d(sys_c, dt, 'zoh');  % Discretize using Zero-Order Hold (ZOH)
Ad = sys_d.A;
Bd = sys_d.B;

%% Define Grid Map  -1 for target 'A', -2 for forbidden 'B')
grid_map = [-2 2 4; 
            4 3  4; 
            2 2 -1];


numberOfSimulations = 1000;
trajectories = cell(numberOfSimulations, 0);

for counterForSimulation = 1 : numberOfSimulations
    
    %% Initial Conditions
    % Initial Position
    x = [0.5; 1.5];
    % Randon initial position
    [x(1), x(2)] = getRandomInitialPosition(grid_map);
    % disp(x);
    
    v = [0; 0];  % Initial velocity
    a = -1;
    b = 1.5;
    r = (b-a).*rand(2,1) + a;
    % random initial velocity
    % v = r;
    
    



    trajectory = [ x ; v];  % Store trajectory
    
    
    %% Simulation Loop
    for t = 0:dt:T
        % Determine the grid cell (round to nearest integer index)
        x_position = ceil(x(1));
        y_position = ceil(x(2));
        
        i = n + 1 - y_position;
        j = x_position;
        
        % if(i > n || j > m)
        %     break;
        % end
    
    if(i > n)
            i = min(i, n);
            % break;
        elseif (i < 1)
            i = 1;
        end
        if(j > m)
            j = min(j, m);
        elseif (j < 1)
            j = 1;
        end
    
        
        % disp([ int2str(x_position), ' , ', int2str(y_position) , '->', int2str(i), ' , ' , int2str(j), ' = ' , int2str(grid_map(i, j))]);
        
        % Desired velocity (lookup from grid_map, mapping 0-7 to angles)
        angles = (0:7) * (pi/4);
        if grid_map(i, j) >= 0
            v_desired = [sin(angles(grid_map(i, j)+1)); cos(angles(grid_map(i, j)+1))];
        else
            v_desired = [0; 0]; % Stop movement in special cells
            % break;
        end
        
        % Discrete-time velocity update
        v = Ad * v + Bd * v_desired;
        
        % Discrete-time position update
        x = x + dt * v;
        
        % Store trajectory
        currentStepTrajectory = [x ; v];
        trajectory = [trajectory, currentStepTrajectory];
        
        % Stop if reaching A or B
        if grid_map(i, j) == -1
            disp('Reached Target (A)');
            break;
        elseif grid_map(i, j) == -2
            disp('Entered Forbidden Zone (B)');
            break;
        end
    end

    trajectories = saveDataInTrejactories(counterForSimulation, trajectories, trajectory);

end
%% 

% --- Plotting the Results ---
figure;
hold on;
grid on;
title('Navigation System Trajectory');
xlabel('x1 (Horizontal Position)');
ylabel('x2 (Vertical Position)');
% set(gca,'FontSize',30);

% Plot grid boundaries
for i = 1:min(n,m)
    plot([0, n], [i-1, i-1], 'k--');
    plot([i-1, i-1], [0, m], 'k--');
end

for i = 1 : size(grid_map,1)
    for j = 1 : size(grid_map, 2)
        if grid_map(i,j) == -1
            text(j-0.5, (n+1)-i-0.5, 'A', 'HorizontalAlignment', 'center');
        elseif grid_map(i,j) == -2
            text(j-0.5, (n+1)-i-0.5, 'B', 'HorizontalAlignment', 'center');
        else
            text(j-0.5, (n+1)-i-0.5, int2str(grid_map(i,j)), 'HorizontalAlignment', 'center');
            quiver(j-0.4, (n+1)-i-0.25,sin( grid_map(i,j)* (pi/4))/4, cos( grid_map(i,j)* (pi/4))/4,'linewidth',0.75);
        end
    end
end
for i = 1 : size(trajectories, 1)
    plot(trajectories{i}(1, :), trajectories{i}(2, :), 'b-', 'LineWidth', 1);
    scatter(trajectories{i}(1, :), trajectories{i}(2, :));
end
hold off;   
%% 



function trajectories = saveDataInTrejactories(index, trajectories,currentStepTrajectory)
    for i = 1 : size(currentStepTrajectory, 1)
        for j = 1 : size(currentStepTrajectory, 2)
            trajectories{index, 1}(i,j) = currentStepTrajectory(i,j);
        end
    end
    
end

% Getting Random Position on Navigation Plane
function [x_position, y_position] = getRandomInitialPosition(grid_map)
    [row, column, dataValues] = find(grid_map > 0.1);
    navigationBlocks_to_strat = [column - 1, size(grid_map,1) - row];
    randomIndex = int32(( size(navigationBlocks_to_strat,1) - 1).*rand(1,1) + 1);
    
    safeSpan = 0.1;

    x_position = navigationBlocks_to_strat(randomIndex,1);
    x_position = ( x_position+1 - x_position).*rand(1,1) + x_position;

    
    y_position = navigationBlocks_to_strat(randomIndex,2);
    y_position = ( y_position+1 - y_position).*rand(1,1) + y_position;

     if (x_position > 2.9)
         x_position = 2.9;
     end
     if y_position < 0.1
         y_position = 0.1;
     end
    
end


%% Exporting the Trejectories
% Starting With one trejectory
DataSet = [];
for i = 1 : size(trajectories,1)
    one_Trejectory = trajectories{i}(1:2,:)';
    size(one_Trejectory);
    
    d_lag = -1;
    lags = [0 d_lag];
    XLag_Trejectory = lagmatrix(one_Trejectory,lags);
    XLag_Trejectory = XLag_Trejectory(1 : end + d_lag, :);
    size(XLag_Trejectory);
     
    % Y_Trejectory = circshift(one_Trejectory,-1);
    % Y_Trejectory = Y_Trejectory(1 : end + d_lag, :);
    % size(Y_Trejectory)
    
    
    DataSet = [DataSet ; XLag_Trejectory]; %[DataSet ; XLag_Trejectory Y_Trejectory]
end

%% 

  csvwrite('./Data_Files/navigation_trejectory_dataset.csv',DataSet)
