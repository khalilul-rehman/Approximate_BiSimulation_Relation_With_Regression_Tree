clc; clear; close all;

%% Parameters
dt = 0.1;  % Time step
T = 12;    % Total simulation time
time_steps = T / dt;

% Continuous System Matrices (from Equation 8 in the paper)
A = [-0.9  0.5  0;
      0.5 -1.3  0.5;
      0    0.5 -0.9];

B = [0.4; 0.3; 0.4]; % External temperature influence
H = diag([6, 7, 8]); % Heater influence matrix

C = eye(3); % Output matrix (identity for temperature tracking)
D = zeros(3,4); % Must match number of inputs (1 for B + 3 for H)

% Convert Continuous System to Discrete using c2d
sys_c = ss(A, [B H], C, D);  % Combine B and H as inputs
sys_d = c2d(sys_c, dt, 'zoh');  % Discretize using Zero-Order Hold (ZOH)
Ad = sys_d.A;
Bd = sys_d.B(:, 1);  % First column corresponds to B (external temp)
Hd = sys_d.B(:, 2:4); % Remaining columns correspond to discretized H

% Outside temperature
u = 4;

numberOfSimulation = 10;
trejectories = cell(numberOfSimulation, 0);
list_of_initial_temprature = getRandomInitialTemprature(numberOfSimulation, 12, 18) %-10, 20

for simulation_counter = 1 : numberOfSimulation
    % Initial conditions
    x = list_of_initial_temprature(simulation_counter, :)'; %[-10; 20; 20];  % Initial temperature of rooms
    h = [1; 1; 0];      % Heater status (1 = on, 0 = off)
    heater_w_r_t_room = [1; 1; 0];

    
    
    % Storage for results
    trajectory = zeros(3, time_steps);
    trajectory(:,1) = x;
    
    % Heater thresholds
    on_thresh = [20; 20; 20];  % Turn heater ON
    off_thresh = [21; 21; 21]; % Turn heater OFF
    
    % Heater movement rules
    get_thresh = [18; 18; 18];  % If room temperature drops below, it may get a heater
    diff_thresh = [1; 1; 1];     % Heater moves if adjacent room is warmer by this amount
    
     
    %% Discrete-Time Simulation Loop
    for k = 2:time_steps
        % Compute the change in temperature using discrete-time model
        x = Ad * x + Bd * u + Hd * h;  % Use Hd instead of H
        
        % Apply thermostat logic
        for i = 1:3
            if x(i) <= on_thresh(i)
                % Only turn on the heater if this room has a heater
                if(heater_w_r_t_room(i))
                    h(i) = 1; % Turn heater ON
                end
            elseif x(i) >= off_thresh(i)
                h(i) = 0; % Turn heater OFF
            end
        end
    
        % Heater movement logic if temperature is low and there is another
        % heater available to move
        
        for i = 1:3
            for j = 1:3
                if i ~= j && heater_w_r_t_room(i) == 0 && heater_w_r_t_room(j) == 1 && x(i) <= get_thresh(i) && (x(j) - x(i)) >= diff_thresh(i)
                    heater_w_r_t_room(i) = 1;
                    h(i) = 1;
                    heater_w_r_t_room(j) = 0;
                    h(j) = 0;
                    
                    % disp([int2str(k), "==  Movement is performed from -> ", int2str(j), " , " , int2str(i) ]);
                    
                end
            end
        end
    
        % Store data
        trajectory(:, k) = x;
        
    end

    trejectories{simulation_counter,1} = trajectory;
end

%% Visualization
figure;
for i = 1 : numberOfSimulation
    plot(0:dt:T-dt, trejectories{i}(1, :), 'r', 'LineWidth', 1.5);
    hold on;
    plot(0:dt:T-dt, trejectories{i}(2, :), 'g', 'LineWidth', 1.5);
    plot(0:dt:T-dt, trejectories{i}(3, :), 'b', 'LineWidth', 1.5);
end
legend('Room 1', 'Room 2', 'Room 3');
xlabel('Time (minutes)');
ylabel('Temperature (°C)');
title('Discrete-Time Room Heating Benchmark Simulation (Corrected)');
grid on;


%% function to get randon initial temprature
function [ list_of_random_tempratures ] = getRandomInitialTemprature(numOfTempratures, minTemprature, maxTemprature)
    list_of_random_tempratures = zeros(numOfTempratures,3);
    counter_for_loop = 1;
    while counter_for_loop <= numOfTempratures
        tempRow = [(maxTemprature-minTemprature).*rand(1,1) + minTemprature ...
            (maxTemprature-minTemprature).*rand(1,1) + minTemprature ...
            (maxTemprature-minTemprature).*rand(1,1) + minTemprature ]
        if ~nnz(ismember(list_of_random_tempratures , tempRow ,'rows'))
            list_of_random_tempratures(counter_for_loop, :) = tempRow;
            counter_for_loop = counter_for_loop + 1;
        end

    end
end

%% Exporting the data

DataSet = [];
time_DataSet = [];
for i = 1 : size(trejectories,1)
    one_Trejectory = trejectories{i}(:,:)';
    size(one_Trejectory);

    
    
    d_lag = -1;
    lags = [0 d_lag];
    XLag_Trejectory = lagmatrix(one_Trejectory,lags);
    XLag_Trejectory = XLag_Trejectory(1 : end + d_lag, :);
    size(XLag_Trejectory);


    % just for storing the time record
    time_steps_record = [dt:dt:T]';
    time_steps_record = time_steps_record(1:end +d_lag, :);
    time_DataSet = [time_DataSet;time_steps_record];
     
    % Y_Trejectory = circshift(one_Trejectory,-1);
    % Y_Trejectory = Y_Trejectory(1 : end + d_lag, :);
    % size(Y_Trejectory)
    
    
    DataSet = [DataSet ; XLag_Trejectory]; %[DataSet ; XLag_Trejectory Y_Trejectory]
end
%% Saving data in file

 % csvwrite('RoomHeating_trejectory_dataset_test_12_18.csv',DataSet)
