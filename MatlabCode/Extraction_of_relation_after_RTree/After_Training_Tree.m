addpath('SupportingFunctions/');
addpath('SupportingFunctions/Polytopes/Plot');
addpath('SupportingFunctions/Polytopes/CustomPolyTope/');
addpath('SupportingFunctions/Optimization/');
addpath("SupportingFunctions/ScenarioBasedApproach/");
%% 
%parpool;

%% 

    parentDirectory = fileparts(cd);

     %dataFiles_folder_path = strcat(parentDirectory, '/CaseStudy_Simulation/NavigationSystem/Data_Files/4D_Data/');% '/Users/khalilulrehman/Academic/Phd Italy 2023_26/University of LAquila/Research Papers tasks/MatlabCodes/RegressionTree/Extraction_of_relation_after_RTree/Data_Files/';
     %dataFiles_folder_path = strcat(parentDirectory, '/CaseStudy_Simulation/NavigationSystem/Data_Files/4D_Data/OneStepSimulation/10Th_MSLeaf250/');
     dataFiles_folder_path = strcat(parentDirectory, '/CaseStudy_Simulation/NavigationSystem/Data_Files/4D_Data/OneStepSimulation/30k_R0_3/');
     
     traning_data_trajectories_file = 'leaf_classified_trejectory_dataset.csv';
     test_data_trajectories_file = 'leaf_classified_test_trejectory_dataset.csv';
     constraint_file_name = '/constrants_array_on_leaves.csv';
     Dimension_of_Attributes = 4;
     Dimension_of_ClassVariables = 4;

     %{
     
     dataFiles_folder_path = strcat(parentDirectory, '/CaseStudy_Simulation/RoomHeatingSystem/DataFiles/'); % '/Users/khalilulrehman/Academic/Phd Italy 2023_26/University of LAquila/Research Papers tasks/MatlabCodes/RegressionTree/CaseStudies/RoomHeatingBenchmark/';
     traning_data_trajectories_file = 'leaf_classified_trejectory_dataset.csv';
     test_data_trajectories_file = 'leaf_classified_test_trejectory_dataset.csv';
     constraint_file_name = 'constrants_array_on_leaves.csv';
     Dimension_of_Attributes = 3;
     Dimension_of_ClassVariables = 3;
     %}

%% 

% This file contains data set 1st column is leaf number, 2nd, 3rd 4th, 
% columns are attributes and 5th, 6th and 7th column is 3D Class
leaf_trajectories_regression = csvread(strcat(dataFiles_folder_path, traning_data_trajectories_file));

%% Split the matrix into cells\
% Find unique values
leaves = unique(leaf_trajectories_regression(:,1));
numLeaves = numel(leaves);

% Separate the matrix
trajectories_on_leaves = cell(size(leaves));
for kk = 1:numel(leaves)
  idx = leaf_trajectories_regression(:,1) == leaves(kk);
  % while saving the trejectories into cells I have deleted the leaf number
  % from data because each cell is a leaf
  trajectories_on_leaves{kk} = leaf_trajectories_regression(idx,2:end);

end

% trajectories_on_leaves

%% Load constraintes on leaves - Each row is for one leaf and in columns 
% we have "X1_min","X1_max", "Y1_min", "Y1_max", "X2_min", "X2_max", "Y2_min", "Y2_max"

constraints_on_leaves = csvread( strcat(dataFiles_folder_path, constraint_file_name));
% constraints_on_leaves

%% 
customPolytope = CustomPolyTope();
vertices_of_polytopes = customPolytope.getVerticesFromBounds(constraints_on_leaves);
%% Create an instance of the class
customPlot = CustomPlotClass();
%% 


if size(constraints_on_leaves,2) == 4
    % Call the method on the instance
    % customPlot.drawConstraintsIn2D(constraints_on_leaves, 'title', 'Constraints as 2D-Rectangles');
    % customPlot.drawConstraintsIn2D(constraints_on_leaves, 'title', 'Constraints as 2D-Rectangles', 'dataPoints2D', trajectories_on_leaves);
    customPlot.draw2DPolytopesWithVertices(vertices_of_polytopes, 'title', 'Constraints as 2D-Rectangles', 'dataPoints2D', trajectories_on_leaves);
    customPlot.draw2DPolytopesWithVertices(vertices_of_polytopes,'color', '#8f9119')
elseif size(constraints_on_leaves,2) == 6
    % customPlot.draw3DPolytopesWithVertices(vertices_of_polytopes, 'color', '#8f9119');
    % customPlot.draw3DPolytopesWithVertices(vertices_of_polytopes, 'title', 'Constraints as 3D-Cubes', 'color', '#8f9119');
    customPlot.draw3DPolytopesWithVertices(vertices_of_polytopes);
else
    customPlot.drawConstraintsIn2D(constraints_on_leaves, 'title', 'More then 3D Constraints in 2D');
    customPlot.draw3DPolytopesWithVertices(vertices_of_polytopes);
end

%% 
minOfAllDimensions = min(constraints_on_leaves);
maxOfAllDimensions = max(constraints_on_leaves);



%% LAMBDAs
LAMBDA = cell(numLeaves,1);
lambda = cell(numLeaves,1);

% We need to know the number of attributes and class variables
for i = 1 : numLeaves
    LAMBDA{i,1} = trajectories_on_leaves{i}(:, 1 : Dimension_of_Attributes);
    lambda{i,1} = trajectories_on_leaves{i}(:, Dimension_of_Attributes + 1 : Dimension_of_Attributes + Dimension_of_ClassVariables);
end

%% Optimization For each leaf
customOptimization = CustomOptimization();
%% 

% h* m* mo*
solutionOptimal_star = cell(numLeaves,1);
% [ ones(numLeaves, 1) zeros(numLeaves, size(trejectories_on_leaves{1}, 2)) ones(numLeaves, 1)];
% Open a parallel pool if it's not already running
if isempty(gcp('nocreate'))
    parpool; % Starts a parallel pool with available workers
end
parfor i = 1 : numLeaves
     disp(['Optimizing h* for leaf = ' , int2str(i)]);
     [h, m , mo] = customOptimization.QuadraticConstraintOptimizer(LAMBDA{i}, lambda{i});
     % [h, m , mo] = QuadraticConstraintOptimizer_MultivariateY( LAMBDA{i}, lambda{i});
      solutionOptimal_star{i,1} = [ones(size(m,1),1)*h m mo];
 end


 %% Temp_ Remove it if do not worked well
 %{
 addpath('SupportingFunctions/temp_/')
 solutionOptimal_star = cell(numLeaves,1);
% [ ones(numLeaves, 1) zeros(numLeaves, size(trejectories_on_leaves{1}, 2)) ones(numLeaves, 1)];
% Open a parallel pool if it's not already running

if isempty(gcp('nocreate'))
    %parpool; % Starts a parallel pool with available workers
    parpool('local', 6); % Start a parallel pool with 4 workers
end


parfor i = 1 : numLeaves
%i = 25;
     disp(['Optimizing h* for leaf = ' , int2str(i)]);
     [h, m , mo] = QuadraticConstraintOptimizer(LAMBDA{i}, lambda{i})
     % [h, m , mo] = QuadraticConstraintOptimizer_MultivariateY( LAMBDA{i}, lambda{i});
      solutionOptimal_star{i,1} = [ones(size(m,1),1)*h m mo];
 end
 %}


 %% Testing h m mo on traning data

% y_hat_with_optimization = (m * trejectories_on_leaves{i}(:,1:4)')' + mo'
NRMSE_Training = customOptimization.NRMSE_Calculation(LAMBDA,lambda,solutionOptimal_star);
disp(["NRMSE on Training data = ", mean(NRMSE_Training)*100]);


%% Loading Testing Data
leaf_test_trajectories_regression = csvread(strcat(dataFiles_folder_path,test_data_trajectories_file));

%% Split the matrix into cells\

% Separate the matrix
test_trajectories_on_leaves = cell(size(leaves));
for kk = 1:numel(leaves)
  idx = leaf_test_trajectories_regression(:,1) == leaves(kk);
  % while saving the trejectories into cells I have deleted the leaf number
  % from data because each cell is a leaf
  test_trajectories_on_leaves{kk} = leaf_test_trajectories_regression(idx,2:end);

end

% test_trajectories_on_leaves

%% LAMBDAs for test
LAMBDA_test = cell(numLeaves,1);
lambda_test = cell(numLeaves,1);

% We need to know the number of attributes and class variables
for i = 1 : numLeaves
    LAMBDA_test{i,1} = test_trajectories_on_leaves{i}(:, 1 : Dimension_of_Attributes);
    lambda_test{i,1} = test_trajectories_on_leaves{i}(:, Dimension_of_Attributes + 1 : Dimension_of_Attributes + Dimension_of_ClassVariables);
end
%% Checking NRMSE on Test data

NRMSE_Testing = customOptimization.NRMSE_Calculation(LAMBDA_test,lambda_test,solutionOptimal_star);
disp(["NRMSE on Testing data = ", mean(NRMSE_Testing)*100]);



%% Boundry Point of each rectangle

% On each cell i will have 2D four boundry points of each rectangle = x1,y1
% x2,y2 x3,y3 ,x4,y4
% vertices_of_rectangles = cell(numLeaves,0);
% 
% for i = 1 : numLeaves
%     vertices_of_rectangles{i,1} = get_vertices_of_hyperrectangle(constraints_on_leaves(i,:));
% end

%% 

% plot_multiple_hyperrectangle(vertices_of_rectangles, 'color' , '#639119');
%% Testing the function
% figure; hold on;
% plot_hyperrectangle(vertices_of_rectangles{i,1});
% plot_hyperrectangle(add_span_to_hyperrectangle(vertices_of_rectangles{i,1}, 10));
% hold off;

%% Transection from main rectangle to elivated rectangles


q_shape_of_x = Dimension_of_Attributes;
p_shape_of_y = Dimension_of_ClassVariables;


% vertices_of_elevated_polytopes = cell(numLeaves,0);
vertices_of_elevated_polytopes_before_span = cell(numLeaves,0);
parfor i = 1 : numLeaves
    h = solutionOptimal_star{i}(1,1);
    m =  solutionOptimal_star{i}( 1 : p_shape_of_y , 2 : 2 - 1 + q_shape_of_x);
    mo = solutionOptimal_star{i}( 1 : p_shape_of_y , end);

    current_elevated_vertices = zeros(size(vertices_of_polytopes{i,1}));
    for j = 1 : size(vertices_of_polytopes{i,1},1)
        
        current_elevated_vertices(j,:) = (m * vertices_of_polytopes{i,1}(j,:)')' + mo';
    end
    
    % current_elevated_vertices
    vertices_of_elevated_polytopes_before_span{i, 1} = current_elevated_vertices;
   % Adding the h span in elevated hyper rectangles
    % vertices_of_elevated_polytopes{i,1} = add_span_to_hyperrectangle(current_elevated_vertices, h);
end





%% show elevated rectangles

if size(vertices_of_elevated_polytopes_before_span{1,1},1) == 4
    % Call the method on the instance
    customPlot.drawDual2DPolytopesWithVertices(vertices_of_polytopes, vertices_of_elevated_polytopes_before_span, 'color1', '#8f9119', 'color2', '#197b91', 'title', 'Evevated 2D-Polytopes');
    
elseif size(vertices_of_elevated_polytopes_before_span{1,1},1) == 8
    
    customPlot.drawDual3DPolytopesWithVertices(vertices_of_polytopes, vertices_of_elevated_polytopes_before_span, 'color1', '#8f9119', 'color2', '#197b91', 'title', 'Evevated 3D-Polytopes');
else
    customPlot.drawDual3DPolytopesWithVertices(vertices_of_polytopes, vertices_of_elevated_polytopes_before_span, 'color1', '#8f9119', 'color2', '#197b91', 'title', 'Evevated 3D-Polytopes');
    customPlot.drawDual2DPolytopesWithVertices(vertices_of_polytopes, vertices_of_elevated_polytopes_before_span, 'title', 'More then 3D Polytopes in 2D');

    
end


%% Checking the Distance between polytopes to check for transition

%{ 
state_graph_transition_matrix = zeros(numLeaves, numLeaves);
for i = 1 : numLeaves
    h = solutionOptimal_star{i}(1,1);
    for j = 1 : numLeaves
        % if  h < 2

            vertices_P1 = vertices_of_elevated_polytopes_before_span{i,1};
            vertices_P2 = vertices_of_polytopes{j,1};
            [P1_opt, P2_opt, min_distance] = customPolytope.minimize_polytope_distance_dual(vertices_P1, vertices_P2);
    

            if min_distance < h   
                %disp([min_distance, "" , h]);
                state_graph_transition_matrix(i,j) = 1;
            end
        % end
    end
end
%}
% Parallel Processing version is below
% Preallocate the state transition matrix
state_graph_transition_matrix = zeros(numLeaves + 1, numLeaves + 1); % this extra 1 is for to save the state goes out of the state space or not

% Open a parallel pool if it's not already running
if isempty(gcp('nocreate'))
    parpool; % Starts a parallel pool with available workers
end

% Parallelize the outer loop
parfor i = 1:numLeaves
    h = solutionOptimal_star{i}(1,1);

    % Initialize local variables (for parfor compatibility)
    local_row = zeros(1, numLeaves); % Store results for row i

    % Iterate over columns of the matrix
    for j = 1:numLeaves
        vertices_P1 = vertices_of_elevated_polytopes_before_span{i,1};
        vertices_P2 = vertices_of_polytopes{j,1};
        [~, ~, min_distance] = customPolytope.minimize_polytope_distance_dual(vertices_P1, vertices_P2);

        if min_distance < h
            local_row(j) = 1; % Assign transition state
        end
    end

    trajectory_go_out = 0;
    x_out_space = ~isempty(find(all(trajectories_on_leaves{i,1}(:,5) < 0 | trajectories_on_leaves{i,1}(:,5) > 3, 2)));
    y_out_space = ~isempty(find(all(trajectories_on_leaves{i,1}(:,6) < 0 | trajectories_on_leaves{i,1}(:,6) > 3, 2)));

    if x_out_space || y_out_space
        trajectory_go_out = 1;
    end

    % Store the result in the global matrix
    state_graph_transition_matrix(i, :) = [local_row trajectory_go_out];
end

% Display or process the resulting state graph
disp('State transition matrix computation completed.');

%% For one time
%{
for i = 1 : numLeaves
    x_out_space = ~isempty(find(all(trajectories_on_leaves{i,1}(:,5) < 0 | trajectories_on_leaves{i,1}(:,5) > 3, 2)));
    y_out_space = ~isempty(find(all(trajectories_on_leaves{i,1}(:,6) < 0 | trajectories_on_leaves{i,1}(:,6) > 3, 2)));
    if x_out_space || y_out_space
        state_graph_transition_matrix(i, numLeaves+1) = 1;
    end
end
%}

%% 
orange_states = find(all(state_graph_transition_matrix(:,numLeaves+1) > 0, 2));
red_states = [numLeaves+1];
showGraphWithHighLightedStates(state_graph_transition_matrix,orange_states, red_states);
 %% Getting stats on # of states found and noramally calculated
 [smallestPolytope, smallestEdgeLength, edgeVertices] = customPolytope.findSmallestPolytopeAndEdge(vertices_of_polytopes);

 
 % Display the results
disp('Smallest Polytope:');
disp(smallestPolytope);

disp('Smallest Edge Length:');
disp(smallestEdgeLength);

disp('Vertices of the Smallest Edge:');
disp(edgeVertices);
 


smallestPolytope = ones(1,size(minOfAllDimensions,2)/2) * smallestEdgeLength;
 
numberOfSmallestpolytopes = calculatePolytopesToFillSpace(minOfAllDimensions(1:2:end), maxOfAllDimensions(2:2:end), smallestPolytope );
disp(['# of Polytopes to fill the state space = ', int2str(numberOfSmallestpolytopes)]);
disp(['# of Polytopes predicted with our model = ', int2str(numLeaves)]);


%% 
%addpath('SupportingFunctions/ScenarioBasedApproach/calculateEpsilon/');
numberOfDataInEachPolytope = zeros(1,numLeaves);
for i = 1 :numLeaves
    numberOfDataInEachPolytope(1,i) = size(LAMBDA{i,1},1);
end

%{
Beta = 1e-5;

epsilon = calculateEpsilon(Beta,4,numberOfDataInEachPolytope);
disp('Max epsilon value:');
disp(epsilon);
%tempiN = find(epsilon == min(epsilon));

betas = calculateBeta(min(epsilon), 4, numberOfDataInEachPolytope);
disp("Beta")
disp(betas)
%}
%% 
%addpath("Extraction_of_relation_after_RTree/SupportingFunctions/ScenarioBasedApproach/");
% Example usage
d = 3;
N_array = [10, 20, 30];
epsilon = 0.1;

% Calculate beta for given epsilon
betas = calculateBeta(epsilon, d, N_array);
disp('Betas:');
disp(betas);

% Calculate epsilon for a given beta
beta_target = 0.5;
epsilons = calculateEpsilon(beta_target, d, N_array);
disp('Epsilons:');
disp(epsilons);

 %% 


for i = 1 : size(state_graph_transition_matrix,1)
    disp(["Transition from State  ", int2str(i), " to "]);
    disp(find(state_graph_transition_matrix(i,:)));
end



%% 



%% Checking max values
 tempH = zeros(size(solutionOptimal_star,1), 3);
 for i = 1 : size(solutionOptimal_star,1)
     % disp(solutionOptimal_star{i,1}(1,1));
     tempH(i,:) = [ i solutionOptimal_star{i,1}(1,1), size(LAMBDA{i,1},1) ];
 end

 sortedTempH = sortrows(tempH,2);

% [maxhValue, index] = max(tempH);
% disp(["maxhValue = " , int2str(maxhValue), " index = " , int2str(index)]);
%% 
%% Working with Scenario Based Approach

%% Start of Support Sub Samples Calculation
 support_subSample_LAMBDA = LAMBDA;
 support_subSample_lambda = lambda;
% 
 %optimizationParametersFor_scenarios = cell(size(LAMBDA,1), 1); %numLeaves

toleranceForHCompare = 0.000001;
 % i_leaf_counter = 1;
 for i_leaf_counter = 1 : size(LAMBDA,1) %numLeaves
    counterForLoop = 1;
     while counterForLoop <= size(support_subSample_LAMBDA{i_leaf_counter}, 1)
         disp(['Leaf-> ', int2str(i_leaf_counter), ' Index-> ', int2str(counterForLoop), 'Total in Leaf: ' , int2str(size(support_subSample_LAMBDA{i_leaf_counter}, 1))] );
        list_of_scenarios = 1 : size(support_subSample_LAMBDA{i_leaf_counter}, 1);
        list_of_scenarios(counterForLoop) = [];
        % disp(list_of_scenarios);
        [h, m , mo] = customOptimization.QuadraticConstraintOptimizer(support_subSample_LAMBDA{i_leaf_counter}(list_of_scenarios,:), support_subSample_lambda{i_leaf_counter}(list_of_scenarios,:));
        %solutionOptimal_star{i,1} = [ones(size(m,1),1)*h m mo];

        %[h, m , mo] = QuadraticConstraintOptimizer(support_subSample_LAMBDA{i_leaf_counter}(list_of_scenarios,:), support_subSample_lambda{i_leaf_counter}(list_of_scenarios));
        % if size(optimizationParametersFor_scenarios{1} , 1) == 0
        %     optimizationParametersFor_scenarios{1}(:,:) = [h m  mo];
        % else
        %optimizedParamAtStep = [h m' mo];
        %optimizationParametersFor_scenarios{i_leaf_counter} = [optimizationParametersFor_scenarios{i_leaf_counter}(:,:) ; optimizedParamAtStep];
        % end
        % optimizationSolutionAtStep = [h m' mo];
        % for tempIndex = 1 : size(optimizationSolutionAtStep, 2)
        %     optimizationParametersFor_scenarios{1}(size(optimizationParametersFor_scenarios{1},1)+1 , optimizationSolutionAtStep(1,tempIndex));
        % end
        % disp(h);
         if (abs(solutionOptimal_star{i_leaf_counter,1}(1,1) -  h) < toleranceForHCompare)
         % if all(abs(solutionOptimal_star(i_leaf_counter, :) - optimizedParamAtStep(1, :)) <= toleranceForHCompare )
            disp('Support Sub Sample found')
            support_subSample_LAMBDA{i_leaf_counter}(counterForLoop, :) = [];
            support_subSample_lambda{i_leaf_counter}(counterForLoop,:) = [];
        else
            counterForLoop = counterForLoop + 1;
            disp('Support Sub Sample not found')
        end
     end
 end
disp('Support Sub Sample Calculation Finised');

%% 
numberOfSupportSubSamples = zeros(size(support_subSample_LAMBDA,1), 3);
for i = 1 : size(numberOfSupportSubSamples,1)
    numberOfSupportSubSamples(i,1) = size(LAMBDA{i},1);
    numberOfSupportSubSamples(i,2) = size(support_subSample_LAMBDA{i},1);
    numberOfSupportSubSamples(i,3) = size(LAMBDA{i},1) - size(support_subSample_LAMBDA{i},1);
end

filename_SSS = strcat(dataFiles_folder_path, 'numberOfSupportSubSamples.xlsx');
writematrix(numberOfSupportSubSamples,filename_SSS,'Sheet',1,'Range','A1');
%% 


















 %% Temp Code for graphical Represeztation
numberOfHighH = size(find(tempH(:,2)<2), 1);
rect_normal = cell(numberOfHighH, 0);
rect_elv = cell(numberOfHighH, 0);
counter_for_rects = 1;
for i = 1 : size(solutionOptimal_star,1)
    h = solutionOptimal_star{i,1}(1,1);
    if h < 2
        rect_normal{counter_for_rects,1} = vertices_of_polytopes{i,1};
        rect_elv{counter_for_rects,1} = vertices_of_elevated_polytopes_before_span;
        counter_for_rects = counter_for_rects + 1;
    end
end

%% show elevated rectangles

if size(vertices_of_elevated_polytopes_before_span{1,1},1) == 4
    % Call the method on the instance
    customPlot.drawDual2DPolytopesWithVertices(vertices_of_polytopes, vertices_of_elevated_polytopes_before_span, 'color1', '#8f9119', 'color2', '#197b91', 'title', 'Evevated 2D-Polytopes');
    
elseif size(vertices_of_elevated_polytopes_before_span{1,1},1) == 8
    
    customPlot.drawDual3DPolytopesWithVertices(vertices_of_polytopes, vertices_of_elevated_polytopes_before_span, 'color1', '#8f9119', 'color2', '#197b91', 'title', 'Evevated 3D-Polytopes');
else
    customPlot.drawDual2DPolytopesWithVertices(vertices_of_polytopes, vertices_of_elevated_polytopes_before_span, 'title', 'More the 3D Polytopes in 2D');
end



%% 

sortedTempH = sortrows(tempH,2);

 [maxhValue, index] = max(tempH);
 disp(["maxhValue = " , int2str(maxhValue), " index = " , int2str(index)]);
 %% 


 for i = 1 : size(sortedTempH)
     disp(["Leaf Number = " , int2str(sortedTempH(i,1)), " H* = " , double(sortedTempH(i,2))])
 end
 
 %% Showing the leaves with high h
list_of_leaf_indexes_above_2h = sortedTempH( find(sortedTempH(:,2) > 2), 1);
rectangles_high_h = cell(size( list_of_leaf_indexes_above_2h, 1), 0);
for i = 1 : size( list_of_leaf_indexes_above_2h, 1)
    rectangles_high_h{i,1} = vertices_of_rectangles{list_of_leaf_indexes_above_2h(i),1};
    disp([ "Leaf #= ", int2str(list_of_leaf_indexes_above_2h(i)), " # of Samples " , size(LAMBDA{list_of_leaf_indexes_above_2h(i),1},1)])
end
plot_multiple_hyperrectangle(rectangles_high_h, 'color', '#197b91');

%% 
% counter_testing_trejectories_graphical = list_of_leaf_indexes_above_2h(1)
counter_testing_trejectories_graphical = 204;
one_LAMBDA_highH = LAMBDA{counter_testing_trejectories_graphical,1};
one_lambda_highH = lambda{counter_testing_trejectories_graphical,1};

one_trajectory_points_highH = [one_LAMBDA_highH one_lambda_highH];
indexs_for_time_data = find(ismember(trajectories_DataSet, one_trajectory_points_highH,'rows'));

%% 
points_for_lines = [ trajectories_time_DataSet(indexs_for_time_data,1) trajectories_DataSet(indexs_for_time_data,1) ...
    trajectories_time_DataSet(indexs_for_time_data,1)+0.1  trajectories_DataSet(indexs_for_time_data,4)];
points_for_lines_two = [ trajectories_time_DataSet(indexs_for_time_data,1) trajectories_DataSet(indexs_for_time_data,2) ...
    trajectories_time_DataSet(indexs_for_time_data,1)+0.1  trajectories_DataSet(indexs_for_time_data,5)];
points_for_lines_three = [ trajectories_time_DataSet(indexs_for_time_data,1) trajectories_DataSet(indexs_for_time_data,3) ...
    trajectories_time_DataSet(indexs_for_time_data,1)+0.1  trajectories_DataSet(indexs_for_time_data,6)];


figure;
hold on;
 % plot([1 3], [2 4]);
 for i = 1 : size(points_for_lines,1)
    % plot(points_for_lines(1,1:2), points_for_lines(1,3:4));
    % scatter(points_for_lines(i,1), points_for_lines(i,2),'MarkerEdgeColor','b');
    % scatter(points_for_lines(i,3), points_for_lines(i,4),'MarkerEdgeColor','b');
    % quiver(points_for_lines(i,1),points_for_lines(i,2),points_for_lines(i,3)-points_for_lines(i,1),points_for_lines(i,4)-points_for_lines(i,2), 'LineWidth',2,  'Color','b');
    drawArrow(points_for_lines(i,1:2), points_for_lines(i,3:4), 'b');
    drawArrow(points_for_lines_two(i,1:2), points_for_lines_two(i,3:4), 'r');
    drawArrow(points_for_lines_three(i,1:2), points_for_lines_three(i,3:4), 'g');
 end
 % scatter(trajectories_time_DataSet(indexs_for_time_data,1), trajectories_DataSet(indexs_for_time_data,1))
 % scatter(trajectories_time_DataSet(indexs_for_time_data,1), trajectories_DataSet(indexs_for_time_data,2))
 % scatter(trajectories_time_DataSet(indexs_for_time_data,1), trajectories_DataSet(indexs_for_time_data,3))
hold off;

function drawArrow(p1,p2,color)
% p1 = [2 3];                         % First Point
% p2 = [9 8];                         % Second Point
dp = p2-p1;                         % Difference
scatter(p1(1),p1(2),'MarkerEdgeColor',color);
scatter(p2(1),p2(2),'MarkerEdgeColor',color);
quiver(p1(1),p1(2),dp(1),dp(2),'LineWidth',2,  'Color',color);
end



%% test one rectangle
 
 test_index_for_sinagle_cube = 212;
 test_rectangles = cell(3,0);
 disp(["Leaf # = " , int2str(test_index_for_sinagle_cube), " H = ", double(solutionOptimal_star{test_index_for_sinagle_cube}(1,1))]);
 test_rectangles{1,1} = vertices_of_rectangles{test_index_for_sinagle_cube,1};
 test_rectangles{2,1} = vertices_of_elevated_rectangles_before_span{test_index_for_sinagle_cube,1};
 % test_rectangles{3,1} = vertices_of_elevated_rectangles{test_index_for_sinagle_cube,1};
 %plot_hyperrectangle(vertices_of_rectangles{test_index_for_sinagle_cube,1});
 % plot_hyperrectangle(vertices_of_elevated_rectangles_before_span{test_index_for_sinagle_cube,1});
 %plot_hyperrectangle(vertices_of_elevated_rectangles{test_index_for_sinagle_cube,1});
 plot_multiple_hyperrectangle(test_rectangles);
 % plot_hyperrectangle(add_span_to_hyperrectangle(vertices_of_rectangles{i,1}, 10));
  figure; hold on;
  scatter(LAMBDA{test_index_for_sinagle_cube}(:,1),LAMBDA{test_index_for_sinagle_cube}(:,2));
  hold off;
%% 

  one_Leaf_LAMBDA = cell(1, 0);
  one_Leaf_LAMBDA{1} = LAMBDA_test{test_index_for_sinagle_cube};
  one_Leaf_lambda = cell(1,0);
  one_Leaf_lambda{1} = lambda_test{test_index_for_sinagle_cube};
 
  NRMSE_Training = NRMSE_Calculation_3D(one_Leaf_LAMBDA,one_Leaf_lambda,solutionOptimal_star);

 disp(["NRMSE on one Leaf with code error = ", NRMSE_Training*100]);





%% showing normal and elivated rectangles
% plot_multiple_hyperrectangle(vertices_of_rectangles);
plot_multiple_hyperrectangle_overlapping(vertices_of_rectangles, vertices_of_elevated_rectangles_before_span, "color1",'#639119', 'color2','#197b91');
%% Calculation of hyperrectangles in complete state space 


space_bounds = zeros(size(minOfAllDimensions));
for i = 1 : size(space_bounds,2)/2
    space_bounds( (i*2) - 1 ) =  floor(minOfAllDimensions( (i*2) - 1 ));
    space_bounds ( i * 2 )  = ceil(maxOfAllDimensions( i*2 )); 
end

% As I have bounds of the hyperrectangles in constraints_on_leaves variable
  % smallest_rect_bounds = get_smallest_hyperrectangle(constraints_on_leaves);
 size_of_smallest_edges = get_smallest_edge_of_rectangle(constraints_on_leaves);
 smallest_uniform_rect = zeros(1, size(space_bounds,2)/2);
 smallest_uniform_rect(1:end) = size_of_smallest_edges;
 count_of_rect = count_hyperrectangles_in_space(space_bounds, smallest_uniform_rect);

 disp(["Complete State Space Bounds = ",  int2str(space_bounds)]);
 disp(["Smallest Hyper Rectangle = ", double(smallest_uniform_rect)]);
 disp(["Total Number of uniform hyper rectangles to fill the state space = ", int2str(count_of_rect) ]);

%% Rectangles Intersection Checking
state_graph_transition_matrix = zeros(numLeaves, numLeaves);
for i = 1 : numLeaves
    for j = 1 : numLeaves
        if check_hyperrectangles_intersection(vertices_of_elevated_rectangles{i,1}, vertices_of_rectangles{j,1})
            state_graph_transition_matrix(i,j) = 1;
        end
    end
end
%% Show Graph

% addpath('/Users/khalilulrehman/Academic/Phd Italy 2023_26/University of LAquila/Research Papers tasks/MatlabCodes/RegressionTree/ApproximateBisimulationRelation/V2/')
% showGraph(state_graph_transition_matrix);



