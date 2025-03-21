addpath('SupportingFunctions/Polytopes/CustomPolyTope/');
addpath('SupportingFunctions/Polytopes/Plot');
addpath('SupportingFunctions/');


customPolytope = CustomPolyTope();
customPlot = CustomPlotClass();





%% 

tempCellArray = cell(2,0);
tempCellArray{1,1} = [1 3; 4 2; 3 -1; 0 -2]; %vertices_of_elevated_polytopes_before_span{1,1};
% customPlot.draw2DPolytopesWithVertices(tempCellArray); 



tempCellArray{2,1} = [-2 -2; -2 0; -3 4; 1 5]; %vertices_of_elevated_polytopes_before_span{1,1};
% temp2CellArray = cell(1,0);
% temp2CellArray = customPolytope.getVerticesFromBounds(customPolytope.getBoundsFromVertices(tempCellArray{1,1}));

customPlot.draw2DPolytopesWithVertices(tempCellArray);  
%% 


function minDist = computePolytopeDistance(vertices1, vertices2)
    % COMPUTEPOLYTOPEDISTANCE: Finds the minimum distance between two n-D polytopes
    % using convex optimization.
    
    % Number of dimensions
    n = size(vertices1, 2);

    % Define optimization variables (points inside polytopes)
    x = optimvar('x', n);
    y = optimvar('y', n);

    % Define the objective function: minimize Euclidean distance ||x - y||
    obj = sum((x - y).^2);
    prob = optimproblem('Objective', obj);

    % Compute convex hulls (H-representation)
    K1 = convhull(vertices1); % Indices of convex hull for polytope 1
    K2 = convhull(vertices2); % Indices of convex hull for polytope 2

    % Convert convex hulls to linear constraints Ax ≤ b
    [A1, b1] = vert2con(vertices1(K1, :));
    [A2, b2] = vert2con(vertices2(K2, :));

    % Constraints: Ensure x and y are inside their respective polytopes
    prob.Constraints.x_inside = A1 * x <= b1;
    prob.Constraints.y_inside = A2 * y <= b2;

    % Initial guess (centroid of convex hulls)
    x0.x = mean(vertices1, 1)';
    x0.y = mean(vertices2, 1)';

    % Solve the optimization problem
    sol = solve(prob, x0);

    % Compute final minimum distance
    minDist = norm(sol.x - sol.y);
end
function [A, b] = vert2con(V)
    % VERT2CON Converts vertices to inequality constraints (A*x ≤ b)
    % Input:
    %   - V: (m x n) matrix of m vertices in n-dimensional space
    % Output:
    %   - A: Matrix representing inequality constraints (A * x ≤ b)
    %   - b: Right-hand side vector for inequality constraints

    % Compute convex hull
    K = convhulln(V);

    % Compute normal vectors
    numFacets = size(K, 1);
    numDims = size(V, 2);
    A = zeros(numFacets, numDims);
    b = zeros(numFacets, 1);

    for i = 1:numFacets
        % Get the facet vertices
        F = V(K(i, :), :);

        % Compute normal vector using cross product
        normal = null(F(2:end, :) - F(1, :))';
        if size(normal, 1) == 1
            normal = normal';
        end

        % Ensure normal is outward-facing
        d = mean(F, 1) * normal;
        if d < 0
            normal = -normal;
            d = -d;
        end

        % Store constraints
        A(i, :) = normal';
        b(i) = d;
    end
end

%% 
computePolytopeDistance(tempCellArray{1,1},tempCellArray{2,1})


%% 



state_graph_transition_matrix = zeros(numLeaves, numLeaves);
for i = 1 : numLeaves
    h = solutionOptimal_star{i}(1,1);
    for j = 1 : numLeaves
        % if  h < 2
            v1 = vertices_of_polytopes{j,1};
            v2 = vertices_of_elevated_polytopes_before_span{i,1};
            [P1_opt, P2_opt, min_distance] = customPolytope.minimize_polytope_distance_dual(v1,v2);
            
            plotPolytopWithDistance(vertices_of_polytopes{j,1}, vertices_of_elevated_polytopes_before_span{i,1}, P1_opt, P2_opt, min_distance);
            if min_distance < h   
                disp([min_distance, "" , h]);
                state_graph_transition_matrix(i,j) = 1;
            end
        % end
    end
end

%% 

showGraph(state_graph_transition_matrix);
%% 



function plotPolytopWithDistance(v1, v2, X_opt, Xp_opt, min_distance)

figure; hold on; axis padded;


 % Iterate over each set of vertices
            % for i = 1:length(verticesCellArray)
                  %verticesCellArray{i};
                
                  if size(v1,2) == 2
                    plotVertives(v1);
                    plotVertives(v2);
                    scatter([X_opt(1) Xp_opt(1)], [X_opt(2) Xp_opt(2)], "filled");
                  else
                      % Define faces of the polytope (assumes rectangular polytope)
                    % faces1 = [1 2 4 3; 1 2 6 5; 1 3 7 5; 2 4 8 6; 3 4 8 7; 5 6 8 7];
                    % patch('Vertices', v1, 'Faces', faces1, ...
                    %    'EdgeColor', 'k', 'FaceAlpha', 0.3);

                    % patch('Vertices', v2, 'Faces', faces1, ...
                    %    'EdgeColor', 'k', 'FaceAlpha', 0.3);

                    scatter3(X_opt(1), X_opt(2) , X_opt(3), 'filled', 'MarkerEdgeColor','k' );
                    scatter3(Xp_opt(1), Xp_opt(2) , Xp_opt(3), 'filled','MarkerEdgeColor','k' );
                    % scatter3([X_opt(1) Xp_opt(1)], [X_opt(2) Xp_opt(2)], [X_opt(3) Xp_opt(3)]);

                    % scatter3(v1(:, 1), v1(:, 2), v1(:, 3), 'filled');
                    trisurf(convhull(v1), v1(:, 1), v1(:, 2), v1(:, 3), 'FaceAlpha', 0.5);

                    % scatter3(v2(:, 1), v2(:, 2), v2(:, 3), 'filled');
                    trisurf(convhull(v2), v2(:, 1), v2(:, 2), v2(:, 3), 'FaceAlpha', 0.5);
                  end
                
                
               
                


            % Set labels and title
            xlabel('X');
            ylabel('Y');
            title( strcat('Distance = ' , num2str( min_distance ) ));

hold off;



end




function plotVertives(polytopeVertices)
                    % Determine the order of vertices using convhull
                    K = convhull(polytopeVertices(:, 1), polytopeVertices(:, 2));
                    
                    plot(polytopeVertices(K, 1), polytopeVertices(K, 2), 'LineWidth', 2); %, 'Color', edgeColor
                end

                
 %%
 %vertices_of_polytopes, vertices_of_elevated_polytopes_before_span

 new_vertices_of_polytopes = cell(size(vertices_of_polytopes,1),0);
 new_vertices_of_elevated_polytopes_before_span = cell(size(vertices_of_elevated_polytopes_before_span,1),0);
 for tempi = 1 : size(vertices_of_polytopes,1)
     new_vertices_of_polytopes{tempi,1} = vertices_of_polytopes{tempi,1}(:,[1 2 4]);
     new_vertices_of_elevated_polytopes_before_span{tempi,1} = vertices_of_elevated_polytopes_before_span{tempi,1}(:,[1 2 4]);
 end

 customPlot.drawDual3DPolytopesWithVertices(new_vertices_of_polytopes,new_vertices_of_elevated_polytopes_before_span);



 %% 
%function nrmse = compareVelocities(trajectories)
    % Input:
    %   - trajectories: Matrix with columns as follows:
    %       Columns 1-2: First position (x1, y1)
    %       Columns 3-4: Velocity at first position (vx1, vy1)
    %       Columns 5-6: Second position (x2, y2)
    %       Columns 7-8: Velocity at second position (vx2, vy2)
    %   - Time difference between positions: 0.1 seconds

    % Time difference
    dt = 0.1; % Time difference between positions (in seconds)

    % Extract positions and velocities
    pos1 = trajectories(:, 1:2); % First position (x1, y1)
    vel1 = trajectories(:, 3:4); % Velocity at first position (vx1, vy1)
    pos2 = trajectories(:, 5:6); % Second position (x2, y2)
    vel2 = trajectories(:, 7:8); % Velocity at second position (vx2, vy2)

    % Calculate velocity using finite differences
    calculated_velocity = (pos2 - pos1) / dt; % Finite difference velocity

    % Calculate average velocity from the given velocities
    average_velocity = (vel1 + vel2) / 2; % Average of vel1 and vel2

    % Calculate NRMSE (Normalized Root Mean Square Error)
    error = calculated_velocity - average_velocity; % Error between calculated and average velocity
    mse = mean(error.^2, 'all'); % Mean Squared Error (MSE)
    rmse = sqrt(mse); % Root Mean Squared Error (RMSE)
    range = max(average_velocity, [], 'all') - min(average_velocity, [], 'all'); % Range of average velocity
    nrmse = rmse / range; % Normalized RMSE

    % Display results
    %{
    fprintf('Calculated Velocity:\n');
    disp(calculated_velocity);
    fprintf('Average Velocity:\n');
    disp(average_velocity);
    %}
    fprintf('NRMSE between calculated and average velocity: %.4f\n', nrmse);
%end

%compareVelocities(trajectories)