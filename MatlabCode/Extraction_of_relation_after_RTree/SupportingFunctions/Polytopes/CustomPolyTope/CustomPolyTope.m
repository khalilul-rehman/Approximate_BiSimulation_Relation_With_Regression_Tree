classdef CustomPolyTope
    methods
        function obj = CustomPolyTope() % Constructor
        end
        
        function vertices = getVerticesFromBounds(obj, bounds)
            % Check if the input bounds is a 2D array and non-empty
            if ~ismatrix(bounds) || isempty(bounds)
                error('Input must be a non-empty 2D array.');
            end
            
            % Number of polytopes and dimensions
            numPolytopes = size(bounds, 1);
            numDimensions = size(bounds, 2) / 2;
            
            % Validate the bounds dimensions
            if mod(size(bounds, 2), 2) ~= 0
                error('Each polytope must have a min and max bound for each dimension.');
            end
            
            % Initialize the cell array to hold vertices of each polytope
            vertices = cell(numPolytopes, 1);
            
            % Iterate over each polytope to compute vertices
            for i = 1:numPolytopes
                % Extract bounds for the current polytope
                polytopeBounds = bounds(i, :);
                
                % Initialize the list of vertices
                vertexList = [];
                
                % Generate all combinations of vertices
                for j = 0:(2^numDimensions - 1)
                    vertex = zeros(1, numDimensions);
                    for k = 1:numDimensions
                        if bitand(j, 2^(k-1))
                            vertex(k) = polytopeBounds(2*k);
                        else
                            vertex(k) = polytopeBounds(2*k-1);
                        end
                    end
                    vertexList = [vertexList; vertex];
                end
                
                % Store vertices in the cell array
                vertices{i} = vertexList;
            end
        end

        function bounds = getBoundsFromVertices(~, vertices)
            % Check if the input is a non-empty matrix
            if isempty(vertices) || ~ismatrix(vertices)
                error('Input must be a non-empty 2D matrix.');
            end

            % Number of dimensions
            numDimensions = size(vertices, 2);

            % Initialize min and max bounds
            minBounds = min(vertices);
            maxBounds = max(vertices);

            % Combine min and max bounds
            bounds = zeros(1, 2 * numDimensions);
            bounds(1:2:end) = minBounds;
            bounds(2:2:end) = maxBounds;
        end


        function [X_opt, Xp_opt, min_distance] = minimize_polytope_distance_dual(~, V1, V2)
            % Inputs:
            % V1: Vertices of Polytope 1 (n x d matrix)
            % V2: Vertices of Polytope 2 (m x d matrix)
            %
            % Outputs:
            % X_opt: Closest point in Polytope 1
            % Xp_opt: Closest point in Polytope 2
            % min_distance: Minimum Euclidean distance
            
            [n, d] = size(V1);  % Polytope 1 has n vertices in d dimensions
            [m, ~] = size(V2);  % Polytope 2 has m vertices in d dimensions
        
            % Initial guess: Uniform distribution of weights
            beta0 = ones(n, 1) / n;
            alpha0 = ones(m, 1) / m;
            z0 = [beta0; alpha0];
        
            % Objective function: Squared Euclidean distance
            objective = @(z) norm((z(1:n)' * V1) - (z(n+1:end)' * V2))^2;
        
            % Lower and upper bounds for beta and alpha
            lb = zeros(n + m, 1);
            ub = ones(n + m, 1);
        
            % Linear equality constraints for convex combinations
            Aeq = [ones(1, n), zeros(1, m);  % sum(beta) = 1
                   zeros(1, n), ones(1, m)]; % sum(alpha) = 1
            beq = [1; 1];
        
            % Optimization settings
            options = optimoptions('fmincon', ...
                'Display', 'off', ...
                'OptimalityTolerance', 1e-8, ...
                'StepTolerance', 1e-8, ...
                'FunctionTolerance', 1e-8, ...
                'MaxIterations', 1000, ...
                'MaxFunctionEvaluations', 50000);
        
            % Solve optimization
            [z_opt, ~] = fmincon(objective, z0, [], [], Aeq, beq, lb, ub, [], options);
        
            % Extract optimal values
            beta_opt = z_opt(1:n);  % Weights for Polytope 1
            alpha_opt = z_opt(n+1:end);  % Weights for Polytope 2
            X_opt = beta_opt' * V1;  % Optimal point in Polytope 1
            Xp_opt = alpha_opt' * V2;  % Optimal point in Polytope 2
        
            % Compute minimum Euclidean distance
            min_distance = norm(X_opt - Xp_opt);
        end


        function [smallestPolytope, smallestEdgeLength, edgeVertices] = findSmallestPolytopeAndEdge(~, polytopes)
            % Find the smallest n-dimensional polytope and its smallest edge.
            %
            % Inputs:
            %   - polytopes: A cell array where each row contains an n-dimensional matrix
            %                representing the vertices of a polytope.
            %
            % Outputs:
            %   - smallestPolytope: The vertices of the smallest polytope.
            %   - smallestEdgeLength: The length of the smallest edge in the smallest polytope.
            %   - edgeVertices: The two vertices defining the smallest edge.
        
            numPolytopes = length(polytopes);
            polytopeVolumes = zeros(numPolytopes, 1);
        
            % Calculate the volume of each polytope
            for i = 1:numPolytopes
                vertices = polytopes{i};
                if size(vertices, 1) > size(vertices, 2)
                    vertices = vertices'; % Ensure vertices are columns for convhulln
                end
                try
                    K = convhulln(vertices'); % Compute convex hull
                    polytopeVolumes(i) = volumeConvexHull(vertices, K); % Compute volume
                catch
                    polytopeVolumes(i) = inf; % Handle degenerate cases
                end
            end
        
            % Find the smallest polytope based on volume
            [~, minIndex] = min(polytopeVolumes);
            smallestPolytope = polytopes{minIndex};
        
            % Compute the smallest edge of the smallest polytope
            numVertices = size(smallestPolytope, 1);
            smallestEdgeLength = inf;
            edgeVertices = [];
        
            for i = 1:numVertices
                for j = i+1:numVertices
                    edgeLength = norm(smallestPolytope(i, :) - smallestPolytope(j, :));
                    if edgeLength < smallestEdgeLength
                        smallestEdgeLength = edgeLength;
                        edgeVertices = [smallestPolytope(i, :); smallestPolytope(j, :)];
                    end
                end
            end
        end
        
        




        %{
        function [X_opt, Xp_opt, min_distance] = minimize_polytope_distance(~,V1, V2)
            % First polytope in the argument will be the uniform polytope,
            % uniform means axis aligned because its min and max values
            % have been used.
            % Inputs:
            % V1: Vertices of Polytope 1 (n x d matrix, where n is the number of vertices and d is the dimension)
            % V2: Vertices of Polytope 2 (m x d matrix, where m is the number of vertices and d is the dimension)
            
            % Get the number of vertices and dimensions
            [m, d] = size(V2);  % Polytope 2 has m vertices in d dimensions
            [n, ~] = size(V1);  % Polytope 1 has n vertices in d dimensions
            
            % Define the objective function
            objective = @(z) norm(z(1:d) - (z(d+1:end)' * V2));
            
            % Initial guess (midpoint of bounds)
            X0 = mean(V1, 1)';  % Midpoint of Polytope 1 (converted to column vector)
            alpha0 = ones(m, 1) / m;  % Uniform weights for Polytope 2 (column vector)
            z0 = [X0; alpha0];  % Initial guess for optimization variables
            
            % Lower and upper bounds for X and alpha
            lb_X = min(V1, [], 1)';  % Lower bounds for X (min of Polytope 1 vertices, column vector)
            ub_X = max(V1, [], 1)';  % Upper bounds for X (max of Polytope 1 vertices, column vector)
            lb_alpha = zeros(m, 1);  % Lower bounds for alpha (0 <= alpha_i <= 1)
            ub_alpha = ones(m, 1);   % Upper bounds for alpha (0 <= alpha_i <= 1)
            lb = [lb_X; lb_alpha];  % Combined lower bounds
            ub = [ub_X; ub_alpha];  % Combined upper bounds
            
            % Linear equality constraint for convex combination (sum(alpha) = 1)
            Aeq = [zeros(1, d), ones(1, m)];  % [0 ... 0 | 1 ... 1]
            beq = 1;  % sum(alpha) = 1
            
            % No inequality constraints
            A = [];
            b = [];
            
            % Perform the optimization using fmincon with tighter tolerances
            options = optimoptions('fmincon', ...
                'Display', 'off', ...
                'OptimalityTolerance', 1e-8, ...  % Tighter optimality tolerance
                'StepTolerance', 1e-2, ...        % 1e-8 Tighter step tolerance
                'FunctionTolerance', 1e-2, ...    % 1e-8 Tighter function tolerance
                'MaxIterations', 10000, ...       %1000 % Increase maximum iterations
                'MaxFunctionEvaluations', 50000);  % Increase maximum function evaluations
            
            [z_opt, ~] = fmincon(objective, z0, A, b, Aeq, beq, lb, ub, [], options);
            
            % Extract the optimized values
            X_opt = z_opt(1:d);  % Optimized point in Polytope 1
            alpha_opt = z_opt(d+1:end);  % Optimized weights for Polytope 2
            Xp_opt = alpha_opt' * V2;  % Optimized point in Polytope 2 (convex combination)
            
            % Calculate the minimum Euclidean distance
            min_distance = norm(X_opt - Xp_opt);
        end
        %}

        % function minDist = findMinDistanceFromBounds(~, vertices1, vertices2)
        %     % FINDMINDISTANCEFROMBOUNDS computes the minimum distance between two polytopes
        %     % using their bounding boxes.
        % 
        %     % Number of dimensions
        %     n = size(vertices1, 2);
        % 
        %     % Compute bounding box min and max for each dimension
        %     minBound1 = min(vertices1, [], 1);
        %     maxBound1 = max(vertices1, [], 1);
        %     minBound2 = min(vertices2, [], 1);
        %     maxBound2 = max(vertices2, [], 1);
        % 
        %     % Define optimization variables (points inside the polytopes)
        %     x = optimvar('x', n);
        %     y = optimvar('y', n);
        % 
        %     % Objective function: minimize Euclidean distance
        %     obj = sum((x - y).^2);
        %     prob = optimproblem('Objective', obj);
        % 
        %     % Constraints: x and y should be inside their respective bounds
        %     prob.Constraints.x_lower = minBound1' <= x;
        %     prob.Constraints.x_upper = x <= maxBound1';
        %     prob.Constraints.y_lower = minBound2' <= y;
        %     prob.Constraints.y_upper = y <= maxBound2';
        % 
        %     % Initial guess (centroids of the bounding boxes)
        %     x0.x = (minBound1 + maxBound1) / 2;
        %     x0.y = (minBound2 + maxBound2) / 2;
        % 
        %     % Solve the optimization problem
        %     sol = solve(prob, x0);
        % 
        %     disp("Point 1 =>");
        %     disp(sol.x);
        %     disp("Point 2");
        %     disp(sol.y);
        % 
        %     % Compute final minimum distance
        %     minDist = norm(sol.x - sol.y);
        % end



    %{
       function [x_opt, y_opt, min_distance] = minimize_polytope_distance(~,bounds_1, bounds_2)
            % Ensure bounds_1 and bounds_2 have the same length
            if length(bounds_1) ~= length(bounds_2)
                error('bounds_1 and bounds_2 must have the same length.');
            end
            
            % Number of dimensions
            n = length(bounds_1) / 2;
            
            % Extract bounds for x and y
            x_min = bounds_1(1:2:end);  % Odd indices: min values for x
            x_max = bounds_1(2:2:end);  % Even indices: max values for x
            y_min = bounds_2(1:2:end);  % Odd indices: min values for y
            y_max = bounds_2(2:2:end);  % Even indices: max values for y
            
            % Initial guess (midpoint of bounds)
            x0 = (x_min + x_max) / 2;
            y0 = (y_min + y_max) / 2;
            z0 = [x0; y0];
            
            % Lower and upper bounds for fmincon
            lb = [x_min; y_min];  % Concatenate lower bounds for x and y
            ub = [x_max; y_max];  % Concatenate upper bounds for x and y
            
            % Objective function: Squared Euclidean distance
            objective = @(z) sum((z(1:n) - z(n+1:end)).^2);
            
            % No linear inequality or equality constraints
            A = [];
            b = [];
            Aeq = [];
            beq = [];
            
            % Perform the optimization using fmincon
            options = optimoptions('fmincon', 'Display', 'off');
            [z_opt, min_squared_distance] = fmincon(objective, z0, A, b, Aeq, beq, lb, ub, [], options);
            
            % Extract the optimized values
            x_opt = z_opt(1:n);
            y_opt = z_opt(n+1:end);
            
            % Calculate the actual Euclidean distance
            min_distance = sqrt(min_squared_distance);

            % Tolerance for considering the distance as zero
            tolerance = 1e-6;
            if min_distance < tolerance
                min_distance = 0;
            end
        end

    %}




    
        % function d = computeDistance(~, vertices1, vertices2)
        %     % Computes the minimum Euclidean distance between two convex polytopes.
        %     % If polytopes overlap, returns 0.
        % 
        %     % Check if polytopes overlap
        %     if any(CustomPolyTope.isInsideConvexHull(vertices1, vertices2)) || ...
        %        any(CustomPolyTope.isInsideConvexHull(vertices2, vertices1))
        %         d = 0;
        %         return;
        %     end
        % 
        %     % Compute minimum distance when polytopes do not overlap
        %     numVertices1 = size(vertices1, 1);
        %     numVertices2 = size(vertices2, 1);
        %     d = inf;
        % 
        %     % Compute minimum vertex-to-polytope distance
        %     for i = 1:numVertices1
        %         d = min(d, CustomPolyTope.pointToPolytopeDistance(vertices1(i, :), vertices2));
        %     end
        % 
        %     for j = 1:numVertices2
        %         d = min(d, CustomPolyTope.pointToPolytopeDistance(vertices2(j, :), vertices1));
        %     end
        % end
        % 
    
   
        
    

    

    end

    methods (Static)
        function volume = volumeConvexHull(vertices, K)
            % Compute the volume of a convex hull using the simplices defined by K
            volume = 0;
            for i = 1:size(K, 1)
                simplex = vertices(K(i, :), :); % Get the vertices of the simplex
                volume = volume + abs(det([simplex, ones(size(simplex, 1), 1)])) / factorial(size(simplex, 1) - 1);
            end
        end
    end

% methods (Static)
% 
%      function inside = isInsideConvexHull(points, polytopeVertices)
%             % Determines if a set of points is inside the convex hull of a polytope
%             % Uses linear programming (linprog)
% 
%             numPoints = size(points, 1);
%             numVertices = size(polytopeVertices, 1);
%             numDims = size(polytopeVertices, 2);
% 
%             % Construct constraint matrices for convex hull formulation
%             A = [polytopeVertices'; ones(1, numVertices)];
%             lb = zeros(numVertices, 1);
%             ub = ones(numVertices, 1);
% 
%             inside = false(numPoints, 1);
% 
%             % Check each point
%             for i = 1:numPoints
%                 b = [points(i, :)'; 1];
%                 f = zeros(numVertices, 1);
%                 options = optimoptions('linprog', 'Display', 'off');
%                 [~, ~, exitflag] = linprog(f, [], [], A, b, lb, ub, options);
% 
%                 if exitflag == 1
%                     inside(i) = true;
%                 end
%             end
%         end
% 
%         function d = pointToPolytopeDistance(point, polytopeVertices)
%             % Computes the minimum distance from a point to a convex polytope
%             numVertices = size(polytopeVertices, 1);
%             d = inf;
% 
%             % Compute Euclidean distance to each vertex
%             for i = 1:numVertices
%                 d = min(d, norm(point - polytopeVertices(i, :)));
%             end
%         end
% 
% 
% end



end
