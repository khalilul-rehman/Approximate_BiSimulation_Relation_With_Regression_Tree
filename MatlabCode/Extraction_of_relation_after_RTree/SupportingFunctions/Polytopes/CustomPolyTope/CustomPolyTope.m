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






        function d = computeDistance(~, vertices1, vertices2)
            % Computes the minimum Euclidean distance between two convex polytopes.
            % If polytopes overlap, returns 0.
            
            % Check if polytopes overlap
            if any(CustomPolyTope.isInsideConvexHull(vertices1, vertices2)) || ...
               any(CustomPolyTope.isInsideConvexHull(vertices2, vertices1))
                d = 0;
                return;
            end
            
            % Compute minimum distance when polytopes do not overlap
            numVertices1 = size(vertices1, 1);
            numVertices2 = size(vertices2, 1);
            d = inf;
            
            % Compute minimum vertex-to-polytope distance
            for i = 1:numVertices1
                d = min(d, CustomPolyTope.pointToPolytopeDistance(vertices1(i, :), vertices2));
            end
            
            for j = 1:numVertices2
                d = min(d, CustomPolyTope.pointToPolytopeDistance(vertices2(j, :), vertices1));
            end
        end
        
       
   
        
    

    

    end


methods (Static)

     function inside = isInsideConvexHull(points, polytopeVertices)
            % Determines if a set of points is inside the convex hull of a polytope
            % Uses linear programming (linprog)
            
            numPoints = size(points, 1);
            numVertices = size(polytopeVertices, 1);
            numDims = size(polytopeVertices, 2);
            
            % Construct constraint matrices for convex hull formulation
            A = [polytopeVertices'; ones(1, numVertices)];
            lb = zeros(numVertices, 1);
            ub = ones(numVertices, 1);
            
            inside = false(numPoints, 1);
            
            % Check each point
            for i = 1:numPoints
                b = [points(i, :)'; 1];
                f = zeros(numVertices, 1);
                options = optimoptions('linprog', 'Display', 'off');
                [~, ~, exitflag] = linprog(f, [], [], A, b, lb, ub, options);
                
                if exitflag == 1
                    inside(i) = true;
                end
            end
        end
        
        function d = pointToPolytopeDistance(point, polytopeVertices)
            % Computes the minimum distance from a point to a convex polytope
            numVertices = size(polytopeVertices, 1);
            d = inf;
            
            % Compute Euclidean distance to each vertex
            for i = 1:numVertices
                d = min(d, norm(point - polytopeVertices(i, :)));
            end
        end


end
    


end
