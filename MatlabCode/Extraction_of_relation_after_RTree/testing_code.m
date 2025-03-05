addpath('SupportingFunctions/Polytopes/CustomPolyTope/');
addpath('SupportingFunctions/Polytopes/Plot');



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
