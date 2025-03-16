
classdef CustomPlotClass
    methods
        function obj = Plot() %Plot(val)
        end
        
        function drawConstraints(obj, constraints_data)
            % Check if the input is a 2D array
            if ismatrix(constraints_data) && ~isempty(constraints_data)
                % Get the number of columns
                [~, numCols] = size(constraints_data);
                
                % Conditional statements based on the number of columns
                if numCols == 4
                    disp('The array has 4 columns.');
                    % Add your code here for 4 columns case
                elseif numCols == 6
                    disp('The array has 6 columns.');
                    % Add your code here for 6 columns case
                else
                    disp(['The array has ', num2str(numCols), ' columns.']);
                end
            else
                error('Input must be a non-empty 2D array.');
            end
        end


        function drawConstraintsIn2D(obj, min_max_bounds, varargin)
            % Parse optional arguments
            p = inputParser;
            addRequired(p, 'constraints_data');
            addParameter(p, 'title', 'Rectangular Representation of Constraints', @ischar);
            addParameter(p, 'dataPoints2D', cell(0,0));

            parse(p, min_max_bounds, varargin{:});
            
            % Retrieve parsed arguments
            min_max_bounds = p.Results.constraints_data;
            plotTitle = p.Results.title;
            dataPoints2D = p.Results.dataPoints2D;

            % Plot the data and rectangles
            numberofelements = size(min_max_bounds, 1);
            colors = lines(numberofelements);
            figure; hold on;
            for i = 1:numberofelements
                rect = min_max_bounds(i, :);
                % Draw rectangle
                rectangle('Position', [rect(1), rect(3), rect(2) - rect(1), rect(4) - rect(3)], ...
                          'EdgeColor', colors(i, :), 'LineWidth', 2);
                % Scatter plot for trajectories if provided
                if ~isempty(dataPoints2D)
                    scatter(dataPoints2D{i}(:, 1), dataPoints2D{i}(:, 2));
                end
            end
            
            % Set plot title and axis labels
            title(plotTitle);
            xlabel('Feature 1');
            ylabel('Feature 2');
            hold off;
        end



        function draw2DPolytopesWithVertices(obj, verticesCellArray, varargin)
            % Parse optional arguments
            p = inputParser;
            addRequired(p, 'verticesCellArray');
            addParameter(p, 'title', '2D Polytopes', @ischar);
            addParameter(p, 'color', [], @(x) isempty(x) || ischar(x) || isnumeric(x));
            addParameter(p, 'dataPoints2D', cell(0, 0));
            
            parse(p, verticesCellArray, varargin{:});
            
            % Retrieve parsed arguments
            verticesCellArray = p.Results.verticesCellArray;
            plotTitle = p.Results.title;
            plotColor = p.Results.color;
            dataPoints2D = p.Results.dataPoints2D;
            
            % Create a new figure
            figure; hold on;
            % axis equal;
            
            % Generate colors using lines function if color is not provided
            if isempty(plotColor)
                colors = lines(length(verticesCellArray));
            end
            
            % Iterate over each set of vertices
            for i = 1:length(verticesCellArray)
                polytopeVertices = verticesCellArray{i};
                
                % Determine the order of vertices using convhull
                K = convhull(polytopeVertices(:, 1), polytopeVertices(:, 2));
                
                % Use colors from lines function if color is not provided
                if isempty(plotColor)
                    edgeColor = colors(i, :);
                else
                    edgeColor = plotColor;
                end
                
                % Plot the polytope using plot function
                plot(polytopeVertices(K, 1), polytopeVertices(K, 2), 'Color', edgeColor, 'LineWidth', 2);
                
                % Scatter plot for data points if provided
                if ~isempty(dataPoints2D)
                    scatter(dataPoints2D{i}(:, 1), dataPoints2D{i}(:, 2));
                end
            end
            
            % Set labels and title
            xlabel('X');
            ylabel('Y');
            title(plotTitle);
            
            hold off;
        end

         function drawDual2DPolytopesWithVertices(obj, verticesCellArray1, verticesCellArray2, varargin)
            % Parse optional arguments
            p = inputParser;
            addRequired(p, 'verticesCellArray1');
            addRequired(p, 'verticesCellArray2');
            addParameter(p, 'title', '2D Polytopes', @ischar);
            addParameter(p, 'color1', [], @(x) isempty(x) || ischar(x) || isnumeric(x));
            addParameter(p, 'color2', [], @(x) isempty(x) || ischar(x) || isnumeric(x));
            
            
            parse(p, verticesCellArray1, verticesCellArray2, varargin{:});
            
            % Retrieve parsed arguments
            verticesCellArray1 = p.Results.verticesCellArray1;
            verticesCellArray2 = p.Results.verticesCellArray2;
            plotTitle = p.Results.title;
            plotColor1 = p.Results.color1;
            plotColor2 = p.Results.color2;
            
            
            % Create a new figure
            figure; hold on;
            % axis equal;
            
            % % Generate random colors if color is not provided
            if isempty(plotColor1)
                plotColor1 = rand(1, 3);
            end
            if isempty(plotColor2)
                plotColor2 = rand(1, 3);
            end
            
            % Iterate over each set of vertices
            for i = 1:length(verticesCellArray1)
                % 1st 
                polytopeVertices1 = verticesCellArray1{i};
                
                % Determine the order of vertices using convhull
                K = convhull(polytopeVertices1(:, 1), polytopeVertices1(:, 2));
                
                % Plot the polytope using plot function
                plot(polytopeVertices1(K, 1), polytopeVertices1(K, 2), 'Color', plotColor1, 'LineWidth', 2);
                
                % 2nd
                polytopeVertices2 = verticesCellArray2{i};
                
                % Determine the order of vertices using convhull
                K2 = convhull(polytopeVertices2(:, 1), polytopeVertices2(:, 2));
                
                % Plot the polytope using plot function
                plot(polytopeVertices2(K2, 1), polytopeVertices2(K2, 2), 'Color', plotColor2, 'LineWidth', 2);
                
                
            end
            
            % Set labels and title
            xlabel('X');
            ylabel('Y');
            title(plotTitle);
            
            hold off;
        end









        function draw3DPolytopesWithVertices(obj, vertices, varargin)
            % Parse optional arguments
            p = inputParser;
            addRequired(p, 'vertices');
            addParameter(p, 'color', [], @(x) isempty(x) || ischar(x) || isnumeric(x));
            addParameter(p, 'title', '3D Polytopes', @ischar);
            
            parse(p, vertices, varargin{:});
            
            % Retrieve parsed arguments
            vertices = p.Results.vertices;
            plotColor = p.Results.color;
            plotTitle = p.Results.title;
            
            % Check if the input is a cell array and non-empty
            if ~iscell(vertices) || isempty(vertices)
                error('Input must be a non-empty cell array.');
            end
            
            % Create a new figure
            figure; hold on;
            view(3); % Set view to 3D
            % axis equal;
            
            % Generate random colors if color is not provided
            if isempty(plotColor)
                colors = lines(length(vertices));
            end
            
            % Iterate over each set of vertices and draw the polytope
            for i = 1:length(vertices)
                polytopeVertices = vertices{i}(1:8,1:3);
                % Define faces of the polytope (assumes rectangular polytope)
                faces = [1 2 4 3; 1 2 6 5; 1 3 7 5; 2 4 8 6; 3 4 8 7; 5 6 8 7];
                
                % Draw the polytope using patch function
                if isempty(plotColor)
                    faceColor = colors(i, :);
                else
                    faceColor = plotColor;
                end
                
                patch('Vertices', polytopeVertices, 'Faces', faces, ...
                      'FaceColor', faceColor, 'EdgeColor', 'k', 'FaceAlpha', 0.3);
            end
            
            % Set labels and title
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title(plotTitle);
            
            hold off;
        end


         function drawDual3DPolytopesWithVertices(obj, vertices1, vertices2, varargin)
            % Parse optional arguments
            p = inputParser;
            addRequired(p, 'vertices1');
            addRequired(p, 'vertices2');
            addParameter(p, 'color1', [], @(x) isempty(x) || ischar(x) || isnumeric(x));
            addParameter(p, 'color2', [], @(x) isempty(x) || ischar(x) || isnumeric(x));
            addParameter(p, 'title', '3D Polytopes', @ischar);
            
            parse(p, vertices1, vertices2, varargin{:});
            
            % Retrieve parsed arguments
            vertices1 = p.Results.vertices1;
            vertices2 = p.Results.vertices2;
            plotColor1 = p.Results.color1;
            plotColor2 = p.Results.color2;
            plotTitle = p.Results.title;
            
            % Check if the input is a cell array and non-empty
            if ~iscell(vertices1) || isempty(vertices1)
                error('Input must be a non-empty cell array.');
            end
            if ~iscell(vertices2) || isempty(vertices2)
                error('Input must be a non-empty cell array.');
            end
            
            % Create a new figure
            figure; hold on;
            view(3); % Set view to 3D
            % axis equal;
            
            % Generate random colors if color is not provided
            if isempty(plotColor1)
                plotColor1 = rand(1,3);
            end
            if isempty(plotColor2)
                plotColor2 = rand(1,3);
            end
            
            % Iterate over each set of vertices and draw the polytope
            for i = 1:length(vertices1)
                polytopeVertices1 = vertices1{i}(:,1:3);
                % Define faces of the polytope (assumes rectangular polytope)
                faces1 = [1 2 4 3; 1 2 6 5; 1 3 7 5; 2 4 8 6; 3 4 8 7; 5 6 8 7];
                
                % Draw the polytope using patch function
                
                
                patch('Vertices', polytopeVertices1, 'Faces', faces1, ...
                      'FaceColor', plotColor1, 'EdgeColor', 'k', 'FaceAlpha', 0.3);
            end

            for i = 1:length(vertices2)
                polytopeVertices2 = vertices2{i}(:,1:3);
                % Define faces of the polytope (assumes rectangular polytope)
                faces2 = [1 2 4 3; 1 2 6 5; 1 3 7 5; 2 4 8 6; 3 4 8 7; 5 6 8 7];
                
                % Draw the polytope using patch function
                
                
                patch('Vertices', polytopeVertices2, 'Faces', faces2, ...
                      'FaceColor', plotColor2, 'EdgeColor', 'k', 'FaceAlpha', 0.3);
            end
            
            % Set labels and title
            xlabel('X');
            ylabel('Y');
            zlabel('Z');
            title(plotTitle);
            
            hold off;
        end




        


    end
end
