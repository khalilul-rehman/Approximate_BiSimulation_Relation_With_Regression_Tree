classdef CustomOptimization
    methods
        function obj = CustomOptimization() % Constructor
        end
        
        function [h, m, mo] = QuadraticConstraintOptimizer(obj, x, y)
            % QuadraticConstraintOptimizer solves a quadratic constraint optimization problem
            % for a dataset (x, y) where y may have multiple dimensions.
            %
            % Inputs:
            %   - x: Input data matrix (size: N × d)
            %   - y: Target data matrix (size: N × p)
            %
            % Outputs:
            %   - h: Optimal objective value
            %   - m: Optimal weight matrix (size: p × d)
            %   - mo: Optimal offset vector (size: p × 1)

            % Get dimensions
            [N, d] = size(x);  % Number of samples (N) and input features (d)
            p = size(y, 2);    % Output dimensions (p)

            % Initial guess for [h, m (p × d), mo (p × 1)]
            x0 = [1; zeros(p * d, 1); zeros(p, 1)]; 

            % Objective function: minimize h
            objective = @(vars) vars(1); % h is the first variable

            % Nonlinear constraint function
            nonlinear_constraint = @(vars) obj.nonlinear_constraints(vars, x, y, p, d);

            % Solve using fmincon
            options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'interior-point');
            [solution, fval] = fmincon(objective, x0, [], [], [], [], [], [], nonlinear_constraint, options);

            % Extract optimized values
            h = solution(1); % h is the first element
            m = reshape(solution(2:p*d+1), p, d); % Reshape to p × d
            mo = solution(p*d+2:end); % Last p elements as vector
        end
        
        function [c, ceq] = nonlinear_constraints(obj, vars, x, y, p, d)
            % Computes the nonlinear constraints for optimization.
            %
            % Inputs:
            %   - vars: Optimization variables (h, m, mo)
            %   - x: Input data (N × d)
            %   - y: Target data (N × p)
            %   - p: Output dimension
            %   - d: Input feature dimension
            %
            % Outputs:
            %   - c: Inequality constraints (each (y - m*x - mo)^2 ≤ h)
            %   - ceq: Empty (no equality constraints)

            h = vars(1); % Extract h
            m = reshape(vars(2:p*d+1), p, d); % Extract weight matrix (p × d)
            mo = vars(p*d+2:end); % Extract offset vector (p × 1)

            % Compute constraint for each (x_i, y_i)
            N = size(x, 1);
            c = zeros(N * p, 1); % Constraints for all samples & outputs

            for i = 1:N
                for j = 1:p
                    % Quadratic constraint: (y_ij - m_j * x_i - mo_j)^2 ≤ h
                    c((i-1)*p + j) = (y(i, j) - m(j, :) * x(i, :)' - mo(j))^2 - h;
                end
            end

            ceq = []; % No equality constraints
        end


        function [NRMSE] = NRMSE_Calculation(obj, LAMBDA, lambda, QCP_Solution_For_Prediction)
            % NRMSE_Calculation computes the Normalized Root Mean Square Error (NRMSE)
            % for multivariate outputs when QCP_Solution_For_Prediction is a cell array.
            %
            % Inputs:
            %   - LAMBDA: Cell array where each cell contains (N_i x d) feature matrix
            %   - lambda: Cell array where each cell contains (N_i x p) true output matrix
            %   - QCP_Solution_For_Prediction: Cell array where each cell contains 
            %     a vector [h, m (p x d), mo (p x 1)] for the corresponding leaf
            %
            % Output:
            %   - NRMSE: Normalized RMSE value averaged over all leaf nodes

            % Get the output dimensions from the first non-empty lambda
            firstNonEmpty = find(~cellfun('isempty', lambda), 1);
            if isempty(firstNonEmpty)
                error('All input cells are empty.');
            end
            p = size(lambda{firstNonEmpty}, 2);
            
            % Initialize errors with correct dimensions
            errors_with_QP = zeros(size(LAMBDA, 1), p); % Initialize errors
            leaf_errors = cell(size(LAMBDA, 1), 1);     % Store errors for each leaf

            for i = 1:size(LAMBDA, 1)
                if isempty(lambda{i}) || isempty(LAMBDA{i}) || isempty(QCP_Solution_For_Prediction{i})
                    continue;  % Skip empty leaves
                end
                
                [num_samples, p] = size(lambda{i}); % Number of samples and output dimensions
                d = size(LAMBDA{i}, 2); % Number of features
                
                leaf_errors{i} = zeros(num_samples, p); % Initialize error matrix
                
                % Extract m (p x d) and mo (p x 1) from cell QCP_Solution_For_Prediction{i}
                solution_vector = QCP_Solution_For_Prediction{i}(1:p, 2:end); % Extract solution vector
                m = solution_vector(1:p, 1:d); % Reshape to (p × d)
                mo = solution_vector(1:p, d+1:end); % Extract mo (p × 1)
                
                % Compute prediction error for each sample in the leaf
                for j = 1:num_samples
                    y_pred = (m * LAMBDA{i}(j, :)')' + mo'; % Predicted y (p × 1)
                    leaf_errors{i}(j, :) = lambda{i}(j, :) - y_pred; % Store error (1 × p)
                end
                
                if size(lambda{i}, 1) >= 1
                    % Compute RMSE for each output dimension
                    RMSE = sqrt(mean(leaf_errors{i}.^2) ./ size(lambda{i}, 1)); % RMSE (1 × p)
                    errors_with_QP(i, :) = RMSE ./ mean(lambda{i}); % Average over output dimensions
                end
            end

            % Compute final NRMSE by averaging over non-empty leaves
            num_non_empty_Leafs = sum(~cellfun(@isempty, lambda));
            if num_non_empty_Leafs > 0
                NRMSE = sum(errors_with_QP) / num_non_empty_Leafs;
            else
                NRMSE = NaN;  % Return NaN if no valid leaves exist
            end
        end



    end
end
