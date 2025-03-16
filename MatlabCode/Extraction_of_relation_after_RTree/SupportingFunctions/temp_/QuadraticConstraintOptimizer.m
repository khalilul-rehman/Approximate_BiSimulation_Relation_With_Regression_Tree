function [h, m, mo] = QuadraticConstraintOptimizer( x, y)


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

    % Step 1: Normalize the data to avoid scaling issues
    x_mean = mean(x);
    x_std = std(x);
    x_norm = (x - x_mean) ./ x_std;  % Normalize x

    y_mean = mean(y);
    y_std = std(y);
    y_norm = (y - y_mean) ./ y_std;  % Normalize y

    % Get dimensions
    [N, d] = size(x_norm);  % Number of samples (N) and input features (d)
    p = size(y_norm, 2);    % Output dimensions (p)

    % Initial guess for [h, m (p × d), mo (p × 1)]
    x0 = [1; zeros(p * d, 1); zeros(p, 1)]; 

    % Objective function: minimize h
    objective = @(vars) vars(1); % h is the first variable

    % Regularization parameter for stability
    lambda = 1e-4;

    % Nonlinear constraint function
    nonlinear_constraint = @(vars) nonlinear_constraints(vars, x_norm, y_norm, p, d, lambda);

    % Solve using fmincon with multiple algorithms
    %'Display', 'iter'
    options = optimoptions('fmincon', 'Display', 'off', 'Algorithm', 'sqp');  % Try 'sqp' as a stable algorithm
    %in below line replased fval with ~
    [solution, ~] = fmincon(objective, x0, [], [], [], [], [], [], nonlinear_constraint, options);

    % Extract optimized values
    h = solution(1); % h is the first element
    m = reshape(solution(2:p*d+1), p, d); % Reshape to p × d
    mo = solution(p*d+2:end); % Last p elements as vector

    % Step 2: Convert the solution back to the original scale
    m = bsxfun(@rdivide, m, x_std');        % Adjust weights based on x scaling
    mo = ((mo .* y_std') + y_mean' - (m * x_mean')); % Ensure m0 remains 4 × 1
    %mo = mo .* y_std + y_mean - m * x_mean'; % Adjust offsets based on scaling and mean shifts
end


function [c, ceq] = nonlinear_constraints( vars, x, y, p, d, lambda)
    % Computes the nonlinear constraints for optimization.
    %
    % Inputs:
    %   - vars: Optimization variables (h, m, mo)
    %   - x: Input data (N × d)
    %   - y: Target data (N × p)
    %   - p: Output dimension
    %   - d: Input feature dimension
    %   - lambda: Regularization parameter
    %
    % Outputs:
    %   - c: Inequality constraints (each (y - m*x - mo)^2 ≤ h + regularization)
    %   - ceq: Empty (no equality constraints)

    h = vars(1); % Extract h
    m = reshape(vars(2:p*d+1), p, d); % Extract weight matrix (p × d)
    mo = vars(p*d+2:end); % Extract offset vector (p × 1)

    % Compute constraint for each (x_i, y_i)
    N = size(x, 1);
    c = zeros(N * p, 1); % Constraints for all samples & outputs

    for i = 1:N
        for j = 1:p
            %disp(["constraint number" , int2str(i)]);
            % Quadratic constraint: (y_ij - m_j * x_i - mo_j)^2 ≤ h + lambda * norm(m, 'fro')^2
            c((i-1)*p + j) = (y(i, j) - m(j, :) * x(i, :)' - mo(j))^2 - h + lambda * (norm(m, 'fro')^2 + norm(mo)^2);
        end
    end

    ceq = []; % No equality constraints
end


%{
function [h, m, m0] = QuadraticConstraintOptimizer( x, y)
    % QuadraticConstraintOptimizer with corrected m0 dimension (p × 1).
    %
    % Inputs:
    %   - x: Input data matrix (size: N × d)
    %   - y: Target data matrix (size: N × p)
    %
    % Outputs:
    %   - h: Optimal objective value
    %   - m: Optimal weight matrix (size: p × d)
    %   - m0: Optimal offset vector (size: p × 1)

    % Get dimensions
    [N, d] = size(x);  % Number of samples (N) and input features (d)
    p = size(y, 2);    % Output dimensions (p)

    % Normalize the input and target data (optional, helps with stability)
    x_mean = mean(x);
    x_std = std(x);
    x = (x - x_mean) ./ x_std;

    y_mean = mean(y);
    y_std = std(y);
    y = (y - y_mean) ./ y_std;

    % Initial guess for [h, m (p × d), m0 (p × 1)]
    x0 = [1; zeros(p * d, 1); zeros(p, 1)];

    % Objective function: minimize h
    objective = @(vars) vars(1); % h is the first variable

    % Nonlinear constraint function
    nonlinear_constraint = @(vars) optimized_constraints(vars, x, y, p, d);

    % Solver options (use lower memory options if available)
    options = optimoptions('fmincon', ...
                           'Display', 'off', ...
                           'Algorithm', 'interior-point', ...
                           'FiniteDifferenceType', 'forward', ... % Reduce memory for gradients
                           'HessianApproximation', 'lbfgs');     % Limited memory for Hessian

    % Solve the optimization problem
    [solution, ~] = fmincon(objective, x0, [], [], [], [], [], [], nonlinear_constraint, options);

    % Extract optimized values
    h = solution(1);                        % h is the first element
    m = reshape(solution(2:p*d+1), p, d);   % Reshape to p × d
    m0 = solution(p*d+2:end);               % Extract m0 as p × 1

    % Convert the solution back to original scale
    m = bsxfun(@rdivide, m, x_std');        % Adjust weights for scaling
    m0 = ((m0 .* y_std') + y_mean' - (m * x_mean')); % Ensure m0 remains 4 × 1

    % m0 = (m0 .* y_std + y_mean' - m * x_mean'); % Ensure consistent dimensions

    %m0 = (m0 .* y_std' + y_mean - m * x_mean'); % Ensure consistent dimensions

    %m0 = (m0 .* y_std + y_mean' - m * x_mean'); % Adjust offsets for scaling
end

function [c, ceq] = optimized_constraints( vars, x, y, p, d)
    % Optimized nonlinear constraints for low memory usage.
    %
    % Inputs:
    %   - vars: Optimization variables (h, m, m0)
    %   - x: Input data (N × d)
    %   - y: Target data (N × p)
    %   - p: Output dimension
    %   - d: Input feature dimension
    %
    % Outputs:
    %   - c: Inequality constraints (each (y - m*x - m0)^2 ≤ h)
    %   - ceq: Empty (no equality constraints)

    h = vars(1);                          % Extract h
    m = reshape(vars(2:p*d+1), p, d);     % Extract weight matrix (p × d)
    m0 = vars(p*d+2:end);                 % Extract offset vector (p × 1)

    % Constraints for all samples
    N = size(x, 1);
    c = zeros(N * p, 1);                  % Preallocate memory for constraints

    % Use a single loop to compute constraints in a memory-efficient way
    for i = 1:N
        % Compute predicted values for all outputs at sample i
        prediction = m * x(i, :)' + m0;   % m0 is already p × 1, no transpose needed
        residuals = y(i, :)' - prediction;
        c((i-1)*p+1:i*p) = residuals.^2 - h; % Quadratic inequality constraint
    end

    ceq = []; % No equality constraints
end

%}