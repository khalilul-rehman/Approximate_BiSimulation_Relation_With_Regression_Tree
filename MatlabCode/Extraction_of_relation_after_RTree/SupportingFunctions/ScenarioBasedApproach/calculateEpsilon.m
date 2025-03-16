%{
function epsilons = calculateEpsilon(beta_target, d, N_array)
    % Calculate epsilon for each N in the array using bisection search
    tol = 1e-10; % Tolerance for bisection search
    epsilons = zeros(size(N_array));
    for idx = 1:length(N_array)
        N = N_array(idx);
        % Bisection search for epsilon in [0, 1]
        epsilons(idx) = bisectionSearch(@(epsilon) calculateBeta(epsilon, d, N) - beta_target, 0, 1, tol);
    end
end

function epsilon = bisectionSearch(fun, a, b, tol)
    % Bisection search to find the root of fun(epsilon) = 0 in [a, b]
    max_iter = 1000; % Maximum number of iterations
    iter = 0;
    
    while (b - a) / 2 > tol && iter < max_iter
        epsilon = (a + b) / 2; % Midpoint
        f_mid = fun(epsilon);
        
        if f_mid == 0
            % Found exact solution
            break;
        elseif fun(a) * f_mid < 0
            % Root is in [a, epsilon]
            b = epsilon;
        else
            % Root is in [epsilon, b]
            a = epsilon;
        end
        
        iter = iter + 1;
    end
    
    epsilon = (a + b) / 2; % Return the midpoint as the solution
end
%}


function epsilons = calculateEpsilon(beta_target, d, N_array)
    epsilons = zeros(size(N_array));
    for idx = 1:length(N_array)
        N = N_array(idx);
        % Define the function to solve for epsilon
        fun = @(epsilon) calculateBeta(epsilon, d, N) - beta_target;
        % Initial guess for epsilon
        epsilon_guess = 0.1;
        % Solve for epsilon
        epsilons(idx) = fzero(fun, epsilon_guess);
    end
end


%{
function epsilons = calculateEpsilon(beta_target, d, N_array)
    epsilons = zeros(size(N_array));
    for idx = 1:length(N_array)
        N = N_array(idx);
        % Define the function to solve for epsilon
        fun = @(epsilon) log(abs(calculateBeta(epsilon, d, N) - beta_target) + 1e-15); % Log scale to handle flat regions
        % Better initial guess based on problem structure
        epsilon_guess = min(0.5, 1 / sqrt(N));
        try
            % Use a reduced search range and better solver settings
            epsilons(idx) = fminbnd(fun, 1e-10, 1 - 1e-5, optimset('TolX', 1e-12));
            % Ensure the result is valid
            if epsilons(idx) >= 1 - 1e-5 || epsilons(idx) <= 1e-10
                epsilons(idx) = NaN; % Mark as invalid if it’s too close to the boundary
            end
        catch
            epsilons(idx) = NaN; % Catch any solver errors and set to NaN
        end
    end
end

%}