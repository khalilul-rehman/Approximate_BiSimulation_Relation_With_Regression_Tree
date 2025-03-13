function epsilon = calculate_epsilon(lj, n, beta)
    epsilon = zeros(length(lj), 1);
    
    options = optimoptions('lsqnonlin', 'Display', 'off');
    
    for j = 1:length(lj)
        L = abs(lj(j));
        
        % Define the equation using arrayfun for size consistency
        equation = @(eps) sum(arrayfun(@(i) ...
            exp(gammaln(L + 1) - gammaln(i + 1) - gammaln(L - i + 1) + ...
                i * log(eps) + (L - i) * log(1 - eps)), 1:n)) - beta;
        
        % Initial guess for epsilon
        initial_guess = 1e-3;
        
        % Solve using lsqnonlin
        epsilon(j) = lsqnonlin(equation, initial_guess, 0, 1, options);
    end
end

%{
function epsilon = calculate_epsilon(lj, n, beta)
    % lj: list of values
    % n: number
    % beta: target probability
    % epsilon: calculated value for each lj

    epsilon = zeros(length(lj), 1);
    
    foptions = optimoptions('fsolve', 'Display', 'off'); % Suppress output
    
    for j = 1:length(lj)
        L = abs(lj(j));
        
        % Define the equation to solve for epsilon
        equation = @(eps) sum(arrayfun(@(i) nchoosek(L, i) * (eps^i) * ((1 - eps)^(L - i)), 0:(n + 1))) - beta;
        
        % Initial guess for epsilon (can adjust if needed)
        initial_guess = 0.5;
        
        % Solve using fsolve
        epsilon(j) = fsolve(equation, initial_guess, foptions);
    end
end
%}