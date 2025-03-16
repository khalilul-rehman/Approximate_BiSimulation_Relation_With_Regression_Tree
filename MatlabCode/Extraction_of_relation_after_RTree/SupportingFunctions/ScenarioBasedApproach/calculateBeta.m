
function betas = calculateBeta(epsilon, d, N_array)
    betas = zeros(size(N_array));
    for idx = 1:length(N_array)
        N = N_array(idx);
        beta = 0;
        for i = 0:(d-1)
            beta = beta + nchoosek(N, i) * (epsilon^i) * ((1 - epsilon)^(N - i));
        end
        betas(idx) = beta;
    end
end


%{
function beta = calculateBeta(epsilon, d, N)
    % Compute beta using log-based binomial coefficient to avoid overflow
    beta = 0;
    for i = 0:(d-1)
        % log(nchoosek) = gammaln(N + 1) - gammaln(i + 1) - gammaln(N - i + 1)
        log_binom = gammaln(N + 1) - gammaln(i + 1) - gammaln(N - i + 1);
        beta = beta + exp(log_binom + i * log(epsilon) + (N - i) * log(1 - epsilon));
    end
end
%}