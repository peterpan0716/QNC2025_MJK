function [isis, mean_hat, var_hat, cv] = simulate_exponential_isi(mean_isi, n_samples)
% mean_isi : vector of mean ISIs (μ, seconds) for each condition
% n_samples: number of ISI draws per condition
% rng_seed : random seed for reproducibility (optional)

K = numel(mean_isi);
isis     = cell(K,1);
mean_hat = zeros(K,1);
var_hat  = zeros(K,1);

for k = 1:K
    mu = mean_isi(k);
    % Draw from Exponential(μ). Use exprnd if available, else inverse-CDF.
    try
        x = exprnd(mu, [n_samples, 1]);    % Statistics Toolbox
    catch
        x = -mu * log(rand(n_samples,1));  % fallback: inverse-CDF
    end
    isis{k}    = x;
    mean_hat(k)= mean(x);
    var_hat(k) = var(x, 0);                % sample variance (N-1)
end

cv = sqrt(var_hat) ./ mean_hat;            % ≈ 1 for exponential
end
