%% QNC HQ 09/05/2025: Probability Distributions I: Concepts and examples 
% Min Jae Kim (minjae.kim@pennmedicine.upenn.edu)
% For my assignment, I chose this paper "The statistical reliability of
% signals in single neurons in cat and monkey visual cortex" that spike counts 
% follows poisson distributions based on λT (expected spike count) in cat and
% monkey.
% URL: https://pubmed.ncbi.nlm.nih.gov/6623937/
%% Initializing Variables
% this code tests different expected spike count (λT) at a given time interval
% (firing rate (λ) x time interval T).
% For each expected spike count tested, it will give:
% 1. counts = values of spike counts across 100 iterations
% 2. mean_hat= mean spike counts in distribution
% 3. var_har = variance of spike counts in distribution
% 4. Fano Factor = var_hat / mean_hat
% variance = variance of disibution

expected_spike_count = [2 4 6 9 12 16 22 30]; 
n_trials=100;
[counts, mean_hat, var_hat, fano] = simulate_poisson_spike_counts(expected_spike_count, n_trials);

%
%% Plotting Variables
% this plot shows how Poisson distribution changes 
% as we test different expected spike count: 
figure('Color','w'); hold on;
for k = 1:numel(expected_spike_count)
    mu = expected_spike_count(k); c = counts{k};
    maxk = max(ceil(mu + 4*sqrt(mu)), max(c));
    ks = 0:maxk;
    plot(ks, poisspdf(ks, mu), 'LineWidth', 1.5);
end
xlabel('Spike count'); ylabel('Probability');
title('Poisson PDFs for different \lambda T'); legend(string(expected_spike_count), 'Location','northeastoutside');
%% Function
function [counts, mean_hat, var_hat, fano] = simulate_poisson_spike_counts(means, ntrials);
% means : vector of expected spike counts per trial for each stimulus level

K = numel(means);
counts  = cell(K,1);
mean_hat = zeros(K,1);
var_hat  = zeros(K,1);

for k = 1:K
    c = poissrnd(means(k), [ntrials, 1]);
    counts{k}  = c;
    mean_hat(k)= mean(c);
    var_hat(k) = var(c, 0);  % sample variance (N-1)
end

fano = var_hat ./ mean_hat;   % ~1 if Poisson-like
end
