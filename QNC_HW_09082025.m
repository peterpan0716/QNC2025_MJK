%% QNC HQ 09/08/2025: Probability Distributions I: Concepts and examples
% Min Jae Kim (minjae.kim@pennmedicine.upenn.edu)
% For my assignment, I chose this paper:
% "The variable discharge of cortical neurons: implications for connectivity, computation, and information coding."
% Shadlen & Newsome (1998) — which describes interspike-interval (ISI) distributions in cortex as approximately exponential.
% URL: https://pubmed.ncbi.nlm.nih.gov/9570816/

% Here, μ_ISI is the mean interspike interval; firing rate λ = 1/μ_ISI.

%% Initializing Variables
% This code tests different mean ISIs (μ_ISI, in seconds).
% For each μ_ISI (mean interspike interval), the output is 
% 1) isis     = simulated ISIs (N samples)
% 2) mean_hat = empirical mean of ISIs
% 3) var_hat  = empirical variance of ISIs
% 4) cv       = coefficient of variation = std(ISI)/mean(ISI)  (≈ 1 for exponential)

mean_isi = [0.05 0.1 0.2 0.3 0.5];    % seconds  (≈ 20, 10, 5, 3.3, 2 Hz)
n_samples = 5000;                     % draws per condition

[isis, mean_hat, var_hat, cv] = simulate_exponential_isi(mean_isi, n_samples);

%% Plotting Variables
% This plot shows how the Exponential PDF changes as we test different mean ISIs.
figure('Color','w'); hold on;
tmax = max(mean_isi) * 6;          % cover most of the mass for all μ
t = linspace(0, tmax, 1000);
for k = 1:numel(mean_isi)
    mu = mean_isi(k);
    % Exponential pdf with mean mu: f(t) = (1/mu) * exp(-t/mu)
    plot(t, (1./mu).*exp(-t./mu), 'LineWidth', 2);
end
xlabel('Interspike interval (s)'); ylabel('Probability Density');
title('Exponential PDFs for different mean ISIs (\mu_{ISI} = 1/\lambda)');
legend(string(mean_isi), 'Location','northeastoutside');
hold off;

%% Function
function [isis, mean_hat, var_hat, cv] = simulate_exponential_isi(mean_isi, n_samples, rng_seed)
% mean_isi : vector of mean ISIs (μ, seconds) for each condition
% n_samples: number of ISI draws per condition
% rng_seed : random seed for reproducibility (optional)

if nargin < 3, rng_seed = 0; end
rng(rng_seed);

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
