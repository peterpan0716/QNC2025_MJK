% LATER Fitting Exercise
%
% Copyright 2023 by Joshua I. Gold, University of Pennsylvania

% The basic idea in fitting a model to data is to find the parameters of
% the model that provide in some sense the best match of the model to the
% data. This match is provided by the "objective function." 
% This exercise is intended to demystify this process by getting you
% to define the initial conditions and objective function for 
% fitting the LATER model to RT data. For a much more thorough, but still 
% very accessible, overview of model fitting (to behavioral data), here
% is a great place to start:
%
% https://elifesciences.org/articles/49547
%
% For this exercise, recall that the point of the LATER model is that 1/RT is
% distributed as a Gaussian, where we can define the parameters
% of the Gaussian (mu and sigma) with respect to the standard parameters
% of the LATER model (muR and deltaS):
%       mu = muR/deltaS
%       sigma = 1/deltaS
%
% So fitting LATER to behavioral data involves finding parameters
% muR and deltaS that provide the best match to the data, according to 
% the appropriate objective function.
%
% Follow along the steps below, some of which will require you to complete
% the code (and therefore hopefully think about how to relate the high-
% level concepts discussed above with the nitty-gritty part of getting 
% everything to actually work.

%%  1. Get the data
%   
%   Use this code to get a data set (array of RTs from a single condition) 
%   to fit, already preprocessed to include correct trials only and remove
%   outliers (including express saccades). See later_getData for details
data = later_getData([], [], 0.2);
RTs = data{1};
clear data

%%  2. Define the objective function
%
% The objective function typically defines the error that you want to 
% minimize between your data and the model predictions. A common objective 
% function is the negative of the sum of the log-likelihoods of the data, 
% given the model parameters. To unpack that for the LATER model:
%
%   1. For each data point (RT from a single trial, in this case) and given
%       set of model parameters, compute the probability of the data, given
%       the model (i.e., the likelihood)
%   2. Take the logarithm
%   3. Sum all these log-likelihoods from all the data points
%   4. Take the negative, because we want to find the minimum (thus
%        corresponding to the maximum likelihood)
%
%   You can define the function simply using an "anonymous function"
%   (https://www.mathworks.com/help/matlab/matlab_prog/anonymous-functions.html), 
%   using this template that assumes that "fits" is a 2x1 vector of
%   [muR, deltaS]:
 
% EXERCISE:

invRT  = 1./RTs;

laterErrFcn = @(fits) -sum( -0.5*log(2*pi*(1/fits(2))^2) - ((invRT - (fits(1)/fits(2))).^2) ./ (2*(1/fits(2))^2) );

%%  3. Define initial conditions
%   
%   For the actual fitting, we will use fmincon
%   (https://www.mathworks.com/help/optim/ug/fmincon.html), which is 
%   "function minimization with constraints." This function allows for 
%   constraints that include upper and lower bounds on the parameters.
%   So here we define those bounds, along with the initial values.
%   We'll use fairly arbitrary values for the lower and upper
%   bounds, but we should pick the initial values more judiciously. HINT: 
%   Recall that the muR and deltaS should be strongly related to 
%   empirical summary statistics of `the (reciprocal) RT distribution.
lowerBounds = [0.001 0.001];
upperBounds = [1000 1000]; 

% EXERCISE:5
m = mean(invRT);
s = std(invRT);
deltaS0 = max(1e-3, 1/s);       % respect lower bound
muR0    = max(1e-3, m/s);

initialValues = [muR0, deltaS0];   % [muR, deltaS]

%%  4. Run the fits
% 
%   We will be using GlobalSearch . The general advantage of this approach 
%   is to avoid local minima; for details, see:
%   https://www.mathworks.com/help/gads/how-globalsearch-and-multistart-work.html
%  
%   These options seem to work well, but I don't have a stronger
%   rationale for using them. See the Matlab documentation if you really
%   want to dive in and understand them, and let me know if you find
%   better settings!
opts = optimoptions(@fmincon, 'Algorithm','active-set', ...
    'MaxIter',3000, 'MaxFunEvals',3000);

problem = createOptimProblem('fmincon', ...
    'objective',   laterErrFcn, ...
    'x0',          initialValues, ...
    'lb',          lowerBounds, ...
    'ub',          upperBounds, ...
    'options',     opts);

gs = GlobalSearch;

[bestFits, nllk] = run(gs, problem);   % bestFits = [muR_hat, deltaS_hat]
muR_hat   = bestFits(1);
deltaS_hat= bestFits(2);
mu_hat    = muR_hat/deltaS_hat;
sigma_hat = 1/deltaS_hat;

%%  5. Evaluate the fits
% 5a) Compare fitted vs. empirical moments of invRT
% If the fitted 𝜇 and  σ closely match the empirical mean and SD of 1 / RT 1/RT, the model captures 
% the data’s central tendency and spread; large gaps flag a poor fit or preprocessing issues.

fprintf('invRT mean (empirical) = %.4f,  invRT std (empirical) = %.4f\n', mean(invRT), std(invRT));
fprintf('invRT mean (fitted)    = %.4f,  invRT std (fitted)    = %.4f\n', mu_hat, sigma_hat);

% 5b) Histogram of invRT with fitted Gaussian overlay
% If the fitted Normal curve tracks the histogram peak and tails, the Normal-in- 1/RT assumption is
% reasonable; systematic mismatches imply skew/heavy tails or mixtures.
figure; hold on
histogram(invRT, 'Normalization','pdf'); 
xg = linspace(min(invRT), max(invRT), 400);
plot(xg, normpdf(xg, mu_hat, sigma_hat), 'LineWidth', 2);
xlabel('1 / RT'); ylabel('Density'); title('LATER fit in 1/RT space'); box on

% 5c) QQ-plot of invRT vs. fitted Gaussian
% Points lying near the 45° line indicate 1/RT is well-approximated by the fitted Normal; 
% curved/S-shaped patterns signal deviations (skew, kurtosis) or outliers.
figure;
qqplot(invRT, makedist('Normal','mu',mu_hat,'sigma',sigma_hat));
title('QQ plot: 1/RT vs fitted Gaussian')

% 5d) Optional GOF: Kolmogorov–Smirnov test in 1/RT space
% A non-significant result (high p) means you can’t reject that 1/RT 
% comes from the fitted Normal distribution.

[h,p] = kstest(invRT, 'CDF', makedist('Normal','mu',mu_hat,'sigma',sigma_hat));
fprintf('KS test vs fitted Gaussian: h=%d, p=%.4g (higher p => better fit)\n', h, p);
