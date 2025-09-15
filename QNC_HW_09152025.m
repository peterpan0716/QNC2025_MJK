%% Post-hoc power for session-wise Spearman r (Poisson LC vs Gaussian pupil)
% Min Jae Kim (minjae.kim@pennmedicine.upenn.edu)
%% Part 0: Initializing Variables

clear; clc; rng(42);

% --- User parameters ---
T_per_session   = 600;                      % timepoints per session
M_null_sessions = 3000;                     % sessions to estimate null SD of r
alpha           = 0.05;                     % two-sided
target_power    = 0.80;
effect_sizes    = [0.05 0.10 0.15 0.20 0.25 0.30];   % mean session-wise r

% LC ~ Poisson, pupil ~ Normal (independent under H0)
lc_lambda   = 5.0;
pupil_mu    = 3.0;
pupil_sigma = 0.5;

%% Part 1 — Estimate null SD of session-wise Spearman r
r_null = nan(M_null_sessions,1);
for i = 1:M_null_sessions
    lc    = poissrnd(lc_lambda,  T_per_session, 1);
    pupil = normrnd(pupil_mu, pupil_sigma, T_per_session, 1);
    r_null(i) = corr(lc, pupil, 'Type','Spearman', 'Rows','complete');
end
null_mean = mean(r_null,'omitnan');
null_sd   = std(r_null,  'omitnan');

fprintf('Null Spearman r (independence): mean=%.4f, SD=%.4f (T=%d, M=%d)\n',...
    null_mean, null_sd, T_per_session, M_null_sessions);

%% Part 2 — Minimal n for 80%% power using exact noncentral t (no helper function)
% One-sample two-sided t-test of mean r against 0 across sessions.
n_max = 2000;
n_req = nan(size(effect_sizes));
pwr_at_n = nan(size(effect_sizes));

for k = 1:numel(effect_sizes)
    eff = effect_sizes(k);
    for n = 2:n_max                                         % need at least 2 for a t-test
        df    = n - 1;
        tcrit = tinv(1 - alpha/2, df);
        ncp   = sqrt(n) * eff / null_sd;                    % noncentrality parameter
        % Two-sided power = P(T < -tcrit) + P(T > tcrit) for T ~ nct(df, ncp)
        p_lo  = nctcdf(-tcrit, df, ncp);
        p_hi  = 1 - nctcdf( tcrit, df, ncp);
        pwr   = p_lo + p_hi;
        if pwr >= target_power
            n_req(k)   = n;
            pwr_at_n(k)= pwr;
            break
        end
    end
end

T = table(effect_sizes(:), n_req(:), pwr_at_n(:), ...
    'VariableNames', {'effect_size_r','n_sessions_80pct','exact_power_at_n'});
disp(T);

%% Part 3 — Plot n (sessions) vs effect size
figure('Color','w'); 
plot(effect_sizes, n_req, '-o','LineWidth',1.5); grid on;
xlabel('Effect size (mean session-wise Spearman r)');
ylabel('n sessions for 80% power (\alpha=0.05, two-sided)');
title({'Sample size vs. effect size', ...
       sprintf('Null SD of r estimated from Poisson–Gaussian independence (SD=%.4f)', null_sd)});

%% Part 4 — Saving Outcome
writetable(T, 'n_by_effectsize_exact_ttest_matlab.csv');
