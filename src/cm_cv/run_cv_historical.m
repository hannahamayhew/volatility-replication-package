%% run_cv_historical.m
% =========================================================================
% CONDITIONAL VOLATILITY BREAKS — HISTORICAL WINDOW (1959:03 – 1999:12)
%
% Uses residuals from the AR(p) conditional-mean SupW model to construct
% the Sensier-style conditional volatility proxy:
%
%   s_t = sqrt(pi/2) * |e_hat_t|
%
% Then tests for a break in the mean of s_t using the same one-break
% SupW + Hansen workflow used for UV.
%
% cv_supw_test() internally applies the sqrt(pi/2) scaling to |e_hat|.
% This script passes the raw residuals from the break-model fit.
%
% PREREQUISITE: run_cm_historical.m must have run first (provides residuals).
%
% INPUTS  (loaded from data_intermediate)
%   cm_hist_results.mat   → cm_supw_resid_break_hist, cm_p_hist, series_meta_cm_hist
%   processed_data.mat    → dates_rep_t, series, include_main_hist
%
% OUTPUTS (saved to data_intermediate/cv_hist_results.mat)
%   series_meta_cv_hist  - table with CV break results
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_historical] Loading CM residuals and processed data...\n');
load(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'cm_supw_resid_break_hist','series_meta_cm_hist');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','series','include_main_hist');

n_series = length(series);
cm_p_hist = series_meta_cm_hist.cm_p_hist;

%% -----------------------------------------------------------------------
% 2. CV ONE-BREAK SUP-W TEST
%
%   For each series with a valid CM fit, pass the CM break-model residuals
%   to cv_supw_test(), which internally constructs s_t = sqrt(pi/2)*|e_hat|.
%
%   HAC lag rule: h = floor(4*(T/100)^(2/9))
%% -----------------------------------------------------------------------
fprintf('[run_cv_historical] Running CV SupW tests (trim=%.0f%%)...\n', TRIM_FRAC*100);

cv_supw_found_hist      = false(n_series,1);
cv_supw_stat_hist       = NaN(n_series,1);
cv_supw_break_date_hist = strings(n_series,1);
cv_supw_hac_lag_hist    = NaN(n_series,1);
cv_ratio_hist           = NaN(n_series,1);
cv_pct_change_hist      = NaN(n_series,1);
cv_mu1_hist             = NaN(n_series,1);
cv_mu2_hist             = NaN(n_series,1);
cv_break_date_hist      = strings(n_series,1);

for j = 1:n_series
    if include_main_hist(j) && series_meta_cm_hist.cm_supw_found_hist(j) && ...
            ~isempty(cm_supw_resid_break_hist{j})

        pj          = cm_p_hist(j);
        dates_reg_j = dates_rep_t((pj+1):end);
        ehat_j      = cm_supw_resid_break_hist{j};

        % Align date vector if residual length differs from expectation
        Tj = length(ehat_j);
        if length(dates_reg_j) > Tj
            dates_reg_j = dates_reg_j(end-Tj+1:end);
        end

        % cv_supw_test constructs s_t = sqrt(pi/2)*|ehat| internally
        res = cv_supw_test(ehat_j, dates_reg_j, TRIM_FRAC, []);

        if res.success
            cv_supw_found_hist(j)      = true;
            cv_supw_stat_hist(j)       = res.supW;
            cv_supw_break_date_hist(j) = string(res.break_date,'dd-MMM-yyyy');
            cv_supw_hac_lag_hist(j)    = res.hac_lag;
            cv_mu1_hist(j)             = res.mu1;
            cv_mu2_hist(j)             = res.mu2;
            cv_ratio_hist(j)           = res.ratio_post_pre;
            cv_pct_change_hist(j)      = res.pct_change;
            cv_break_date_hist(j)      = string(res.break_date,'dd-MMM-yyyy');
        end
    end
end

%% -----------------------------------------------------------------------
% 3. HANSEN (1997) P-VALUES  (m=1, pi0=TRIM_FRAC)
%% -----------------------------------------------------------------------
cv_hansen_pval_hist = NaN(n_series,1);
cv_hansen_sig5_hist = false(n_series,1);

for j = 1:n_series
    if cv_supw_found_hist(j)
        cv_hansen_pval_hist(j) = hansen_supf_pval(cv_supw_stat_hist(j), 1, TRIM_FRAC);
        cv_hansen_sig5_hist(j) = cv_hansen_pval_hist(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cv_hist = table( ...
    string(series(:)), include_main_hist(:), ...
    cv_supw_found_hist, cv_supw_stat_hist, cv_supw_break_date_hist, ...
    cv_supw_hac_lag_hist, ...
    cv_mu1_hist, cv_mu2_hist, cv_ratio_hist, cv_pct_change_hist, ...
    cv_break_date_hist, ...
    cv_hansen_pval_hist, cv_hansen_sig5_hist, ...
    'VariableNames',{ ...
        'series_name','include_main_hist', ...
        'cv_supw_found_hist','cv_supw_stat_hist','cv_supw_break_date_hist', ...
        'cv_supw_hac_lag_hist', ...
        'cv_mu1_hist','cv_mu2_hist','cv_ratio_hist','cv_pct_change_hist', ...
        'cv_break_date_hist', ...
        'cv_hansen_pval_hist','cv_hansen_sig5_hist'});

%% -----------------------------------------------------------------------
% 5. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'cv_hist_results.mat'), 'series_meta_cv_hist');
writetable(series_meta_cv_hist, fullfile(OUT_TABLES,'cv_hist_results.csv'));

fprintf('\n[run_cv_historical] SUMMARY\n');
fprintf('  Series in main sample with CM fit: %d\n', ...
    sum(include_main_hist(:) & series_meta_cm_hist.cm_supw_found_hist(:)));
fprintf('  Significant CV breaks (p<5%%):    %d\n', sum(cv_hansen_sig5_hist));
fprintf('[run_cv_historical] Done.\n\n');
