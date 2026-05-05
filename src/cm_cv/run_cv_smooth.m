%% run_cv_smooth.m
% =========================================================================
% CV SMOOTH-BREAK ESTIMATION — ALL THREE WINDOWS
%
% Fits a logistic smooth-transition mean-shift model to each CV series:
%
%   s_t = a + b * F(t; gamma, tau) + u_t
%
% where s_t = sqrt(pi/2)*|e_hat_t| is the Sensier-style conditional
% volatility proxy from the CM break-model residuals, F is the logistic CDF,
% and (gamma, tau) are chosen by minimising SSR over a grid search.
%
% ESTIMATOR: smooth_break_search_mean (no change to logic)
% GAMMA GRID: GAMMA_GRID from config.m  ([1 2 5 10 20 50 100])
% TRIMMING:   TRIM_FRAC from config.m   (0.15)
%
% PREREQUISITES:
%   run_cm_historical, run_cm_2007, run_cm_fullsample
%   (need cm_supw_resid_break_* saved in cm_*_results.mat)
%
% OUTPUTS (data_intermediate/)
%   cv_smooth_hist_results.mat  → series_meta_cv_smooth_hist
%   cv_smooth_2007_results.mat  → series_meta_cv_smooth_2007
%   cv_smooth_full_results.mat  → series_meta_cv_smooth_full
%
% OUTPUTS (output/tables/)
%   cv_smooth_hist_results.csv
%   cv_smooth_2007_results.csv
%   cv_smooth_full_results.csv
%
% KEY COLUMNS (per window, naming matches termpaper610.m and presentation layer)
%   cv_smooth_found_*    - logical: smooth-break search succeeded
%   cv_smooth_gamma_*    - best-fit logistic speed parameter
%   cv_smooth_midpoint_* - midpoint date (date string, dd-MMM-yyyy)
%   cv_smooth_ssr_gain_* - SSR reduction (%) vs. constant mean
%   cv_smooth_ratio_*    - post/pre mean ratio at midpoint
%   cv_smooth_pct_change_* - percent change in mean at midpoint
%   cv_smooth_abrupt_*   - logical: gamma_best >= 20
% =========================================================================

run('config.m');

if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end
if ~exist(DATA_INT,'dir');   mkdir(DATA_INT);   end

load(fullfile(DATA_INT,'processed_data.mat'), 'series');
n_series = length(series);

%% -----------------------------------------------------------------------
% SECTION 1 — HISTORICAL WINDOW (1959:03–1999:12)
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_smooth] Historical window...\n');

load(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'cm_supw_resid_break_hist','series_meta_cm_hist');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','include_main_hist');

cm_p_hist = series_meta_cm_hist.cm_p_hist;

cv_smooth_found_hist      = false(n_series,1);
cv_smooth_gamma_hist      = NaN(n_series,1);
cv_smooth_midpoint_hist   = strings(n_series,1);
cv_smooth_ssr_gain_hist   = NaN(n_series,1);
cv_smooth_ratio_hist      = NaN(n_series,1);
cv_smooth_pct_change_hist = NaN(n_series,1);
cv_smooth_abrupt_hist     = false(n_series,1);

for j = 1:n_series
    if ~(include_main_hist(j) && series_meta_cm_hist.cm_supw_found_hist(j) && ...
            ~isempty(cm_supw_resid_break_hist{j}))
        continue;
    end

    pj          = cm_p_hist(j);
    dates_reg_j = dates_rep_t((pj+1):end);
    ehat_j      = cm_supw_resid_break_hist{j};
    s_j         = sqrt(pi/2) * abs(ehat_j);

    Tj = length(s_j);
    if length(dates_reg_j) > Tj
        dates_reg_j = dates_reg_j(end-Tj+1:end);
    end

    res = smooth_break_search_mean(s_j, dates_reg_j, TRIM_FRAC, GAMMA_GRID);
    if res.success
        cv_smooth_found_hist(j)      = true;
        cv_smooth_gamma_hist(j)      = res.gamma_best;
        cv_smooth_midpoint_hist(j)   = string(res.midpoint_date,'dd-MMM-yyyy');
        cv_smooth_ssr_gain_hist(j)   = res.ssr_gain_pct;
        cv_smooth_ratio_hist(j)      = res.ratio_post_pre;
        cv_smooth_pct_change_hist(j) = res.pct_change;
        cv_smooth_abrupt_hist(j)     = res.gamma_best >= 20;
    end
end

series_meta_cv_smooth_hist = table( ...
    string(series(:)), include_main_hist(:), ...
    cv_smooth_found_hist, cv_smooth_gamma_hist, cv_smooth_midpoint_hist, ...
    cv_smooth_ssr_gain_hist, cv_smooth_ratio_hist, ...
    cv_smooth_pct_change_hist, cv_smooth_abrupt_hist, ...
    'VariableNames', { ...
        'series_name','include_main_hist', ...
        'cv_smooth_found_hist','cv_smooth_gamma_hist','cv_smooth_midpoint_hist', ...
        'cv_smooth_ssr_gain_hist','cv_smooth_ratio_hist', ...
        'cv_smooth_pct_change_hist','cv_smooth_abrupt_hist'});

save(fullfile(DATA_INT,'cv_smooth_hist_results.mat'), 'series_meta_cv_smooth_hist');
writetable(series_meta_cv_smooth_hist, fullfile(OUT_TABLES,'cv_smooth_hist_results.csv'));

fprintf('[run_cv_smooth] Hist: %d/%d series estimated, %d abrupt (gamma>=20)\n', ...
    sum(cv_smooth_found_hist), sum(include_main_hist), ...
    sum(cv_smooth_abrupt_hist));

%% -----------------------------------------------------------------------
% SECTION 2 — 2007 WINDOW (1959:03–2007:12)
%% -----------------------------------------------------------------------
fprintf('[run_cv_smooth] 2007 window...\n');

load(fullfile(DATA_INT,'cm_2007_results.mat'), ...
    'cm_supw_resid_break_2007','series_meta_cm_2007');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_2007','include_main_2007');

cm2007_p = series_meta_cm_2007.cm2007_p;

cv2007_smooth_found      = false(n_series,1);
cv2007_smooth_gamma      = NaN(n_series,1);
cv2007_smooth_midpoint   = strings(n_series,1);
cv2007_smooth_ssr_gain   = NaN(n_series,1);
cv2007_smooth_ratio      = NaN(n_series,1);
cv2007_smooth_pct_change = NaN(n_series,1);
cv2007_smooth_abrupt     = false(n_series,1);

for j = 1:n_series
    if ~(include_main_2007(j) && series_meta_cm_2007.cm2007_supw_found(j) && ...
            ~isempty(cm_supw_resid_break_2007{j}))
        continue;
    end

    pj          = cm2007_p(j);
    dates_reg_j = dates_ext_2007((pj+1):end);
    ehat_j      = cm_supw_resid_break_2007{j};
    s_j         = sqrt(pi/2) * abs(ehat_j);

    Tj = length(s_j);
    if length(dates_reg_j) > Tj
        dates_reg_j = dates_reg_j(end-Tj+1:end);
    end

    res = smooth_break_search_mean(s_j, dates_reg_j, TRIM_FRAC, GAMMA_GRID);
    if res.success
        cv2007_smooth_found(j)      = true;
        cv2007_smooth_gamma(j)      = res.gamma_best;
        cv2007_smooth_midpoint(j)   = string(res.midpoint_date,'dd-MMM-yyyy');
        cv2007_smooth_ssr_gain(j)   = res.ssr_gain_pct;
        cv2007_smooth_ratio(j)      = res.ratio_post_pre;
        cv2007_smooth_pct_change(j) = res.pct_change;
        cv2007_smooth_abrupt(j)     = res.gamma_best >= 20;
    end
end

series_meta_cv_smooth_2007 = table( ...
    string(series(:)), include_main_2007(:), ...
    cv2007_smooth_found, cv2007_smooth_gamma, cv2007_smooth_midpoint, ...
    cv2007_smooth_ssr_gain, cv2007_smooth_ratio, ...
    cv2007_smooth_pct_change, cv2007_smooth_abrupt, ...
    'VariableNames', { ...
        'series_name','include_main_2007', ...
        'cv2007_smooth_found','cv2007_smooth_gamma','cv2007_smooth_midpoint', ...
        'cv2007_smooth_ssr_gain','cv2007_smooth_ratio', ...
        'cv2007_smooth_pct_change','cv2007_smooth_abrupt'});

save(fullfile(DATA_INT,'cv_smooth_2007_results.mat'), 'series_meta_cv_smooth_2007');
writetable(series_meta_cv_smooth_2007, fullfile(OUT_TABLES,'cv_smooth_2007_results.csv'));

fprintf('[run_cv_smooth] 2007: %d/%d series estimated, %d abrupt (gamma>=20)\n', ...
    sum(cv2007_smooth_found), sum(include_main_2007), ...
    sum(cv2007_smooth_abrupt));

%% -----------------------------------------------------------------------
% SECTION 3 — FULL SAMPLE (1959:03–latest)
%% -----------------------------------------------------------------------
fprintf('[run_cv_smooth] Full sample...\n');

load(fullfile(DATA_INT,'cm_full_results.mat'), ...
    'cm_supw_resid_break_full','series_meta_cm_full');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_full','include_main_full_ext');

cmfull_p = series_meta_cm_full.cmfull_p;

cvfull_smooth_found      = false(n_series,1);
cvfull_smooth_gamma      = NaN(n_series,1);
cvfull_smooth_midpoint   = strings(n_series,1);
cvfull_smooth_ssr_gain   = NaN(n_series,1);
cvfull_smooth_ratio      = NaN(n_series,1);
cvfull_smooth_pct_change = NaN(n_series,1);
cvfull_smooth_abrupt     = false(n_series,1);

for j = 1:n_series
    if ~(include_main_full_ext(j) && series_meta_cm_full.cmfull_supw_found(j) && ...
            ~isempty(cm_supw_resid_break_full{j}))
        continue;
    end

    pj          = cmfull_p(j);
    dates_reg_j = dates_ext_full((pj+1):end);
    ehat_j      = cm_supw_resid_break_full{j};
    s_j         = sqrt(pi/2) * abs(ehat_j);

    Tj = length(s_j);
    if length(dates_reg_j) > Tj
        dates_reg_j = dates_reg_j(end-Tj+1:end);
    end

    res = smooth_break_search_mean(s_j, dates_reg_j, TRIM_FRAC, GAMMA_GRID);
    if res.success
        cvfull_smooth_found(j)      = true;
        cvfull_smooth_gamma(j)      = res.gamma_best;
        cvfull_smooth_midpoint(j)   = string(res.midpoint_date,'dd-MMM-yyyy');
        cvfull_smooth_ssr_gain(j)   = res.ssr_gain_pct;
        cvfull_smooth_ratio(j)      = res.ratio_post_pre;
        cvfull_smooth_pct_change(j) = res.pct_change;
        cvfull_smooth_abrupt(j)     = res.gamma_best >= 20;
    end
end

series_meta_cv_smooth_full = table( ...
    string(series(:)), include_main_full_ext(:), ...
    cvfull_smooth_found, cvfull_smooth_gamma, cvfull_smooth_midpoint, ...
    cvfull_smooth_ssr_gain, cvfull_smooth_ratio, ...
    cvfull_smooth_pct_change, cvfull_smooth_abrupt, ...
    'VariableNames', { ...
        'series_name','include_main_full_ext', ...
        'cvfull_smooth_found','cvfull_smooth_gamma','cvfull_smooth_midpoint', ...
        'cvfull_smooth_ssr_gain','cvfull_smooth_ratio', ...
        'cvfull_smooth_pct_change','cvfull_smooth_abrupt'});

save(fullfile(DATA_INT,'cv_smooth_full_results.mat'), 'series_meta_cv_smooth_full');
writetable(series_meta_cv_smooth_full, fullfile(OUT_TABLES,'cv_smooth_full_results.csv'));

fprintf('[run_cv_smooth] Full: %d/%d series estimated, %d abrupt (gamma>=20)\n', ...
    sum(cvfull_smooth_found), sum(include_main_full_ext), ...
    sum(cvfull_smooth_abrupt));

fprintf('[run_cv_smooth] Done.\n\n');
