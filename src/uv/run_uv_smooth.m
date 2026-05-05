%% run_uv_smooth.m
% =========================================================================
% UV SMOOTH-BREAK ESTIMATION — ALL THREE WINDOWS
%
% Fits a logistic smooth-transition mean-shift model to each UV series:
%
%   z_t = a + b * F(t; gamma, tau) + u_t
%
% where F is the logistic CDF and (gamma, tau) are chosen by minimising SSR
% over a grid search.  A series is classified as "abrupt" if gamma_best >= 20
% (pragmatic rule from Sensier & van Dijk 2004 / termpaper610.m).
%
% ESTIMATOR: smooth_break_search_mean (no change to logic)
% GAMMA GRID: GAMMA_GRID from config.m  ([1 2 5 10 20 50 100])
% TRIMMING:   TRIM_FRAC from config.m   (0.15)
%
% PREREQUISITES: run_uv_historical, run_uv_2007, run_uv_fullsample
%   (need vol_obj* saved in uv_*_results.mat)
%
% OUTPUTS (data_intermediate/)
%   uv_smooth_hist_results.mat  → series_meta_uv_smooth_hist
%   uv_smooth_2007_results.mat  → series_meta_uv_smooth_2007
%   uv_smooth_full_results.mat  → series_meta_uv_smooth_full
%
% OUTPUTS (output/tables/)
%   uv_smooth_hist_results.csv
%   uv_smooth_2007_results.csv
%   uv_smooth_full_results.csv
%
% KEY COLUMNS (per window, naming matches termpaper610.m and presentation layer)
%   uv_smooth_found_*   - logical: smooth-break search succeeded
%   uv_smooth_gamma_*   - best-fit logistic speed parameter
%   uv_smooth_midpoint_* - midpoint date (date string, dd-MMM-yyyy)
%   uv_smooth_ssr_gain_* - SSR reduction (%) vs. constant mean
%   uv_smooth_ratio_*   - post/pre mean ratio at midpoint
%   uv_smooth_pct_change_* - percent change in mean at midpoint
%   uv_smooth_abrupt_*  - logical: gamma_best >= 20
% =========================================================================

run('config.m');

if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end
if ~exist(DATA_INT,'dir');   mkdir(DATA_INT);   end

%% -----------------------------------------------------------------------
% SECTION 1 — HISTORICAL WINDOW (1959:03–1999:12)
%% -----------------------------------------------------------------------
fprintf('\n[run_uv_smooth] Historical window...\n');

load(fullfile(DATA_INT,'uv_hist_results.mat'), 'vol_obj');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','series','include_main_hist');

n_series = length(series);

uv_smooth_found_hist      = false(n_series,1);
uv_smooth_gamma_hist      = NaN(n_series,1);
uv_smooth_midpoint_hist   = strings(n_series,1);
uv_smooth_ssr_gain_hist   = NaN(n_series,1);
uv_smooth_ratio_hist      = NaN(n_series,1);
uv_smooth_pct_change_hist = NaN(n_series,1);
uv_smooth_abrupt_hist     = false(n_series,1);

for j = 1:n_series
    if ~include_main_hist(j); continue; end
    vj  = vol_obj(:,j);
    res = smooth_break_search_mean(vj, dates_rep_t, TRIM_FRAC, GAMMA_GRID);
    if res.success
        uv_smooth_found_hist(j)      = true;
        uv_smooth_gamma_hist(j)      = res.gamma_best;
        uv_smooth_midpoint_hist(j)   = string(res.midpoint_date,'dd-MMM-yyyy');
        uv_smooth_ssr_gain_hist(j)   = res.ssr_gain_pct;
        uv_smooth_ratio_hist(j)      = res.ratio_post_pre;
        uv_smooth_pct_change_hist(j) = res.pct_change;
        uv_smooth_abrupt_hist(j)     = res.gamma_best >= 20;
    end
end

series_meta_uv_smooth_hist = table( ...
    string(series(:)), include_main_hist(:), ...
    uv_smooth_found_hist, uv_smooth_gamma_hist, uv_smooth_midpoint_hist, ...
    uv_smooth_ssr_gain_hist, uv_smooth_ratio_hist, ...
    uv_smooth_pct_change_hist, uv_smooth_abrupt_hist, ...
    'VariableNames', { ...
        'series_name','include_main_hist', ...
        'uv_smooth_found_hist','uv_smooth_gamma_hist','uv_smooth_midpoint_hist', ...
        'uv_smooth_ssr_gain_hist','uv_smooth_ratio_hist', ...
        'uv_smooth_pct_change_hist','uv_smooth_abrupt_hist'});

save(fullfile(DATA_INT,'uv_smooth_hist_results.mat'), 'series_meta_uv_smooth_hist');
writetable(series_meta_uv_smooth_hist, fullfile(OUT_TABLES,'uv_smooth_hist_results.csv'));

fprintf('[run_uv_smooth] Hist: %d/%d series estimated, %d abrupt (gamma>=20)\n', ...
    sum(uv_smooth_found_hist), sum(include_main_hist), ...
    sum(uv_smooth_abrupt_hist));

%% -----------------------------------------------------------------------
% SECTION 2 — 2007 WINDOW (1959:03–2007:12)
%% -----------------------------------------------------------------------
fprintf('[run_uv_smooth] 2007 window...\n');

load(fullfile(DATA_INT,'uv_2007_results.mat'), 'vol_obj_2007');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_2007','include_main_2007');

uv2007_smooth_found      = false(n_series,1);
uv2007_smooth_gamma      = NaN(n_series,1);
uv2007_smooth_midpoint   = strings(n_series,1);
uv2007_smooth_ssr_gain   = NaN(n_series,1);
uv2007_smooth_ratio      = NaN(n_series,1);
uv2007_smooth_pct_change = NaN(n_series,1);
uv2007_smooth_abrupt     = false(n_series,1);

for j = 1:n_series
    if ~include_main_2007(j); continue; end
    vj  = vol_obj_2007(:,j);
    res = smooth_break_search_mean(vj, dates_ext_2007, TRIM_FRAC, GAMMA_GRID);
    if res.success
        uv2007_smooth_found(j)      = true;
        uv2007_smooth_gamma(j)      = res.gamma_best;
        uv2007_smooth_midpoint(j)   = string(res.midpoint_date,'dd-MMM-yyyy');
        uv2007_smooth_ssr_gain(j)   = res.ssr_gain_pct;
        uv2007_smooth_ratio(j)      = res.ratio_post_pre;
        uv2007_smooth_pct_change(j) = res.pct_change;
        uv2007_smooth_abrupt(j)     = res.gamma_best >= 20;
    end
end

series_meta_uv_smooth_2007 = table( ...
    string(series(:)), include_main_2007(:), ...
    uv2007_smooth_found, uv2007_smooth_gamma, uv2007_smooth_midpoint, ...
    uv2007_smooth_ssr_gain, uv2007_smooth_ratio, ...
    uv2007_smooth_pct_change, uv2007_smooth_abrupt, ...
    'VariableNames', { ...
        'series_name','include_main_2007', ...
        'uv2007_smooth_found','uv2007_smooth_gamma','uv2007_smooth_midpoint', ...
        'uv2007_smooth_ssr_gain','uv2007_smooth_ratio', ...
        'uv2007_smooth_pct_change','uv2007_smooth_abrupt'});

save(fullfile(DATA_INT,'uv_smooth_2007_results.mat'), 'series_meta_uv_smooth_2007');
writetable(series_meta_uv_smooth_2007, fullfile(OUT_TABLES,'uv_smooth_2007_results.csv'));

fprintf('[run_uv_smooth] 2007: %d/%d series estimated, %d abrupt (gamma>=20)\n', ...
    sum(uv2007_smooth_found), sum(include_main_2007), ...
    sum(uv2007_smooth_abrupt));

%% -----------------------------------------------------------------------
% SECTION 3 — FULL SAMPLE (1959:03–latest)
%% -----------------------------------------------------------------------
fprintf('[run_uv_smooth] Full sample...\n');

load(fullfile(DATA_INT,'uv_full_results.mat'), 'vol_obj_full');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_full','include_main_full_ext');

uvfull_smooth_found      = false(n_series,1);
uvfull_smooth_gamma      = NaN(n_series,1);
uvfull_smooth_midpoint   = strings(n_series,1);
uvfull_smooth_ssr_gain   = NaN(n_series,1);
uvfull_smooth_ratio      = NaN(n_series,1);
uvfull_smooth_pct_change = NaN(n_series,1);
uvfull_smooth_abrupt     = false(n_series,1);

for j = 1:n_series
    if ~include_main_full_ext(j); continue; end
    vj  = vol_obj_full(:,j);
    res = smooth_break_search_mean(vj, dates_ext_full, TRIM_FRAC, GAMMA_GRID);
    if res.success
        uvfull_smooth_found(j)      = true;
        uvfull_smooth_gamma(j)      = res.gamma_best;
        uvfull_smooth_midpoint(j)   = string(res.midpoint_date,'dd-MMM-yyyy');
        uvfull_smooth_ssr_gain(j)   = res.ssr_gain_pct;
        uvfull_smooth_ratio(j)      = res.ratio_post_pre;
        uvfull_smooth_pct_change(j) = res.pct_change;
        uvfull_smooth_abrupt(j)     = res.gamma_best >= 20;
    end
end

series_meta_uv_smooth_full = table( ...
    string(series(:)), include_main_full_ext(:), ...
    uvfull_smooth_found, uvfull_smooth_gamma, uvfull_smooth_midpoint, ...
    uvfull_smooth_ssr_gain, uvfull_smooth_ratio, ...
    uvfull_smooth_pct_change, uvfull_smooth_abrupt, ...
    'VariableNames', { ...
        'series_name','include_main_full_ext', ...
        'uvfull_smooth_found','uvfull_smooth_gamma','uvfull_smooth_midpoint', ...
        'uvfull_smooth_ssr_gain','uvfull_smooth_ratio', ...
        'uvfull_smooth_pct_change','uvfull_smooth_abrupt'});

save(fullfile(DATA_INT,'uv_smooth_full_results.mat'), 'series_meta_uv_smooth_full');
writetable(series_meta_uv_smooth_full, fullfile(OUT_TABLES,'uv_smooth_full_results.csv'));

fprintf('[run_uv_smooth] Full: %d/%d series estimated, %d abrupt (gamma>=20)\n', ...
    sum(uvfull_smooth_found), sum(include_main_full_ext), ...
    sum(uvfull_smooth_abrupt));

fprintf('[run_uv_smooth] Done.\n\n');
