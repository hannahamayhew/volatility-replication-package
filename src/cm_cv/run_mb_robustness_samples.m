%% run_mb_robustness_samples.m
% =========================================================================
% ROBUSTNESS — RESTRICTED-SAMPLE CHECKS
%
% Check 1 — Balanced panel
%   Requires ≥99% non-missing observations AND first valid date ≤ 1959:03.
%   These flags (include_bal_hist, include_bal_2007) are computed in
%   run_data_prep.m using BALANCED_THRESH=0.99 and BALANCED_START_CUTOFF.
%   Expected sample size: ~116 series.
%
% Check 2 — Exclude McCracken-Ng (2016) flagged series
%   Eight series flagged for explicit adjustments or splices:
%     Bilateral exchange rates : EXSZUSx, EXJPUSx, EXUSUKx, EXCAUSx
%     Broad exchange rate index: TWEXAFEGSMTHx
%     Commercial paper         : CP3Mx
%     Oil price                : OILPRICEx
%     Consumer sentiment       : UMCSENTx
%   Source: McCracken & Ng (2016), Appendix.
%
% Both checks use identical parameters to the main analysis:
%   Historical window : 15%*T rule (naturally = 74 months for T≈490)
%   2007 window       : MIN_SEG_EXT = 74 months (fixed)
%   Full window       : MIN_SEG_EXT = 74 months (fixed)
%   MAX_BREAKS        = 5
%
% PREREQUISITES
%   run_data_prep.m, run_uv_historical.m, run_uv_2007.m, run_uv_fullsample.m,
%   run_cm_historical.m (and optionally run_cm_2007.m, run_cm_fullsample.m)
%
% OUTPUTS (saved to output/tables/)
%   mb_robust_balanced_results.csv     — per-series, balanced subset
%   mb_robust_noflag_results.csv       — per-series, no-flag subset
%   mb_robust_samples_summary.csv      — cross-sample break-distribution summary
% =========================================================================

run('config.m');

fprintf('\n[run_mb_robustness_samples] Starting restricted-sample robustness checks...\n');

% -------------------------------------------------------------------------
% Flagged series (McCracken & Ng 2016 — explicit adjustments / splices)
% -------------------------------------------------------------------------
FLAGGED_SERIES = {'EXSZUSx','EXJPUSx','EXUSUKx','EXCAUSx', ...
                  'TWEXAFEGSMTHx','CP3Mx','OILPRICEx','UMCSENTx'};

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('[run_mb_robustness_samples] Loading inputs...\n');

load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','dates_ext_2007','dates_ext_full', ...
    'series','fred_group', ...
    'include_main_hist','include_bal_hist', ...
    'include_main_2007','include_bal_2007', ...
    'include_main_full_ext');
load(fullfile(DATA_INT,'uv_hist_results.mat'),  'vol_obj');
load(fullfile(DATA_INT,'uv_2007_results.mat'),  'vol_obj_2007');
load(fullfile(DATA_INT,'uv_full_results.mat'),  'vol_obj_full');
load(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'cm_supw_resid_break_hist','series_meta_cm_hist');

n_series      = length(series);
cm_p_hist     = series_meta_cm_hist.cm_p_hist;
cm_supw_found = series_meta_cm_hist.cm_supw_found_hist;

% -------------------------------------------------------------------------
% Build no-flag mask (exclude 8 McCracken-Ng flagged series from main mask)
% -------------------------------------------------------------------------
series_str    = string(series(:));
flag_mask     = ismember(series_str, string(FLAGGED_SERIES));
n_flagged_found = sum(flag_mask);
fprintf('[run_mb_robustness_samples] Flagged series found in dataset: %d / %d requested.\n', ...
    n_flagged_found, numel(FLAGGED_SERIES));
if n_flagged_found < numel(FLAGGED_SERIES)
    missing = FLAGGED_SERIES(~ismember(string(FLAGGED_SERIES), series_str));
    fprintf('  WARNING: Not found: %s\n', strjoin(missing,', '));
end

include_noflag_hist     = include_main_hist(:)     & ~flag_mask;
include_noflag_2007     = include_main_2007(:)     & ~flag_mask;
include_noflag_full_ext = include_main_full_ext(:) & ~flag_mask;

% Report sample sizes
fprintf('[run_mb_robustness_samples] Historical window sample sizes:\n');
fprintf('  Main:     %d series\n', sum(include_main_hist));
fprintf('  Balanced: %d series\n', sum(include_bal_hist));
fprintf('  No-flag:  %d series\n', sum(include_noflag_hist));

fprintf('[run_mb_robustness_samples] 2007 window sample sizes:\n');
fprintf('  Main:     %d series\n', sum(include_main_2007));
fprintf('  Balanced: %d series\n', sum(include_bal_2007));
fprintf('  No-flag:  %d series\n', sum(include_noflag_2007));

fprintf('[run_mb_robustness_samples] Full window sample sizes:\n');
fprintf('  Main:     %d series\n', sum(include_main_full_ext));
fprintf('  No-flag:  %d series\n', sum(include_noflag_full_ext));

%% -----------------------------------------------------------------------
% 2. DEFINE SAMPLES AND RUN LOOP
%% -----------------------------------------------------------------------
% For each sample: [label, include_hist, include_2007, include_full]
samples = { ...
    'balanced', include_bal_hist(:),      include_bal_2007(:),      include_main_full_ext(:) & ~flag_mask(:); ...
    'noflag',   include_noflag_hist(:),   include_noflag_2007(:),   include_noflag_full_ext(:)};

% Load optional CM for 2007 and full windows
cm2007_path  = fullfile(DATA_INT,'cm_2007_results.mat');
cm_full_path = fullfile(DATA_INT,'cm_full_results.mat');

if exist(cm2007_path,'file')
    load(cm2007_path,'cm_supw_resid_break_2007','series_meta_cm_2007');
    cm_p_2007    = series_meta_cm_2007.cm2007_p;
    cm_found_2007 = series_meta_cm_2007.cm2007_supw_found;
    have_cm2007  = true;
else
    have_cm2007  = false;
end

if exist(cm_full_path,'file')
    load(cm_full_path,'cm_supw_resid_break_full','series_meta_cm_full');
    cm_p_full    = series_meta_cm_full.cmfull_p;
    cm_found_full = series_meta_cm_full.cmfull_supw_found;
    have_cm_full = true;
else
    have_cm_full = false;
end

% Storage for summary rows
sum_lbl  = strings(0,1); sum_win = strings(0,1); sum_meas = strings(0,1);
sum_n    = []; sum_0b = []; sum_1b = []; sum_2p = []; sum_3p = [];
sum_med  = []; sum_max = [];

all_results = struct();

for si = 1:size(samples,1)
    samp_lbl   = samples{si,1};
    inc_hist_s = samples{si,2};
    inc_2007_s = samples{si,3};
    inc_full_s = samples{si,4};

    fprintf('\n[run_mb_robustness_samples] === Sample: %s ===\n', upper(samp_lbl));

    % --- Preallocate ---
    uv_nb_hist = NaN(n_series,1); uv_dt_hist = strings(n_series,1);
    cv_nb_hist = NaN(n_series,1); cv_dt_hist = strings(n_series,1);
    uv_nb_2007 = NaN(n_series,1); uv_dt_2007 = strings(n_series,1);
    cv_nb_2007 = NaN(n_series,1); cv_dt_2007 = strings(n_series,1);
    uv_nb_full = NaN(n_series,1); uv_dt_full = strings(n_series,1);
    cv_nb_full = NaN(n_series,1); cv_dt_full = strings(n_series,1);

    % --- UV historical ---
    fprintf('[run_mb_robustness_samples] UV historical (%s)...\n', samp_lbl);
    for j = 1:n_series
        if ~inc_hist_s(j); continue; end
        res = uv_recursive_breaks(vol_obj(:,j), dates_rep_t, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN);
        uv_nb_hist(j) = res.n_breaks;
        if res.n_breaks > 0
            uv_dt_hist(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'),' | ');
        end
    end
    tmp = uv_nb_hist(inc_hist_s);
    fprintf('  0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f (n=%d)\n', ...
        mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
        mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
        median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(~isnan(tmp)));
    [sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max] = ...
        append_row(sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max, ...
        samp_lbl,'hist','UV',tmp);

    % --- CV historical ---
    fprintf('[run_mb_robustness_samples] CV historical (%s)...\n', samp_lbl);
    for j = 1:n_series
        if ~(inc_hist_s(j) && cm_supw_found(j) && ~isempty(cm_supw_resid_break_hist{j})); continue; end
        pj = cm_p_hist(j); ehat_j = cm_supw_resid_break_hist{j};
        s_j = sqrt(pi/2)*abs(ehat_j); Tj = length(s_j);
        dj  = dates_rep_t((pj+1):end); dj = dj(end-Tj+1:end);
        res = cv_recursive_breaks(s_j, dj, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN);
        cv_nb_hist(j) = res.n_breaks;
        if res.n_breaks > 0
            cv_dt_hist(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'),' | ');
        end
    end
    tmp = cv_nb_hist(inc_hist_s);
    fprintf('  0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f (n=%d)\n', ...
        mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
        mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
        median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(~isnan(tmp)));
    [sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max] = ...
        append_row(sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max, ...
        samp_lbl,'hist','CV',tmp);

    % --- UV 2007 ---
    fprintf('[run_mb_robustness_samples] UV 2007 (%s)...\n', samp_lbl);
    for j = 1:n_series
        if ~inc_2007_s(j); continue; end
        res = uv_recursive_breaks(vol_obj_2007(:,j), dates_ext_2007, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN, MIN_SEG_EXT);
        uv_nb_2007(j) = res.n_breaks;
        if res.n_breaks > 0
            uv_dt_2007(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'),' | ');
        end
    end
    tmp = uv_nb_2007(inc_2007_s);
    fprintf('  0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f (n=%d)\n', ...
        mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
        mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
        median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(~isnan(tmp)));
    [sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max] = ...
        append_row(sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max, ...
        samp_lbl,'2007','UV',tmp);

    % --- CV 2007 ---
    if have_cm2007
        fprintf('[run_mb_robustness_samples] CV 2007 (%s)...\n', samp_lbl);
        for j = 1:n_series
            if ~(inc_2007_s(j) && cm_found_2007(j) && ~isempty(cm_supw_resid_break_2007{j})); continue; end
            pj = cm_p_2007(j); ehat_j = cm_supw_resid_break_2007{j};
            s_j = sqrt(pi/2)*abs(ehat_j); Tj = length(s_j);
            dj  = dates_ext_2007((pj+1):end); dj = dj(end-Tj+1:end);
            res = cv_recursive_breaks(s_j, dj, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN, MIN_SEG_EXT);
            cv_nb_2007(j) = res.n_breaks;
            if res.n_breaks > 0
                cv_dt_2007(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'),' | ');
            end
        end
        tmp = cv_nb_2007(inc_2007_s);
        fprintf('  0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f (n=%d)\n', ...
            mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
            mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
            median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(~isnan(tmp)));
        [sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max] = ...
            append_row(sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max, ...
            samp_lbl,'2007','CV',tmp);
    end

    % --- UV full ---
    fprintf('[run_mb_robustness_samples] UV full (%s)...\n', samp_lbl);
    for j = 1:n_series
        if ~inc_full_s(j); continue; end
        res = uv_recursive_breaks(vol_obj_full(:,j), dates_ext_full, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN, MIN_SEG_EXT);
        uv_nb_full(j) = res.n_breaks;
        if res.n_breaks > 0
            uv_dt_full(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'),' | ');
        end
    end
    tmp = uv_nb_full(inc_full_s);
    fprintf('  0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f (n=%d)\n', ...
        mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
        mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
        median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(~isnan(tmp)));
    [sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max] = ...
        append_row(sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max, ...
        samp_lbl,'full','UV',tmp);

    % --- CV full ---
    if have_cm_full
        fprintf('[run_mb_robustness_samples] CV full (%s)...\n', samp_lbl);
        for j = 1:n_series
            if ~(inc_full_s(j) && cm_found_full(j) && ~isempty(cm_supw_resid_break_full{j})); continue; end
            pj = cm_p_full(j); ehat_j = cm_supw_resid_break_full{j};
            s_j = sqrt(pi/2)*abs(ehat_j); Tj = length(s_j);
            dj  = dates_ext_full((pj+1):end);
            if length(dj) > Tj; dj = dj(end-Tj+1:end); end
            res = cv_recursive_breaks(s_j, dj, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN, MIN_SEG_EXT);
            cv_nb_full(j) = res.n_breaks;
            if res.n_breaks > 0
                cv_dt_full(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'),' | ');
            end
        end
        tmp = cv_nb_full(inc_full_s);
        fprintf('  0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f (n=%d)\n', ...
            mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
            mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
            median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(~isnan(tmp)));
        [sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max] = ...
            append_row(sum_lbl,sum_win,sum_meas,sum_n,sum_0b,sum_1b,sum_2p,sum_3p,sum_med,sum_max, ...
            samp_lbl,'full','CV',tmp);
    end

    % --- Save per-series table ---
    out_tbl = table( ...
        series_str, string(fred_group(:)), ...
        inc_hist_s, inc_2007_s, inc_full_s, ...
        uv_nb_hist, uv_dt_hist, cv_nb_hist, cv_dt_hist, ...
        uv_nb_2007, uv_dt_2007, cv_nb_2007, cv_dt_2007, ...
        uv_nb_full, uv_dt_full, cv_nb_full, cv_dt_full, ...
        'VariableNames',{ ...
            'series_name','fred_group', ...
            'include_hist','include_2007','include_full', ...
            'uv_nb_hist','uv_break_dates_hist','cv_nb_hist','cv_break_dates_hist', ...
            'uv_nb_2007','uv_break_dates_2007','cv_nb_2007','cv_break_dates_2007', ...
            'uv_nb_full','uv_break_dates_full','cv_nb_full','cv_break_dates_full'});

    fname = sprintf('mb_robust_%s_results.csv', samp_lbl);
    writetable(out_tbl, fullfile(OUT_TABLES, fname));
    fprintf('[run_mb_robustness_samples] Saved %s\n', fname);
    all_results.(samp_lbl) = out_tbl;
end

%% -----------------------------------------------------------------------
% 3. CROSS-SAMPLE SUMMARY TABLE
%% -----------------------------------------------------------------------
summary_tbl = table(sum_lbl, sum_win, sum_meas, sum_n, ...
    sum_0b, sum_1b, sum_2p, sum_3p, sum_med, sum_max, ...
    'VariableNames',{ ...
        'sample','window','measure','n_series', ...
        'share_0b','share_1b','share_2plus','share_3plus', ...
        'median_breaks','max_breaks'});

writetable(summary_tbl, fullfile(OUT_TABLES,'mb_robust_samples_summary.csv'));
fprintf('[run_mb_robustness_samples] Saved mb_robust_samples_summary.csv\n');
fprintf('[run_mb_robustness_samples] Done.\n\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [sl,sw,sm,sn,s0,s1,s2,s3,smed,smax] = append_row( ...
        sl,sw,sm,sn,s0,s1,s2,s3,smed,smax, lbl,win,meas,tmp)
    tmp = tmp(~isnan(tmp));
    sl(end+1,1)   = lbl;
    sw(end+1,1)   = win;
    sm(end+1,1)   = meas;
    sn(end+1,1)   = numel(tmp);
    s0(end+1,1)   = mean(tmp==0);
    s1(end+1,1)   = mean(tmp==1);
    s2(end+1,1)   = mean(tmp>=2);
    s3(end+1,1)   = mean(tmp>=3);
    smed(end+1,1) = median(tmp);
    smax(end+1,1) = max(tmp);
end
