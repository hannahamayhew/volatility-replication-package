%% run_mb_fixed74.m
% =========================================================================
% MAIN SPECIFICATION — RECURSIVE MULTIPLE-BREAK, FIXED 74-MONTH MIN SEG
%
%   min_seg_len_eff = 74 months (FIXED, uniform across all three windows)
%   max_breaks      = 5  (consistent with main analysis)
%
% This is the main specification for the paper. The 74-month floor is
% defensible because it matches the natural h-rule value for the historical
% window (ceil(0.15 * 490) = 74) and is therefore consistent with
% Sensier & van Dijk (2004) in practice. Fixing it uniformly across all
% windows ensures that minimum-segment geometry is comparable across windows.
%
% RATIONALE
%   The standard h-rule (ceil(0.15 * T)) gives:
%     Historical (T≈490):   74 months  ← natural value for this window
%     Pre-GFC   (T≈586):   88 months
%     Full      (T≈803):  121 months
%   Locking all windows at 74 makes the geometry COMPARABLE across windows
%   and eliminates the mechanical suppression of late-sample breaks that
%   inflated 0-break counts for Labor Market and Output & Income in the
%   full-sample window.
%
% OUTPUTS (new files — nothing existing is overwritten):
%   data_intermediate/mb_fixed74_results.mat
%   output/tables/mb_fixed74_results.csv        (per-series)
%   output/tables/mb_fixed74_group_summary.csv  (sector × window × measure)
%
% PREREQUISITES: same as run_mb_historical.m
%   run_uv_historical, run_uv_2007, run_uv_fullsample
%   run_cm_historical, run_cm_2007 (optional), run_cm_fullsample (optional)
% =========================================================================

run('config.m');

MAX_BREAKS_F74   = 5;   % consistent with main analysis
MIN_SEG_EFF_F74  = 74;  % fixed override — bypasses sample-relative h-rule

fprintf('\n[run_mb_fixed74] Fixed-74 specification: max_breaks=%d, min_seg_eff=%d\n', ...
    MAX_BREAKS_F74, MIN_SEG_EFF_F74);

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('[run_mb_fixed74] Loading inputs...\n');

load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','dates_ext_2007','dates_ext_full', ...
    'series','include_main_hist','include_main_2007','include_main_full_ext', ...
    'fred_group');
load(fullfile(DATA_INT,'uv_hist_results.mat'),  'vol_obj');
load(fullfile(DATA_INT,'uv_2007_results.mat'),  'vol_obj_2007');
load(fullfile(DATA_INT,'uv_full_results.mat'),  'vol_obj_full');
load(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'cm_supw_resid_break_hist','series_meta_cm_hist');

n_series      = length(series);
cm_p_hist     = series_meta_cm_hist.cm_p_hist;
cm_supw_found = series_meta_cm_hist.cm_supw_found_hist;

%% -----------------------------------------------------------------------
% 2. UV HISTORICAL
%% -----------------------------------------------------------------------
fprintf('[run_mb_fixed74] UV historical (max=%d, min_seg_eff=%d)...\n', ...
    MAX_BREAKS_F74, MIN_SEG_EFF_F74);

uv_mb_nbreaks_hist     = NaN(n_series,1);
uv_mb_break_dates_hist = strings(n_series,1);
uv_mb_break_stats_hist = strings(n_series,1);
uv_mb_break_pvals_hist = strings(n_series,1);
uv_hist_rp_iters     = NaN(n_series,1);
uv_hist_rp_converged = false(n_series,1);
uv_hist_rp_maxdelta  = NaN(n_series,1);

for j = 1:n_series
    if ~include_main_hist(j); continue; end
    res = uv_recursive_breaks(vol_obj(:,j), dates_rep_t, TRIM_FRAC, ...
        MAX_BREAKS_F74, MIN_SEG_LEN, MIN_SEG_EFF_F74);
    uv_mb_nbreaks_hist(j) = res.n_breaks;
    if res.n_breaks > 0
        uv_mb_break_dates_hist(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
        uv_mb_break_stats_hist(j) = strjoin(string(round(res.break_stats,3)), ' | ');
        uv_mb_break_pvals_hist(j) = strjoin(string(round(res.break_pvals,6)), ' | ');
    else
        uv_mb_break_dates_hist(j) = "";
        uv_mb_break_stats_hist(j) = "";
        uv_mb_break_pvals_hist(j) = "";
    end
    uv_hist_rp_iters(j)     = res.repartition_iters;
    uv_hist_rp_converged(j) = res.repartition_converged;
    uv_hist_rp_maxdelta(j)  = res.repartition_max_delta;
    if mod(j,20)==0; fprintf('  UV hist: %d / %d\n', j, n_series); end
end

tmp = uv_mb_nbreaks_hist(include_main_hist);
fprintf('  UV hist: 0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f cap5=%d\n', ...
    mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
    mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
    median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(tmp==MAX_BREAKS_F74,'omitnan'));
rp_mask = include_main_hist(:) & ~isnan(uv_hist_rp_iters(:));
fprintf('  Repartition: mean_iters=%.1f max_iters=%.0f non-converged=%d/%d\n', ...
    mean(uv_hist_rp_iters(rp_mask)), max(uv_hist_rp_iters(rp_mask)), ...
    sum(~uv_hist_rp_converged(rp_mask)), sum(rp_mask));

%% -----------------------------------------------------------------------
% 3. UV 2007 (PRE-GFC)
%% -----------------------------------------------------------------------
fprintf('[run_mb_fixed74] UV 2007 (pre-GFC)...\n');

uv2007_mb_nbreaks     = NaN(n_series,1);
uv2007_mb_break_dates = strings(n_series,1);
uv2007_mb_break_stats = strings(n_series,1);
uv2007_mb_break_pvals = strings(n_series,1);

for j = 1:n_series
    if ~include_main_2007(j); continue; end
    res = uv_recursive_breaks(vol_obj_2007(:,j), dates_ext_2007, TRIM_FRAC, ...
        MAX_BREAKS_F74, MIN_SEG_LEN, MIN_SEG_EFF_F74);
    uv2007_mb_nbreaks(j) = res.n_breaks;
    if res.n_breaks > 0
        uv2007_mb_break_dates(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
        uv2007_mb_break_stats(j) = strjoin(string(round(res.break_stats,3)), ' | ');
        uv2007_mb_break_pvals(j) = strjoin(string(round(res.break_pvals,6)), ' | ');
    else
        uv2007_mb_break_dates(j) = "";
        uv2007_mb_break_stats(j) = "";
        uv2007_mb_break_pvals(j) = "";
    end
    if mod(j,20)==0; fprintf('  UV 2007: %d / %d\n', j, n_series); end
end

tmp = uv2007_mb_nbreaks(include_main_2007);
fprintf('  UV 2007: 0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f cap5=%d\n', ...
    mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
    mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
    median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(tmp==MAX_BREAKS_F74,'omitnan'));

%% -----------------------------------------------------------------------
% 4. UV FULL SAMPLE
%% -----------------------------------------------------------------------
fprintf('[run_mb_fixed74] UV full sample...\n');

uvfull_mb_nbreaks     = NaN(n_series,1);
uvfull_mb_break_dates = strings(n_series,1);
uvfull_mb_break_stats = strings(n_series,1);
uvfull_mb_break_pvals = strings(n_series,1);

for j = 1:n_series
    if ~include_main_full_ext(j); continue; end
    res = uv_recursive_breaks(vol_obj_full(:,j), dates_ext_full, TRIM_FRAC, ...
        MAX_BREAKS_F74, MIN_SEG_LEN, MIN_SEG_EFF_F74);
    uvfull_mb_nbreaks(j) = res.n_breaks;
    if res.n_breaks > 0
        uvfull_mb_break_dates(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
        uvfull_mb_break_stats(j) = strjoin(string(round(res.break_stats,3)), ' | ');
        uvfull_mb_break_pvals(j) = strjoin(string(round(res.break_pvals,6)), ' | ');
    else
        uvfull_mb_break_dates(j) = "";
        uvfull_mb_break_stats(j) = "";
        uvfull_mb_break_pvals(j) = "";
    end
    if mod(j,20)==0; fprintf('  UV full: %d / %d\n', j, n_series); end
end

tmp = uvfull_mb_nbreaks(include_main_full_ext);
fprintf('  UV full: 0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f cap5=%d\n', ...
    mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
    mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
    median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(tmp==MAX_BREAKS_F74,'omitnan'));

%% -----------------------------------------------------------------------
% 5. CV HISTORICAL
%% -----------------------------------------------------------------------
fprintf('[run_mb_fixed74] CV historical...\n');

cv_mb_nbreaks_hist     = NaN(n_series,1);
cv_mb_break_dates_hist = strings(n_series,1);
cv_mb_break_stats_hist = strings(n_series,1);
cv_mb_break_pvals_hist = strings(n_series,1);
cv_hist_rp_iters     = NaN(n_series,1);
cv_hist_rp_converged = false(n_series,1);
cv_hist_rp_maxdelta  = NaN(n_series,1);

for j = 1:n_series
    if ~(include_main_hist(j) && cm_supw_found(j) && ~isempty(cm_supw_resid_break_hist{j}))
        continue;
    end
    pj          = cm_p_hist(j);
    ehat_j      = cm_supw_resid_break_hist{j};
    s_j         = sqrt(pi/2) * abs(ehat_j);
    Tj          = length(s_j);
    dates_reg_j = dates_rep_t((pj+1):end);
    dates_reg_j = dates_reg_j(end-Tj+1:end);

    res = cv_recursive_breaks(s_j, dates_reg_j, TRIM_FRAC, ...
        MAX_BREAKS_F74, MIN_SEG_LEN, MIN_SEG_EFF_F74);
    cv_mb_nbreaks_hist(j) = res.n_breaks;
    if res.n_breaks > 0
        cv_mb_break_dates_hist(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
        cv_mb_break_stats_hist(j) = strjoin(string(round(res.break_stats,3)), ' | ');
        cv_mb_break_pvals_hist(j) = strjoin(string(round(res.break_pvals,6)), ' | ');
    else
        cv_mb_break_dates_hist(j) = "";
        cv_mb_break_stats_hist(j) = "";
        cv_mb_break_pvals_hist(j) = "";
    end
    cv_hist_rp_iters(j)     = res.repartition_iters;
    cv_hist_rp_converged(j) = res.repartition_converged;
    cv_hist_rp_maxdelta(j)  = res.repartition_max_delta;
    if mod(j,20)==0; fprintf('  CV hist: %d / %d\n', j, n_series); end
end

tmp = cv_mb_nbreaks_hist(include_main_hist);
fprintf('  CV hist: 0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f cap5=%d\n', ...
    mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
    mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
    median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(tmp==MAX_BREAKS_F74,'omitnan'));
rp_mask_cv = include_main_hist(:) & ~isnan(cv_hist_rp_iters(:));
fprintf('  Repartition: mean_iters=%.1f max_iters=%.0f non-converged=%d/%d\n', ...
    mean(cv_hist_rp_iters(rp_mask_cv)), max(cv_hist_rp_iters(rp_mask_cv)), ...
    sum(~cv_hist_rp_converged(rp_mask_cv)), sum(rp_mask_cv));

%% -----------------------------------------------------------------------
% 6. CV 2007 (conditional)
%% -----------------------------------------------------------------------
cv2007_mb_nbreaks     = NaN(n_series,1);
cv2007_mb_break_dates = strings(n_series,1);
cv2007_mb_break_stats = strings(n_series,1);
cv2007_mb_break_pvals = strings(n_series,1);

cm2007_path = fullfile(DATA_INT,'cm_2007_results.mat');
if exist(cm2007_path,'file')
    fprintf('[run_mb_fixed74] CV 2007...\n');
    load(cm2007_path, 'cm_supw_resid_break_2007','series_meta_cm_2007');

    for j = 1:n_series
        if ~(include_main_2007(j) && series_meta_cm_2007.cm2007_supw_found(j) && ...
                ~isempty(cm_supw_resid_break_2007{j}))
            continue;
        end
        pj          = series_meta_cm_2007.cm2007_p(j);
        ehat_j      = cm_supw_resid_break_2007{j};
        s_j         = sqrt(pi/2) * abs(ehat_j);
        Tj          = length(s_j);
        dates_reg_j = dates_ext_2007((pj+1):end);
        dates_reg_j = dates_reg_j(end-Tj+1:end);

        res = cv_recursive_breaks(s_j, dates_reg_j, TRIM_FRAC, ...
            MAX_BREAKS_F74, MIN_SEG_LEN, MIN_SEG_EFF_F74);
        cv2007_mb_nbreaks(j) = res.n_breaks;
        if res.n_breaks > 0
            cv2007_mb_break_dates(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
            cv2007_mb_break_stats(j) = strjoin(string(round(res.break_stats,3)), ' | ');
            cv2007_mb_break_pvals(j) = strjoin(string(round(res.break_pvals,6)), ' | ');
        else
            cv2007_mb_break_dates(j) = "";
            cv2007_mb_break_stats(j) = "";
            cv2007_mb_break_pvals(j) = "";
        end
        if mod(j,20)==0; fprintf('  CV 2007: %d / %d\n', j, n_series); end
    end

    tmp = cv2007_mb_nbreaks(include_main_2007);
    fprintf('  CV 2007: 0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f cap5=%d\n', ...
        mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
        mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
        median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(tmp==MAX_BREAKS_F74,'omitnan'));
else
    fprintf('[run_mb_fixed74] cm_2007_results.mat not found — CV 2007 skipped.\n');
end

%% -----------------------------------------------------------------------
% 7. CV FULL SAMPLE (conditional)
%% -----------------------------------------------------------------------
cvfull_mb_nbreaks     = NaN(n_series,1);
cvfull_mb_break_dates = strings(n_series,1);
cvfull_mb_break_stats = strings(n_series,1);
cvfull_mb_break_pvals = strings(n_series,1);

cm_full_path = fullfile(DATA_INT,'cm_full_results.mat');
cmfull_found   = false(n_series,1);
cmfull_p_arr   = NaN(n_series,1);
cm_supw_resid_break_full = [];

if exist(cm_full_path,'file')
    fprintf('[run_mb_fixed74] CV full sample...\n');
    load(cm_full_path, 'cm_supw_resid_break_full','series_meta_cm_full');
    cmfull_found = series_meta_cm_full.cmfull_supw_found;
    cmfull_p_arr = series_meta_cm_full.cmfull_p;

    for j = 1:n_series
        if ~(include_main_full_ext(j) && cmfull_found(j) && ...
                ~isempty(cm_supw_resid_break_full{j}))
            continue;
        end
        pj          = cmfull_p_arr(j);
        ehat_j      = cm_supw_resid_break_full{j};
        s_j         = sqrt(pi/2) * abs(ehat_j);
        Tj          = length(s_j);
        dates_reg_j = dates_ext_full((pj+1):end);
        if length(dates_reg_j) > Tj
            dates_reg_j = dates_reg_j(end-Tj+1:end);
        end

        res = cv_recursive_breaks(s_j, dates_reg_j, TRIM_FRAC, ...
            MAX_BREAKS_F74, MIN_SEG_LEN, MIN_SEG_EFF_F74);
        cvfull_mb_nbreaks(j) = res.n_breaks;
        if res.n_breaks > 0
            cvfull_mb_break_dates(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
            cvfull_mb_break_stats(j) = strjoin(string(round(res.break_stats,3)), ' | ');
            cvfull_mb_break_pvals(j) = strjoin(string(round(res.break_pvals,6)), ' | ');
        else
            cvfull_mb_break_dates(j) = "";
            cvfull_mb_break_stats(j) = "";
            cvfull_mb_break_pvals(j) = "";
        end
        if mod(j,20)==0; fprintf('  CV full: %d / %d\n', j, n_series); end
    end

    tmp = cvfull_mb_nbreaks(include_main_full_ext);
    fprintf('  CV full: 0b=%.3f 1b=%.3f 2+=%.3f 3+=%.3f med=%.1f max=%.0f cap5=%d\n', ...
        mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
        mean(tmp>=2,'omitnan'), mean(tmp>=3,'omitnan'), ...
        median(tmp,'omitnan'), max(tmp,[],'omitnan'), sum(tmp==MAX_BREAKS_F74,'omitnan'));
else
    fprintf('[run_mb_fixed74] cm_full_results.mat not found — CV full skipped.\n');
end

%% -----------------------------------------------------------------------
% 8. ASSEMBLE PER-SERIES RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_mb_f74 = table( ...
    string(series(:)), string(fred_group(:)), ...
    include_main_hist(:), include_main_2007(:), include_main_full_ext(:), ...
    uv_mb_nbreaks_hist,     uv_mb_break_dates_hist, ...
    uv_mb_break_stats_hist, uv_mb_break_pvals_hist, ...
    cv_mb_nbreaks_hist,     cv_mb_break_dates_hist, ...
    cv_mb_break_stats_hist, cv_mb_break_pvals_hist, ...
    uv2007_mb_nbreaks,      uv2007_mb_break_dates, ...
    uv2007_mb_break_stats,  uv2007_mb_break_pvals, ...
    cv2007_mb_nbreaks,      cv2007_mb_break_dates, ...
    cv2007_mb_break_stats,  cv2007_mb_break_pvals, ...
    uvfull_mb_nbreaks,      uvfull_mb_break_dates, ...
    uvfull_mb_break_stats,  uvfull_mb_break_pvals, ...
    cvfull_mb_nbreaks,      cvfull_mb_break_dates, ...
    cvfull_mb_break_stats,  cvfull_mb_break_pvals, ...
    'VariableNames',{ ...
        'series_name','fred_group', ...
        'include_main_hist','include_main_2007','include_main_full_ext', ...
        'uv_mb_nbreaks_hist',     'uv_mb_break_dates_hist', ...
        'uv_mb_break_stats_hist', 'uv_mb_break_pvals_hist', ...
        'cv_mb_nbreaks_hist',     'cv_mb_break_dates_hist', ...
        'cv_mb_break_stats_hist', 'cv_mb_break_pvals_hist', ...
        'uv2007_mb_nbreaks',      'uv2007_mb_break_dates', ...
        'uv2007_mb_break_stats',  'uv2007_mb_break_pvals', ...
        'cv2007_mb_nbreaks',      'cv2007_mb_break_dates', ...
        'cv2007_mb_break_stats',  'cv2007_mb_break_pvals', ...
        'uvfull_mb_nbreaks',      'uvfull_mb_break_dates', ...
        'uvfull_mb_break_stats',  'uvfull_mb_break_pvals', ...
        'cvfull_mb_nbreaks',      'cvfull_mb_break_dates', ...
        'cvfull_mb_break_stats',  'cvfull_mb_break_pvals'});

%% -----------------------------------------------------------------------
% 9. SAVE PER-SERIES OUTPUTS
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'mb_fixed74_results.mat'), 'series_meta_mb_f74', ...
    'uv_mb_nbreaks_hist','cv_mb_nbreaks_hist', ...
    'uv2007_mb_nbreaks','cv2007_mb_nbreaks', ...
    'uvfull_mb_nbreaks','cvfull_mb_nbreaks', ...
    'uv_mb_break_dates_hist','cv_mb_break_dates_hist', ...
    'uv2007_mb_break_dates','cv2007_mb_break_dates', ...
    'uvfull_mb_break_dates','cvfull_mb_break_dates', ...
    'MAX_BREAKS_F74','MIN_SEG_EFF_F74');
writetable(series_meta_mb_f74, fullfile(OUT_TABLES,'mb_fixed74_results.csv'));
fprintf('[run_mb_fixed74] Saved mb_fixed74_results.mat and mb_fixed74_results.csv\n');

%% -----------------------------------------------------------------------
% 10. SECTOR-LEVEL GROUP SUMMARY
%     Mirrors format of grouped_recursive_mb_summary.csv, adding n_at_cap5
%% -----------------------------------------------------------------------
fprintf('[run_mb_fixed74] Building sector group summary...\n');

groups = unique(string(fred_group(:)));
groups = sort(groups);
groups = groups(~ismember(groups, ["Stock Market","Unassigned",""]));

% Window definitions  [nbreaks_vec,  include_mask,  dates_vec,  label_string]
win_defs = { ...
    '1959:03-1999:12', include_main_hist(:), ...
    uv_mb_nbreaks_hist, uv_mb_break_dates_hist, ...
    cv_mb_nbreaks_hist, cv_mb_break_dates_hist, ...
    vol_obj, dates_rep_t, cmfull_found, cmfull_p_arr; ...
    '1959:03-2007:12', include_main_2007(:), ...
    uv2007_mb_nbreaks, uv2007_mb_break_dates, ...
    cv2007_mb_nbreaks, cv2007_mb_break_dates, ...
    vol_obj_2007, dates_ext_2007, cmfull_found, cmfull_p_arr; ...
    '1959:03-2026:01', include_main_full_ext(:), ...
    uvfull_mb_nbreaks, uvfull_mb_break_dates, ...
    cvfull_mb_nbreaks, cvfull_mb_break_dates, ...
    vol_obj_full, dates_ext_full, cmfull_found, cmfull_p_arr};

% Pre-load CM residuals vectors needed for CV sigma ratio
% (already loaded above for full; hist is in workspace; 2007 may be loaded)
cv_resid_by_win = {cm_supw_resid_break_hist};
if exist('cm_supw_resid_break_2007','var')
    cv_resid_by_win{2} = cm_supw_resid_break_2007;
else
    cv_resid_by_win{2} = {};
end
if exist('cm_supw_resid_break_full','var') && ~isempty(cm_supw_resid_break_full)
    cv_resid_by_win{3} = cm_supw_resid_break_full;
else
    cv_resid_by_win{3} = {};
end
cv_p_by_win = {cm_p_hist};
if exist('series_meta_cm_2007','var')
    cv_p_by_win{2} = series_meta_cm_2007.cm2007_p;
else
    cv_p_by_win{2} = NaN(n_series,1);
end
if exist('series_meta_cm_full','var')
    cv_p_by_win{3} = series_meta_cm_full.cmfull_p;
else
    cv_p_by_win{3} = NaN(n_series,1);
end

row_win    = strings(0,1);
row_grp    = strings(0,1);
row_meas   = strings(0,1);
row_n      = [];
row_0b     = [];
row_1b     = [];
row_2p     = [];
row_3p     = [];
row_med    = [];
row_maxb   = [];
row_cap8   = [];
row_pmod   = [];
row_lsig   = [];

for wi = 1:size(win_defs,1)
    win_lbl     = win_defs{wi,1};
    inc_mask    = win_defs{wi,2};
    uv_nb_w     = win_defs{wi,3};
    uv_date_w   = win_defs{wi,4};
    cv_nb_w     = win_defs{wi,5};
    cv_date_w   = win_defs{wi,6};
    vobj_w      = win_defs{wi,7};
    dates_w     = win_defs{wi,8};
    % cv_found/p not needed here — use per-window cv residuals below

    cv_resid_w  = cv_resid_by_win{wi};
    cv_p_w      = cv_p_by_win{wi};

    for gi = 1:numel(groups)
        grp   = groups(gi);
        gmask = string(fred_group(:)) == grp & inc_mask;

        for meas = ["UV","CV"]
            if meas == "UV"
                nb_v   = uv_nb_w(gmask);
                date_v = uv_date_w(gmask);
            else
                nb_v   = cv_nb_w(gmask);
                date_v = cv_date_w(gmask);
            end

            % Sigma ratio computation
            all_lsig = [];
            idx_g = find(gmask);
            for jj = 1:numel(idx_g)
                j = idx_g(jj);  % global series index
                if nb_v(jj) < 1; continue; end
                if meas == "UV"
                    [~, br_ls] = parse_break_sigma_local( ...
                        vobj_w(:,j), dates_w, date_v(jj));
                else
                    if isempty(cv_resid_w) || j > numel(cv_resid_w) || ...
                            isempty(cv_resid_w{j}); continue; end
                    pj  = cv_p_w(j);
                    ej  = cv_resid_w{j};
                    sj  = sqrt(pi/2)*abs(ej);
                    Tj2 = length(sj);
                    dj  = dates_w((pj+1):end);
                    if length(dj) > Tj2; dj = dj(end-Tj2+1:end); end
                    [~, br_ls] = parse_break_sigma_local(sj, dj, date_v(jj));
                end
                all_lsig = [all_lsig; br_ls(:)]; %#ok<AGROW>
            end

            row_win(end+1,1)  = win_lbl;
            row_grp(end+1,1)  = grp;
            row_meas(end+1,1) = meas;
            row_n(end+1,1)    = sum(~isnan(nb_v));
            row_0b(end+1,1)   = mean(nb_v == 0, 'omitnan');
            row_1b(end+1,1)   = mean(nb_v == 1, 'omitnan');
            row_2p(end+1,1)   = mean(nb_v >= 2, 'omitnan');
            row_3p(end+1,1)   = mean(nb_v >= 3, 'omitnan');
            row_med(end+1,1)  = median(nb_v, 'omitnan');
            row_maxb(end+1,1) = max(nb_v, [], 'omitnan');
            row_cap8(end+1,1) = sum(nb_v == MAX_BREAKS_F74, 'omitnan');
            if ~isempty(all_lsig)
                row_lsig(end+1,1) = median(all_lsig,'omitnan');
                row_pmod(end+1,1) = mean(all_lsig < 0, 'omitnan');
            else
                row_lsig(end+1,1) = NaN;
                row_pmod(end+1,1) = NaN;
            end
        end
    end
end

grp_tbl = table(row_win, row_grp, row_meas, row_n, ...
    row_0b, row_1b, row_2p, row_3p, row_med, row_maxb, row_cap8, ...
    row_lsig, row_pmod, ...
    'VariableNames',{ ...
        'window','fred_group','measure','n_series', ...
        'share_0b','share_1b','share_2plus','share_3plus', ...
        'median_breaks','max_breaks','n_at_cap5', ...
        'median_log_sigma_ratio','share_moderating'});

writetable(grp_tbl, fullfile(OUT_TABLES,'mb_fixed74_group_summary.csv'));
fprintf('[run_mb_fixed74] Saved mb_fixed74_group_summary.csv\n');

fprintf('\n[run_mb_fixed74] Done. Main spec: min_seg_eff=%d, max_breaks=%d\n', ...
    MIN_SEG_EFF_F74, MAX_BREAKS_F74);
fprintf('  Files written:\n');
fprintf('    data_intermediate/mb_fixed74_results.mat\n');
fprintf('    output/tables/mb_fixed74_results.csv\n');
fprintf('    output/tables/mb_fixed74_group_summary.csv\n');

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function [br_yrs, br_lsig] = parse_break_sigma_local(zt, dates, bstr)
    br_yrs  = [];
    br_lsig = [];
    if bstr == "" || strtrim(bstr) == ""; return; end
    parts = strtrim(strsplit(bstr, '|'));
    parts = parts(~cellfun(@(x) isempty(strtrim(x)), cellstr(parts)));
    if isempty(parts); return; end
    try
        break_dates = datetime(parts, 'InputFormat','dd-MMM-yyyy');
    catch; return; end
    break_dates = sort(break_dates(:));
    sentinel    = dates(end) + caldays(32);
    boundaries  = [dates(1); break_dates; sentinel];
    for k = 1:numel(break_dates)
        pre  = dates >= boundaries(k)   & dates < boundaries(k+1);
        post = dates >= boundaries(k+1) & dates < boundaries(k+2);
        mp   = nanmean(zt(pre));
        mn   = nanmean(zt(post));
        if isnan(mp)||isnan(mn)||mp<=0||mn<=0; continue; end
        br_yrs(end+1,1)  = year(break_dates(k)); %#ok<AGROW>
        br_lsig(end+1,1) = log(mn/mp);           %#ok<AGROW>
    end
end
