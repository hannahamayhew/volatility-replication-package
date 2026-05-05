%% run_mb_lm_stationarity_robust.m
% =========================================================================
% APPENDIX ROBUSTNESS CHECK — Labor Market: Stationarity-Based Exclusions
%
% Motivation: Full-sample zero-break rates are elevated for Labor Market.
% Two series fail both ADF and KPSS (nonstationary); several others pass
% ADF but are rejected by KPSS (ambiguous), indicating persistent level or
% variance shifts that are difficult to distinguish from unit roots in long
% samples. This check drops the flagged Labor Market series and re-runs the
% fixed-74 recursive MB, asking whether zero-break dominance survives.
%
% EXCLUSION CRITERION: based on full-sample stationarity tests
%   (stationarity_full_results.csv, status IN {'nonstationary','ambiguous'})
%   applied to the Labor Market group only.
%
% OUTPUTS (new files — nothing existing overwritten):
%   output/tables/mb_lm_stat_excl_table.csv    — excluded series + stats
%   output/tables/mb_lm_stat_robust_results.csv — per-series MB results
%   output/tables/mb_lm_stat_robust_compare.csv — break dist comparison
%
% PREREQUISITES:
%   run_stationarity_tests.m (writes stationarity_full_results.csv)
%   run_uv_historical / _2007 / _fullsample
%   run_cm_historical (required); run_cm_2007, run_cm_fullsample (optional)
%   run_mb_fixed74.m (writes mb_fixed74_results.mat for full-LM baseline)
% =========================================================================

run('config.m');

MAX_BREAKS_LMSR = 5;
MIN_SEG_LMSR    = 74;

fprintf('\n[run_mb_lm_stationarity_robust] Labor Market stationarity robustness check\n');
fprintf('  Spec: max_breaks=%d, min_seg_eff=%d\n', MAX_BREAKS_LMSR, MIN_SEG_LMSR);
fprintf('  Exclusion basis: full-sample ADF+KPSS (nonstationary or ambiguous)\n');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','dates_ext_2007','dates_ext_full', ...
    'series','include_main_hist','include_main_2007','include_main_full_ext', ...
    'fred_group');
load(fullfile(DATA_INT,'uv_hist_results.mat'),  'vol_obj');
load(fullfile(DATA_INT,'uv_2007_results.mat'),  'vol_obj_2007');
load(fullfile(DATA_INT,'uv_full_results.mat'),  'vol_obj_full');
load(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'cm_supw_resid_break_hist','series_meta_cm_hist');

n_series           = length(series);
cm_p_hist          = series_meta_cm_hist.cm_p_hist;
cm_supw_found_hist = series_meta_cm_hist.cm_supw_found_hist;

% Pre-compute group mask for efficiency
is_lm = string(fred_group(:)) == "Labor Market";

%% -----------------------------------------------------------------------
% 2. BUILD EXCLUSION MASK FROM FULL-SAMPLE STATIONARITY RESULTS
%% -----------------------------------------------------------------------
stat_csv = fullfile(OUT_TABLES,'stationarity_full_results.csv');
if ~exist(stat_csv,'file')
    error('[run_mb_lm_stationarity_robust] stationarity_full_results.csv not found. Run run_stationarity_tests.m first.');
end
stat_tbl = readtable(stat_csv, 'TextType','string');

is_excl     = false(n_series,1);  % true = excluded (LM + nonstationary or ambiguous)
excl_names  = strings(0,1);
excl_tcode  = NaN(0,1);
excl_adf_h  = false(0,1);
excl_adf_stat = NaN(0,1);
excl_adf_p  = NaN(0,1);
excl_kpss_stat = NaN(0,1);
excl_kpss_h = false(0,1);
excl_status = strings(0,1);

for j = 1:n_series
    if ~is_lm(j); continue; end
    sname = string(series{j});
    row   = stat_tbl(stat_tbl.series_name == sname, :);
    if isempty(row); continue; end
    st = string(row.status(1));
    if ismember(st, ["nonstationary","ambiguous"])
        is_excl(j) = true;
        excl_names(end+1,1)    = sname;
        excl_tcode(end+1,1)    = row.tcode(1);
        excl_adf_h(end+1,1)    = logical(row.adf_reject_unit_root(1));
        excl_adf_stat(end+1,1) = row.adf_stat(1);
        excl_adf_p(end+1,1)    = row.adf_pval(1);
        excl_kpss_stat(end+1,1)= row.kpss_stat(1);
        excl_kpss_h(end+1,1)   = logical(row.kpss_reject_stationarity(1));
        excl_status(end+1,1)   = st;
    end
end

n_lm   = sum(is_lm & include_main_hist(:));
n_excl = sum(is_excl);
n_ns   = sum(excl_status == "nonstationary");
n_am   = sum(excl_status == "ambiguous");
fprintf('[run_mb_lm_stationarity_robust] Labor Market in hist sample: %d\n', n_lm);
fprintf('  Excluded: %d total (%d nonstationary, %d ambiguous)\n', n_excl, n_ns, n_am);

% Sort exclusion table: nonstationary first, then ambiguous, alpha within
sort_key = double(excl_status == "ambiguous");
[~, ord] = sortrows([sort_key, double(sort(excl_names))]);  % secondary: name
[~, ord] = sort(sort_key);
excl_tbl = table(excl_names(ord), excl_tcode(ord), ...
    excl_adf_h(ord), excl_adf_stat(ord), excl_adf_p(ord), ...
    excl_kpss_stat(ord), excl_kpss_h(ord), excl_status(ord), ...
    'VariableNames',{'series_name','tcode', ...
        'adf_reject_H0_unit_root','adf_stat','adf_pval', ...
        'kpss_stat','kpss_reject_H0_stationary','status'});
writetable(excl_tbl, fullfile(OUT_TABLES,'mb_lm_stat_excl_table.csv'));
fprintf('[run_mb_lm_stationarity_robust] Saved mb_lm_stat_excl_table.csv\n');
disp(excl_tbl(:,{'series_name','tcode','adf_pval','kpss_stat','status'}));

%% -----------------------------------------------------------------------
% 3. UV BREAKS — RESTRICTED LABOR MARKET
%% -----------------------------------------------------------------------
uv_lmr_nbreaks_hist = NaN(n_series,1);
uv_lmr_nbreaks_2007 = NaN(n_series,1);
uv_lmr_nbreaks_full = NaN(n_series,1);

% Historical
fprintf('[run_mb_lm_stationarity_robust] UV historical...\n');
for j = 1:n_series
    if ~(include_main_hist(j) && is_lm(j) && ~is_excl(j)); continue; end
    res = uv_recursive_breaks(vol_obj(:,j), dates_rep_t, TRIM_FRAC, ...
        MAX_BREAKS_LMSR, MIN_SEG_LEN, MIN_SEG_LMSR);
    uv_lmr_nbreaks_hist(j) = res.n_breaks;
end
tmp = uv_lmr_nbreaks_hist(include_main_hist(:) & is_lm & ~is_excl);
fprintf('  UV hist (n=%d): 0b=%.3f 1b=%.3f 2+=%.3f med=%.1f\n', ...
    sum(~isnan(tmp)), mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
    mean(tmp>=2,'omitnan'), median(tmp,'omitnan'));

% 2007
fprintf('[run_mb_lm_stationarity_robust] UV 2007...\n');
for j = 1:n_series
    if ~(include_main_2007(j) && is_lm(j) && ~is_excl(j)); continue; end
    res = uv_recursive_breaks(vol_obj_2007(:,j), dates_ext_2007, TRIM_FRAC, ...
        MAX_BREAKS_LMSR, MIN_SEG_LEN, MIN_SEG_LMSR);
    uv_lmr_nbreaks_2007(j) = res.n_breaks;
end
tmp = uv_lmr_nbreaks_2007(include_main_2007(:) & is_lm & ~is_excl);
fprintf('  UV 2007 (n=%d): 0b=%.3f 1b=%.3f 2+=%.3f med=%.1f\n', ...
    sum(~isnan(tmp)), mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
    mean(tmp>=2,'omitnan'), median(tmp,'omitnan'));

% Full sample
fprintf('[run_mb_lm_stationarity_robust] UV full sample...\n');
for j = 1:n_series
    if ~(include_main_full_ext(j) && is_lm(j) && ~is_excl(j)); continue; end
    res = uv_recursive_breaks(vol_obj_full(:,j), dates_ext_full, TRIM_FRAC, ...
        MAX_BREAKS_LMSR, MIN_SEG_LEN, MIN_SEG_LMSR);
    uv_lmr_nbreaks_full(j) = res.n_breaks;
end
tmp = uv_lmr_nbreaks_full(include_main_full_ext(:) & is_lm & ~is_excl);
fprintf('  UV full (n=%d): 0b=%.3f 1b=%.3f 2+=%.3f med=%.1f\n', ...
    sum(~isnan(tmp)), mean(tmp==0,'omitnan'), mean(tmp==1,'omitnan'), ...
    mean(tmp>=2,'omitnan'), median(tmp,'omitnan'));

%% -----------------------------------------------------------------------
% 4. CV BREAKS — RESTRICTED LABOR MARKET
%% -----------------------------------------------------------------------
cv_lmr_nbreaks_hist = NaN(n_series,1);
cv_lmr_nbreaks_2007 = NaN(n_series,1);
cv_lmr_nbreaks_full = NaN(n_series,1);

% Historical
fprintf('[run_mb_lm_stationarity_robust] CV historical...\n');
for j = 1:n_series
    if ~(include_main_hist(j) && is_lm(j) && ~is_excl(j)); continue; end
    if ~(cm_supw_found_hist(j) && ~isempty(cm_supw_resid_break_hist{j})); continue; end
    pj          = cm_p_hist(j);
    ehat_j      = cm_supw_resid_break_hist{j};
    s_j         = sqrt(pi/2) * abs(ehat_j);
    Tj          = length(s_j);
    dates_reg_j = dates_rep_t((pj+1):end);
    dates_reg_j = dates_reg_j(end-Tj+1:end);
    res = cv_recursive_breaks(s_j, dates_reg_j, TRIM_FRAC, ...
        MAX_BREAKS_LMSR, MIN_SEG_LEN, MIN_SEG_LMSR);
    cv_lmr_nbreaks_hist(j) = res.n_breaks;
end
tmp = cv_lmr_nbreaks_hist(include_main_hist(:) & is_lm & ~is_excl);
tmp = tmp(~isnan(tmp));
fprintf('  CV hist (n=%d): 0b=%.3f 1b=%.3f 2+=%.3f med=%.1f\n', ...
    numel(tmp), mean(tmp==0), mean(tmp==1), mean(tmp>=2), median(tmp));

% 2007 (optional — skipped if cm_2007_results.mat absent)
cm2007_path = fullfile(DATA_INT,'cm_2007_results.mat');
if exist(cm2007_path,'file')
    fprintf('[run_mb_lm_stationarity_robust] CV 2007...\n');
    load(cm2007_path, 'cm_supw_resid_break_2007','series_meta_cm_2007');
    for j = 1:n_series
        if ~(include_main_2007(j) && is_lm(j) && ~is_excl(j)); continue; end
        if ~(series_meta_cm_2007.cm2007_supw_found(j) && ~isempty(cm_supw_resid_break_2007{j})); continue; end
        pj          = series_meta_cm_2007.cm2007_p(j);
        ehat_j      = cm_supw_resid_break_2007{j};
        s_j         = sqrt(pi/2) * abs(ehat_j);
        Tj          = length(s_j);
        dates_reg_j = dates_ext_2007((pj+1):end);
        dates_reg_j = dates_reg_j(end-Tj+1:end);
        res = cv_recursive_breaks(s_j, dates_reg_j, TRIM_FRAC, ...
            MAX_BREAKS_LMSR, MIN_SEG_LEN, MIN_SEG_LMSR);
        cv_lmr_nbreaks_2007(j) = res.n_breaks;
    end
    tmp = cv_lmr_nbreaks_2007(include_main_2007(:) & is_lm & ~is_excl);
    tmp = tmp(~isnan(tmp));
    fprintf('  CV 2007 (n=%d): 0b=%.3f 1b=%.3f 2+=%.3f med=%.1f\n', ...
        numel(tmp), mean(tmp==0), mean(tmp==1), mean(tmp>=2), median(tmp));
else
    fprintf('[run_mb_lm_stationarity_robust] cm_2007_results.mat not found — CV 2007 skipped.\n');
end

% Full sample (optional)
cm_full_path = fullfile(DATA_INT,'cm_full_results.mat');
if exist(cm_full_path,'file')
    fprintf('[run_mb_lm_stationarity_robust] CV full sample...\n');
    load(cm_full_path, 'cm_supw_resid_break_full','series_meta_cm_full');
    for j = 1:n_series
        if ~(include_main_full_ext(j) && is_lm(j) && ~is_excl(j)); continue; end
        if ~(series_meta_cm_full.cmfull_supw_found(j) && ~isempty(cm_supw_resid_break_full{j})); continue; end
        pj          = series_meta_cm_full.cmfull_p(j);
        ehat_j      = cm_supw_resid_break_full{j};
        s_j         = sqrt(pi/2) * abs(ehat_j);
        Tj          = length(s_j);
        dates_reg_j = dates_ext_full((pj+1):end);
        if length(dates_reg_j) > Tj; dates_reg_j = dates_reg_j(end-Tj+1:end); end
        res = cv_recursive_breaks(s_j, dates_reg_j, TRIM_FRAC, ...
            MAX_BREAKS_LMSR, MIN_SEG_LEN, MIN_SEG_LMSR);
        cv_lmr_nbreaks_full(j) = res.n_breaks;
    end
    tmp = cv_lmr_nbreaks_full(include_main_full_ext(:) & is_lm & ~is_excl);
    tmp = tmp(~isnan(tmp));
    fprintf('  CV full (n=%d): 0b=%.3f 1b=%.3f 2+=%.3f med=%.1f\n', ...
        numel(tmp), mean(tmp==0), mean(tmp==1), mean(tmp>=2), median(tmp));
else
    fprintf('[run_mb_lm_stationarity_robust] cm_full_results.mat not found — CV full skipped.\n');
end

%% -----------------------------------------------------------------------
% 5. PER-SERIES RESULTS TABLE (Labor Market only)
%% -----------------------------------------------------------------------
lm_idx = find(is_lm);
res_names  = strings(numel(lm_idx),1);
res_excl   = false(numel(lm_idx),1);
res_status = strings(numel(lm_idx),1);
res_inc_h  = false(numel(lm_idx),1);
res_inc_07 = false(numel(lm_idx),1);
res_inc_f  = false(numel(lm_idx),1);
res_uv_h   = NaN(numel(lm_idx),1);
res_uv_07  = NaN(numel(lm_idx),1);
res_uv_f   = NaN(numel(lm_idx),1);
res_cv_h   = NaN(numel(lm_idx),1);
res_cv_07  = NaN(numel(lm_idx),1);
res_cv_f   = NaN(numel(lm_idx),1);

for ii = 1:numel(lm_idx)
    j = lm_idx(ii);
    res_names(ii)  = string(series{j});
    res_excl(ii)   = is_excl(j);
    res_inc_h(ii)  = logical(include_main_hist(j));
    res_inc_07(ii) = logical(include_main_2007(j));
    res_inc_f(ii)  = logical(include_main_full_ext(j));
    res_uv_h(ii)   = uv_lmr_nbreaks_hist(j);
    res_uv_07(ii)  = uv_lmr_nbreaks_2007(j);
    res_uv_f(ii)   = uv_lmr_nbreaks_full(j);
    res_cv_h(ii)   = cv_lmr_nbreaks_hist(j);
    res_cv_07(ii)  = cv_lmr_nbreaks_2007(j);
    res_cv_f(ii)   = cv_lmr_nbreaks_full(j);
    % Retrieve stationarity status for retained series
    sname = string(series{j});
    row   = stat_tbl(stat_tbl.series_name == sname, :);
    if ~isempty(row)
        res_status(ii) = string(row.status(1));
    else
        res_status(ii) = "not_tested";
    end
end

res_tbl = table(res_names, res_excl, res_status, ...
    res_inc_h, res_inc_07, res_inc_f, ...
    res_uv_h, res_uv_07, res_uv_f, ...
    res_cv_h, res_cv_07, res_cv_f, ...
    'VariableNames',{'series_name','excluded','stat_status_full', ...
        'include_hist','include_2007','include_full', ...
        'uv_nbreaks_hist','uv_nbreaks_2007','uv_nbreaks_full', ...
        'cv_nbreaks_hist','cv_nbreaks_2007','cv_nbreaks_full'});

writetable(res_tbl, fullfile(OUT_TABLES,'mb_lm_stat_robust_results.csv'));
fprintf('[run_mb_lm_stationarity_robust] Saved mb_lm_stat_robust_results.csv (%d LM series)\n', height(res_tbl));

%% -----------------------------------------------------------------------
% 6. COMPARISON TABLE: full LM vs restricted LM × window × measure
%    Loads full-LM break counts from mb_fixed74_results.mat
%% -----------------------------------------------------------------------
f74_path = fullfile(DATA_INT,'mb_fixed74_results.mat');
if ~exist(f74_path,'file')
    warning('[run_mb_lm_stationarity_robust] mb_fixed74_results.mat not found — comparison table will have NaN for full-LM columns.');
    uv_f74_hist = NaN(n_series,1); cv_f74_hist = NaN(n_series,1);
    uv_f74_2007 = NaN(n_series,1); cv_f74_2007 = NaN(n_series,1);
    uv_f74_full = NaN(n_series,1); cv_f74_full = NaN(n_series,1);
else
    load(f74_path, ...
        'uv_mb_nbreaks_hist','cv_mb_nbreaks_hist', ...
        'uv2007_mb_nbreaks','cv2007_mb_nbreaks', ...
        'uvfull_mb_nbreaks','cvfull_mb_nbreaks');
    uv_f74_hist = uv_mb_nbreaks_hist; cv_f74_hist = cv_mb_nbreaks_hist;
    uv_f74_2007 = uv2007_mb_nbreaks;  cv_f74_2007 = cv2007_mb_nbreaks;
    uv_f74_full = uvfull_mb_nbreaks;  cv_f74_full = cvfull_mb_nbreaks;
end

% Define windows and corresponding break-count vectors
win_labels  = {'1959:03-1999:12','1959:03-2007:12','1959:03-2026:01'};
inc_masks   = {include_main_hist(:), include_main_2007(:), include_main_full_ext(:)};
uv_full_nb  = {uv_f74_hist,         uv_f74_2007,          uv_f74_full};
cv_full_nb  = {cv_f74_hist,         cv_f74_2007,          cv_f74_full};
uv_restr_nb = {uv_lmr_nbreaks_hist, uv_lmr_nbreaks_2007,  uv_lmr_nbreaks_full};
cv_restr_nb = {cv_lmr_nbreaks_hist, cv_lmr_nbreaks_2007,  cv_lmr_nbreaks_full};

cmp_win  = strings(0,1);
cmp_meas = strings(0,1);
cmp_samp = strings(0,1);
cmp_n    = [];
cmp_0b   = [];
cmp_1b   = [];
cmp_2p   = [];
cmp_3p   = [];
cmp_med  = [];

for wi = 1:3
    inc_full  = is_lm & logical(inc_masks{wi});
    inc_restr = inc_full & ~is_excl;
    for meas = ["UV","CV"]
        if meas == "UV"
            nb_full_v  = uv_full_nb{wi};
            nb_restr_v = uv_restr_nb{wi};
        else
            nb_full_v  = cv_full_nb{wi};
            nb_restr_v = cv_restr_nb{wi};
        end
        for samp = ["full_lm","restricted_lm"]
            if samp == "full_lm"
                nb = nb_full_v(inc_full);
            else
                nb = nb_restr_v(inc_restr);
            end
            nb = nb(~isnan(nb));
            cmp_win(end+1,1)  = win_labels{wi};
            cmp_meas(end+1,1) = meas;
            cmp_samp(end+1,1) = samp;
            cmp_n(end+1,1)    = numel(nb);
            if isempty(nb)
                cmp_0b(end+1,1) = NaN; cmp_1b(end+1,1) = NaN;
                cmp_2p(end+1,1) = NaN; cmp_3p(end+1,1) = NaN;
                cmp_med(end+1,1)= NaN;
            else
                cmp_0b(end+1,1) = mean(nb == 0);
                cmp_1b(end+1,1) = mean(nb == 1);
                cmp_2p(end+1,1) = mean(nb >= 2);
                cmp_3p(end+1,1) = mean(nb >= 3);
                cmp_med(end+1,1)= median(nb);
            end
        end
    end
end

cmp_tbl = table(cmp_win, cmp_meas, cmp_samp, cmp_n, ...
    cmp_0b, cmp_1b, cmp_2p, cmp_3p, cmp_med, ...
    'VariableNames',{'window','measure','sample','n_series', ...
        'share_0b','share_1b','share_2plus','share_3plus','median_breaks'});
writetable(cmp_tbl, fullfile(OUT_TABLES,'mb_lm_stat_robust_compare.csv'));
fprintf('[run_mb_lm_stationarity_robust] Saved mb_lm_stat_robust_compare.csv\n');

%% -----------------------------------------------------------------------
% 7. PRINT SUMMARY
%% -----------------------------------------------------------------------
fprintf('\n[run_mb_lm_stationarity_robust] Break distribution: full LM vs stationary-only LM\n');
fprintf('  %-20s %-4s %-15s %4s %6s %6s %6s %5s\n', ...
    'Window','Meas','Sample','N','0-brk','1-brk','2+brk','Med');
fprintf('  %s\n', repmat('-',1,75));
for i = 1:height(cmp_tbl)
    fprintf('  %-20s %-4s %-15s %4d %6.3f %6.3f %6.3f %5.1f\n', ...
        cmp_tbl.window(i), cmp_tbl.measure(i), cmp_tbl.sample(i), ...
        cmp_tbl.n_series(i), cmp_tbl.share_0b(i), ...
        cmp_tbl.share_1b(i), cmp_tbl.share_2plus(i), cmp_tbl.median_breaks(i));
end

fprintf('\n[run_mb_lm_stationarity_robust] Done.\n');
fprintf('  output/tables/mb_lm_stat_excl_table.csv\n');
fprintf('  output/tables/mb_lm_stat_robust_results.csv\n');
fprintf('  output/tables/mb_lm_stat_robust_compare.csv\n\n');
