%% run_mb_modern.m
% =========================================================================
% STAGE — RECURSIVE MULTIPLE-BREAK ESTIMATION: MODERN-ERA APPENDIX WINDOW
%
% Applies recursive SupW multiple-break detection on the modern-era subsample
% (1999:01 – latest available month) for:
%   UV: z_t = sqrt(pi/2) * |y_t - mean(y)|   (Sensier & van Dijk 2004)
%   CV: s_t = sqrt(pi/2) * |e_hat_t|          (CM break-model residuals)
%
% MOTIVATION (Appendix exercise):
%   The full 1959–2026 recursive window may underdetect post-1999 regime
%   changes because many earlier volatility regimes are embedded in the
%   sample.  Restricting to 1999:01–latest creates a cleaner window for
%   studying modern-era dynamics.
%
% METHODOLOGY:
%   Identical to the historical/2007/full windows in run_mb_historical.m:
%   - Same recursive SupW logic (uv_recursive_breaks / cv_recursive_breaks)
%   - Same Hansen (1997) significance criterion (p < 0.05)
%   - Same TRIM_FRAC = 0.15, MAX_BREAKS = 5
%   - Same Bai (1997) iterative repartition to convergence
%   - Uses MIN_SEG_EXT = 74 (fixed override) rather than 15%*T rule, so
%     segment geometry is comparable to the 2007 and full-sample windows.
%
% WINDOW-LENGTH NOTE:
%   T ≈ 325 months (1999:01–2026:02).
%   min_seg_len_eff = MIN_SEG_EXT = 74 (fixed; 15%*T = 49 < 74).
%   Achieving MAX_BREAKS = 5 requires (5+1)×74 = 444 obs > T ≈ 325.
%   The algorithm handles this gracefully; series at n_breaks=5 are flagged.
%   MAX_BREAKS is NOT lowered (additive-only policy).
%
% PREREQUISITES:
%   run_data_prep.m    → processed_data.mat  (yt_mod, dates_mod, include_main_mod)
%   run_uv_modern.m    → uv_modern_results.mat  (vol_obj_mod)
%   run_cm_modern.m    → cm_modern_results.mat  (cm_supw_resid_break_mod)
%   run_mb_historical.m (optional) → mb_hist_results.mat  (for comparison table)
%
% OUTPUTS
%   data_intermediate/mb_modern_results.mat
%   output/tables/mb_modern_results.csv       — per-series UV + CV MB results
%   output/tables/mb_window_comparison.csv    — cross-window summary (4 windows)
%   output/tables/mb_modern_by_group.csv      — FRED-MD group summaries (modern era)
%
% KEY COLUMNS
%   uvmod_mb_nbreaks, uvmod_mb_break_dates, uvmod_mb_break_stats, uvmod_mb_break_pvals
%   cvmod_mb_nbreaks, cvmod_mb_break_dates, cvmod_mb_break_stats, cvmod_mb_break_pvals
%
% PARAMETERS (from config.m)
%   MAX_BREAKS  = 5
%   MIN_SEG_LEN = 60  (hard floor)
%   MIN_SEG_EXT = 74  (fixed effective override used for this window)
%   MAX_BREAKS  = 5
%   MIN_SEG_LEN = 60
%   TRIM_FRAC   = 0.15
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_mb_modern] Loading inputs...\n');

load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_mod','series','include_main_mod','fred_group');
load(fullfile(DATA_INT,'uv_modern_results.mat'),  'vol_obj_mod');
load(fullfile(DATA_INT,'cm_modern_results.mat'), ...
    'cm_supw_resid_break_mod','series_meta_cm_mod');

n_series         = length(series);
cmmod_p          = series_meta_cm_mod.cmmod_p;
cmmod_supw_found = series_meta_cm_mod.cmmod_supw_found;

T_mod           = length(dates_mod);
min_seg_eff_mod = max(ceil(TRIM_FRAC * T_mod), MIN_SEG_LEN);

fprintf('[run_mb_modern] Window: %s – %s  (T=%d, min_seg_eff=%d)\n', ...
    string(dates_mod(1),'yyyy-MM'), string(dates_mod(end),'yyyy-MM'), ...
    T_mod, min_seg_eff_mod);
fprintf('[run_mb_modern] MAX_BREAKS=%d  MIN_SEG_LEN=%d  TRIM_FRAC=%.2f\n', ...
    MAX_BREAKS, MIN_SEG_LEN, TRIM_FRAC);

%% -----------------------------------------------------------------------
% 2. UV MODERN-ERA MULTIPLE-BREAK LOOP
%% -----------------------------------------------------------------------
fprintf('[run_mb_modern] UV modern-era (max=%d, min_seg=%d months fixed)...\n', ...
    MAX_BREAKS, MIN_SEG_EXT);

uvmod_mb_nbreaks     = NaN(n_series,1);
uvmod_mb_break_dates = strings(n_series,1);
uvmod_mb_break_stats = strings(n_series,1);
uvmod_mb_break_pvals = strings(n_series,1);

uvmod_rp_iters     = NaN(n_series,1);
uvmod_rp_converged = false(n_series,1);
uvmod_rp_maxdelta  = NaN(n_series,1);

for j = 1:n_series
    if ~include_main_mod(j); continue; end
    res = uv_recursive_breaks(vol_obj_mod(:,j), dates_mod, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN, MIN_SEG_EXT);
    uvmod_mb_nbreaks(j) = res.n_breaks;
    if res.n_breaks > 0
        uvmod_mb_break_dates(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
        uvmod_mb_break_stats(j) = strjoin(string(round(res.break_stats,3)),      ' | ');
        uvmod_mb_break_pvals(j) = strjoin(string(round(res.break_pvals,6)),      ' | ');
    else
        uvmod_mb_break_dates(j) = "";
        uvmod_mb_break_stats(j) = "";
        uvmod_mb_break_pvals(j) = "";
    end
    uvmod_rp_iters(j)     = res.repartition_iters;
    uvmod_rp_converged(j) = res.repartition_converged;
    uvmod_rp_maxdelta(j)  = res.repartition_max_delta;
    if mod(j,20)==0; fprintf('  UV mod: %d / %d\n', j, n_series); end
end

tmp_uv      = uvmod_mb_nbreaks(include_main_mod);
n_at_cap_uv = sum(tmp_uv == MAX_BREAKS, 'omitnan');
fprintf('  0b:%.3f  1b:%.3f  2+:%.3f  3+:%.3f  med:%.1f  max:%.0f  n_at_cap(%d)=%d\n', ...
    mean(tmp_uv==0,'omitnan'), mean(tmp_uv==1,'omitnan'), ...
    mean(tmp_uv>=2,'omitnan'), mean(tmp_uv>=3,'omitnan'), ...
    median(tmp_uv,'omitnan'), max(tmp_uv,[],'omitnan'), MAX_BREAKS, n_at_cap_uv);
if n_at_cap_uv > 0
    fprintf('  FLAG (UV): %d series hit MAX_BREAKS=%d cap. True break count may exceed cap.\n', ...
        n_at_cap_uv, MAX_BREAKS);
end

rp_mask_uv = include_main_mod(:) & ~isnan(uvmod_rp_iters(:));
if any(rp_mask_uv)
    fprintf('  Repartition: iters mean=%.1f max=%.0f  non-converged=%d/%d  max_delta_final=%.0f\n', ...
        mean(uvmod_rp_iters(rp_mask_uv)), max(uvmod_rp_iters(rp_mask_uv)), ...
        sum(~uvmod_rp_converged(rp_mask_uv)), sum(rp_mask_uv), ...
        max(uvmod_rp_maxdelta(rp_mask_uv)));
end

%% -----------------------------------------------------------------------
% 3. CV MODERN-ERA MULTIPLE-BREAK LOOP
%% -----------------------------------------------------------------------
fprintf('[run_mb_modern] CV modern-era (max=%d, min_seg=%d months fixed)...\n', ...
    MAX_BREAKS, MIN_SEG_EXT);

cvmod_mb_nbreaks     = NaN(n_series,1);
cvmod_mb_break_dates = strings(n_series,1);
cvmod_mb_break_stats = strings(n_series,1);
cvmod_mb_break_pvals = strings(n_series,1);

cvmod_rp_iters     = NaN(n_series,1);
cvmod_rp_converged = false(n_series,1);
cvmod_rp_maxdelta  = NaN(n_series,1);

for j = 1:n_series
    if ~(include_main_mod(j) && cmmod_supw_found(j) && ~isempty(cm_supw_resid_break_mod{j}))
        continue;
    end
    pj          = cmmod_p(j);
    ehat_j      = cm_supw_resid_break_mod{j};
    s_j         = sqrt(pi/2) * abs(ehat_j);
    Tj          = length(s_j);
    dates_reg_j = dates_mod((pj+1):end);
    dates_reg_j = dates_reg_j(end-Tj+1:end);

    res = cv_recursive_breaks(s_j, dates_reg_j, TRIM_FRAC, MAX_BREAKS, MIN_SEG_LEN, MIN_SEG_EXT);
    cvmod_mb_nbreaks(j) = res.n_breaks;
    if res.n_breaks > 0
        cvmod_mb_break_dates(j) = strjoin(string(res.break_dates,'dd-MMM-yyyy'), ' | ');
        cvmod_mb_break_stats(j) = strjoin(string(round(res.break_stats,3)),      ' | ');
        cvmod_mb_break_pvals(j) = strjoin(string(round(res.break_pvals,6)),      ' | ');
    else
        cvmod_mb_break_dates(j) = "";
        cvmod_mb_break_stats(j) = "";
        cvmod_mb_break_pvals(j) = "";
    end
    cvmod_rp_iters(j)     = res.repartition_iters;
    cvmod_rp_converged(j) = res.repartition_converged;
    cvmod_rp_maxdelta(j)  = res.repartition_max_delta;
    if mod(j,20)==0; fprintf('  CV mod: %d / %d\n', j, n_series); end
end

cv_include_mod = include_main_mod(:) & cmmod_supw_found(:) & ...
    ~cellfun(@isempty, cm_supw_resid_break_mod);
tmp_cv      = cvmod_mb_nbreaks(cv_include_mod);
n_at_cap_cv = sum(tmp_cv == MAX_BREAKS, 'omitnan');
fprintf('  0b:%.3f  1b:%.3f  2+:%.3f  3+:%.3f  med:%.1f  max:%.0f  n_at_cap(%d)=%d\n', ...
    mean(tmp_cv==0,'omitnan'), mean(tmp_cv==1,'omitnan'), ...
    mean(tmp_cv>=2,'omitnan'), mean(tmp_cv>=3,'omitnan'), ...
    median(tmp_cv,'omitnan'), max(tmp_cv,[],'omitnan'), MAX_BREAKS, n_at_cap_cv);
if n_at_cap_cv > 0
    fprintf('  FLAG (CV): %d series hit MAX_BREAKS=%d cap. True break count may exceed cap.\n', ...
        n_at_cap_cv, MAX_BREAKS);
end

rp_mask_cv = cv_include_mod(:) & ~isnan(cvmod_rp_iters(:));
if any(rp_mask_cv)
    fprintf('  Repartition: iters mean=%.1f max=%.0f  non-converged=%d/%d  max_delta_final=%.0f\n', ...
        mean(cvmod_rp_iters(rp_mask_cv)), max(cvmod_rp_iters(rp_mask_cv)), ...
        sum(~cvmod_rp_converged(rp_mask_cv)), sum(rp_mask_cv), ...
        max(cvmod_rp_maxdelta(rp_mask_cv)));
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE PER-SERIES RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_mb_modern = table( ...
    string(series(:)), include_main_mod(:), cv_include_mod(:), ...
    uvmod_mb_nbreaks,     uvmod_mb_break_dates, ...
    uvmod_mb_break_stats, uvmod_mb_break_pvals, ...
    cvmod_mb_nbreaks,     cvmod_mb_break_dates, ...
    cvmod_mb_break_stats, cvmod_mb_break_pvals, ...
    'VariableNames',{ ...
        'series_name','include_main_mod','include_cv_mod', ...
        'uvmod_mb_nbreaks',     'uvmod_mb_break_dates', ...
        'uvmod_mb_break_stats', 'uvmod_mb_break_pvals', ...
        'cvmod_mb_nbreaks',     'cvmod_mb_break_dates', ...
        'cvmod_mb_break_stats', 'cvmod_mb_break_pvals'});

%% -----------------------------------------------------------------------
% 5. SAVE PER-SERIES RESULTS
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'mb_modern_results.mat'), 'series_meta_mb_modern');
writetable(series_meta_mb_modern, fullfile(OUT_TABLES,'mb_modern_results.csv'));
fprintf('[run_mb_modern] Saved mb_modern_results.mat and mb_modern_results.csv\n');

%% -----------------------------------------------------------------------
% 6. CROSS-WINDOW COMPARISON TABLE
%   Loads existing window results from mb_hist_results.mat (if available)
%   and assembles a single comparison table.
%   If mb_hist_results.mat is absent, the hist/2007/full rows are NaN.
%% -----------------------------------------------------------------------
fprintf('[run_mb_modern] Building cross-window comparison table...\n');

% -- compute_win_stats: safe stats for a break-count vector after masking --
% Defined as a local function at the end of this file.
% Returns struct with fields: n, s0, s1, s2, s3, med, mx, ncap.
% All fields are NaN if the filtered vector is empty.
% ----- (called below) -----

% Summary stats for modern-era rows (always available at this point)
% tmp_uv and tmp_cv already exclude NaN (computed in sections 2 and 3)
sw_uvmod = compute_win_stats(tmp_uv, true(length(tmp_uv),1), MAX_BREAKS);
sw_cvmod = compute_win_stats(tmp_cv, true(length(tmp_cv),1), MAX_BREAKS);
% Override ncap with the already-computed value that also handles NaN
sw_uvmod.ncap = n_at_cap_uv;
sw_cvmod.ncap = n_at_cap_cv;

% Row layout: 1=hist-UV, 2=hist-CV, 3=2007-UV, 4=2007-CV,
%             5=full-UV,  6=full-CV,  7=mod-UV,  8=mod-CV
comp_windows  = {'hist';'hist';'2007';'2007';'full';'full';'modern';'modern'};
comp_measures = {'UV';'CV';'UV';'CV';'UV';'CV';'UV';'CV'};
comp_n        = NaN(8,1);  comp_0b   = NaN(8,1);  comp_1b   = NaN(8,1);
comp_2bp      = NaN(8,1);  comp_3bp  = NaN(8,1);  comp_med  = NaN(8,1);
comp_max      = NaN(8,1);  comp_ncap = NaN(8,1);  comp_capflag = repmat("",8,1);

% Fill modern-era row (always available)
sw_list = [{struct()},{struct()},{struct()},{struct()}, ...
           {struct()},{struct()},{sw_uvmod},{sw_cvmod}];
for ri_fill = 7:8
    sw = sw_list{ri_fill};
    comp_n(ri_fill)=sw.n; comp_0b(ri_fill)=sw.s0; comp_1b(ri_fill)=sw.s1;
    comp_2bp(ri_fill)=sw.s2; comp_3bp(ri_fill)=sw.s3;
    comp_med(ri_fill)=sw.med; comp_max(ri_fill)=sw.mx;
    comp_ncap(ri_fill)=sw.ncap;
    if sw.ncap > 0
        comp_capflag(ri_fill) = sprintf('FLAGGED:%d@cap%d', sw.ncap, MAX_BREAKS);
    end
end

mb_hist_path = fullfile(DATA_INT,'mb_hist_results.mat');
if exist(mb_hist_path,'file')
    load(mb_hist_path, 'series_meta_mb_hist');

    pairs = { ...
        series_meta_mb_hist.include_main_hist(:), ...
            series_meta_mb_hist.uv_mb_nbreaks_hist, 1; ...
        series_meta_mb_hist.include_main_hist(:), ...
            series_meta_mb_hist.cv_mb_nbreaks_hist, 2; ...
        series_meta_mb_hist.include_main_2007(:), ...
            series_meta_mb_hist.uv2007_mb_nbreaks, 3; ...
        series_meta_mb_hist.include_main_2007(:), ...
            series_meta_mb_hist.cv2007_mb_nbreaks, 4; ...
        series_meta_mb_hist.include_main_full_ext(:), ...
            series_meta_mb_hist.uvfull_mb_nbreaks, 5; ...
        series_meta_mb_hist.include_main_full_ext(:), ...
            series_meta_mb_hist.cvfull_mb_nbreaks, 6 };

    for pk = 1:size(pairs,1)
        inc_k  = pairs{pk,1};
        nb_k   = pairs{pk,2};
        ri_k   = pairs{pk,3};
        sw     = compute_win_stats(nb_k, inc_k, MAX_BREAKS);
        comp_n(ri_k)=sw.n; comp_0b(ri_k)=sw.s0; comp_1b(ri_k)=sw.s1;
        comp_2bp(ri_k)=sw.s2; comp_3bp(ri_k)=sw.s3;
        comp_med(ri_k)=sw.med; comp_max(ri_k)=sw.mx;
        comp_ncap(ri_k)=sw.ncap;
        if sw.ncap > 0
            comp_capflag(ri_k) = sprintf('FLAGGED:%d@cap%d', sw.ncap, MAX_BREAKS);
        end
    end

    fprintf('[run_mb_modern] Loaded mb_hist_results.mat — hist/2007/full comparison rows filled.\n');
else
    fprintf('[run_mb_modern] WARNING: mb_hist_results.mat not found.\n');
    fprintf('[run_mb_modern]   hist/2007/full comparison rows will be NaN.\n');
    fprintf('[run_mb_modern]   Run run_mb_historical.m first to populate those rows.\n');
end

comp_tbl = table( ...
    string(comp_windows), string(comp_measures), ...
    comp_n, comp_0b, comp_1b, comp_2bp, comp_3bp, ...
    comp_med, comp_max, comp_ncap, comp_capflag, ...
    'VariableNames',{ ...
        'window','measure','n_series', ...
        'share_0b','share_1b','share_2b_plus','share_3b_plus', ...
        'med_breaks','max_breaks_obs','n_at_cap','cap_flag'});

writetable(comp_tbl, fullfile(OUT_TABLES,'mb_window_comparison.csv'));
fprintf('[run_mb_modern] Saved mb_window_comparison.csv\n');

fprintf('\n  === CROSS-WINDOW COMPARISON (UV) ===\n');
fprintf('  %-8s  %4s  %5s  %5s  %5s  %5s  %3s  %3s  %5s\n', ...
    'Window','N','0b','1b','2+b','3+b','Med','Max','nCap');
for ri = 1:2:8
    n_=comp_n(ri); s0=comp_0b(ri); s1=comp_1b(ri);
    s2=comp_2bp(ri); s3=comp_3bp(ri); md=comp_med(ri); mx=comp_max(ri); nc=comp_ncap(ri);
    if isnan(n_)
        fprintf('  %-8s  %4s  %5s  %5s  %5s  %5s  %3s  %3s  %5s\n', ...
            comp_windows{ri},'NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN');
    else
        fprintf('  %-8s  %4d  %5.3f  %5.3f  %5.3f  %5.3f  %3.1f  %3.0f  %5d\n', ...
            comp_windows{ri}, n_, s0, s1, s2, s3, md, mx, nc);
    end
end

fprintf('\n  === CROSS-WINDOW COMPARISON (CV) ===\n');
fprintf('  %-8s  %4s  %5s  %5s  %5s  %5s  %3s  %3s  %5s\n', ...
    'Window','N','0b','1b','2+b','3+b','Med','Max','nCap');
for ri = 2:2:8
    n_=comp_n(ri); s0=comp_0b(ri); s1=comp_1b(ri);
    s2=comp_2bp(ri); s3=comp_3bp(ri); md=comp_med(ri); mx=comp_max(ri); nc=comp_ncap(ri);
    if isnan(n_)
        fprintf('  %-8s  %4s  %5s  %5s  %5s  %5s  %3s  %3s  %5s\n', ...
            comp_windows{ri},'NaN','NaN','NaN','NaN','NaN','NaN','NaN','NaN');
    else
        fprintf('  %-8s  %4d  %5.3f  %5.3f  %5.3f  %5.3f  %3.1f  %3.0f  %5d\n', ...
            comp_windows{ri}, n_, s0, s1, s2, s3, md, mx, nc);
    end
end

%% -----------------------------------------------------------------------
% 7. FRED-MD GROUP SUMMARIES (MODERN-ERA WINDOW)
%
%   Priority sectors for the modern-era exercise:
%     Housing / Construction, Interest and Exchange Rates,
%     Prices, Labor Market
%   All groups are reported; interpretation should be cautious for N < 5.
%% -----------------------------------------------------------------------
fprintf('[run_mb_modern] Computing FRED-MD group summaries...\n');

group_names_all = unique(string(fred_group(:)));
group_names_all = group_names_all(group_names_all ~= "");

priority_groups = ["Housing / Construction", "Interest and Exchange Rates", ...
                   "Prices", "Labor Market"];

if isempty(group_names_all)
    fprintf('[run_mb_modern] WARNING: fred_group not populated — group summary skipped.\n');
    fprintf('[run_mb_modern]   Ensure series_metadata_original_order.csv is in %s\n', ...
        fileparts(DATA_RAW));
else
    n_groups = length(group_names_all);

    grp_label      = strings(n_groups,1);
    grp_priority   = false(n_groups,1);
    grp_n_uv       = zeros(n_groups,1);
    grp_uv_0b      = NaN(n_groups,1);
    grp_uv_1b      = NaN(n_groups,1);
    grp_uv_2bplus  = NaN(n_groups,1);
    grp_uv_med     = NaN(n_groups,1);
    grp_uv_max     = NaN(n_groups,1);
    grp_n_cv       = zeros(n_groups,1);
    grp_cv_0b      = NaN(n_groups,1);
    grp_cv_1b      = NaN(n_groups,1);
    grp_cv_2bplus  = NaN(n_groups,1);
    grp_cv_med     = NaN(n_groups,1);
    grp_cv_max     = NaN(n_groups,1);

    fg_col = string(fred_group(:));

    for gi = 1:n_groups
        gname  = group_names_all(gi);
        g_mask = fg_col == gname;

        uv_mask = g_mask & include_main_mod(:);
        uv_nb   = uvmod_mb_nbreaks(uv_mask);
        uv_nb   = uv_nb(~isnan(uv_nb));

        cv_mask = g_mask & cv_include_mod(:);
        cv_nb   = cvmod_mb_nbreaks(cv_mask);
        cv_nb   = cv_nb(~isnan(cv_nb));

        grp_label(gi)    = gname;
        grp_priority(gi) = any(priority_groups == gname);
        grp_n_uv(gi)     = length(uv_nb);
        grp_n_cv(gi)     = length(cv_nb);
        if ~isempty(uv_nb)
            grp_uv_0b(gi)     = mean(uv_nb==0);
            grp_uv_1b(gi)     = mean(uv_nb==1);
            grp_uv_2bplus(gi) = mean(uv_nb>=2);
            grp_uv_med(gi)    = median(uv_nb);
            grp_uv_max(gi)    = max(uv_nb);
        end
        if ~isempty(cv_nb)
            grp_cv_0b(gi)     = mean(cv_nb==0);
            grp_cv_1b(gi)     = mean(cv_nb==1);
            grp_cv_2bplus(gi) = mean(cv_nb>=2);
            grp_cv_med(gi)    = median(cv_nb);
            grp_cv_max(gi)    = max(cv_nb);
        end
    end

    grp_tbl = table( ...
        grp_label, grp_priority, ...
        grp_n_uv, grp_uv_0b, grp_uv_1b, grp_uv_2bplus, grp_uv_med, grp_uv_max, ...
        grp_n_cv, grp_cv_0b, grp_cv_1b, grp_cv_2bplus, grp_cv_med, grp_cv_max, ...
        'VariableNames',{ ...
            'group','priority_group', ...
            'n_uv','uv_share_0b','uv_share_1b','uv_share_2bplus', ...
            'uv_med_breaks','uv_max_breaks', ...
            'n_cv','cv_share_0b','cv_share_1b','cv_share_2bplus', ...
            'cv_med_breaks','cv_max_breaks'});

    % Sort: priority groups first, then alphabetically within priority
    [~, sort_ord] = sortrows(table(~grp_priority, grp_label));
    grp_tbl = grp_tbl(sort_ord,:);

    writetable(grp_tbl, fullfile(OUT_TABLES,'mb_modern_by_group.csv'));
    fprintf('[run_mb_modern] Saved mb_modern_by_group.csv\n');

    fprintf('\n  === GROUP SUMMARIES — MODERN ERA (%s – %s) ===\n', ...
        string(dates_mod(1),'yyyy-MM'), string(dates_mod(end),'yyyy-MM'));
    fprintf('  (* = priority sector for modern-era interpretation)\n');
    fprintf('  %-38s  N_UV  UV:0b  1b  2+b  med  |  N_CV  CV:0b  1b  2+b  med\n', 'Group');
    for gi = 1:height(grp_tbl)
        g   = grp_tbl(gi,:);
        pfx = ' '; if g.priority_group; pfx = '*'; end
        nflag = ''; if g.n_uv < 5; nflag = ' (n_UV<5)'; end
        fprintf('  %s %-37s %4d  %4.2f %4.2f %4.2f %4.1f  |  %4d  %4.2f %4.2f %4.2f %4.1f%s\n', ...
            pfx, char(g.group), g.n_uv, ...
            g.uv_share_0b, g.uv_share_1b, g.uv_share_2bplus, g.uv_med_breaks, ...
            g.n_cv, ...
            g.cv_share_0b, g.cv_share_1b, g.cv_share_2bplus, g.cv_med_breaks, nflag);
    end
end

%% -----------------------------------------------------------------------
% 8. FINAL SUMMARY
%% -----------------------------------------------------------------------
fprintf('\n[run_mb_modern] n_mod_uv=%d, n_mod_cv=%d\n', ...
    sum(include_main_mod), sum(cv_include_mod));
fprintf('[run_mb_modern] Outputs written:\n');
fprintf('  data_intermediate/mb_modern_results.mat\n');
fprintf('  output/tables/mb_modern_results.csv\n');
fprintf('  output/tables/mb_window_comparison.csv\n');
if ~isempty(group_names_all)
    fprintf('  output/tables/mb_modern_by_group.csv\n');
end
fprintf('[run_mb_modern] Done.\n\n');

% =========================================================================
% LOCAL FUNCTION — must appear after all executable script code (R2016b+)
% =========================================================================

function sw = compute_win_stats(nbreaks_vec, inc_mask, cap)
% COMPUTE_WIN_STATS  Safe summary statistics for a multiple-break count vector.
%
%   sw = compute_win_stats(nbreaks_vec, inc_mask, cap)
%
%   Applies the logical mask inc_mask, drops NaN values, then computes
%   summary statistics.  Returns NaN for all statistics if the resulting
%   vector is empty (handles the case where a window's CV results were
%   never populated, leaving all NaN).
%
%   Fields of output struct sw:
%     n     - number of included non-NaN series
%     s0    - share with 0 breaks
%     s1    - share with 1 break
%     s2    - share with 2+ breaks
%     s3    - share with 3+ breaks
%     med   - median breaks
%     mx    - maximum breaks observed (NaN if empty — NOT [] to avoid
%             subscripted-assignment errors in the caller)
%     ncap  - number of series at the cap (== cap)

    v = nbreaks_vec(inc_mask(:));
    v = v(~isnan(v));

    if isempty(v)
        sw = struct('n',0,'s0',NaN,'s1',NaN,'s2',NaN,'s3',NaN, ...
                    'med',NaN,'mx',NaN,'ncap',0);
    else
        sw = struct( ...
            'n',    length(v), ...
            's0',   mean(v == 0), ...
            's1',   mean(v == 1), ...
            's2',   mean(v >= 2), ...
            's3',   mean(v >= 3), ...
            'med',  median(v), ...
            'mx',   max(v), ...
            'ncap', sum(v == cap));
    end
end
