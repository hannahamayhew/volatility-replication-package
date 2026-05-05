%% run_mb_fixed74_smooth_char.m
% =========================================================================
% POST-HOC SMOOTH-BREAK CHARACTERIZATION OF RECURSIVE BREAK DATES
%
% For each series and each recursively-detected break date τ_i, fits a
% logistic smooth-transition model within the surrounding segment
% [τ_{i-1}+1, τ_{i+1}] and records:
%
%   gamma_best  — logistic speed parameter (higher = more abrupt)
%   is_abrupt   — gamma_best >= 20 (Sensier & van Dijk threshold)
%   smooth_mid  — smooth-model midpoint date (comparison to recursive τ_i)
%   mid_lag_mo  — signed lag in months: smooth_mid − τ_i (positive = smooth
%                 midpoint is LATER than recursive break date)
%
% This adds interpretive richness WITHOUT changing the detection methodology.
% The recursive break dates remain authoritative; smooth characterization is
% supplementary, describing the SHAPE of each transition.
%
% INPUTS (loaded from data_intermediate/)
%   mb_fixed74_results.mat   — break dates from run_mb_fixed74.m
%   uv_*_results.mat         — UV vol objects
%   cm_*_results.mat         — CM residuals (for CV)
%   processed_data.mat       — date vectors, include masks
%
% OUTPUTS (appended to mb_fixed74 outputs):
%   mb_fixed74_smooth_char.csv   — per-series smooth characterization
%   mb_fixed74_smooth_summary.csv — sector × window × measure summary
%       Columns: window, fred_group, measure, n_series, n_with_breaks,
%                pct_abrupt_all_breaks, pct_series_all_abrupt,
%                pct_series_any_abrupt, median_gamma, mean_mid_lag_months
%
% PREREQUISITES:
%   run_mb_fixed74.m must have completed successfully.
% =========================================================================

run('config.m');

ABRUPT_THRESH = 20;   % gamma >= 20 classified as abrupt (S&vD 2004)

fprintf('\n[smooth_char] Loading fixed-74 results...\n');

load(fullfile(DATA_INT,'mb_fixed74_results.mat'));
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','dates_ext_2007','dates_ext_full', ...
    'series','include_main_hist','include_main_2007','include_main_full_ext', ...
    'fred_group');
load(fullfile(DATA_INT,'uv_hist_results.mat'),  'vol_obj');
load(fullfile(DATA_INT,'uv_2007_results.mat'),  'vol_obj_2007');
load(fullfile(DATA_INT,'uv_full_results.mat'),  'vol_obj_full');
load(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'cm_supw_resid_break_hist','series_meta_cm_hist');

n_series  = length(series);
cm_p_hist = series_meta_cm_hist.cm_p_hist;

% Optional CM 2007 / full
cm2007_path  = fullfile(DATA_INT,'cm_2007_results.mat');
cm_full_path = fullfile(DATA_INT,'cm_full_results.mat');

have_cm_2007 = exist(cm2007_path,'file');
have_cm_full = exist(cm_full_path,'file');

if have_cm_2007
    load(cm2007_path, 'cm_supw_resid_break_2007','series_meta_cm_2007');
end
if have_cm_full
    load(cm_full_path, 'cm_supw_resid_break_full','series_meta_cm_full');
end

% =========================================================================
% HELPER: parse pipe-separated break dates → datetime array
% =========================================================================
function bd = parse_break_dates(bstr)
    bd = NaT(0,1);
    if bstr == "" || strtrim(bstr) == ""; return; end
    parts = strtrim(strsplit(bstr,'|'));
    parts = parts(~cellfun(@(x) isempty(strtrim(x)), cellstr(parts)));
    if isempty(parts); return; end
    try
        bd = sort(datetime(parts,'InputFormat','dd-MMM-yyyy'))';
    catch; end
end

% =========================================================================
% HELPER: smooth-characterise one series' break dates
%   zt        = T×1 series (UV or CV proxy)
%   dates     = T×1 datetime vector aligned with zt
%   break_str = pipe-separated break date string
%   Returns: gamma_str, abrupt_str, smid_str, lag_str (all pipe-separated)
% =========================================================================
function [gamma_str, abrupt_str, smid_str, lag_str, all_gamma, all_abrupt, all_lag] = ...
        char_breaks(zt, dates, break_str, gamma_grid, trim_frac, abrupt_thresh)

    gamma_str  = "";
    abrupt_str = "";
    smid_str   = "";
    lag_str    = "";
    all_gamma  = [];
    all_abrupt = [];
    all_lag    = [];

    bd = parse_break_dates(break_str);
    nb = numel(bd);
    if nb == 0; return; end

    % Drop NaNs from zt/dates (smooth_break_search_mean handles them but
    % let's use consistent indexing)
    keep = ~isnan(zt);
    zt_c    = zt(keep);
    dates_c = dates(keep);
    T = length(zt_c);
    if T < 20; return; end

    % Map break dates to indices in dates_c
    [~, tau_idx] = ismember(bd, dates_c);
    % If a break date isn't exact (floating dates), find nearest
    for ii = 1:nb
        if tau_idx(ii) == 0
            [~, tau_idx(ii)] = min(abs(dates_c - bd(ii)));
        end
    end
    tau_idx = sort(tau_idx);

    gamma_vals  = NaN(nb,1);
    abrupt_vals = false(nb,1);
    smid_vals   = NaT(nb,1);
    lag_vals    = NaN(nb,1);

    % Boundary indices: 0 = before series start, T+1 = after series end
    tau_L_arr = [0;         tau_idx(1:end-1)];
    tau_R_arr = [tau_idx(2:end); T          ];

    for ii = 1:nb
        tau_L = tau_L_arr(ii);
        tau_R = tau_R_arr(ii);

        i1 = tau_L + 1;
        i2 = tau_R;
        if i2 - i1 + 1 < 20; continue; end

        seg_z = zt_c(i1:i2);
        seg_d = dates_c(i1:i2);

        res = smooth_break_search_mean(seg_z, seg_d, trim_frac, gamma_grid);
        if ~res.success; continue; end

        gamma_vals(ii)  = res.gamma_best;
        abrupt_vals(ii) = res.gamma_best >= abrupt_thresh;
        smid_vals(ii)   = res.midpoint_date;

        % Lag = smooth midpoint − recursive break date (in months, approx)
        if ~isnat(res.midpoint_date)
            lag_vals(ii) = calmonths(between(bd(ii), res.midpoint_date, 'months'));
        end
    end

    valid = ~isnan(gamma_vals);
    if ~any(valid); return; end

    gamma_str  = strjoin(string(round(gamma_vals,1)),  ' | ');
    abrupt_str = strjoin(string(double(abrupt_vals)),   ' | ');
    smid_str   = strjoin(string(smid_vals,'dd-MMM-yyyy'), ' | ');
    lag_str    = strjoin(string(round(lag_vals,0)),     ' | ');

    all_gamma  = gamma_vals(valid);
    all_abrupt = abrupt_vals(valid);
    all_lag    = lag_vals(valid);
end

% =========================================================================
% WINDOW DEFINITIONS: struct array with all needed info per window
% =========================================================================
%   win_label, uv_nb, uv_dates, cv_nb, cv_dates,
%   vobj, dates_vec, cv_resid, cv_p, include_mask
% =========================================================================

W(1).label   = '1959:03-1999:12';
W(1).inc     = include_main_hist;
W(1).uv_nb   = uv_mb_nbreaks_hist;
W(1).uv_bd   = uv_mb_break_dates_hist;
W(1).cv_nb   = cv_mb_nbreaks_hist;
W(1).cv_bd   = cv_mb_break_dates_hist;
W(1).vobj    = vol_obj;
W(1).dates   = dates_rep_t;
W(1).cv_res  = cm_supw_resid_break_hist;
W(1).cv_p    = cm_p_hist;
W(1).has_cv  = true;

W(2).label   = '1959:03-2007:12';
W(2).inc     = include_main_2007;
W(2).uv_nb   = uv2007_mb_nbreaks;
W(2).uv_bd   = uv2007_mb_break_dates;
W(2).cv_nb   = cv2007_mb_nbreaks;
W(2).cv_bd   = cv2007_mb_break_dates;
W(2).vobj    = vol_obj_2007;
W(2).dates   = dates_ext_2007;
W(2).has_cv  = have_cm_2007;
if have_cm_2007
    W(2).cv_res = cm_supw_resid_break_2007;
    W(2).cv_p   = series_meta_cm_2007.cm2007_p;
else
    W(2).cv_res = {};
    W(2).cv_p   = NaN(n_series,1);
end

W(3).label   = '1959:03-2026:01';
W(3).inc     = include_main_full_ext;
W(3).uv_nb   = uvfull_mb_nbreaks;
W(3).uv_bd   = uvfull_mb_break_dates;
W(3).cv_nb   = cvfull_mb_nbreaks;
W(3).cv_bd   = cvfull_mb_break_dates;
W(3).vobj    = vol_obj_full;
W(3).dates   = dates_ext_full;
W(3).has_cv  = have_cm_full;
if have_cm_full
    W(3).cv_res = cm_supw_resid_break_full;
    W(3).cv_p   = series_meta_cm_full.cmfull_p;
else
    W(3).cv_res = {};
    W(3).cv_p   = NaN(n_series,1);
end

% =========================================================================
% MAIN LOOP: smooth-characterise all series × windows
% =========================================================================

% Output columns (one per window × measure combination)
fld_names = {'hist_uv','hist_cv','p2007_uv','p2007_cv','full_uv','full_cv'};
wi_map    = {1,'UV'; 1,'CV'; 2,'UV'; 2,'CV'; 3,'UV'; 3,'CV'};

% Preallocate string arrays for gamma/abrupt/smid/lag per field
for fi = 1:numel(fld_names)
    fn = fld_names{fi};
    gamma_out.(fn)  = strings(n_series,1);
    abrupt_out.(fn) = strings(n_series,1);
    smid_out.(fn)   = strings(n_series,1);
    lag_out.(fn)    = strings(n_series,1);
end

for fi = 1:numel(fld_names)
    wi   = wi_map{fi,1};
    meas = wi_map{fi,2};
    fn   = fld_names{fi};

    fprintf('[smooth_char] %s %s...\n', W(wi).label, meas);

    if strcmp(meas,'UV')
        nb_v  = W(wi).uv_nb;
        bd_v  = W(wi).uv_bd;
    else
        nb_v  = W(wi).cv_nb;
        bd_v  = W(wi).cv_bd;
    end

    for j = 1:n_series
        if ~W(wi).inc(j); continue; end
        if isnan(nb_v(j)) || nb_v(j) < 1; continue; end

        % Build series
        if strcmp(meas,'UV')
            zt_j = W(wi).vobj(:,j);
            dt_j = W(wi).dates;
        else
            if ~W(wi).has_cv || isempty(W(wi).cv_res) || ...
                    j > numel(W(wi).cv_res) || isempty(W(wi).cv_res{j})
                continue;
            end
            pj      = W(wi).cv_p(j);
            ej      = W(wi).cv_res{j};
            zt_j    = sqrt(pi/2) * abs(ej);
            Tj      = length(zt_j);
            dt_j    = W(wi).dates((pj+1):end);
            if length(dt_j) > Tj; dt_j = dt_j(end-Tj+1:end); end
        end

        [gstr, astr, sstr, lstr, ~, ~, ~] = char_breaks( ...
            zt_j, dt_j, bd_v(j), GAMMA_GRID, TRIM_FRAC, ABRUPT_THRESH);

        gamma_out.(fn)(j)  = gstr;
        abrupt_out.(fn)(j) = astr;
        smid_out.(fn)(j)   = sstr;
        lag_out.(fn)(j)    = lstr;

        if mod(j,30)==0
            fprintf('  %d / %d\n', j, n_series);
        end
    end
end

% =========================================================================
% SAVE PER-SERIES TABLE
% =========================================================================
sc_tbl = table( ...
    string(series(:)), string(fred_group(:)), ...
    include_main_hist(:), include_main_2007(:), include_main_full_ext(:), ...
    gamma_out.hist_uv,  abrupt_out.hist_uv,  smid_out.hist_uv,  lag_out.hist_uv, ...
    gamma_out.hist_cv,  abrupt_out.hist_cv,  smid_out.hist_cv,  lag_out.hist_cv, ...
    gamma_out.p2007_uv, abrupt_out.p2007_uv, smid_out.p2007_uv, lag_out.p2007_uv, ...
    gamma_out.p2007_cv, abrupt_out.p2007_cv, smid_out.p2007_cv, lag_out.p2007_cv, ...
    gamma_out.full_uv,  abrupt_out.full_uv,  smid_out.full_uv,  lag_out.full_uv, ...
    gamma_out.full_cv,  abrupt_out.full_cv,  smid_out.full_cv,  lag_out.full_cv, ...
    'VariableNames',{ ...
        'series_name','fred_group', ...
        'include_main_hist','include_main_2007','include_main_full_ext', ...
        'hist_uv_gamma','hist_uv_abrupt','hist_uv_smooth_mid','hist_uv_lag_mo', ...
        'hist_cv_gamma','hist_cv_abrupt','hist_cv_smooth_mid','hist_cv_lag_mo', ...
        'p2007_uv_gamma','p2007_uv_abrupt','p2007_uv_smooth_mid','p2007_uv_lag_mo', ...
        'p2007_cv_gamma','p2007_cv_abrupt','p2007_cv_smooth_mid','p2007_cv_lag_mo', ...
        'full_uv_gamma','full_uv_abrupt','full_uv_smooth_mid','full_uv_lag_mo', ...
        'full_cv_gamma','full_cv_abrupt','full_cv_smooth_mid','full_cv_lag_mo'});

writetable(sc_tbl, fullfile(OUT_TABLES,'mb_fixed74_smooth_char.csv'));
save(fullfile(DATA_INT,'mb_fixed74_smooth_char.mat'), 'sc_tbl');
fprintf('[smooth_char] Saved mb_fixed74_smooth_char.csv\n');

% =========================================================================
% SECTOR-LEVEL SMOOTH SUMMARY
% =========================================================================
fprintf('[smooth_char] Building sector smooth summary...\n');

groups = unique(string(fred_group(:)));
groups = sort(groups);
groups = groups(~ismember(groups, ["Stock Market","Unassigned",""]));

win_labels = {'1959:03-1999:12','1959:03-2007:12','1959:03-2026:01'};
fn_uv = {'hist_uv','p2007_uv','full_uv'};
fn_cv = {'hist_cv','p2007_cv','full_cv'};
inc_masks = {include_main_hist, include_main_2007, include_main_full_ext};
nb_uv_arr = {uv_mb_nbreaks_hist, uv2007_mb_nbreaks, uvfull_mb_nbreaks};
nb_cv_arr = {cv_mb_nbreaks_hist, cv2007_mb_nbreaks, cvfull_mb_nbreaks};

row_win  = strings(0,1); row_grp = strings(0,1); row_meas = strings(0,1);
row_n    = []; row_nbr = []; row_pct_ab = []; row_ser_all = []; row_ser_any = [];
row_med_g = []; row_mean_lag = [];

for wi = 1:3
    inc_m = inc_masks{wi};
    for gi = 1:numel(groups)
        grp   = groups(gi);
        gmask = string(fred_group(:)) == grp & inc_m(:);
        idx_g = find(gmask);

        for mi = 1:2
            if mi == 1
                meas = "UV"; fn = fn_uv{wi}; nb_v = nb_uv_arr{wi};
            else
                meas = "CV"; fn = fn_cv{wi}; nb_v = nb_cv_arr{wi};
            end

            % Restrict to indices within the range of saved arrays
            % (fred_group in processed_data.mat may span all raw FRED-MD series
            % while nb_v and gamma_out were initialized with n_series)
            valid_idx_g = idx_g(idx_g <= n_series);

            all_g = []; all_ab = []; all_lag = [];
            n_with_breaks = 0;
            n_all_abrupt  = 0;
            n_any_abrupt  = 0;

            for jj = 1:numel(valid_idx_g)
                j = valid_idx_g(jj);
                if isnan(nb_v(j)) || nb_v(j) < 1; continue; end
                n_with_breaks = n_with_breaks + 1;

                gstr = gamma_out.(fn)(j);
                astr = abrupt_out.(fn)(j);
                lstr = lag_out.(fn)(j);
                if gstr == ""; continue; end

                g_parts = strtrim(strsplit(gstr,'|'));
                a_parts = strtrim(strsplit(astr,'|'));
                l_parts = strtrim(strsplit(lstr,'|'));

                g_vals = cellfun(@str2double, cellstr(g_parts));
                a_vals = cellfun(@str2double, cellstr(a_parts));
                l_vals = cellfun(@str2double, cellstr(l_parts));

                g_vals = g_vals(~isnan(g_vals));
                a_vals = a_vals(~isnan(a_vals));
                l_vals = l_vals(~isnan(l_vals));

                all_g   = [all_g;   g_vals(:)]; %#ok<AGROW>
                all_ab  = [all_ab;  a_vals(:)]; %#ok<AGROW>
                all_lag = [all_lag; l_vals(:)]; %#ok<AGROW>

                if ~isempty(a_vals)
                    if all(a_vals == 1); n_all_abrupt = n_all_abrupt + 1; end
                    if any(a_vals == 1); n_any_abrupt = n_any_abrupt + 1; end
                end
            end

            row_win(end+1,1)    = win_labels{wi};
            row_grp(end+1,1)    = grp;
            row_meas(end+1,1)   = meas;
            row_n(end+1,1)      = sum(gmask);
            row_nbr(end+1,1)    = n_with_breaks;
            row_pct_ab(end+1,1) = mean(all_ab,'omitnan');
            row_ser_all(end+1,1)= n_all_abrupt / max(n_with_breaks,1);
            row_ser_any(end+1,1)= n_any_abrupt  / max(n_with_breaks,1);
            row_med_g(end+1,1)  = median(all_g,'omitnan');
            row_mean_lag(end+1,1)= mean(all_lag,'omitnan');
        end
    end
end

sm_tbl = table(row_win, row_grp, row_meas, row_n, row_nbr, ...
    row_pct_ab, row_ser_all, row_ser_any, row_med_g, row_mean_lag, ...
    'VariableNames',{ ...
        'window','fred_group','measure','n_series','n_with_breaks', ...
        'pct_abrupt_all_breaks','pct_series_all_abrupt', ...
        'pct_series_any_abrupt','median_gamma','mean_lag_months'});

writetable(sm_tbl, fullfile(OUT_TABLES,'mb_fixed74_smooth_summary.csv'));
fprintf('[smooth_char] Saved mb_fixed74_smooth_summary.csv\n');

% =========================================================================
% PRINT SUMMARY
% =========================================================================
fprintf('\n[smooth_char] === SMOOTH CHARACTERIZATION SUMMARY ===\n');
fprintf('  Threshold: gamma >= %d = abrupt\n\n', ABRUPT_THRESH);

for wi = 1:3
    fprintf('  Window: %s\n', win_labels{wi});
    for meas = ["UV","CV"]
        rows = sm_tbl(sm_tbl.window == win_labels{wi} & sm_tbl.measure == meas,:);
        fprintf('    %s:\n', meas);
        for ri = 1:height(rows)
            if rows.n_with_breaks(ri) == 0; continue; end
            fprintf('      %-32s  med_gamma=%4.0f  pct_abrupt=%3.0f%%  any_abrupt=%3.0f%%\n', ...
                rows.fred_group{ri}, rows.median_gamma(ri), ...
                100*rows.pct_abrupt_all_breaks(ri), ...
                100*rows.pct_series_any_abrupt(ri));
        end
    end
    fprintf('\n');
end

fprintf('[smooth_char] Done.\n');
fprintf('  Files written:\n');
fprintf('    output/tables/mb_fixed74_smooth_char.csv   (per-series gamma/abrupt/lag per break)\n');
fprintf('    output/tables/mb_fixed74_smooth_summary.csv (sector summary)\n');
fprintf('    data_intermediate/mb_fixed74_smooth_char.mat\n');
