%% run_stationarity_tests.m
% =========================================================================
% STATIONARITY CHECKS — ADF + KPSS
%
% Tests each series in the main historical sample for stationarity using:
%   ADF  (Augmented Dickey-Fuller): H0 = unit root. Reject → stationary.
%   KPSS (Kwiatkowski-Phillips-Schmidt-Shin): H0 = stationary. Reject → non-stationary.
%
% The two tests are complementary:
%   ADF reject  + KPSS not-reject  → OK           (both agree: stationary)
%   ADF reject  + KPSS reject      → ambiguous     (conflict; inspect individually)
%   ADF fail    + KPSS not-reject  → adf_only      (ADF weak; KPSS OK — borderline)
%   ADF fail    + KPSS reject      → nonstationary (both agree: unit root present)
%
% Nominal FRED-MD groups (Prices, Money and Credit, Interest and Exchange
% Rates, Stock Market) are highlighted because the McCracken-Ng codes vary
% more across these series (e.g. some rate series use first difference rather
% than log-difference).
%
% TESTS RUN ON:
%   Historical window (yt_rep, 1959:03–1999:12) — primary
%   Full sample (yt_ext_full, 1959:03–2026:01)  — supplementary
%
% MODEL SPECIFICATION:
%   ADF: constant, no trend ('AR' model). Lag selected by BIC over 0:PMAX.
%        adftest returns results for all lags; we pick the BIC-minimising lag.
%   KPSS: level stationarity (trend=false). Lag = Newey-West rule ceil(4*(n/100)^0.25).
%   Rationale: transformed FRED-MD series should be mean-stationary without a
%   deterministic trend after the McCracken-Ng codes are applied.
%
% PREREQUISITE: run_data_prep.m
%
% OUTPUTS (saved to output/tables/)
%   stationarity_hist_results.csv   — per-series, historical window
%   stationarity_full_results.csv   — per-series, full sample
%   stationarity_summary.csv        — counts by fred_group and status
% =========================================================================

run('config.m');

fprintf('\n[run_stationarity_tests] Starting ADF + KPSS stationarity checks...\n');

NOMINAL_GROUPS = ["Prices","Money and Credit","Interest and Exchange Rates","Stock Market"];

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_rep','dates_rep_t', ...
    'yt_ext_full','dates_ext_full', ...
    'series','tcode','fred_group', ...
    'include_main_hist','include_main_full_ext');

n_series = length(series);

% Initialise cross-window group-summary accumulators
sum_win  = strings(0,1); sum_grp = strings(0,1); sum_nom = false(0,1);
sum_n    = []; sum_ok = []; sum_am = []; sum_af = []; sum_ns = [];

%% -----------------------------------------------------------------------
% 2. RUN TESTS ON BOTH WINDOWS
%% -----------------------------------------------------------------------
windows = { ...
    'hist', yt_rep,      dates_rep_t,    include_main_hist; ...
    'full', yt_ext_full, dates_ext_full, include_main_full_ext};

for wi = 1:size(windows,1)
    win_lbl  = windows{wi,1};
    yt_win   = windows{wi,2};
    inc_mask = logical(windows{wi,4}(:));

    fprintf('\n[run_stationarity_tests] Window: %s  (T=%d, n_include=%d)\n', ...
        upper(win_lbl), size(yt_win,1), sum(inc_mask));

    % Pre-allocate output vectors
    out_series  = strings(n_series,1);
    out_group   = strings(n_series,1);
    out_nominal = false(n_series,1);
    out_tcode   = NaN(n_series,1);
    out_n       = NaN(n_series,1);
    out_adf_lag = NaN(n_series,1);
    out_adf_stat= NaN(n_series,1);
    out_adf_p   = NaN(n_series,1);
    out_adf_h   = false(n_series,1);
    out_kpss_lag= NaN(n_series,1);
    out_kpss_stat=NaN(n_series,1);
    out_kpss_h  = false(n_series,1);
    out_status  = strings(n_series,1);
    out_include = false(n_series,1);

    for j = 1:n_series
        out_series(j)  = string(series{j});
        out_group(j)   = string(fred_group(j));
        out_nominal(j) = ismember(string(fred_group(j)), NOMINAL_GROUPS);
        out_tcode(j)   = tcode(j);
        out_include(j) = inc_mask(j);

        if ~inc_mask(j); out_status(j) = "excluded"; continue; end

        % Strip NaNs — use all non-NaN observations (typically just leading NaNs)
        y = yt_win(:,j);
        y_clean = y(~isnan(y));
        n_obs = numel(y_clean);
        out_n(j) = n_obs;

        if n_obs < 20
            out_status(j) = "too_short";
            fprintf('  WARNING: %s — only %d valid obs, skipping.\n', series{j}, n_obs);
            continue;
        end

        % --- ADF: constant, BIC-selected lag over 0:PMAX ---
        max_lag = min(PMAX, floor((n_obs - 2) / 3));
        try
            [h_vec, p_vec, stat_vec, ~, reg_vec] = adftest(y_clean, ...
                'model','AR', 'lags', 0:max_lag);
            % Pick BIC-minimising lag
            bic_vals = [reg_vec.BIC];
            [~, best] = min(bic_vals);
            out_adf_h(j)    = h_vec(best);
            out_adf_p(j)    = p_vec(best);
            out_adf_stat(j) = stat_vec(best);
            % Lag = number of lagged-difference terms = length(coeff) - 2
            % (coeff = [intercept, rho, delta_1 ... delta_k])
            out_adf_lag(j)  = numel(reg_vec(best).coeff) - 2;
        catch ME
            fprintf('  WARNING: ADF failed for %s: %s\n', series{j}, ME.message);
            out_status(j) = "adf_error";
            continue;
        end

        % --- KPSS: level stationarity, Newey-West bandwidth ---
        kpss_lag = ceil(4 * (n_obs/100)^0.25);
        out_kpss_lag(j) = kpss_lag;
        try
            [h_kpss, ~, stat_kpss] = kpsstest(y_clean, 'trend', false, 'lags', kpss_lag);
            out_kpss_h(j)    = h_kpss;
            out_kpss_stat(j) = stat_kpss;
        catch ME
            fprintf('  WARNING: KPSS failed for %s: %s\n', series{j}, ME.message);
            h_kpss = NaN;
        end

        % --- Classify ---
        h_adf  = out_adf_h(j);
        h_kpss = out_kpss_h(j);
        if h_adf && ~h_kpss
            out_status(j) = "OK";
        elseif h_adf && h_kpss
            out_status(j) = "ambiguous";
        elseif ~h_adf && ~h_kpss
            out_status(j) = "adf_only";
        elseif ~h_adf && h_kpss
            out_status(j) = "nonstationary";
        else
            out_status(j) = "error";
        end
    end

    % --- Print overall counts ---
    inc = inc_mask;
    fprintf('\n[run_stationarity_tests] Results (%s window, n=%d):\n', win_lbl, sum(inc));
    for s = ["OK","ambiguous","adf_only","nonstationary","too_short","error","adf_error"]
        n_s = sum(out_status(inc) == s);
        if n_s > 0
            fprintf('  %-16s: %d\n', s, n_s);
        end
    end

    % --- Print flagged series ---
    flag_mask = inc & ~ismember(out_status, ["OK","excluded","too_short"]);
    if any(flag_mask)
        fprintf('\n  Flagged series (%s window):\n', win_lbl);
        fprintf('  %-20s %-35s %5s %5s %7s %6s  %s\n', ...
            'Series','Group','tcode','ADF_h','ADF_p','KPSS_h','Status');
        fprintf('  %s\n', repmat('-',1,90));
        idx_flag = find(flag_mask);
        % Sort: nonstationary → ambiguous → adf_only
        [~,ord] = sort(out_status(idx_flag));
        for ii = 1:numel(idx_flag)
            j = idx_flag(ord(ii));
            nom_tag = ''; if out_nominal(j); nom_tag = ' *'; end
            fprintf('  %-20s %-35s %5d %5d %7.4f %6d  %s%s\n', ...
                out_series(j), out_group(j), out_tcode(j), ...
                out_adf_h(j), out_adf_p(j), out_kpss_h(j), ...
                out_status(j), nom_tag);
        end
        fprintf('  (* = nominal group)\n');
    else
        fprintf('  No flagged series — all pass both tests.\n');
    end

    % --- Group-level summary ---
    fprintf('\n  Group-level summary (%s window):\n', win_lbl);
    fprintf('  %-40s %4s %4s %4s %5s %5s\n', 'Group','n','OK','ambig','adf_f','nons');
    groups_u = unique(out_group(inc));
    groups_u = sort(groups_u(groups_u ~= ""));
    for gi = 1:numel(groups_u)
        grp   = groups_u(gi);
        gmask = inc & out_group == grp;
        n_g   = sum(gmask);
        n_ok  = sum(out_status(gmask) == "OK");
        n_am  = sum(out_status(gmask) == "ambiguous");
        n_af  = sum(out_status(gmask) == "adf_only");
        n_ns  = sum(out_status(gmask) == "nonstationary");
        nom_tag = ''; if ismember(grp, NOMINAL_GROUPS); nom_tag = ' *'; end
        fprintf('  %-40s %4d %4d %4d %5d %5d%s\n', grp, n_g, n_ok, n_am, n_af, n_ns, nom_tag);

        % Accumulate cross-window group summary
        sum_win(end+1,1)  = win_lbl;
        sum_grp(end+1,1)  = grp;
        sum_nom(end+1,1)  = ismember(grp, NOMINAL_GROUPS);
        sum_n(end+1,1)    = n_g;
        sum_ok(end+1,1)   = n_ok;
        sum_am(end+1,1)   = n_am;
        sum_af(end+1,1)   = n_af;
        sum_ns(end+1,1)   = n_ns;
    end
    fprintf('  (* = nominal group)\n');

    % --- Save per-series CSV ---
    out_tbl = table(out_series, out_group, out_nominal, out_tcode, out_include, out_n, ...
        out_adf_lag, out_adf_stat, out_adf_p, out_adf_h, ...
        out_kpss_lag, out_kpss_stat, out_kpss_h, out_status, ...
        'VariableNames',{ ...
            'series_name','fred_group','is_nominal','tcode','include_main','n_obs', ...
            'adf_lag_bic','adf_stat','adf_pval','adf_reject_unit_root', ...
            'kpss_lag_nw','kpss_stat','kpss_reject_stationarity','status'});

    fname = sprintf('stationarity_%s_results.csv', win_lbl);
    writetable(out_tbl, fullfile(OUT_TABLES, fname));
    fprintf('\n[run_stationarity_tests] Saved %s\n', fname);
end

%% -----------------------------------------------------------------------
% 3. CROSS-WINDOW GROUP SUMMARY
%% -----------------------------------------------------------------------
grp_tbl = table(sum_win, sum_grp, sum_nom, sum_n, sum_ok, sum_am, sum_af, sum_ns, ...
    'VariableNames',{'window','fred_group','is_nominal','n_series', ...
        'n_OK','n_ambiguous','n_adf_only','n_nonstationary'});
writetable(grp_tbl, fullfile(OUT_TABLES,'stationarity_summary.csv'));
fprintf('[run_stationarity_tests] Saved stationarity_summary.csv\n');
fprintf('[run_stationarity_tests] Done.\n\n');
