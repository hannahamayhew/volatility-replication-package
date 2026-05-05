%% run_uv_historical.m
% =========================================================================
% UNCONDITIONAL VOLATILITY — HISTORICAL WINDOW (1959:03 – 1999:12)
%
% Implements the one-break UV workflow following Sensier & van Dijk (2004):
%
%   Step 1.  Build z_t = sqrt(pi/2) * |y_t - mean(y)|  via build_uv_object()
%   Step 2.  Search candidate break dates with 15% trimming
%   Step 3.  uv_onebreak_search   →  SSR-minimising break date
%   Step 4.  uv_supw_test         →  SupW statistic with HAC covariance
%   Step 5.  hansen_supf_pval     →  Hansen (1997) asymptotic p-values
%
% NOTE: The reported one-break date is the SSR-minimising date from
% uv_onebreak_search.  The SupW statistic and its significance (from
% uv_supw_test + Hansen p-values) are the formal test of whether a break
% exists; the SupW-maximising date may differ from the SSR-minimising date.
%
% HAC lag rule: h = floor(4*(T/100)^(2/9))  [Newey-West type]
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_rep, dates_rep_t, series, include_main_hist
%
% OUTPUTS (saved to data_intermediate/uv_hist_results.mat)
%   series_meta_uv_hist  - table with UV break results
%   vol_obj              - T_rep x N scaled volatility matrix
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_uv_historical] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_rep','dates_rep_t','series','include_main_hist');

n_series = length(series);

%% -----------------------------------------------------------------------
% 2. BUILD UV OBJECT  [Sensier & van Dijk (2004)]
%
%   z_t = sqrt(pi/2) * |y_t - mean(y)|
%
%   Demeaning uses the full historical window mean (nanmean column-wise).
%   The sqrt(pi/2) factor is applied by build_uv_object().
%% -----------------------------------------------------------------------
vol_obj = build_uv_object(yt_rep);   % T_rep x N

%% -----------------------------------------------------------------------
% 3. PRE/POST 1984 DESCRIPTIVE SPLIT
%% -----------------------------------------------------------------------
pre84_idx  = dates_rep_t < datetime(1984,1,1);
post84_idx = dates_rep_t >= datetime(1984,1,1);

uv_mean_pre84     = NaN(n_series,1);
uv_mean_post84    = NaN(n_series,1);
uv_ratio_post_pre = NaN(n_series,1);
uv_pct_change     = NaN(n_series,1);
uv_nobs           = NaN(n_series,1);

for j = 1:n_series
    if include_main_hist(j)
        vj = vol_obj(:,j);
        uv_nobs(j)           = sum(~isnan(vj));
        uv_mean_pre84(j)     = nanmean(vj(pre84_idx));
        uv_mean_post84(j)    = nanmean(vj(post84_idx));
        uv_ratio_post_pre(j) = uv_mean_post84(j) / uv_mean_pre84(j);
        uv_pct_change(j)     = 100*(uv_mean_post84(j)-uv_mean_pre84(j)) / uv_mean_pre84(j);
    end
end

%% -----------------------------------------------------------------------
% 4. ONE-BREAK SSR SEARCH
%
%   Searches all candidate break dates in [trim_frac*T, (1-trim_frac)*T].
%   Returns the date that minimises the within-regime sum of squared
%   residuals.  This is the reported break date.
%% -----------------------------------------------------------------------
fprintf('[run_uv_historical] Running one-break SSR search (trim=%.0f%%)...\n', ...
    TRIM_FRAC*100);

uv_break_found      = false(n_series,1);
uv_break_nobs       = NaN(n_series,1);
uv_break_date_str   = strings(n_series,1);
uv_mu1              = NaN(n_series,1);
uv_mu2              = NaN(n_series,1);
uv_delta            = NaN(n_series,1);
uv_ratio_break      = NaN(n_series,1);
uv_pct_change_break = NaN(n_series,1);
uv_ssr_nobreak      = NaN(n_series,1);
uv_ssr_best         = NaN(n_series,1);
uv_ssr_gain_pct     = NaN(n_series,1);

for j = 1:n_series
    if include_main_hist(j)
        vj  = vol_obj(:,j);
        res = uv_onebreak_search(vj, dates_rep_t, TRIM_FRAC);

        if res.success
            uv_break_found(j)      = true;
            uv_break_nobs(j)       = res.nobs;
            uv_break_date_str(j)   = string(res.break_date, 'dd-MMM-yyyy');
            uv_mu1(j)              = res.mu1;
            uv_mu2(j)              = res.mu2;
            uv_delta(j)            = res.delta;
            uv_ratio_break(j)      = res.ratio_post_pre;
            uv_pct_change_break(j) = res.pct_change;
            uv_ssr_nobreak(j)      = res.ssr_nobreak;
            uv_ssr_best(j)         = res.ssr_best;
            if res.ssr_nobreak > 0
                uv_ssr_gain_pct(j) = 100*(res.ssr_nobreak-res.ssr_best)/res.ssr_nobreak;
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 5. SUP-W TEST WITH HAC COVARIANCE
%
%   For each candidate break date tau, estimates:
%       v_t = alpha + delta * I(t > tau) + u_t
%   using HAC (Newey-West) variance with h = floor(4*(T/100)^(2/9)) lags.
%   SupW = max over trimmed candidates of the Wald statistic for delta=0.
%
%   The SupW-maximising break date may differ from the SSR-minimising date
%   because the HAC correction weights the test differently across dates.
%% -----------------------------------------------------------------------
fprintf('[run_uv_historical] Running SupW test (trim=%.0f%%, HAC auto-lag)...\n', ...
    TRIM_FRAC*100);

uv_supw_found      = false(n_series,1);
uv_supw_stat       = NaN(n_series,1);
uv_supw_hac_lag    = NaN(n_series,1);
uv_supw_break_date = strings(n_series,1);
uv_supw_mu1        = NaN(n_series,1);
uv_supw_mu2        = NaN(n_series,1);
uv_supw_ratio      = NaN(n_series,1);
uv_supw_pct_change = NaN(n_series,1);

for j = 1:n_series
    if include_main_hist(j)
        vj  = vol_obj(:,j);
        res = uv_supw_test(vj, dates_rep_t, TRIM_FRAC, []);

        if res.success
            uv_supw_found(j)      = true;
            uv_supw_stat(j)       = res.supW;
            uv_supw_hac_lag(j)    = res.hac_lag;
            uv_supw_break_date(j) = string(res.break_date, 'dd-MMM-yyyy');
            uv_supw_mu1(j)        = res.mu1;
            uv_supw_mu2(j)        = res.mu2;
            uv_supw_ratio(j)      = res.ratio_post_pre;
            uv_supw_pct_change(j) = res.pct_change;
        end
    end
end

%% -----------------------------------------------------------------------
% 6. HANSEN (1997) P-VALUES
%
%   Asymptotic p-value for the SupW statistic under H0: no break.
%   Approximation valid for m=1 restriction (break in mean), pi0=0.15.
%   Significance at 5% level: p < 0.05.
%% -----------------------------------------------------------------------
uv_hansen_pval = NaN(n_series,1);
uv_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if uv_supw_found(j)
        uv_hansen_pval(j) = hansen_supf_pval(uv_supw_stat(j), 1, TRIM_FRAC);
        uv_hansen_sig5(j) = uv_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 7. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_uv_hist = table( ...
    string(series(:)), ...
    include_main_hist(:), ...
    uv_nobs, ...
    uv_mean_pre84, uv_mean_post84, uv_ratio_post_pre, uv_pct_change, ...
    uv_break_found, uv_break_date_str, ...
    uv_mu1, uv_mu2, uv_delta, uv_ratio_break, uv_pct_change_break, ...
    uv_ssr_nobreak, uv_ssr_best, uv_ssr_gain_pct, ...
    uv_supw_found, uv_supw_stat, uv_supw_hac_lag, uv_supw_break_date, ...
    uv_supw_mu1, uv_supw_mu2, uv_supw_ratio, uv_supw_pct_change, ...
    uv_hansen_pval, uv_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_hist', ...
        'uv_nobs_hist', ...
        'uv_mean_pre84','uv_mean_post84','uv_ratio_post_pre84','uv_pct_change_post_pre84', ...
        'uv_break_found_hist','uv_break_date_hist', ...
        'uv_mu1_hist','uv_mu2_hist','uv_delta_hist', ...
        'uv_ratio_break_hist','uv_pct_change_break_hist', ...
        'uv_ssr_nobreak_hist','uv_ssr_best_hist','uv_ssr_gain_pct_hist', ...
        'uv_supw_found_hist','uv_supw_stat_hist','uv_supw_hac_lag_hist', ...
        'uv_supw_break_date_hist', ...
        'uv_supw_mu1_hist','uv_supw_mu2_hist', ...
        'uv_supw_ratio_hist','uv_supw_pct_change_hist', ...
        'uv_hansen_pval_hist','uv_hansen_sig5_hist'});

%% -----------------------------------------------------------------------
% 8. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'uv_hist_results.mat'), ...
    'series_meta_uv_hist','vol_obj');

writetable(series_meta_uv_hist, fullfile(OUT_TABLES,'uv_hist_results.csv'));

%% -----------------------------------------------------------------------
% 9. SUMMARY PRINT
%% -----------------------------------------------------------------------
% Use (:) to ensure scalar result regardless of include_main_hist row/col shape.
n_main = sum(include_main_hist(:));
n_sig  = sum(uv_hansen_sig5(:) & include_main_hist(:));
fprintf('\n[run_uv_historical] SUMMARY\n');
fprintf('  Historical window: %s – %s\n', ...
    string(HIST_START,'yyyy-MM'), string(HIST_END,'yyyy-MM'));
fprintf('  Series in main sample:           %d\n', n_main);
fprintf('  Significant UV breaks (p<5%%):   %d of %d\n', n_sig, n_main);
fprintf('  Median SupW statistic:           %.3f\n', nanmedian(uv_supw_stat(include_main_hist)));
fprintf('  Median Hansen p-value:           %.4f\n', nanmedian(uv_hansen_pval(include_main_hist)));
fprintf('  Median post/pre ratio (SSR date):%.3f\n', nanmedian(uv_ratio_break(include_main_hist)));
fprintf('[run_uv_historical] Done.\n\n');
