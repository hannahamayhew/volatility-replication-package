%% run_uv_modern.m
% =========================================================================
% UNCONDITIONAL VOLATILITY — MODERN-ERA APPENDIX WINDOW (1999:01 – latest)
%
% Appendix robustness exercise: re-runs the UV one-break SupW workflow
% on the modern-era subsample only.  Motivated by the hypothesis that the
% full 1959–2026 window contains too many embedded regimes to cleanly
% detect post-1999 volatility changes.
%
% Implements the same Sensier & van Dijk (2004) workflow as run_uv_fullsample.m:
%
%   Step 1.  Build z_t = sqrt(pi/2) * |y_t - mean(y)|  via build_uv_object()
%            Mean is computed over the MODERN-ERA window (1999:01–latest).
%   Step 2.  Search candidate break dates with 15% trimming
%   Step 3.  uv_onebreak_search   →  SSR-minimising break date
%   Step 4.  uv_supw_test         →  SupW statistic with HAC covariance
%   Step 5.  hansen_supf_pval_general → Hansen (1997) p-values (m=1)
%
% WINDOW NOTE:
%   T ≈ 325 months.  The descriptive pre/post-1984 split used in the
%   historical and full-sample UV scripts is not applicable here (1984 is
%   outside the window).  No split is reported.
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_mod, dates_mod, series, include_main_mod
%
% OUTPUTS (saved to data_intermediate/uv_modern_results.mat)
%   series_meta_uv_mod  - table with modern-era UV one-break results
%   vol_obj_mod         - T_mod x N scaled volatility matrix
%                         (required by run_mb_modern.m)
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_uv_modern] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_mod','dates_mod','series','include_main_mod');

n_series = length(series);

fprintf('[run_uv_modern] Window: %s – %s  (T=%d)\n', ...
    string(dates_mod(1),'yyyy-MM'), string(dates_mod(end),'yyyy-MM'), length(dates_mod));

%% -----------------------------------------------------------------------
% 2. BUILD UV OBJECT  [Sensier & van Dijk (2004)]
%
%   z_t = sqrt(pi/2) * |y_t - mean(y)|
%   Mean is computed over the modern-era window only.
%% -----------------------------------------------------------------------
vol_obj_mod = build_uv_object(yt_mod);   % T_mod x N

%% -----------------------------------------------------------------------
% 3. ONE-BREAK SSR SEARCH
%% -----------------------------------------------------------------------
fprintf('[run_uv_modern] Running one-break SSR search (trim=%.0f%%)...\n', ...
    TRIM_FRAC*100);

uvmod_break_found      = false(n_series,1);
uvmod_break_nobs       = NaN(n_series,1);
uvmod_break_date       = strings(n_series,1);
uvmod_mu1              = NaN(n_series,1);
uvmod_mu2              = NaN(n_series,1);
uvmod_delta            = NaN(n_series,1);
uvmod_ratio_break      = NaN(n_series,1);
uvmod_pct_change_break = NaN(n_series,1);
uvmod_ssr_nobreak      = NaN(n_series,1);
uvmod_ssr_best         = NaN(n_series,1);
uvmod_ssr_gain_pct     = NaN(n_series,1);

for j = 1:n_series
    if include_main_mod(j)
        vj  = vol_obj_mod(:,j);
        res = uv_onebreak_search(vj, dates_mod, TRIM_FRAC);

        if res.success
            uvmod_break_found(j)      = true;
            uvmod_break_nobs(j)       = res.nobs;
            uvmod_break_date(j)       = string(res.break_date,'dd-MMM-yyyy');
            uvmod_mu1(j)              = res.mu1;
            uvmod_mu2(j)              = res.mu2;
            uvmod_delta(j)            = res.delta;
            uvmod_ratio_break(j)      = res.ratio_post_pre;
            uvmod_pct_change_break(j) = res.pct_change;
            uvmod_ssr_nobreak(j)      = res.ssr_nobreak;
            uvmod_ssr_best(j)         = res.ssr_best;
            if res.ssr_nobreak > 0
                uvmod_ssr_gain_pct(j) = 100*(res.ssr_nobreak-res.ssr_best)/res.ssr_nobreak;
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 4. SUP-W TEST WITH HAC COVARIANCE
%% -----------------------------------------------------------------------
fprintf('[run_uv_modern] Running SupW test (trim=%.0f%%, HAC auto-lag)...\n', ...
    TRIM_FRAC*100);

uvmod_supw_found      = false(n_series,1);
uvmod_supw_stat       = NaN(n_series,1);
uvmod_supw_hac_lag    = NaN(n_series,1);
uvmod_supw_break_date = strings(n_series,1);
uvmod_supw_mu1        = NaN(n_series,1);
uvmod_supw_mu2        = NaN(n_series,1);
uvmod_supw_ratio      = NaN(n_series,1);
uvmod_supw_pct_change = NaN(n_series,1);

for j = 1:n_series
    if include_main_mod(j)
        vj  = vol_obj_mod(:,j);
        res = uv_supw_test(vj, dates_mod, TRIM_FRAC, []);

        if res.success
            uvmod_supw_found(j)      = true;
            uvmod_supw_stat(j)       = res.supW;
            uvmod_supw_hac_lag(j)    = res.hac_lag;
            uvmod_supw_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
            uvmod_supw_mu1(j)        = res.mu1;
            uvmod_supw_mu2(j)        = res.mu2;
            uvmod_supw_ratio(j)      = res.ratio_post_pre;
            uvmod_supw_pct_change(j) = res.pct_change;
        end
    end
end

%% -----------------------------------------------------------------------
% 5. HANSEN (1997) P-VALUES
%% -----------------------------------------------------------------------
uvmod_hansen_pval = NaN(n_series,1);
uvmod_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if uvmod_supw_found(j)
        uvmod_hansen_pval(j) = hansen_supf_pval_general(uvmod_supw_stat(j), 1, TRIM_FRAC);
        uvmod_hansen_sig5(j) = uvmod_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 6. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_uv_mod = table( ...
    string(series(:)), include_main_mod(:), ...
    uvmod_break_found, uvmod_break_date, ...
    uvmod_mu1, uvmod_mu2, uvmod_delta, ...
    uvmod_ratio_break, uvmod_pct_change_break, ...
    uvmod_ssr_nobreak, uvmod_ssr_best, uvmod_ssr_gain_pct, ...
    uvmod_supw_found, uvmod_supw_stat, uvmod_supw_hac_lag, ...
    uvmod_supw_break_date, ...
    uvmod_supw_mu1, uvmod_supw_mu2, ...
    uvmod_supw_ratio, uvmod_supw_pct_change, ...
    uvmod_hansen_pval, uvmod_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_mod', ...
        'uvmod_break_found','uvmod_break_date', ...
        'uvmod_mu1','uvmod_mu2','uvmod_delta', ...
        'uvmod_ratio_break','uvmod_pct_change_break', ...
        'uvmod_ssr_nobreak','uvmod_ssr_best','uvmod_ssr_gain_pct', ...
        'uvmod_supw_found','uvmod_supw_stat','uvmod_supw_hac_lag', ...
        'uvmod_supw_break_date', ...
        'uvmod_supw_mu1','uvmod_supw_mu2', ...
        'uvmod_supw_ratio','uvmod_supw_pct_change', ...
        'uvmod_hansen_pval','uvmod_hansen_sig5'});

%% -----------------------------------------------------------------------
% 7. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'uv_modern_results.mat'), ...
    'series_meta_uv_mod','vol_obj_mod');

writetable(series_meta_uv_mod, fullfile(OUT_TABLES,'uv_modern_results.csv'));

%% -----------------------------------------------------------------------
% 8. SUMMARY PRINT
%% -----------------------------------------------------------------------
n_main = sum(include_main_mod(:));
n_sig  = sum(uvmod_hansen_sig5(:) & include_main_mod(:));
fprintf('\n[run_uv_modern] SUMMARY\n');
fprintf('  Window: %s – %s  (T=%d)\n', ...
    string(dates_mod(1),'yyyy-MM'), string(dates_mod(end),'yyyy-MM'), length(dates_mod));
fprintf('  Series in main sample:          %d\n', n_main);
fprintf('  Significant UV breaks (p<5%%):  %d of %d\n', n_sig, n_main);
fprintf('  Median SupW statistic:          %.3f\n', nanmedian(uvmod_supw_stat(include_main_mod)));
fprintf('  Median Hansen p-value:          %.4f\n', nanmedian(uvmod_hansen_pval(include_main_mod)));
fprintf('[run_uv_modern] Done.\n\n');
