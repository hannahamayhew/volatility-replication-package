%% run_uv_2007.m
% =========================================================================
% UNCONDITIONAL VOLATILITY — 2007 WINDOW (1959:03 – 2007:12)
%
% Implements the one-break UV workflow following Sensier & van Dijk (2004):
%
%   Step 1.  Build z_t = sqrt(pi/2) * |y_t - mean(y)|  via build_uv_object()
%            Mean is computed over the 2007 window (1959:03–2007:12).
%   Step 2.  Search candidate break dates with 15% trimming
%   Step 3.  uv_onebreak_search   →  SSR-minimising break date
%   Step 4.  uv_supw_test         →  SupW statistic with HAC covariance
%   Step 5.  hansen_supf_pval_general → Hansen (1997) p-values (m=1)
%
% HAC lag rule: h = floor(4*(T/100)^(2/9))  [Newey-West type]
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_ext_2007, dates_ext_2007, series, include_main_2007
%
% OUTPUTS (saved to data_intermediate/uv_2007_results.mat)
%   series_meta_uv_2007  - table with UV 2007-window break results
%   vol_obj_2007         - T_2007 x N scaled volatility matrix
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_uv_2007] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_ext_2007','dates_ext_2007','series','include_main_2007');

n_series = length(series);

%% -----------------------------------------------------------------------
% 2. BUILD UV OBJECT  [Sensier & van Dijk (2004)]
%
%   z_t = sqrt(pi/2) * |y_t - mean(y)|
%   Mean is computed over the 2007 window (1959:03–2007:12).
%% -----------------------------------------------------------------------
vol_obj_2007 = build_uv_object(yt_ext_2007);   % T_2007 x N

%% -----------------------------------------------------------------------
% 3. PRE/POST 1984 DESCRIPTIVE SPLIT
%% -----------------------------------------------------------------------
pre84_idx_2007  = dates_ext_2007 < datetime(1984,1,1);
post84_idx_2007 = dates_ext_2007 >= datetime(1984,1,1);

uv2007_nobs             = NaN(n_series,1);
uv2007_mean_pre84       = NaN(n_series,1);
uv2007_mean_post84      = NaN(n_series,1);
uv2007_ratio_post_pre84 = NaN(n_series,1);
uv2007_pct_change_post84 = NaN(n_series,1);

for j = 1:n_series
    if include_main_2007(j)
        vj = vol_obj_2007(:,j);
        uv2007_nobs(j)              = sum(~isnan(vj));
        uv2007_mean_pre84(j)        = nanmean(vj(pre84_idx_2007));
        uv2007_mean_post84(j)       = nanmean(vj(post84_idx_2007));
        uv2007_ratio_post_pre84(j)  = uv2007_mean_post84(j) / uv2007_mean_pre84(j);
        uv2007_pct_change_post84(j) = 100*(uv2007_mean_post84(j)-uv2007_mean_pre84(j)) / uv2007_mean_pre84(j);
    end
end

%% -----------------------------------------------------------------------
% 4. ONE-BREAK SSR SEARCH
%% -----------------------------------------------------------------------
fprintf('[run_uv_2007] Running one-break SSR search (trim=%.0f%%)...\n', TRIM_FRAC*100);

uv2007_break_found      = false(n_series,1);
uv2007_break_nobs       = NaN(n_series,1);
uv2007_break_date       = strings(n_series,1);
uv2007_mu1              = NaN(n_series,1);
uv2007_mu2              = NaN(n_series,1);
uv2007_delta            = NaN(n_series,1);
uv2007_ratio_break      = NaN(n_series,1);
uv2007_pct_change_break = NaN(n_series,1);
uv2007_ssr_nobreak      = NaN(n_series,1);
uv2007_ssr_best         = NaN(n_series,1);
uv2007_ssr_gain_pct     = NaN(n_series,1);

for j = 1:n_series
    if include_main_2007(j)
        vj  = vol_obj_2007(:,j);
        res = uv_onebreak_search(vj, dates_ext_2007, TRIM_FRAC);

        if res.success
            uv2007_break_found(j)      = true;
            uv2007_break_nobs(j)       = res.nobs;
            uv2007_break_date(j)       = string(res.break_date,'dd-MMM-yyyy');
            uv2007_mu1(j)              = res.mu1;
            uv2007_mu2(j)              = res.mu2;
            uv2007_delta(j)            = res.delta;
            uv2007_ratio_break(j)      = res.ratio_post_pre;
            uv2007_pct_change_break(j) = res.pct_change;
            uv2007_ssr_nobreak(j)      = res.ssr_nobreak;
            uv2007_ssr_best(j)         = res.ssr_best;
            if res.ssr_nobreak > 0
                uv2007_ssr_gain_pct(j) = 100*(res.ssr_nobreak-res.ssr_best)/res.ssr_nobreak;
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 5. SUP-W TEST WITH HAC COVARIANCE
%% -----------------------------------------------------------------------
fprintf('[run_uv_2007] Running SupW test (trim=%.0f%%, HAC auto-lag)...\n', TRIM_FRAC*100);

uv2007_supw_found      = false(n_series,1);
uv2007_supw_stat       = NaN(n_series,1);
uv2007_supw_hac_lag    = NaN(n_series,1);
uv2007_supw_break_date = strings(n_series,1);
uv2007_supw_mu1        = NaN(n_series,1);
uv2007_supw_mu2        = NaN(n_series,1);
uv2007_supw_ratio      = NaN(n_series,1);
uv2007_supw_pct_change = NaN(n_series,1);

for j = 1:n_series
    if include_main_2007(j)
        vj  = vol_obj_2007(:,j);
        res = uv_supw_test(vj, dates_ext_2007, TRIM_FRAC, []);

        if res.success
            uv2007_supw_found(j)      = true;
            uv2007_supw_stat(j)       = res.supW;
            uv2007_supw_hac_lag(j)    = res.hac_lag;
            uv2007_supw_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
            uv2007_supw_mu1(j)        = res.mu1;
            uv2007_supw_mu2(j)        = res.mu2;
            uv2007_supw_ratio(j)      = res.ratio_post_pre;
            uv2007_supw_pct_change(j) = res.pct_change;
        end
    end
end

%% -----------------------------------------------------------------------
% 6. HANSEN (1997) P-VALUES
%% -----------------------------------------------------------------------
uv2007_hansen_pval = NaN(n_series,1);
uv2007_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if uv2007_supw_found(j)
        uv2007_hansen_pval(j) = hansen_supf_pval_general(uv2007_supw_stat(j), 1, TRIM_FRAC);
        uv2007_hansen_sig5(j) = uv2007_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 7. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_uv_2007 = table( ...
    string(series(:)), include_main_2007(:), ...
    uv2007_nobs, ...
    uv2007_mean_pre84, uv2007_mean_post84, ...
    uv2007_ratio_post_pre84, uv2007_pct_change_post84, ...
    uv2007_break_found, uv2007_break_date, ...
    uv2007_mu1, uv2007_mu2, uv2007_delta, ...
    uv2007_ratio_break, uv2007_pct_change_break, ...
    uv2007_ssr_nobreak, uv2007_ssr_best, uv2007_ssr_gain_pct, ...
    uv2007_supw_found, uv2007_supw_stat, uv2007_supw_hac_lag, ...
    uv2007_supw_break_date, ...
    uv2007_supw_mu1, uv2007_supw_mu2, ...
    uv2007_supw_ratio, uv2007_supw_pct_change, ...
    uv2007_hansen_pval, uv2007_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_2007', ...
        'uv2007_nobs', ...
        'uv2007_mean_pre84','uv2007_mean_post84', ...
        'uv2007_ratio_post_pre84','uv2007_pct_change_post84', ...
        'uv2007_break_found','uv2007_break_date', ...
        'uv2007_mu1','uv2007_mu2','uv2007_delta', ...
        'uv2007_ratio_break','uv2007_pct_change_break', ...
        'uv2007_ssr_nobreak','uv2007_ssr_best','uv2007_ssr_gain_pct', ...
        'uv2007_supw_found','uv2007_supw_stat','uv2007_supw_hac_lag', ...
        'uv2007_supw_break_date', ...
        'uv2007_supw_mu1','uv2007_supw_mu2', ...
        'uv2007_supw_ratio','uv2007_supw_pct_change', ...
        'uv2007_hansen_pval','uv2007_hansen_sig5'});

%% -----------------------------------------------------------------------
% 8. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'uv_2007_results.mat'), ...
    'series_meta_uv_2007','vol_obj_2007');

writetable(series_meta_uv_2007, fullfile(OUT_TABLES,'uv_2007_results.csv'));

%% -----------------------------------------------------------------------
% 9. SUMMARY PRINT
%% -----------------------------------------------------------------------
n_main = sum(include_main_2007(:));
n_sig  = sum(uv2007_hansen_sig5(:) & include_main_2007(:));
fprintf('\n[run_uv_2007] SUMMARY\n');
fprintf('  Window: %s – %s\n', ...
    string(HIST_START,'yyyy-MM'), string(WIN_2007_END,'yyyy-MM'));
fprintf('  Series in main sample:          %d\n', n_main);
fprintf('  Significant UV breaks (p<5%%):  %d of %d\n', n_sig, n_main);
fprintf('  Median SupW statistic:          %.3f\n', nanmedian(uv2007_supw_stat(include_main_2007)));
fprintf('  Median Hansen p-value:          %.4f\n', nanmedian(uv2007_hansen_pval(include_main_2007)));
fprintf('[run_uv_2007] Done.\n\n');
