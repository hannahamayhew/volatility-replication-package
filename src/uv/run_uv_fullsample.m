%% run_uv_fullsample.m
% =========================================================================
% UNCONDITIONAL VOLATILITY — FULL SAMPLE (1959:03 – latest vintage)
%
% Implements the one-break UV workflow following Sensier & van Dijk (2004):
%
%   Step 1.  Build z_t = sqrt(pi/2) * |y_t - mean(y)|  via build_uv_object()
%            Mean is computed over the full available sample.
%   Step 2.  Search candidate break dates with 15% trimming
%   Step 3.  uv_onebreak_search   →  SSR-minimising break date
%   Step 4.  uv_supw_test         →  SupW statistic with HAC covariance
%   Step 5.  hansen_supf_pval_general → Hansen (1997) p-values (m=1)
%
% HAC lag rule: h = floor(4*(T/100)^(2/9))  [Newey-West type]
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_ext_full, dates_ext_full, series, include_main_full_ext
%
% OUTPUTS (saved to data_intermediate/uv_full_results.mat)
%   series_meta_uv_full  - table with full-sample UV break results
%   vol_obj_full         - T_full x N scaled volatility matrix
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_uv_fullsample] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_ext_full','dates_ext_full','series','include_main_full_ext');

n_series = length(series);

%% -----------------------------------------------------------------------
% 2. BUILD UV OBJECT  [Sensier & van Dijk (2004)]
%
%   z_t = sqrt(pi/2) * |y_t - mean(y)|
%   Mean is computed over the full available sample.
%% -----------------------------------------------------------------------
vol_obj_full = build_uv_object(yt_ext_full);   % T_full x N

%% -----------------------------------------------------------------------
% 3. PRE/POST 1984 DESCRIPTIVE SPLIT
%% -----------------------------------------------------------------------
pre84_idx_full  = dates_ext_full < datetime(1984,1,1);
post84_idx_full = dates_ext_full >= datetime(1984,1,1);

uvfull_nobs             = NaN(n_series,1);
uvfull_mean_pre84       = NaN(n_series,1);
uvfull_mean_post84      = NaN(n_series,1);
uvfull_ratio_post_pre84 = NaN(n_series,1);
uvfull_pct_change_post84 = NaN(n_series,1);

for j = 1:n_series
    if include_main_full_ext(j)
        vj = vol_obj_full(:,j);
        uvfull_nobs(j)              = sum(~isnan(vj));
        uvfull_mean_pre84(j)        = nanmean(vj(pre84_idx_full));
        uvfull_mean_post84(j)       = nanmean(vj(post84_idx_full));
        uvfull_ratio_post_pre84(j)  = uvfull_mean_post84(j) / uvfull_mean_pre84(j);
        uvfull_pct_change_post84(j) = 100*(uvfull_mean_post84(j)-uvfull_mean_pre84(j)) / uvfull_mean_pre84(j);
    end
end

%% -----------------------------------------------------------------------
% 4. ONE-BREAK SSR SEARCH
%% -----------------------------------------------------------------------
fprintf('[run_uv_fullsample] Running one-break SSR search (trim=%.0f%%)...\n', TRIM_FRAC*100);

uvfull_break_found      = false(n_series,1);
uvfull_break_nobs       = NaN(n_series,1);
uvfull_break_date       = strings(n_series,1);
uvfull_mu1              = NaN(n_series,1);
uvfull_mu2              = NaN(n_series,1);
uvfull_delta            = NaN(n_series,1);
uvfull_ratio_break      = NaN(n_series,1);
uvfull_pct_change_break = NaN(n_series,1);
uvfull_ssr_nobreak      = NaN(n_series,1);
uvfull_ssr_best         = NaN(n_series,1);
uvfull_ssr_gain_pct     = NaN(n_series,1);

for j = 1:n_series
    if include_main_full_ext(j)
        vj  = vol_obj_full(:,j);
        res = uv_onebreak_search(vj, dates_ext_full, TRIM_FRAC);

        if res.success
            uvfull_break_found(j)      = true;
            uvfull_break_nobs(j)       = res.nobs;
            uvfull_break_date(j)       = string(res.break_date,'dd-MMM-yyyy');
            uvfull_mu1(j)              = res.mu1;
            uvfull_mu2(j)              = res.mu2;
            uvfull_delta(j)            = res.delta;
            uvfull_ratio_break(j)      = res.ratio_post_pre;
            uvfull_pct_change_break(j) = res.pct_change;
            uvfull_ssr_nobreak(j)      = res.ssr_nobreak;
            uvfull_ssr_best(j)         = res.ssr_best;
            if res.ssr_nobreak > 0
                uvfull_ssr_gain_pct(j) = 100*(res.ssr_nobreak-res.ssr_best)/res.ssr_nobreak;
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 5. SUP-W TEST WITH HAC COVARIANCE
%% -----------------------------------------------------------------------
fprintf('[run_uv_fullsample] Running SupW test (trim=%.0f%%, HAC auto-lag)...\n', TRIM_FRAC*100);

uvfull_supw_found      = false(n_series,1);
uvfull_supw_stat       = NaN(n_series,1);
uvfull_supw_hac_lag    = NaN(n_series,1);
uvfull_supw_break_date = strings(n_series,1);
uvfull_supw_mu1        = NaN(n_series,1);
uvfull_supw_mu2        = NaN(n_series,1);
uvfull_supw_ratio      = NaN(n_series,1);
uvfull_supw_pct_change = NaN(n_series,1);

for j = 1:n_series
    if include_main_full_ext(j)
        vj  = vol_obj_full(:,j);
        res = uv_supw_test(vj, dates_ext_full, TRIM_FRAC, []);

        if res.success
            uvfull_supw_found(j)      = true;
            uvfull_supw_stat(j)       = res.supW;
            uvfull_supw_hac_lag(j)    = res.hac_lag;
            uvfull_supw_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
            uvfull_supw_mu1(j)        = res.mu1;
            uvfull_supw_mu2(j)        = res.mu2;
            uvfull_supw_ratio(j)      = res.ratio_post_pre;
            uvfull_supw_pct_change(j) = res.pct_change;
        end
    end
end

%% -----------------------------------------------------------------------
% 6. HANSEN (1997) P-VALUES
%% -----------------------------------------------------------------------
uvfull_hansen_pval = NaN(n_series,1);
uvfull_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if uvfull_supw_found(j)
        uvfull_hansen_pval(j) = hansen_supf_pval_general(uvfull_supw_stat(j), 1, TRIM_FRAC);
        uvfull_hansen_sig5(j) = uvfull_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 7. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_uv_full = table( ...
    string(series(:)), include_main_full_ext(:), ...
    uvfull_nobs, ...
    uvfull_mean_pre84, uvfull_mean_post84, ...
    uvfull_ratio_post_pre84, uvfull_pct_change_post84, ...
    uvfull_break_found, uvfull_break_date, ...
    uvfull_mu1, uvfull_mu2, uvfull_delta, ...
    uvfull_ratio_break, uvfull_pct_change_break, ...
    uvfull_ssr_nobreak, uvfull_ssr_best, uvfull_ssr_gain_pct, ...
    uvfull_supw_found, uvfull_supw_stat, uvfull_supw_hac_lag, ...
    uvfull_supw_break_date, ...
    uvfull_supw_mu1, uvfull_supw_mu2, ...
    uvfull_supw_ratio, uvfull_supw_pct_change, ...
    uvfull_hansen_pval, uvfull_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_full_ext', ...
        'uvfull_nobs', ...
        'uvfull_mean_pre84','uvfull_mean_post84', ...
        'uvfull_ratio_post_pre84','uvfull_pct_change_post84', ...
        'uvfull_break_found','uvfull_break_date', ...
        'uvfull_mu1','uvfull_mu2','uvfull_delta', ...
        'uvfull_ratio_break','uvfull_pct_change_break', ...
        'uvfull_ssr_nobreak','uvfull_ssr_best','uvfull_ssr_gain_pct', ...
        'uvfull_supw_found','uvfull_supw_stat','uvfull_supw_hac_lag', ...
        'uvfull_supw_break_date', ...
        'uvfull_supw_mu1','uvfull_supw_mu2', ...
        'uvfull_supw_ratio','uvfull_supw_pct_change', ...
        'uvfull_hansen_pval','uvfull_hansen_sig5'});

%% -----------------------------------------------------------------------
% 8. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'uv_full_results.mat'), ...
    'series_meta_uv_full','vol_obj_full');

writetable(series_meta_uv_full, fullfile(OUT_TABLES,'uv_full_results.csv'));

%% -----------------------------------------------------------------------
% 9. SUMMARY PRINT
%% -----------------------------------------------------------------------
n_main = sum(include_main_full_ext(:));
n_sig  = sum(uvfull_hansen_sig5(:) & include_main_full_ext(:));
fprintf('\n[run_uv_fullsample] SUMMARY\n');
fprintf('  Window: %s – %s\n', ...
    string(HIST_START,'yyyy-MM'), string(dates_ext_full(end),'yyyy-MM'));
fprintf('  Series in main sample:          %d\n', n_main);
fprintf('  Significant UV breaks (p<5%%):  %d of %d\n', n_sig, n_main);
fprintf('  Median SupW statistic:          %.3f\n', nanmedian(uvfull_supw_stat(include_main_full_ext)));
fprintf('  Median Hansen p-value:          %.4f\n', nanmedian(uvfull_hansen_pval(include_main_full_ext)));
fprintf('[run_uv_fullsample] Done.\n\n');
