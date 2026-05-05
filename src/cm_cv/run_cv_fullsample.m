%% run_cv_fullsample.m
% =========================================================================
% CONDITIONAL VOLATILITY BREAKS — FULL SAMPLE (1959:03 – latest vintage)
%
% Uses residuals from the AR(p) conditional-mean SupW model (full sample)
% to construct the Sensier-style conditional volatility proxy:
%
%   s_t = sqrt(pi/2) * |e_hat_t|
%
% Then tests for a break in the mean of s_t using the same one-break
% SupW + Hansen workflow used for the historical CV.
%
% PREREQUISITE: run_cm_fullsample.m must have run first (provides residuals).
%
% INPUTS  (loaded from data_intermediate)
%   cm_full_results.mat   → cm_supw_resid_break_full, series_meta_cm_full
%   processed_data.mat    → dates_ext_full, series, include_main_full_ext
%
% OUTPUTS (saved to data_intermediate/cv_full_results.mat)
%   series_meta_cv_full  - table with CV full-sample break results
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_fullsample] Loading CM residuals and processed data...\n');
load(fullfile(DATA_INT,'cm_full_results.mat'), ...
    'cm_supw_resid_break_full','series_meta_cm_full');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_full','series','include_main_full_ext');

n_series  = length(series);
cmfull_p  = series_meta_cm_full.cmfull_p;

%% -----------------------------------------------------------------------
% 2. CV ONE-BREAK SUP-W TEST
%% -----------------------------------------------------------------------
fprintf('[run_cv_fullsample] Running CV SupW tests (trim=%.0f%%)...\n', TRIM_FRAC*100);

cvfull_found      = false(n_series,1);
cvfull_supw_stat  = NaN(n_series,1);
cvfull_break_date = strings(n_series,1);
cvfull_hac_lag    = NaN(n_series,1);
cvfull_mu1        = NaN(n_series,1);
cvfull_mu2        = NaN(n_series,1);
cvfull_ratio      = NaN(n_series,1);
cvfull_pct_change = NaN(n_series,1);

for j = 1:n_series
    if include_main_full_ext(j) && series_meta_cm_full.cmfull_supw_found(j) && ...
            ~isempty(cm_supw_resid_break_full{j})

        pj          = cmfull_p(j);
        dates_reg_j = dates_ext_full((pj+1):end);
        ehat_j      = cm_supw_resid_break_full{j};

        Tj = length(ehat_j);
        if length(dates_reg_j) > Tj
            dates_reg_j = dates_reg_j(end-Tj+1:end);
        end

        res = cv_supw_test(ehat_j, dates_reg_j, TRIM_FRAC, []);

        if res.success
            cvfull_found(j)      = true;
            cvfull_supw_stat(j)  = res.supW;
            cvfull_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
            cvfull_hac_lag(j)    = res.hac_lag;
            cvfull_mu1(j)        = res.mu1;
            cvfull_mu2(j)        = res.mu2;
            cvfull_ratio(j)      = res.ratio_post_pre;
            cvfull_pct_change(j) = res.pct_change;
        end
    end
end

%% -----------------------------------------------------------------------
% 3. HANSEN (1997) P-VALUES  (m=1, pi0=TRIM_FRAC)
%% -----------------------------------------------------------------------
cvfull_hansen_pval = NaN(n_series,1);
cvfull_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if cvfull_found(j)
        cvfull_hansen_pval(j) = hansen_supf_pval(cvfull_supw_stat(j), 1, TRIM_FRAC);
        cvfull_hansen_sig5(j) = cvfull_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cv_full = table( ...
    string(series(:)), include_main_full_ext(:), ...
    cvfull_found, cvfull_supw_stat, cvfull_break_date, ...
    cvfull_hac_lag, ...
    cvfull_mu1, cvfull_mu2, cvfull_ratio, cvfull_pct_change, ...
    cvfull_hansen_pval, cvfull_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_full_ext', ...
        'cvfull_found','cvfull_supw_stat','cvfull_break_date', ...
        'cvfull_hac_lag', ...
        'cvfull_mu1','cvfull_mu2','cvfull_ratio','cvfull_pct_change', ...
        'cvfull_hansen_pval','cvfull_hansen_sig5'});

%% -----------------------------------------------------------------------
% 5. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'cv_full_results.mat'), 'series_meta_cv_full');
writetable(series_meta_cv_full, fullfile(OUT_TABLES,'cv_full_results.csv'));

fprintf('\n[run_cv_fullsample] SUMMARY\n');
fprintf('  Series with CM fit (main full):    %d\n', ...
    sum(include_main_full_ext(:) & series_meta_cm_full.cmfull_supw_found(:)));
fprintf('  Significant CV breaks (p<5%%):     %d\n', sum(cvfull_hansen_sig5));
fprintf('[run_cv_fullsample] Done.\n\n');
