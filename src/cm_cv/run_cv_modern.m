%% run_cv_modern.m
% =========================================================================
% CONDITIONAL VOLATILITY BREAKS — MODERN-ERA APPENDIX WINDOW (1999:01 – latest)
%
% Uses residuals from the AR(p) conditional-mean SupW model (modern-era window)
% to construct the Sensier-style conditional volatility proxy:
%
%   s_t = sqrt(pi/2) * |e_hat_t|
%
% Then tests for a break in the mean of s_t using the same one-break
% SupW + Hansen workflow used for other windows.
%
% PREREQUISITE: run_cm_modern.m must have run first (provides residuals).
%
% INPUTS  (loaded from data_intermediate)
%   cm_modern_results.mat  → cm_supw_resid_break_mod, series_meta_cm_mod
%   processed_data.mat     → dates_mod, series, include_main_mod
%
% OUTPUTS (saved to data_intermediate/cv_modern_results.mat)
%   series_meta_cv_mod  - table with CV modern-era one-break results
%
% NOTE: This script tests for a SINGLE break in CV.  For the recursive
% multiple-break analysis, see run_mb_modern.m, which uses the CM residuals
% directly from cm_modern_results.mat.
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_modern] Loading CM residuals and processed data...\n');
load(fullfile(DATA_INT,'cm_modern_results.mat'), ...
    'cm_supw_resid_break_mod','series_meta_cm_mod');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_mod','series','include_main_mod');

n_series = length(series);
cmmod_p  = series_meta_cm_mod.cmmod_p;

fprintf('[run_cv_modern] Window: %s – %s\n', ...
    string(dates_mod(1),'yyyy-MM'), string(dates_mod(end),'yyyy-MM'));

%% -----------------------------------------------------------------------
% 2. CV ONE-BREAK SUP-W TEST
%% -----------------------------------------------------------------------
fprintf('[run_cv_modern] Running CV SupW tests (trim=%.0f%%)...\n', TRIM_FRAC*100);

cvmod_found      = false(n_series,1);
cvmod_supw_stat  = NaN(n_series,1);
cvmod_break_date = strings(n_series,1);
cvmod_hac_lag    = NaN(n_series,1);
cvmod_mu1        = NaN(n_series,1);
cvmod_mu2        = NaN(n_series,1);
cvmod_ratio      = NaN(n_series,1);
cvmod_pct_change = NaN(n_series,1);

for j = 1:n_series
    if include_main_mod(j) && series_meta_cm_mod.cmmod_supw_found(j) && ...
            ~isempty(cm_supw_resid_break_mod{j})

        pj          = cmmod_p(j);
        ehat_j      = cm_supw_resid_break_mod{j};
        Tj          = length(ehat_j);
        dates_reg_j = dates_mod((pj+1):end);
        if length(dates_reg_j) > Tj
            dates_reg_j = dates_reg_j(end-Tj+1:end);
        end

        res = cv_supw_test(ehat_j, dates_reg_j, TRIM_FRAC, []);

        if res.success
            cvmod_found(j)      = true;
            cvmod_supw_stat(j)  = res.supW;
            cvmod_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
            cvmod_hac_lag(j)    = res.hac_lag;
            cvmod_mu1(j)        = res.mu1;
            cvmod_mu2(j)        = res.mu2;
            cvmod_ratio(j)      = res.ratio_post_pre;
            cvmod_pct_change(j) = res.pct_change;
        end
    end
end

%% -----------------------------------------------------------------------
% 3. HANSEN (1997) P-VALUES  (m=1, pi0=TRIM_FRAC)
%% -----------------------------------------------------------------------
cvmod_hansen_pval = NaN(n_series,1);
cvmod_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if cvmod_found(j)
        cvmod_hansen_pval(j) = hansen_supf_pval(cvmod_supw_stat(j), 1, TRIM_FRAC);
        cvmod_hansen_sig5(j) = cvmod_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cv_mod = table( ...
    string(series(:)), include_main_mod(:), ...
    cvmod_found, cvmod_supw_stat, cvmod_break_date, ...
    cvmod_hac_lag, ...
    cvmod_mu1, cvmod_mu2, cvmod_ratio, cvmod_pct_change, ...
    cvmod_hansen_pval, cvmod_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_mod', ...
        'cvmod_found','cvmod_supw_stat','cvmod_break_date', ...
        'cvmod_hac_lag', ...
        'cvmod_mu1','cvmod_mu2','cvmod_ratio','cvmod_pct_change', ...
        'cvmod_hansen_pval','cvmod_hansen_sig5'});

%% -----------------------------------------------------------------------
% 5. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'cv_modern_results.mat'), 'series_meta_cv_mod');
writetable(series_meta_cv_mod, fullfile(OUT_TABLES,'cv_modern_results.csv'));

fprintf('\n[run_cv_modern] SUMMARY\n');
fprintf('  Series with CM fit (main mod):  %d\n', ...
    sum(include_main_mod(:) & series_meta_cm_mod.cmmod_supw_found(:)));
fprintf('  Significant CV breaks (p<5%%):  %d\n', sum(cvmod_hansen_sig5));
fprintf('[run_cv_modern] Done.\n\n');
