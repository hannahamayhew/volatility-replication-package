%% run_cv_expansions_fullsample.m
% =========================================================================
% ROBUSTNESS: EXPANSIONS-ONLY CONDITIONAL VOLATILITY (FULL SAMPLE)
%
% Constructs the Sensier-style CV volatility proxy from CM break-model
% residuals (full sample), then restricts to NBER expansion months before
% running the same one-break SupW + Hansen workflow as run_cv_fullsample.m.
%
%   s_t = sqrt(pi/2) * |e_hat_t|   (expansion months only)
%
% Econometric logic is identical to run_cv_expansions_only.m;
% only the estimation window changes (1959:03–latest vintage).
%
% PREREQUISITE: run_cm_fullsample.m must have run first.
%
% INPUTS  (loaded from data_intermediate/)
%   cm_full_results.mat  → cm_supw_resid_break_full, series_meta_cm_full
%   processed_data.mat   → dates_ext_full, series, include_main_full_ext
%
% OUTPUTS (saved to data_intermediate/cv_expansions_full_results.mat)
%   series_meta_cvexpfull  - table with expansion-only CV break results
%
% KEY COLUMNS
%   cvexpfull_found      - logical: SupW test converged
%   cvexpfull_stat       - SupW statistic
%   cvexpfull_break_date - date string of SupW-maximising break
%   cvexpfull_ratio      - post/pre mean ratio among expansion months
%   cvexpfull_pct_change - percent change in s_t at break (expansion months)
%   cvexpfull_nobs_exp   - number of expansion months used
%   cvexpfull_pval       - Hansen (1997) p-value (m=1, pi0=TRIM_FRAC)
%   cvexpfull_sig5       - logical: p < 0.05
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_expansions_fullsample] Loading CM residuals and processed data...\n');
load(fullfile(DATA_INT,'cm_full_results.mat'), ...
    'cm_supw_resid_break_full','series_meta_cm_full');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_full','series','include_main_full_ext');

n_series   = length(series);
cmfull_p   = series_meta_cm_full.cmfull_p;

%% -----------------------------------------------------------------------
% 2. NBER RECESSION INDICATOR
%% -----------------------------------------------------------------------
rec_full = make_nber_recession_indicator(dates_ext_full);

%% -----------------------------------------------------------------------
% 3. EXPANSION-ONLY SUP-W LOOP
%% -----------------------------------------------------------------------
fprintf('[run_cv_expansions_fullsample] Running expansion-only SupW tests (trim=%.0f%%)...\n', ...
    TRIM_FRAC*100);

cvexpfull_found      = false(n_series,1);
cvexpfull_stat       = NaN(n_series,1);
cvexpfull_break_date = strings(n_series,1);
cvexpfull_hac_lag    = NaN(n_series,1);
cvexpfull_mu1        = NaN(n_series,1);
cvexpfull_mu2        = NaN(n_series,1);
cvexpfull_ratio      = NaN(n_series,1);
cvexpfull_pct_change = NaN(n_series,1);
cvexpfull_nobs_exp   = NaN(n_series,1);
cvexpfull_pval       = NaN(n_series,1);
cvexpfull_sig5       = false(n_series,1);

for j = 1:n_series
    if ~(include_main_full_ext(j) && series_meta_cm_full.cmfull_supw_found(j) && ...
            ~isempty(cm_supw_resid_break_full{j}))
        continue;
    end

    pj               = cmfull_p(j);
    dates_reg_j_full = dates_ext_full((pj+1):end);
    ehat_j           = cm_supw_resid_break_full{j};
    s_j              = sqrt(pi/2) * abs(ehat_j);

    Tj = length(s_j);
    if Tj < 10 || length(dates_reg_j_full) < Tj
        continue;
    end

    % Align date and recession vectors to residual length
    dates_reg_j = dates_reg_j_full(end-Tj+1:end);
    rec_reg_j   = rec_full((pj+1):end);
    rec_reg_j   = rec_reg_j(end-Tj+1:end);

    % Restrict to expansion months
    s_exp = s_j(~rec_reg_j);
    d_exp = dates_reg_j(~rec_reg_j);

    cvexpfull_nobs_exp(j) = length(s_exp);

    if length(s_exp) < 20
        continue;
    end

    % SupW test on expansion-only s_t
    res = uv_supw_test(s_exp, d_exp, TRIM_FRAC, []);

    if res.success
        cvexpfull_found(j)      = true;
        cvexpfull_stat(j)       = res.supW;
        cvexpfull_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
        cvexpfull_hac_lag(j)    = res.hac_lag;
        cvexpfull_mu1(j)        = res.mu1;
        cvexpfull_mu2(j)        = res.mu2;
        cvexpfull_ratio(j)      = res.ratio_post_pre;
        cvexpfull_pct_change(j) = res.pct_change;
        cvexpfull_pval(j)       = hansen_supf_pval(res.supW, 1, TRIM_FRAC);
        cvexpfull_sig5(j)       = cvexpfull_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cvexpfull = table( ...
    string(series(:)), include_main_full_ext(:), ...
    cvexpfull_found, cvexpfull_stat, cvexpfull_break_date, ...
    cvexpfull_hac_lag, cvexpfull_mu1, cvexpfull_mu2, ...
    cvexpfull_ratio, cvexpfull_pct_change, cvexpfull_nobs_exp, ...
    cvexpfull_pval, cvexpfull_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_full_ext', ...
        'cvexpfull_found','cvexpfull_stat','cvexpfull_break_date', ...
        'cvexpfull_hac_lag','cvexpfull_mu1','cvexpfull_mu2', ...
        'cvexpfull_ratio','cvexpfull_pct_change','cvexpfull_nobs_exp', ...
        'cvexpfull_pval','cvexpfull_sig5'});

%% -----------------------------------------------------------------------
% 5. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'cv_expansions_full_results.mat'), 'series_meta_cvexpfull');
writetable(series_meta_cvexpfull, fullfile(OUT_TABLES,'cv_expansions_full_results.csv'));

fprintf('\n[run_cv_expansions_fullsample] SUMMARY\n');
n_elig = sum(include_main_full_ext(:) & series_meta_cm_full.cmfull_supw_found(:));
fprintf('  Series eligible (main sample + CM fit):  %d\n', n_elig);
fprintf('  Tests run (sufficient expansion obs):    %d\n', sum(cvexpfull_found));
fprintf('  Significant breaks (p<5%%):              %d\n', sum(cvexpfull_sig5));
if n_elig > 0
    fprintf('  Share significant (of eligible):         %.3f\n', ...
        sum(cvexpfull_sig5) / n_elig);
end
fprintf('[run_cv_expansions_fullsample] Done.\n\n');
