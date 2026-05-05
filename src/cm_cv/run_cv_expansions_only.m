%% run_cv_expansions_only.m
% =========================================================================
% ROBUSTNESS: EXPANSIONS-ONLY CONDITIONAL VOLATILITY (HISTORICAL WINDOW)
%
% Constructs the Sensier-style CV volatility proxy from CM break-model
% residuals, then restricts to NBER expansion months before running the
% same one-break SupW + Hansen workflow as run_cv_historical.m.
%
%   s_t = sqrt(pi/2) * |e_hat_t|   (expansion months only)
%
% Econometric logic is identical to run_cv_historical.m / termpaper610.m;
% only the sample mask changes (recession months dropped).
%
% PREREQUISITE: run_cm_historical.m must have run first.
%
% INPUTS  (loaded from data_intermediate/)
%   cm_hist_results.mat  → cm_supw_resid_break_hist, series_meta_cm_hist
%   processed_data.mat   → dates_rep_t, series, include_main_hist
%
% OUTPUTS (saved to data_intermediate/cv_expansions_results.mat)
%   series_meta_cvexp    - table with expansion-only CV break results
%
% KEY COLUMNS
%   cvexp_found      - logical: SupW test converged
%   cvexp_stat       - SupW statistic
%   cvexp_break_date - date string of SupW-maximising break
%   cvexp_ratio      - post/pre mean ratio among expansion months
%   cvexp_pct_change - percent change in s_t at break (expansion months)
%   cvexp_nobs_exp   - number of expansion months used
%   cvexp_pval       - Hansen (1997) p-value (m=1, pi0=TRIM_FRAC)
%   cvexp_sig5       - logical: p < 0.05
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_expansions_only] Loading CM residuals and processed data...\n');
load(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'cm_supw_resid_break_hist','series_meta_cm_hist');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','series','include_main_hist');

n_series  = length(series);
cm_p_hist = series_meta_cm_hist.cm_p_hist;

%% -----------------------------------------------------------------------
% 2. NBER RECESSION INDICATOR
%% -----------------------------------------------------------------------
rec_hist = make_nber_recession_indicator(dates_rep_t);

%% -----------------------------------------------------------------------
% 3. EXPANSION-ONLY SUP-W LOOP
%% -----------------------------------------------------------------------
fprintf('[run_cv_expansions_only] Running expansion-only SupW tests (trim=%.0f%%)...\n', ...
    TRIM_FRAC*100);

cvexp_found      = false(n_series,1);
cvexp_stat       = NaN(n_series,1);
cvexp_break_date = strings(n_series,1);
cvexp_hac_lag    = NaN(n_series,1);
cvexp_mu1        = NaN(n_series,1);
cvexp_mu2        = NaN(n_series,1);
cvexp_ratio      = NaN(n_series,1);
cvexp_pct_change = NaN(n_series,1);
cvexp_nobs_exp   = NaN(n_series,1);
cvexp_pval       = NaN(n_series,1);
cvexp_sig5       = false(n_series,1);

for j = 1:n_series
    if ~(include_main_hist(j) && series_meta_cm_hist.cm_supw_found_hist(j) && ...
            ~isempty(cm_supw_resid_break_hist{j}))
        continue;
    end

    pj               = cm_p_hist(j);
    dates_reg_j_full = dates_rep_t((pj+1):end);
    ehat_j           = cm_supw_resid_break_hist{j};
    s_j              = sqrt(pi/2) * abs(ehat_j);

    Tj = length(s_j);
    if Tj < 10 || length(dates_reg_j_full) < Tj
        continue;
    end

    % Align date and recession vectors to residual length
    dates_reg_j = dates_reg_j_full(end-Tj+1:end);
    rec_reg_j   = rec_hist((pj+1):end);
    rec_reg_j   = rec_reg_j(end-Tj+1:end);

    % Restrict to expansion months
    s_exp = s_j(~rec_reg_j);
    d_exp = dates_reg_j(~rec_reg_j);

    cvexp_nobs_exp(j) = length(s_exp);

    if length(s_exp) < 20
        continue;
    end

    % SupW test on expansion-only s_t
    res = uv_supw_test(s_exp, d_exp, TRIM_FRAC, []);

    if res.success
        cvexp_found(j)      = true;
        cvexp_stat(j)       = res.supW;
        cvexp_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
        cvexp_hac_lag(j)    = res.hac_lag;
        cvexp_mu1(j)        = res.mu1;
        cvexp_mu2(j)        = res.mu2;
        cvexp_ratio(j)      = res.ratio_post_pre;
        cvexp_pct_change(j) = res.pct_change;
        cvexp_pval(j)       = hansen_supf_pval(res.supW, 1, TRIM_FRAC);
        cvexp_sig5(j)       = cvexp_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cvexp = table( ...
    string(series(:)), include_main_hist(:), ...
    cvexp_found, cvexp_stat, cvexp_break_date, ...
    cvexp_hac_lag, cvexp_mu1, cvexp_mu2, ...
    cvexp_ratio, cvexp_pct_change, cvexp_nobs_exp, ...
    cvexp_pval, cvexp_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_hist', ...
        'cvexp_found','cvexp_stat','cvexp_break_date', ...
        'cvexp_hac_lag','cvexp_mu1','cvexp_mu2', ...
        'cvexp_ratio','cvexp_pct_change','cvexp_nobs_exp', ...
        'cvexp_pval','cvexp_sig5'});

%% -----------------------------------------------------------------------
% 5. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'cv_expansions_results.mat'), 'series_meta_cvexp');
writetable(series_meta_cvexp, fullfile(OUT_TABLES,'cv_expansions_results.csv'));

fprintf('\n[run_cv_expansions_only] SUMMARY\n');
n_elig = sum(include_main_hist(:) & series_meta_cm_hist.cm_supw_found_hist(:));
fprintf('  Series eligible (main sample + CM fit):  %d\n', n_elig);
fprintf('  Tests run (sufficient expansion obs):    %d\n', sum(cvexp_found));
fprintf('  Significant breaks (p<5%%):              %d\n', sum(cvexp_sig5));
if n_elig > 0
    fprintf('  Share significant (of eligible):         %.3f\n', ...
        sum(cvexp_sig5) / n_elig);
end
fprintf('[run_cv_expansions_only] Done.\n\n');
