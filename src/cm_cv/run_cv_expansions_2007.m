%% run_cv_expansions_2007.m
% =========================================================================
% ROBUSTNESS: EXPANSIONS-ONLY CONDITIONAL VOLATILITY (2007 WINDOW)
%
% Constructs the Sensier-style CV volatility proxy from CM break-model
% residuals (2007 window), then restricts to NBER expansion months before
% running the same one-break SupW + Hansen workflow as run_cv_2007.m.
%
%   s_t = sqrt(pi/2) * |e_hat_t|   (expansion months only)
%
% Econometric logic is identical to run_cv_expansions_only.m;
% only the estimation window changes (1959:03–2007:12).
%
% PREREQUISITE: run_cm_2007.m must have run first.
%
% INPUTS  (loaded from data_intermediate/)
%   cm_2007_results.mat  → cm_supw_resid_break_2007, series_meta_cm_2007
%   processed_data.mat   → dates_ext_2007, series, include_main_2007
%
% OUTPUTS (saved to data_intermediate/cv_expansions_2007_results.mat)
%   series_meta_cvexp2007  - table with expansion-only CV break results
%
% KEY COLUMNS
%   cvexp2007_found      - logical: SupW test converged
%   cvexp2007_stat       - SupW statistic
%   cvexp2007_break_date - date string of SupW-maximising break
%   cvexp2007_ratio      - post/pre mean ratio among expansion months
%   cvexp2007_pct_change - percent change in s_t at break (expansion months)
%   cvexp2007_nobs_exp   - number of expansion months used
%   cvexp2007_pval       - Hansen (1997) p-value (m=1, pi0=TRIM_FRAC)
%   cvexp2007_sig5       - logical: p < 0.05
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_expansions_2007] Loading CM residuals and processed data...\n');
load(fullfile(DATA_INT,'cm_2007_results.mat'), ...
    'cm_supw_resid_break_2007','series_meta_cm_2007');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_2007','series','include_main_2007');

n_series   = length(series);
cm2007_p   = series_meta_cm_2007.cm2007_p;

%% -----------------------------------------------------------------------
% 2. NBER RECESSION INDICATOR
%% -----------------------------------------------------------------------
rec_2007 = make_nber_recession_indicator(dates_ext_2007);

%% -----------------------------------------------------------------------
% 3. EXPANSION-ONLY SUP-W LOOP
%% -----------------------------------------------------------------------
fprintf('[run_cv_expansions_2007] Running expansion-only SupW tests (trim=%.0f%%)...\n', ...
    TRIM_FRAC*100);

cvexp2007_found      = false(n_series,1);
cvexp2007_stat       = NaN(n_series,1);
cvexp2007_break_date = strings(n_series,1);
cvexp2007_hac_lag    = NaN(n_series,1);
cvexp2007_mu1        = NaN(n_series,1);
cvexp2007_mu2        = NaN(n_series,1);
cvexp2007_ratio      = NaN(n_series,1);
cvexp2007_pct_change = NaN(n_series,1);
cvexp2007_nobs_exp   = NaN(n_series,1);
cvexp2007_pval       = NaN(n_series,1);
cvexp2007_sig5       = false(n_series,1);

for j = 1:n_series
    if ~(include_main_2007(j) && series_meta_cm_2007.cm2007_supw_found(j) && ...
            ~isempty(cm_supw_resid_break_2007{j}))
        continue;
    end

    pj               = cm2007_p(j);
    dates_reg_j_full = dates_ext_2007((pj+1):end);
    ehat_j           = cm_supw_resid_break_2007{j};
    s_j              = sqrt(pi/2) * abs(ehat_j);

    Tj = length(s_j);
    if Tj < 10 || length(dates_reg_j_full) < Tj
        continue;
    end

    % Align date and recession vectors to residual length
    dates_reg_j = dates_reg_j_full(end-Tj+1:end);
    rec_reg_j   = rec_2007((pj+1):end);
    rec_reg_j   = rec_reg_j(end-Tj+1:end);

    % Restrict to expansion months
    s_exp = s_j(~rec_reg_j);
    d_exp = dates_reg_j(~rec_reg_j);

    cvexp2007_nobs_exp(j) = length(s_exp);

    if length(s_exp) < 20
        continue;
    end

    % SupW test on expansion-only s_t
    res = uv_supw_test(s_exp, d_exp, TRIM_FRAC, []);

    if res.success
        cvexp2007_found(j)      = true;
        cvexp2007_stat(j)       = res.supW;
        cvexp2007_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
        cvexp2007_hac_lag(j)    = res.hac_lag;
        cvexp2007_mu1(j)        = res.mu1;
        cvexp2007_mu2(j)        = res.mu2;
        cvexp2007_ratio(j)      = res.ratio_post_pre;
        cvexp2007_pct_change(j) = res.pct_change;
        cvexp2007_pval(j)       = hansen_supf_pval(res.supW, 1, TRIM_FRAC);
        cvexp2007_sig5(j)       = cvexp2007_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cvexp2007 = table( ...
    string(series(:)), include_main_2007(:), ...
    cvexp2007_found, cvexp2007_stat, cvexp2007_break_date, ...
    cvexp2007_hac_lag, cvexp2007_mu1, cvexp2007_mu2, ...
    cvexp2007_ratio, cvexp2007_pct_change, cvexp2007_nobs_exp, ...
    cvexp2007_pval, cvexp2007_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_2007', ...
        'cvexp2007_found','cvexp2007_stat','cvexp2007_break_date', ...
        'cvexp2007_hac_lag','cvexp2007_mu1','cvexp2007_mu2', ...
        'cvexp2007_ratio','cvexp2007_pct_change','cvexp2007_nobs_exp', ...
        'cvexp2007_pval','cvexp2007_sig5'});

%% -----------------------------------------------------------------------
% 5. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'cv_expansions_2007_results.mat'), 'series_meta_cvexp2007');
writetable(series_meta_cvexp2007, fullfile(OUT_TABLES,'cv_expansions_2007_results.csv'));

fprintf('\n[run_cv_expansions_2007] SUMMARY\n');
n_elig = sum(include_main_2007(:) & series_meta_cm_2007.cm2007_supw_found(:));
fprintf('  Series eligible (main sample + CM fit):  %d\n', n_elig);
fprintf('  Tests run (sufficient expansion obs):    %d\n', sum(cvexp2007_found));
fprintf('  Significant breaks (p<5%%):              %d\n', sum(cvexp2007_sig5));
if n_elig > 0
    fprintf('  Share significant (of eligible):         %.3f\n', ...
        sum(cvexp2007_sig5) / n_elig);
end
fprintf('[run_cv_expansions_2007] Done.\n\n');
