%% run_cv_2007.m
% =========================================================================
% CONDITIONAL VOLATILITY BREAKS — 2007 WINDOW (1959:03 – 2007:12)
%
% Uses residuals from the AR(p) conditional-mean SupW model (2007 window)
% to construct the Sensier-style conditional volatility proxy:
%
%   s_t = sqrt(pi/2) * |e_hat_t|
%
% Then tests for a break in the mean of s_t using the same one-break
% SupW + Hansen workflow used for the historical CV.
%
% PREREQUISITE: run_cm_2007.m must have run first (provides residuals).
%
% INPUTS  (loaded from data_intermediate)
%   cm_2007_results.mat   → cm_supw_resid_break_2007, series_meta_cm_2007
%   processed_data.mat    → dates_ext_2007, series, include_main_2007
%
% OUTPUTS (saved to data_intermediate/cv_2007_results.mat)
%   series_meta_cv_2007  - table with CV 2007-window break results
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cv_2007] Loading CM residuals and processed data...\n');
load(fullfile(DATA_INT,'cm_2007_results.mat'), ...
    'cm_supw_resid_break_2007','series_meta_cm_2007');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_ext_2007','series','include_main_2007');

n_series   = length(series);
cm2007_p   = series_meta_cm_2007.cm2007_p;

%% -----------------------------------------------------------------------
% 2. CV ONE-BREAK SUP-W TEST
%% -----------------------------------------------------------------------
fprintf('[run_cv_2007] Running CV SupW tests (trim=%.0f%%)...\n', TRIM_FRAC*100);

cv2007_found      = false(n_series,1);
cv2007_supw_stat  = NaN(n_series,1);
cv2007_break_date = strings(n_series,1);
cv2007_hac_lag    = NaN(n_series,1);
cv2007_mu1        = NaN(n_series,1);
cv2007_mu2        = NaN(n_series,1);
cv2007_ratio      = NaN(n_series,1);
cv2007_pct_change = NaN(n_series,1);

for j = 1:n_series
    if include_main_2007(j) && series_meta_cm_2007.cm2007_supw_found(j) && ...
            ~isempty(cm_supw_resid_break_2007{j})

        pj          = cm2007_p(j);
        dates_reg_j = dates_ext_2007((pj+1):end);
        ehat_j      = cm_supw_resid_break_2007{j};

        Tj = length(ehat_j);
        if length(dates_reg_j) > Tj
            dates_reg_j = dates_reg_j(end-Tj+1:end);
        end

        res = cv_supw_test(ehat_j, dates_reg_j, TRIM_FRAC, []);

        if res.success
            cv2007_found(j)      = true;
            cv2007_supw_stat(j)  = res.supW;
            cv2007_break_date(j) = string(res.break_date,'dd-MMM-yyyy');
            cv2007_hac_lag(j)    = res.hac_lag;
            cv2007_mu1(j)        = res.mu1;
            cv2007_mu2(j)        = res.mu2;
            cv2007_ratio(j)      = res.ratio_post_pre;
            cv2007_pct_change(j) = res.pct_change;
        end
    end
end

%% -----------------------------------------------------------------------
% 3. HANSEN (1997) P-VALUES  (m=1, pi0=TRIM_FRAC)
%% -----------------------------------------------------------------------
cv2007_hansen_pval = NaN(n_series,1);
cv2007_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if cv2007_found(j)
        cv2007_hansen_pval(j) = hansen_supf_pval(cv2007_supw_stat(j), 1, TRIM_FRAC);
        cv2007_hansen_sig5(j) = cv2007_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 4. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cv_2007 = table( ...
    string(series(:)), include_main_2007(:), ...
    cv2007_found, cv2007_supw_stat, cv2007_break_date, ...
    cv2007_hac_lag, ...
    cv2007_mu1, cv2007_mu2, cv2007_ratio, cv2007_pct_change, ...
    cv2007_hansen_pval, cv2007_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_2007', ...
        'cv2007_found','cv2007_supw_stat','cv2007_break_date', ...
        'cv2007_hac_lag', ...
        'cv2007_mu1','cv2007_mu2','cv2007_ratio','cv2007_pct_change', ...
        'cv2007_hansen_pval','cv2007_hansen_sig5'});

%% -----------------------------------------------------------------------
% 5. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'cv_2007_results.mat'), 'series_meta_cv_2007');
writetable(series_meta_cv_2007, fullfile(OUT_TABLES,'cv_2007_results.csv'));

fprintf('\n[run_cv_2007] SUMMARY\n');
fprintf('  Series with CM fit (main 2007):    %d\n', ...
    sum(include_main_2007(:) & series_meta_cm_2007.cm2007_supw_found(:)));
fprintf('  Significant CV breaks (p<5%%):     %d\n', sum(cv2007_hansen_sig5));
fprintf('[run_cv_2007] Done.\n\n');
