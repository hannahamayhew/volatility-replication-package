%% run_cm_historical.m
% =========================================================================
% CONDITIONAL MEAN BREAKS — HISTORICAL WINDOW (1959:03 – 1999:12)
%
% Estimates AR(p) conditional-mean break model for each FRED-MD series:
%   1. BIC-select AR lag order p via select_ar_lag_bic()
%   2. cm_onebreak_search  →  SSR-minimising break date in AR(p)
%   3. cm_supw_test        →  SupW with m=p+1 restrictions, HAC covariance
%   4. hansen_supf_pval_general  →  Hansen (1997) p-values for m restrictions
%
% HAC lag rule: h = floor(4*(T/100)^(2/9))
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_rep, dates_rep_t, series, include_main_hist
%
% OUTPUTS (saved to data_intermediate/cm_hist_results.mat)
%   series_meta_cm_hist        - table with CM break results
%   cm_resid_nobreak_hist      - {N x 1} cell of no-break residuals
%   cm_resid_break_hist        - {N x 1} cell of break-model residuals
%   cm_supw_resid_break_hist   - {N x 1} cell of SupW break-model residuals
%                                (needed by run_cv_historical.m)
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cm_historical] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_rep','dates_rep_t','series','include_main_hist');

n_series = length(series);

%% -----------------------------------------------------------------------
% 2. LAG SELECTION + ONE-BREAK SEARCH
%% -----------------------------------------------------------------------
fprintf('[run_cm_historical] BIC lag selection (pmax=%d) + one-break search...\n', PMAX);

cm_p_hist             = NaN(n_series,1);
cm_bic_hist           = NaN(n_series,1);
cm_Tlag_hist          = NaN(n_series,1);
cm_found_hist         = false(n_series,1);
cm_break_date_hist    = strings(n_series,1);
cm_nobs_raw_hist      = NaN(n_series,1);
cm_nobs_reg_hist      = NaN(n_series,1);
cm_ssr_nobreak_hist   = NaN(n_series,1);
cm_ssr_best_hist      = NaN(n_series,1);
cm_ssr_gain_pct_hist  = NaN(n_series,1);
cm_trim_n_hist        = NaN(n_series,1);
cm_ncand_hist         = NaN(n_series,1);

% Cell arrays for residuals (needed by CV scripts)
cm_resid_nobreak_hist   = cell(n_series,1);
cm_resid_break_hist     = cell(n_series,1);

for j = 1:n_series
    if include_main_hist(j)
        yj     = yt_rep(:,j);
        lagres = select_ar_lag_bic(yj, PMAX);

        if ~isnan(lagres.p_best)
            cm_p_hist(j)    = lagres.p_best;
            cm_bic_hist(j)  = lagres.bic_best;
            cm_Tlag_hist(j) = lagres.T_used;

            cmres = cm_onebreak_search(yj, dates_rep_t, lagres.p_best, TRIM_FRAC);

            if cmres.success
                cm_found_hist(j)         = true;
                cm_break_date_hist(j)    = string(cmres.break_date,'dd-MMM-yyyy');
                cm_nobs_raw_hist(j)      = cmres.nobs_raw;
                cm_nobs_reg_hist(j)      = cmres.nobs_reg;
                cm_ssr_nobreak_hist(j)   = cmres.ssr_nobreak;
                cm_ssr_best_hist(j)      = cmres.ssr_best;
                cm_ssr_gain_pct_hist(j)  = cmres.ssr_gain_pct;
                cm_trim_n_hist(j)        = cmres.trim_n;
                cm_ncand_hist(j)         = cmres.n_candidates;
                cm_resid_nobreak_hist{j} = cmres.resid_nobreak;
                cm_resid_break_hist{j}   = cmres.resid_break;
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 3. SUP-W TEST WITH HAC COVARIANCE  (m = p+1 restrictions)
%% -----------------------------------------------------------------------
fprintf('[run_cm_historical] Running CM SupW tests...\n');

cm_supw_found_hist       = false(n_series,1);
cm_supw_stat_hist        = NaN(n_series,1);
cm_supw_break_date_hist  = strings(n_series,1);
cm_supw_hac_lag_hist     = NaN(n_series,1);
cm_supw_m_hist           = NaN(n_series,1);
cm_supw_nobs_raw_hist    = NaN(n_series,1);
cm_supw_nobs_reg_hist    = NaN(n_series,1);
cm_supw_trim_n_hist      = NaN(n_series,1);
cm_supw_ncand_hist       = NaN(n_series,1);
cm_supw_ssr_nobreak_hist = NaN(n_series,1);
cm_supw_ssr_best_hist    = NaN(n_series,1);
cm_supw_ssr_gain_pct_hist = NaN(n_series,1);

cm_supw_resid_break_hist = cell(n_series,1);

for j = 1:n_series
    if include_main_hist(j) && ~isnan(cm_p_hist(j))
        yj  = yt_rep(:,j);
        pj  = cm_p_hist(j);
        res = cm_supw_test(yj, dates_rep_t, pj, TRIM_FRAC, []);

        if res.success
            cm_supw_found_hist(j)       = true;
            cm_supw_stat_hist(j)        = res.supW;
            cm_supw_break_date_hist(j)  = string(res.break_date,'dd-MMM-yyyy');
            cm_supw_hac_lag_hist(j)     = res.hac_lag;
            cm_supw_m_hist(j)           = res.m_restr;  % cm_supw_test stores m as m_restr
            cm_supw_nobs_raw_hist(j)    = res.nobs_raw;
            cm_supw_nobs_reg_hist(j)    = res.nobs_reg;
            cm_supw_trim_n_hist(j)      = res.trim_n;
            cm_supw_ncand_hist(j)       = res.n_candidates;
            cm_supw_ssr_nobreak_hist(j) = res.ssr_nobreak;
            cm_supw_ssr_best_hist(j)    = res.ssr_best;
            if res.ssr_nobreak > 0
                cm_supw_ssr_gain_pct_hist(j) = ...
                    100*(res.ssr_nobreak-res.ssr_best)/res.ssr_nobreak;
            end
            cm_supw_resid_break_hist{j} = res.resid_break;
        end
    end
end

%% -----------------------------------------------------------------------
% 4. HANSEN P-VALUES  (m = p+1, pi0 = 0.15)
%% -----------------------------------------------------------------------
cm_hansen_pval_hist = NaN(n_series,1);
cm_hansen_sig5_hist = false(n_series,1);

for j = 1:n_series
    if cm_supw_found_hist(j)
        m_j = cm_supw_m_hist(j);
        cm_hansen_pval_hist(j) = hansen_supf_pval_general(cm_supw_stat_hist(j), m_j, TRIM_FRAC);
        cm_hansen_sig5_hist(j) = cm_hansen_pval_hist(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 5. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cm_hist = table( ...
    string(series(:)), include_main_hist(:), ...
    cm_p_hist, cm_bic_hist, cm_Tlag_hist, ...
    cm_found_hist, cm_break_date_hist, ...
    cm_nobs_raw_hist, cm_nobs_reg_hist, ...
    cm_ssr_nobreak_hist, cm_ssr_best_hist, cm_ssr_gain_pct_hist, ...
    cm_supw_found_hist, cm_supw_stat_hist, cm_supw_break_date_hist, ...
    cm_supw_hac_lag_hist, cm_supw_m_hist, ...
    cm_supw_ssr_nobreak_hist, cm_supw_ssr_best_hist, cm_supw_ssr_gain_pct_hist, ...
    cm_hansen_pval_hist, cm_hansen_sig5_hist, ...
    'VariableNames',{ ...
        'series_name','include_main_hist', ...
        'cm_p_hist','cm_bic_hist','cm_Tlag_hist', ...
        'cm_found_hist','cm_break_date_hist', ...
        'cm_nobs_raw_hist','cm_nobs_reg_hist', ...
        'cm_ssr_nobreak_hist','cm_ssr_best_hist','cm_ssr_gain_pct_hist', ...
        'cm_supw_found_hist','cm_supw_stat_hist','cm_supw_break_date_hist', ...
        'cm_supw_hac_lag_hist','cm_supw_m_hist', ...
        'cm_supw_ssr_nobreak_hist','cm_supw_ssr_best_hist','cm_supw_ssr_gain_pct_hist', ...
        'cm_hansen_pval_hist','cm_hansen_sig5_hist'});

%% -----------------------------------------------------------------------
% 6. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'cm_hist_results.mat'), ...
    'series_meta_cm_hist', ...
    'cm_resid_nobreak_hist', 'cm_resid_break_hist', ...
    'cm_supw_resid_break_hist');

writetable(series_meta_cm_hist, fullfile(OUT_TABLES,'cm_hist_results.csv'));

fprintf('\n[run_cm_historical] SUMMARY\n');
fprintf('  Series in main sample:              %d\n', sum(include_main_hist));
fprintf('  Successful one-break searches:      %d\n', sum(cm_found_hist));
fprintf('  Significant CM breaks (p<5%%):      %d\n', sum(cm_hansen_sig5_hist));
fprintf('[run_cm_historical] Done.\n\n');
