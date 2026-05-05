%% run_cm_fullsample.m
% =========================================================================
% CONDITIONAL MEAN BREAKS — FULL SAMPLE (1959:03 – latest vintage)
%
% Estimates AR(p) conditional-mean break model for each FRED-MD series:
%   1. BIC-select AR lag order p via select_ar_lag_bic()
%   2. cm_onebreak_search  →  SSR-minimising break date in AR(p)
%   3. cm_supw_test        →  SupW with m=p+1 restrictions, HAC covariance
%   4. hansen_supf_pval_general  →  Hansen (1997) p-values for m restrictions
%
% Identical econometric logic to run_cm_historical.m; only the sample window
% changes (1959:03–latest instead of 1959:03–1999:12).
%
% HAC lag rule: h = floor(4*(T/100)^(2/9))
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_ext_full, dates_ext_full, series, include_main_full_ext
%
% OUTPUTS (saved to data_intermediate/cm_full_results.mat)
%   series_meta_cm_full        - table with CM full-sample break results
%   cm_resid_nobreak_full      - {N x 1} cell of no-break residuals
%   cm_resid_break_full        - {N x 1} cell of break-model residuals
%   cm_supw_resid_break_full   - {N x 1} cell of SupW break-model residuals
%                                (needed by run_cv_fullsample.m and run_mb_historical.m)
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cm_fullsample] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_ext_full','dates_ext_full','series','include_main_full_ext');

n_series = length(series);

%% -----------------------------------------------------------------------
% 2. LAG SELECTION + ONE-BREAK SEARCH
%% -----------------------------------------------------------------------
fprintf('[run_cm_fullsample] BIC lag selection (pmax=%d) + one-break search...\n', PMAX);

cmfull_p             = NaN(n_series,1);
cmfull_bic           = NaN(n_series,1);
cmfull_Tlag          = NaN(n_series,1);
cmfull_found         = false(n_series,1);
cmfull_break_date    = strings(n_series,1);
cmfull_nobs_raw      = NaN(n_series,1);
cmfull_nobs_reg      = NaN(n_series,1);
cmfull_ssr_nobreak   = NaN(n_series,1);
cmfull_ssr_best      = NaN(n_series,1);
cmfull_ssr_gain_pct  = NaN(n_series,1);

cm_resid_nobreak_full = cell(n_series,1);
cm_resid_break_full   = cell(n_series,1);

for j = 1:n_series
    if include_main_full_ext(j)
        yj     = yt_ext_full(:,j);
        lagres = select_ar_lag_bic(yj, PMAX);

        if ~isnan(lagres.p_best)
            cmfull_p(j)   = lagres.p_best;
            cmfull_bic(j) = lagres.bic_best;
            cmfull_Tlag(j) = lagres.T_used;

            cmres = cm_onebreak_search(yj, dates_ext_full, lagres.p_best, TRIM_FRAC);

            if cmres.success
                cmfull_found(j)        = true;
                cmfull_break_date(j)   = string(cmres.break_date,'dd-MMM-yyyy');
                cmfull_nobs_raw(j)     = cmres.nobs_raw;
                cmfull_nobs_reg(j)     = cmres.nobs_reg;
                cmfull_ssr_nobreak(j)  = cmres.ssr_nobreak;
                cmfull_ssr_best(j)     = cmres.ssr_best;
                cmfull_ssr_gain_pct(j) = cmres.ssr_gain_pct;
                cm_resid_nobreak_full{j} = cmres.resid_nobreak;
                cm_resid_break_full{j}   = cmres.resid_break;
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 3. SUP-W TEST WITH HAC COVARIANCE  (m = p+1 restrictions)
%% -----------------------------------------------------------------------
fprintf('[run_cm_fullsample] Running CM SupW tests...\n');

cmfull_supw_found        = false(n_series,1);
cmfull_supw_stat         = NaN(n_series,1);
cmfull_supw_break_date   = strings(n_series,1);
cmfull_supw_hac_lag      = NaN(n_series,1);
cmfull_supw_m            = NaN(n_series,1);
cmfull_supw_ssr_nobreak  = NaN(n_series,1);
cmfull_supw_ssr_best     = NaN(n_series,1);
cmfull_supw_ssr_gain_pct = NaN(n_series,1);

cm_supw_resid_break_full = cell(n_series,1);

for j = 1:n_series
    if include_main_full_ext(j) && ~isnan(cmfull_p(j))
        yj  = yt_ext_full(:,j);
        pj  = cmfull_p(j);
        res = cm_supw_test(yj, dates_ext_full, pj, TRIM_FRAC, []);

        if res.success
            cmfull_supw_found(j)        = true;
            cmfull_supw_stat(j)         = res.supW;
            cmfull_supw_break_date(j)   = string(res.break_date,'dd-MMM-yyyy');
            cmfull_supw_hac_lag(j)      = res.hac_lag;
            cmfull_supw_m(j)            = res.m_restr;
            cmfull_supw_ssr_nobreak(j)  = res.ssr_nobreak;
            cmfull_supw_ssr_best(j)     = res.ssr_best;
            if res.ssr_nobreak > 0
                cmfull_supw_ssr_gain_pct(j) = ...
                    100*(res.ssr_nobreak - res.ssr_best) / res.ssr_nobreak;
            end
            cm_supw_resid_break_full{j} = res.resid_break;
        end
    end
end

%% -----------------------------------------------------------------------
% 4. HANSEN P-VALUES  (m = p+1, pi0 = TRIM_FRAC)
%% -----------------------------------------------------------------------
cmfull_hansen_pval = NaN(n_series,1);
cmfull_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if cmfull_supw_found(j)
        m_j = cmfull_supw_m(j);
        cmfull_hansen_pval(j) = hansen_supf_pval_general(cmfull_supw_stat(j), m_j, TRIM_FRAC);
        cmfull_hansen_sig5(j) = cmfull_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 5. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cm_full = table( ...
    string(series(:)), include_main_full_ext(:), ...
    cmfull_p, cmfull_bic, cmfull_Tlag, ...
    cmfull_found, cmfull_break_date, ...
    cmfull_nobs_raw, cmfull_nobs_reg, ...
    cmfull_ssr_nobreak, cmfull_ssr_best, cmfull_ssr_gain_pct, ...
    cmfull_supw_found, cmfull_supw_stat, cmfull_supw_break_date, ...
    cmfull_supw_hac_lag, cmfull_supw_m, ...
    cmfull_supw_ssr_nobreak, cmfull_supw_ssr_best, cmfull_supw_ssr_gain_pct, ...
    cmfull_hansen_pval, cmfull_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_full_ext', ...
        'cmfull_p','cmfull_bic','cmfull_Tlag', ...
        'cmfull_found','cmfull_break_date', ...
        'cmfull_nobs_raw','cmfull_nobs_reg', ...
        'cmfull_ssr_nobreak','cmfull_ssr_best','cmfull_ssr_gain_pct', ...
        'cmfull_supw_found','cmfull_supw_stat','cmfull_supw_break_date', ...
        'cmfull_supw_hac_lag','cmfull_supw_m', ...
        'cmfull_supw_ssr_nobreak','cmfull_supw_ssr_best','cmfull_supw_ssr_gain_pct', ...
        'cmfull_hansen_pval','cmfull_hansen_sig5'});

%% -----------------------------------------------------------------------
% 6. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'cm_full_results.mat'), ...
    'series_meta_cm_full', ...
    'cm_resid_nobreak_full', 'cm_resid_break_full', ...
    'cm_supw_resid_break_full');

writetable(series_meta_cm_full, fullfile(OUT_TABLES,'cm_full_results.csv'));

fprintf('\n[run_cm_fullsample] SUMMARY\n');
fprintf('  Series in main sample:              %d\n', sum(include_main_full_ext));
fprintf('  Successful one-break searches:      %d\n', sum(cmfull_found));
fprintf('  Significant CM breaks (p<5%%):      %d\n', sum(cmfull_hansen_sig5));
fprintf('[run_cm_fullsample] Done.\n\n');
