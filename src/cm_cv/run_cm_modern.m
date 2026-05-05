%% run_cm_modern.m
% =========================================================================
% CONDITIONAL MEAN BREAKS — MODERN-ERA APPENDIX WINDOW (1999:01 – latest)
%
% Appendix exercise: estimates the AR(p) conditional-mean break model for
% each FRED-MD series over the modern-era subsample.
%
% Identical econometric logic to run_cm_fullsample.m; only the sample window
% changes (1999:01–latest instead of 1959:03–latest).
%
%   1. BIC-select AR lag order p via select_ar_lag_bic()
%   2. cm_onebreak_search  →  SSR-minimising break date in AR(p)
%   3. cm_supw_test        →  SupW with m=p+1 restrictions, HAC covariance
%   4. hansen_supf_pval_general  →  Hansen (1997) p-values for m restrictions
%
% CRITICAL OUTPUT:
%   cm_supw_resid_break_mod  — residuals from the SupW CM model at the
%   best break date.  These are the inputs for CV in run_mb_modern.m.
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_mod, dates_mod, series, include_main_mod
%
% OUTPUTS (saved to data_intermediate/cm_modern_results.mat)
%   series_meta_cm_mod        - table with CM modern-era break results
%   cm_resid_nobreak_mod      - {N x 1} cell of no-break residuals
%   cm_resid_break_mod        - {N x 1} cell of break-model residuals
%   cm_supw_resid_break_mod   - {N x 1} cell of SupW break-model residuals
%                               (used by run_cv_modern.m and run_mb_modern.m)
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cm_modern] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_mod','dates_mod','series','include_main_mod');

n_series = length(series);

fprintf('[run_cm_modern] Window: %s – %s  (T=%d)\n', ...
    string(dates_mod(1),'yyyy-MM'), string(dates_mod(end),'yyyy-MM'), length(dates_mod));

%% -----------------------------------------------------------------------
% 2. LAG SELECTION + ONE-BREAK SEARCH
%% -----------------------------------------------------------------------
fprintf('[run_cm_modern] BIC lag selection (pmax=%d) + one-break search...\n', PMAX);

cmmod_p             = NaN(n_series,1);
cmmod_bic           = NaN(n_series,1);
cmmod_Tlag          = NaN(n_series,1);
cmmod_found         = false(n_series,1);
cmmod_break_date    = strings(n_series,1);
cmmod_nobs_raw      = NaN(n_series,1);
cmmod_nobs_reg      = NaN(n_series,1);
cmmod_ssr_nobreak   = NaN(n_series,1);
cmmod_ssr_best      = NaN(n_series,1);
cmmod_ssr_gain_pct  = NaN(n_series,1);

cm_resid_nobreak_mod = cell(n_series,1);
cm_resid_break_mod   = cell(n_series,1);

for j = 1:n_series
    if include_main_mod(j)
        yj     = yt_mod(:,j);
        lagres = select_ar_lag_bic(yj, PMAX);

        if ~isnan(lagres.p_best)
            cmmod_p(j)    = lagres.p_best;
            cmmod_bic(j)  = lagres.bic_best;
            cmmod_Tlag(j) = lagres.T_used;

            cmres = cm_onebreak_search(yj, dates_mod, lagres.p_best, TRIM_FRAC);

            if cmres.success
                cmmod_found(j)        = true;
                cmmod_break_date(j)   = string(cmres.break_date,'dd-MMM-yyyy');
                cmmod_nobs_raw(j)     = cmres.nobs_raw;
                cmmod_nobs_reg(j)     = cmres.nobs_reg;
                cmmod_ssr_nobreak(j)  = cmres.ssr_nobreak;
                cmmod_ssr_best(j)     = cmres.ssr_best;
                cmmod_ssr_gain_pct(j) = cmres.ssr_gain_pct;
                cm_resid_nobreak_mod{j} = cmres.resid_nobreak;
                cm_resid_break_mod{j}   = cmres.resid_break;
            end
        end
    end
    if mod(j,20)==0; fprintf('  CM mod: %d / %d\n', j, n_series); end
end

%% -----------------------------------------------------------------------
% 3. SUP-W TEST WITH HAC COVARIANCE  (m = p+1 restrictions)
%% -----------------------------------------------------------------------
fprintf('[run_cm_modern] Running CM SupW tests...\n');

cmmod_supw_found        = false(n_series,1);
cmmod_supw_stat         = NaN(n_series,1);
cmmod_supw_break_date   = strings(n_series,1);
cmmod_supw_hac_lag      = NaN(n_series,1);
cmmod_supw_m            = NaN(n_series,1);
cmmod_supw_ssr_nobreak  = NaN(n_series,1);
cmmod_supw_ssr_best     = NaN(n_series,1);
cmmod_supw_ssr_gain_pct = NaN(n_series,1);

cm_supw_resid_break_mod = cell(n_series,1);

for j = 1:n_series
    if include_main_mod(j) && ~isnan(cmmod_p(j))
        yj  = yt_mod(:,j);
        pj  = cmmod_p(j);
        res = cm_supw_test(yj, dates_mod, pj, TRIM_FRAC, []);

        if res.success
            cmmod_supw_found(j)        = true;
            cmmod_supw_stat(j)         = res.supW;
            cmmod_supw_break_date(j)   = string(res.break_date,'dd-MMM-yyyy');
            cmmod_supw_hac_lag(j)      = res.hac_lag;
            cmmod_supw_m(j)            = res.m_restr;
            cmmod_supw_ssr_nobreak(j)  = res.ssr_nobreak;
            cmmod_supw_ssr_best(j)     = res.ssr_best;
            if res.ssr_nobreak > 0
                cmmod_supw_ssr_gain_pct(j) = ...
                    100*(res.ssr_nobreak - res.ssr_best) / res.ssr_nobreak;
            end
            cm_supw_resid_break_mod{j} = res.resid_break;
        end
    end
end

%% -----------------------------------------------------------------------
% 4. HANSEN P-VALUES  (m = p+1, pi0 = TRIM_FRAC)
%% -----------------------------------------------------------------------
cmmod_hansen_pval = NaN(n_series,1);
cmmod_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if cmmod_supw_found(j)
        m_j = cmmod_supw_m(j);
        cmmod_hansen_pval(j) = hansen_supf_pval_general(cmmod_supw_stat(j), m_j, TRIM_FRAC);
        cmmod_hansen_sig5(j) = cmmod_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 5. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cm_mod = table( ...
    string(series(:)), include_main_mod(:), ...
    cmmod_p, cmmod_bic, cmmod_Tlag, ...
    cmmod_found, cmmod_break_date, ...
    cmmod_nobs_raw, cmmod_nobs_reg, ...
    cmmod_ssr_nobreak, cmmod_ssr_best, cmmod_ssr_gain_pct, ...
    cmmod_supw_found, cmmod_supw_stat, cmmod_supw_break_date, ...
    cmmod_supw_hac_lag, cmmod_supw_m, ...
    cmmod_supw_ssr_nobreak, cmmod_supw_ssr_best, cmmod_supw_ssr_gain_pct, ...
    cmmod_hansen_pval, cmmod_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_mod', ...
        'cmmod_p','cmmod_bic','cmmod_Tlag', ...
        'cmmod_found','cmmod_break_date', ...
        'cmmod_nobs_raw','cmmod_nobs_reg', ...
        'cmmod_ssr_nobreak','cmmod_ssr_best','cmmod_ssr_gain_pct', ...
        'cmmod_supw_found','cmmod_supw_stat','cmmod_supw_break_date', ...
        'cmmod_supw_hac_lag','cmmod_supw_m', ...
        'cmmod_supw_ssr_nobreak','cmmod_supw_ssr_best','cmmod_supw_ssr_gain_pct', ...
        'cmmod_hansen_pval','cmmod_hansen_sig5'});

%% -----------------------------------------------------------------------
% 6. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'cm_modern_results.mat'), ...
    'series_meta_cm_mod', ...
    'cm_resid_nobreak_mod', 'cm_resid_break_mod', ...
    'cm_supw_resid_break_mod');

writetable(series_meta_cm_mod, fullfile(OUT_TABLES,'cm_modern_results.csv'));

fprintf('\n[run_cm_modern] SUMMARY\n');
fprintf('  Series in main modern-era sample:  %d\n', sum(include_main_mod));
fprintf('  Successful one-break searches:     %d\n', sum(cmmod_found));
fprintf('  Significant CM breaks (p<5%%):     %d\n', sum(cmmod_hansen_sig5));
fprintf('[run_cm_modern] Done.\n\n');
