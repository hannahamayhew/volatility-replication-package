%% run_cm_2007.m
% =========================================================================
% CONDITIONAL MEAN BREAKS — 2007 WINDOW (1959:03 – 2007:12)
%
% Estimates AR(p) conditional-mean break model for each FRED-MD series:
%   1. BIC-select AR lag order p via select_ar_lag_bic()
%   2. cm_onebreak_search  →  SSR-minimising break date in AR(p)
%   3. cm_supw_test        →  SupW with m=p+1 restrictions, HAC covariance
%   4. hansen_supf_pval_general  →  Hansen (1997) p-values for m restrictions
%
% Identical econometric logic to run_cm_historical.m; only the sample window
% changes (1959:03–2007:12 instead of 1959:03–1999:12).
%
% HAC lag rule: h = floor(4*(T/100)^(2/9))
%
% INPUTS  (loaded from data_intermediate/processed_data.mat)
%   yt_ext_2007, dates_ext_2007, series, include_main_2007
%
% OUTPUTS (saved to data_intermediate/cm_2007_results.mat)
%   series_meta_cm_2007        - table with CM 2007-window break results
%   cm_resid_nobreak_2007      - {N x 1} cell of no-break residuals
%   cm_resid_break_2007        - {N x 1} cell of break-model residuals
%   cm_supw_resid_break_2007   - {N x 1} cell of SupW break-model residuals
%                                (needed by run_cv_2007.m and run_mb_historical.m)
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_cm_2007] Loading processed data...\n');
load(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt_ext_2007','dates_ext_2007','series','include_main_2007');

n_series = length(series);

%% -----------------------------------------------------------------------
% 2. LAG SELECTION + ONE-BREAK SEARCH
%% -----------------------------------------------------------------------
fprintf('[run_cm_2007] BIC lag selection (pmax=%d) + one-break search...\n', PMAX);

cm2007_p             = NaN(n_series,1);
cm2007_bic           = NaN(n_series,1);
cm2007_Tlag          = NaN(n_series,1);
cm2007_found         = false(n_series,1);
cm2007_break_date    = strings(n_series,1);
cm2007_nobs_raw      = NaN(n_series,1);
cm2007_nobs_reg      = NaN(n_series,1);
cm2007_ssr_nobreak   = NaN(n_series,1);
cm2007_ssr_best      = NaN(n_series,1);
cm2007_ssr_gain_pct  = NaN(n_series,1);

cm_resid_nobreak_2007 = cell(n_series,1);
cm_resid_break_2007   = cell(n_series,1);

for j = 1:n_series
    if include_main_2007(j)
        yj     = yt_ext_2007(:,j);
        lagres = select_ar_lag_bic(yj, PMAX);

        if ~isnan(lagres.p_best)
            cm2007_p(j)   = lagres.p_best;
            cm2007_bic(j) = lagres.bic_best;
            cm2007_Tlag(j) = lagres.T_used;

            cmres = cm_onebreak_search(yj, dates_ext_2007, lagres.p_best, TRIM_FRAC);

            if cmres.success
                cm2007_found(j)        = true;
                cm2007_break_date(j)   = string(cmres.break_date,'dd-MMM-yyyy');
                cm2007_nobs_raw(j)     = cmres.nobs_raw;
                cm2007_nobs_reg(j)     = cmres.nobs_reg;
                cm2007_ssr_nobreak(j)  = cmres.ssr_nobreak;
                cm2007_ssr_best(j)     = cmres.ssr_best;
                cm2007_ssr_gain_pct(j) = cmres.ssr_gain_pct;
                cm_resid_nobreak_2007{j} = cmres.resid_nobreak;
                cm_resid_break_2007{j}   = cmres.resid_break;
            end
        end
    end
end

%% -----------------------------------------------------------------------
% 3. SUP-W TEST WITH HAC COVARIANCE  (m = p+1 restrictions)
%% -----------------------------------------------------------------------
fprintf('[run_cm_2007] Running CM SupW tests...\n');

cm2007_supw_found        = false(n_series,1);
cm2007_supw_stat         = NaN(n_series,1);
cm2007_supw_break_date   = strings(n_series,1);
cm2007_supw_hac_lag      = NaN(n_series,1);
cm2007_supw_m            = NaN(n_series,1);
cm2007_supw_ssr_nobreak  = NaN(n_series,1);
cm2007_supw_ssr_best     = NaN(n_series,1);
cm2007_supw_ssr_gain_pct = NaN(n_series,1);

cm_supw_resid_break_2007 = cell(n_series,1);

for j = 1:n_series
    if include_main_2007(j) && ~isnan(cm2007_p(j))
        yj  = yt_ext_2007(:,j);
        pj  = cm2007_p(j);
        res = cm_supw_test(yj, dates_ext_2007, pj, TRIM_FRAC, []);

        if res.success
            cm2007_supw_found(j)        = true;
            cm2007_supw_stat(j)         = res.supW;
            cm2007_supw_break_date(j)   = string(res.break_date,'dd-MMM-yyyy');
            cm2007_supw_hac_lag(j)      = res.hac_lag;
            cm2007_supw_m(j)            = res.m_restr;
            cm2007_supw_ssr_nobreak(j)  = res.ssr_nobreak;
            cm2007_supw_ssr_best(j)     = res.ssr_best;
            if res.ssr_nobreak > 0
                cm2007_supw_ssr_gain_pct(j) = ...
                    100*(res.ssr_nobreak - res.ssr_best) / res.ssr_nobreak;
            end
            cm_supw_resid_break_2007{j} = res.resid_break;
        end
    end
end

%% -----------------------------------------------------------------------
% 4. HANSEN P-VALUES  (m = p+1, pi0 = TRIM_FRAC)
%% -----------------------------------------------------------------------
cm2007_hansen_pval = NaN(n_series,1);
cm2007_hansen_sig5 = false(n_series,1);

for j = 1:n_series
    if cm2007_supw_found(j)
        m_j = cm2007_supw_m(j);
        cm2007_hansen_pval(j) = hansen_supf_pval_general(cm2007_supw_stat(j), m_j, TRIM_FRAC);
        cm2007_hansen_sig5(j) = cm2007_hansen_pval(j) < 0.05;
    end
end

%% -----------------------------------------------------------------------
% 5. ASSEMBLE RESULTS TABLE
%% -----------------------------------------------------------------------
series_meta_cm_2007 = table( ...
    string(series(:)), include_main_2007(:), ...
    cm2007_p, cm2007_bic, cm2007_Tlag, ...
    cm2007_found, cm2007_break_date, ...
    cm2007_nobs_raw, cm2007_nobs_reg, ...
    cm2007_ssr_nobreak, cm2007_ssr_best, cm2007_ssr_gain_pct, ...
    cm2007_supw_found, cm2007_supw_stat, cm2007_supw_break_date, ...
    cm2007_supw_hac_lag, cm2007_supw_m, ...
    cm2007_supw_ssr_nobreak, cm2007_supw_ssr_best, cm2007_supw_ssr_gain_pct, ...
    cm2007_hansen_pval, cm2007_hansen_sig5, ...
    'VariableNames',{ ...
        'series_name','include_main_2007', ...
        'cm2007_p','cm2007_bic','cm2007_Tlag', ...
        'cm2007_found','cm2007_break_date', ...
        'cm2007_nobs_raw','cm2007_nobs_reg', ...
        'cm2007_ssr_nobreak','cm2007_ssr_best','cm2007_ssr_gain_pct', ...
        'cm2007_supw_found','cm2007_supw_stat','cm2007_supw_break_date', ...
        'cm2007_supw_hac_lag','cm2007_supw_m', ...
        'cm2007_supw_ssr_nobreak','cm2007_supw_ssr_best','cm2007_supw_ssr_gain_pct', ...
        'cm2007_hansen_pval','cm2007_hansen_sig5'});

%% -----------------------------------------------------------------------
% 6. SAVE
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end
if ~exist(OUT_TABLES,'dir'); mkdir(OUT_TABLES); end

save(fullfile(DATA_INT,'cm_2007_results.mat'), ...
    'series_meta_cm_2007', ...
    'cm_resid_nobreak_2007', 'cm_resid_break_2007', ...
    'cm_supw_resid_break_2007');

writetable(series_meta_cm_2007, fullfile(OUT_TABLES,'cm_2007_results.csv'));

fprintf('\n[run_cm_2007] SUMMARY\n');
fprintf('  Series in main sample:              %d\n', sum(include_main_2007));
fprintf('  Successful one-break searches:      %d\n', sum(cm2007_found));
fprintf('  Significant CM breaks (p<5%%):      %d\n', sum(cm2007_hansen_sig5));
fprintf('[run_cm_2007] Done.\n\n');
