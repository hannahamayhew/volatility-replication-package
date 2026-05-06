%% master_run.m
% =========================================================================
% MASTER ORCHESTRATOR — FRED-MD Volatility Break Replication Package
%
% Runs the full replication pipeline in order:
%
%   Stage 0  config.m            — paths and parameters
%   Stage 1  run_data_prep       — load, transform, save panel data
%   Stage 2a run_uv_historical   — UV breaks: historical window (–1999:12)
%   Stage 2b run_uv_2007         — UV breaks: 2007 window (–2007:12)
%   Stage 2c run_uv_fullsample   — UV breaks: full vintage (–latest)
%   Stage 3a run_cm_historical             — CM breaks + residuals: historical window
%   Stage 3b run_cv_historical             — CV breaks: historical window
%   Stage 3c run_cv_expansions_only        — CV robustness: expansion months only (hist)
%   Stage 3d run_cm_2007                   — CM breaks + residuals: 2007 window
%   Stage 3e run_cv_2007                   — CV breaks: 2007 window
%   Stage 3f run_cv_expansions_2007        — CV robustness: expansion months only (2007)
%   Stage 3g run_cm_fullsample             — CM breaks + residuals: full sample
%   Stage 3h run_cv_fullsample             — CV breaks: full sample
%   Stage 3i run_cv_expansions_fullsample  — CV robustness: expansion months only (full)
%   Stage 3j run_mb_fixed74               — UV + CV multiple-break estimation: all windows (74-month min-seg, cap=8)
%   Stage 3k run_uv_smooth                 — UV smooth-break estimation: all windows
%   Stage 3l run_cv_smooth                 — CV smooth-break estimation: all windows
%   Stage 4  run_tables_figures            — merge results and produce outputs
%
% ECONOMETRIC METHODS
%   Unconditional volatility (UV) follows Sensier & van Dijk (2004):
%       z_t = sqrt(pi/2) * |y_t - mean(y)|
%   Break search uses 15% trimming; significance uses Hansen (1997) p-values.
%   The reported one-break date is the SSR-minimising break date from
%   uv_onebreak_search; the SupW statistic is the formal test of existence.
%
% USAGE
%   cd to replication_package/ (or add it to the path), then:
%       master_run
%
% REQUIREMENTS
%   - 2026-02-MD.csv must be in data_raw/ (included)
%   - series_metadata_original_order.csv must be in data_raw/ (included)
%   - All helper functions are in helpers/ (config.m adds them to path)
%   - No additional toolboxes required
% =========================================================================

clear; close all; clc;

% Anchor CWD to replication_package/ so run('config.m') resolves correctly
% and mfilename('fullpath') inside config.m returns the right path.
cd(fileparts(mfilename('fullpath')));

t_start = tic;
fprintf('==========================================================\n');
fprintf('  FRED-MD Volatility Break Replication Package\n');
fprintf('  %s\n', datestr(now));
fprintf('==========================================================\n\n');

%% ------------------------------------------------------------------
% Stage 0: Configuration
%% ------------------------------------------------------------------
run('config.m');

%% ------------------------------------------------------------------
% Stage 1: Data preparation
%% ------------------------------------------------------------------
fprintf('\n--- Stage 1: Data Preparation ---\n');
run(fullfile(SRC_DATA, 'run_data_prep.m'));

%% ------------------------------------------------------------------
% Stage 1b: Stationarity checks — ADF + KPSS
%   Verifies that McCracken-Ng transformation codes render each series I(0).
%   Flags any failures, with nominal-side groups highlighted.
%   Results saved to output/tables/stationarity_*.csv.
%% ------------------------------------------------------------------
fprintf('\n--- Stage 1b: Stationarity Checks (ADF + KPSS) ---\n');
run(fullfile(SRC_DATA, 'run_stationarity_tests.m'));

%% ------------------------------------------------------------------
% Stage 2: Unconditional volatility
%% ------------------------------------------------------------------
fprintf('\n--- Stage 2a: UV — Historical Window ---\n');
run(fullfile(SRC_UV, 'run_uv_historical.m'));

fprintf('\n--- Stage 2b: UV — 2007 Window ---\n');
run(fullfile(SRC_UV, 'run_uv_2007.m'));

fprintf('\n--- Stage 2c: UV — Full Sample ---\n');
run(fullfile(SRC_UV, 'run_uv_fullsample.m'));

%% ------------------------------------------------------------------
% Stage 3: Conditional mean / conditional volatility
%% ------------------------------------------------------------------
fprintf('\n--- Stage 3a: CM — Historical Window ---\n');
run(fullfile(SRC_CMCV, 'run_cm_historical.m'));

fprintf('\n--- Stage 3b: CV — Historical Window ---\n');
run(fullfile(SRC_CMCV, 'run_cv_historical.m'));

fprintf('\n--- Stage 3c: CV — Expansions-Only Robustness ---\n');
run(fullfile(SRC_CMCV, 'run_cv_expansions_only.m'));

fprintf('\n--- Stage 3d: CM — 2007 Window ---\n');
run(fullfile(SRC_CMCV, 'run_cm_2007.m'));

fprintf('\n--- Stage 3e: CV — 2007 Window ---\n');
run(fullfile(SRC_CMCV, 'run_cv_2007.m'));

fprintf('\n--- Stage 3f: CV — Expansions-Only Robustness (2007 Window) ---\n');
run(fullfile(SRC_CMCV, 'run_cv_expansions_2007.m'));

fprintf('\n--- Stage 3g: CM — Full Sample ---\n');
run(fullfile(SRC_CMCV, 'run_cm_fullsample.m'));

fprintf('\n--- Stage 3h: CV — Full Sample ---\n');
run(fullfile(SRC_CMCV, 'run_cv_fullsample.m'));

fprintf('\n--- Stage 3i: CV — Expansions-Only Robustness (Full Sample) ---\n');
run(fullfile(SRC_CMCV, 'run_cv_expansions_fullsample.m'));

fprintf('\n--- Stage 3j: Multiple Breaks — All Windows (74-month min-seg, cap=8) ---\n');
run(fullfile(SRC_CMCV, 'run_mb_fixed74.m'));

fprintf('\n--- Stage 3k: UV Smooth-Break — All Windows ---\n');
run(fullfile(SRC_UV, 'run_uv_smooth.m'));

fprintf('\n--- Stage 3l: CV Smooth-Break — All Windows ---\n');
run(fullfile(SRC_CMCV, 'run_cv_smooth.m'));

%% ------------------------------------------------------------------
% Stages 3m–3p: APPENDIX — Modern-era robustness exercise (1999–latest)
%   These stages are ADDITIVE and do not alter any main results.
%   They implement the supplementary recursive multiple-break analysis
%   restricted to the modern era (1999:01 – latest available month).
%% ------------------------------------------------------------------
fprintf('\n--- Stage 3m (Appendix): UV — Modern-Era Window ---\n');
run(fullfile(SRC_UV, 'run_uv_modern.m'));

fprintf('\n--- Stage 3n (Appendix): CM — Modern-Era Window ---\n');
run(fullfile(SRC_CMCV, 'run_cm_modern.m'));

fprintf('\n--- Stage 3o (Appendix): CV — Modern-Era Window ---\n');
run(fullfile(SRC_CMCV, 'run_cv_modern.m'));

fprintf('\n--- Stage 3p (Appendix): Multiple Breaks — Modern-Era Window + Comparison ---\n');
run(fullfile(SRC_CMCV, 'run_mb_modern.m'));

%% ------------------------------------------------------------------
% Stage 3r: APPENDIX — Robustness: restricted samples
%   (1) Balanced panel  (2) Excluding 8 McCracken-Ng flagged series
%   Run manually when needed: run(fullfile(SRC_CMCV,'run_mb_robustness_samples.m'))
%% ------------------------------------------------------------------
% run(fullfile(SRC_CMCV, 'run_mb_robustness_samples.m'));

%% ------------------------------------------------------------------
% Stage 3s: APPENDIX — Robustness: Labor Market stationarity exclusions
%   Drops nonstationary + ambiguous LM series (full-sample ADF+KPSS flags)
%   and re-runs fixed-74 recursive MB to test whether zero-break dominance
%   survives. Requires: run_stationarity_tests, all UV/CM stages, run_mb_fixed74.
%   Run manually when needed: run(fullfile(SRC_CMCV,'run_mb_lm_stationarity_robust.m'))
%% ------------------------------------------------------------------
% run(fullfile(SRC_CMCV, 'run_mb_lm_stationarity_robust.m'));

%% ------------------------------------------------------------------
% Stage 4: Tables and figures
%% ------------------------------------------------------------------
fprintf('\n--- Stage 4a: Merged Series-Meta and Presentation Layer ---\n');
run(fullfile(SRC_FIGS, 'run_tables_figures.m'));

fprintf('\n--- Stage 4b: Paper Figures (Figures 1–4b and A1–A7) ---\n');
run(fullfile(SRC_FIGS, 'produce_paper_figures.m'));

%% ------------------------------------------------------------------
% Done
%% ------------------------------------------------------------------
run('config.m');
fprintf('\n==========================================================\n');
try
    fprintf('  Replication complete.  Elapsed: %.1f minutes\n', toc(t_start)/60);
catch
    fprintf('  Replication complete.\n');
end
fprintf('  Output tables → %s\n', OUT_TABLES);
fprintf('  Figures       → %s\n', OUT_FIGS);
fprintf('==========================================================\n');
