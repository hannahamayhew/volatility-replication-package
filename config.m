%% config.m
% Central configuration for the FRED-MD volatility-break replication package.
%
% Run this file (or call it from master_run.m) to populate all path and
% parameter variables into the caller's workspace.
%
% Usage:
%   run('config.m')          % from within the replication_package folder
%   config                   % if replication_package/ is on the path

% =========================================================================
% ROOT PATHS
% =========================================================================

% Absolute root of this replication package
% mfilename('fullpath') is unreliable when called via run() — fall back to pwd,
% which master_run.m anchors to replication_package/ with cd() before loading config.
pkg_root_try = fileparts(mfilename('fullpath'));
if isempty(pkg_root_try)
    pkg_root_try = pwd;
end
PKG_ROOT = pkg_root_try;

% Raw input data (place FRED-MD .csv here)
DATA_RAW  = fullfile(PKG_ROOT, 'data_raw');

% Intermediate products (processed panels, saved .mat files)
DATA_INT  = fullfile(PKG_ROOT, 'data_intermediate');

% Source sub-directories
SRC_DATA   = fullfile(PKG_ROOT, 'src', 'data_prep');
SRC_UV     = fullfile(PKG_ROOT, 'src', 'uv');
SRC_CMCV   = fullfile(PKG_ROOT, 'src', 'cm_cv');
SRC_FIGS   = fullfile(PKG_ROOT, 'src', 'tables_figures');

% Output sub-directories
OUT_TABLES = fullfile(PKG_ROOT, 'output', 'tables');
OUT_FIGS   = fullfile(PKG_ROOT, 'output', 'figures');
OUT_LOGS   = fullfile(PKG_ROOT, 'output', 'logs');

% =========================================================================
% INPUT DATA
% =========================================================================

% FRED-MD vintage CSV — lives in data_raw/ inside this package
CSV_IN = fullfile(DATA_RAW, '2026-02-MD.csv');

% =========================================================================
% BREAK-TEST PARAMETERS
% =========================================================================

% Trimming fraction for candidate break-date search (Sensier & van Dijk 2004 use 15%)
% Trimming excludes the first and last trim_frac*T observations as candidates,
% guaranteeing a minimum number of observations in each regime.
TRIM_FRAC = 0.15;

% Smooth-break logistic grid (gamma controls transition speed; larger = more abrupt)
GAMMA_GRID = [1 2 5 10 20 50 100];

% Recursive multiple-break settings
MAX_BREAKS  = 5;   % maximum number of breaks searched recursively
MIN_SEG_LEN = 60;  % hard floor: segment length never allowed below this
% Historical window uses the sample-relative 15%*T rule (naturally yields 74 months).
% Extended windows (2007, full-sample, modern-era) use a fixed 74-month floor so
% that minimum segment geometry is comparable across windows.
MIN_SEG_EXT = 74;  % fixed effective min-seg for 2007, full-sample, and modern-era windows

% Maximum AR lag order for BIC selection (conditional-mean models)
PMAX = 12;

% =========================================================================
% FACTOR ESTIMATION (FRED-MD)
% =========================================================================

DEMEAN = 2;  % 0=none, 1=demean, 2=demean+standardize, 3=recursive demean+standardize
JJ     = 2;  % information criterion: 1=PC_p1, 2=PC_p2, 3=PC_p3
KMAX   = 8;  % maximum number of factors (or 99 to force 8)

% =========================================================================
% SAMPLE WINDOWS
% =========================================================================

% Historical replication window (Sensier & van Dijk 2004)
HIST_START = datetime(1959,3,1);   % first usable date after differencing
HIST_END   = datetime(1999,12,1);  % end of historical window

% Extended window (post-Great Moderation)
WIN_2007_END = datetime(2007,12,1);  % through end of 2007

% Full sample: HIST_START through the last date available in the vintage

% Modern-era appendix window (Appendix robustness exercise)
%   Motivated by the embedded-regime hypothesis: the long 1959–2026 window
%   contains many historical regimes that may make post-1999 volatility
%   changes harder to detect cleanly.
%
% FEASIBILITY NOTE:
%   T ≈ 325 months (1999:01–2026:02).  min_seg_len_eff = max(ceil(0.15×325),60)
%   = max(49,60) = 60.  The 60-month hard floor (MIN_SEG_LEN) is always
%   binding for this window.  Achieving 5 breaks requires that all 6 resulting
%   segments each contain ≥60 obs → 360 obs needed vs T ≈ 325 available.
%   The algorithm handles this gracefully (no crash, no code change needed),
%   but the reported cap of MAX_BREAKS=5 is rarely reachable in practice for
%   this window.  Series reporting 5 breaks should be examined individually.
MOD_ERA_START = datetime(1999,1,1);  % 1999:01 through latest available month

% =========================================================================
% SERIES-INCLUSION THRESHOLDS
% =========================================================================

MAIN_THRESH     = 0.80;  % minimum share non-missing for main sample
BALANCED_THRESH = 0.99;  % minimum share non-missing for balanced sample
BALANCED_START_CUTOFF = datetime(1959,3,1);  % latest allowed first-valid date in balanced sample

% =========================================================================
% ADD ALL SOURCE FOLDERS TO PATH
% =========================================================================

addpath(SRC_DATA);
addpath(SRC_UV);
addpath(SRC_CMCV);
addpath(SRC_FIGS);

% Add helper functions (all shared .m functions bundled in helpers/)
addpath(PKG_ROOT);
addpath(fullfile(PKG_ROOT, 'helpers'));

fprintf('config.m loaded. Package root: %s\n', PKG_ROOT);
