%% run_tables_figures.m
% =========================================================================
% TABLES AND FIGURES
%
% Merges UV, CM, and CV result tables from data_intermediate/ and
% delegates to the original tables_figures_sensier_fredmd.m presentation
% layer, which requires a consolidated series_meta table.
%
% This script:
%   1. Loads all per-window result tables
%   2. Merges them into a single series_meta table (keyed on series_name)
%   3. Calls the original tables_figures_sensier_fredmd.m
%
% PREREQUISITE: All run_* scripts must have completed successfully.
%
% INPUTS  (loaded from data_intermediate/)
%   uv_hist_results.mat, uv_2007_results.mat, uv_full_results.mat
%   cm_hist_results.mat, cv_hist_results.mat
%   processed_data.mat   (for series names / groups)
%
% OUTPUTS (forwarded to tables_figures_sensier_fredmd.m)
%   Saves series_meta.mat in DATA_INT for use by the presentation layer.
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD ALL RESULT TABLES
%% -----------------------------------------------------------------------
fprintf('\n[run_tables_figures] Loading result tables...\n');

load(fullfile(DATA_INT,'processed_data.mat'), 'series', 'fred_group');

load(fullfile(DATA_INT,'uv_hist_results.mat'),  'series_meta_uv_hist');
load(fullfile(DATA_INT,'uv_2007_results.mat'),  'series_meta_uv_2007');
load(fullfile(DATA_INT,'uv_full_results.mat'),  'series_meta_uv_full');
load(fullfile(DATA_INT,'cm_hist_results.mat'),  'series_meta_cm_hist');
load(fullfile(DATA_INT,'cv_hist_results.mat'),  'series_meta_cv_hist');

%% -----------------------------------------------------------------------
% 2. MERGE INTO SINGLE series_meta TABLE
% (keyed on series_name, joining all per-window results)
%% -----------------------------------------------------------------------
% Start from UV hist (canonical series ordering)
series_meta = series_meta_uv_hist;

% Join 2007 UV
vars_2007 = setdiff(series_meta_uv_2007.Properties.VariableNames, ...
                    {'series_name','include_main_2007'});
for k = 1:numel(vars_2007)
    series_meta.(vars_2007{k}) = series_meta_uv_2007.(vars_2007{k});
end
series_meta.include_main_2007 = series_meta_uv_2007.include_main_2007;

% Join full UV
vars_full = setdiff(series_meta_uv_full.Properties.VariableNames, ...
                    {'series_name','include_main_full_ext'});
for k = 1:numel(vars_full)
    series_meta.(vars_full{k}) = series_meta_uv_full.(vars_full{k});
end
series_meta.include_main_full_ext = series_meta_uv_full.include_main_full_ext;

% Join CM hist
vars_cm = setdiff(series_meta_cm_hist.Properties.VariableNames, ...
                  {'series_name','include_main_hist'});
for k = 1:numel(vars_cm)
    series_meta.(vars_cm{k}) = series_meta_cm_hist.(vars_cm{k});
end

% Join CV hist
vars_cv = setdiff(series_meta_cv_hist.Properties.VariableNames, ...
                  {'series_name','include_main_hist'});
for k = 1:numel(vars_cv)
    series_meta.(vars_cv{k}) = series_meta_cv_hist.(vars_cv{k});
end

% Join CM 2007 (if Stage 3d has run)
cm2007_path = fullfile(DATA_INT,'cm_2007_results.mat');
if exist(cm2007_path,'file')
    load(cm2007_path, 'series_meta_cm_2007');
    vars_cm2007 = setdiff(series_meta_cm_2007.Properties.VariableNames, ...
                          {'series_name','include_main_2007'});
    for k = 1:numel(vars_cm2007)
        series_meta.(vars_cm2007{k}) = series_meta_cm_2007.(vars_cm2007{k});
    end
    fprintf('[run_tables_figures] Loaded cm_2007_results.mat (%d fields)\n', numel(vars_cm2007));
end

% Join CV 2007 (if Stage 3e has run)
cv2007_path = fullfile(DATA_INT,'cv_2007_results.mat');
if exist(cv2007_path,'file')
    load(cv2007_path, 'series_meta_cv_2007');
    vars_cv2007 = setdiff(series_meta_cv_2007.Properties.VariableNames, ...
                          {'series_name','include_main_2007'});
    for k = 1:numel(vars_cv2007)
        series_meta.(vars_cv2007{k}) = series_meta_cv_2007.(vars_cv2007{k});
    end
    fprintf('[run_tables_figures] Loaded cv_2007_results.mat (%d fields)\n', numel(vars_cv2007));
end

% Join CM full (if Stage 3f has run)
cm_full_path = fullfile(DATA_INT,'cm_full_results.mat');
if exist(cm_full_path,'file')
    load(cm_full_path, 'series_meta_cm_full');
    vars_cmfull = setdiff(series_meta_cm_full.Properties.VariableNames, ...
                          {'series_name','include_main_full_ext'});
    for k = 1:numel(vars_cmfull)
        series_meta.(vars_cmfull{k}) = series_meta_cm_full.(vars_cmfull{k});
    end
    fprintf('[run_tables_figures] Loaded cm_full_results.mat (%d fields)\n', numel(vars_cmfull));
end

% Join CV full (if Stage 3g has run)
cv_full_path = fullfile(DATA_INT,'cv_full_results.mat');
if exist(cv_full_path,'file')
    load(cv_full_path, 'series_meta_cv_full');
    vars_cvfull = setdiff(series_meta_cv_full.Properties.VariableNames, ...
                          {'series_name','include_main_full_ext'});
    for k = 1:numel(vars_cvfull)
        series_meta.(vars_cvfull{k}) = series_meta_cv_full.(vars_cvfull{k});
    end
    fprintf('[run_tables_figures] Loaded cv_full_results.mat (%d fields)\n', numel(vars_cvfull));
end

fprintf('[run_tables_figures] Merged table: %d rows x %d columns\n', ...
    height(series_meta), width(series_meta));

%% -----------------------------------------------------------------------
% 2b. JOIN SMOOTH-BREAK RESULTS (Stages 3k, 3l)
%% -----------------------------------------------------------------------
smooth_files = { ...
    'uv_smooth_hist_results.mat',  'series_meta_uv_smooth_hist',  'include_main_hist'; ...
    'uv_smooth_2007_results.mat',  'series_meta_uv_smooth_2007',  'include_main_2007'; ...
    'uv_smooth_full_results.mat',  'series_meta_uv_smooth_full',  'include_main_full_ext'; ...
    'cv_smooth_hist_results.mat',  'series_meta_cv_smooth_hist',  'include_main_hist'; ...
    'cv_smooth_2007_results.mat',  'series_meta_cv_smooth_2007',  'include_main_2007'; ...
    'cv_smooth_full_results.mat',  'series_meta_cv_smooth_full',  'include_main_full_ext'};

for sf = 1:size(smooth_files,1)
    spath = fullfile(DATA_INT, smooth_files{sf,1});
    svar  = smooth_files{sf,2};
    excl  = {'series_name', smooth_files{sf,3}};
    if exist(spath,'file')
        tmp = load(spath, svar);
        sm  = tmp.(svar);
        vars_sm = setdiff(sm.Properties.VariableNames, excl);
        for k = 1:numel(vars_sm)
            series_meta.(vars_sm{k}) = sm.(vars_sm{k});
        end
        fprintf('[run_tables_figures] Loaded %s (%d fields)\n', smooth_files{sf,1}, numel(vars_sm));
    else
        fprintf('[run_tables_figures] %s not found — skipping.\n', smooth_files{sf,1});
    end
end

%% -----------------------------------------------------------------------
% 2f. DERIVE PER-SERIES DIRECTIONALITY COLUMNS
%
%   UV direction: based on SupW break ratio (post/pre).
%   Assigned only for series with a significant break (p<0.05).
%   CV direction: based on cv_ratio (post/pre mean of s_t).
%   Logic matches termpaper610.m (ratio < 1 → "down", > 1 → "up").
%% -----------------------------------------------------------------------
n = height(series_meta);
vars = series_meta.Properties.VariableNames;

derive_direction = @(sig, ratio) ...
    deal(sig & ratio < 1, sig & ratio > 1, sig & abs(ratio - 1) < 1e-12);

% UV historical
if ismember('uv_hansen_sig5_hist',vars) && ismember('uv_supw_ratio_hist',vars)
    sig   = series_meta.uv_hansen_sig5_hist;
    ratio = series_meta.uv_supw_ratio_hist;
    dir   = strings(n,1);
    [dn,up,fl] = derive_direction(sig, ratio);
    dir(dn) = "down"; dir(up) = "up"; dir(fl) = "flat";
    series_meta.uv_direction_hist = dir;
end

% UV 2007
if ismember('uv2007_hansen_sig5',vars) && ismember('uv2007_supw_ratio',vars)
    sig   = series_meta.uv2007_hansen_sig5;
    ratio = series_meta.uv2007_supw_ratio;
    dir   = strings(n,1);
    [dn,up,fl] = derive_direction(sig, ratio);
    dir(dn) = "down"; dir(up) = "up"; dir(fl) = "flat";
    series_meta.uv2007_direction = dir;
end

% UV full
if ismember('uvfull_hansen_sig5',vars) && ismember('uvfull_supw_ratio',vars)
    sig   = series_meta.uvfull_hansen_sig5;
    ratio = series_meta.uvfull_supw_ratio;
    dir   = strings(n,1);
    [dn,up,fl] = derive_direction(sig, ratio);
    dir(dn) = "down"; dir(up) = "up"; dir(fl) = "flat";
    series_meta.uvfull_direction = dir;
end

% CV historical
if ismember('cv_hansen_sig5_hist',vars) && ismember('cv_ratio_hist',vars)
    sig   = series_meta.cv_hansen_sig5_hist;
    ratio = series_meta.cv_ratio_hist;
    dir   = strings(n,1);
    [dn,up,fl] = derive_direction(sig, ratio);
    dir(dn) = "down"; dir(up) = "up"; dir(fl) = "flat";
    series_meta.cv_direction_hist = dir;
end

% CV 2007
if ismember('cv2007_hansen_sig5',vars) && ismember('cv2007_ratio',vars)
    sig   = series_meta.cv2007_hansen_sig5;
    ratio = series_meta.cv2007_ratio;
    dir   = strings(n,1);
    [dn,up,fl] = derive_direction(sig, ratio);
    dir(dn) = "down"; dir(up) = "up"; dir(fl) = "flat";
    series_meta.cv2007_direction = dir;
end

% CV full
if ismember('cvfull_hansen_sig5',vars) && ismember('cvfull_ratio',vars)
    sig   = series_meta.cvfull_hansen_sig5;
    ratio = series_meta.cvfull_ratio;
    dir   = strings(n,1);
    [dn,up,fl] = derive_direction(sig, ratio);
    dir(dn) = "down"; dir(up) = "up"; dir(fl) = "flat";
    series_meta.cvfull_direction = dir;
end

fprintf('[run_tables_figures] Derived directionality columns (uv/cv _hist/2007/full)\n');

%% -----------------------------------------------------------------------
% 2e. ATTACH fred_group
%% -----------------------------------------------------------------------
% fred_group is a 1×N string array saved by run_data_prep.m.
% Attach it as a column so tables_figures_sensier_fredmd.m can group results
% by FRED-MD category (e.g. "Labor Market", "Housing / Construction").
if exist('fred_group','var') && isstring(fred_group)
    if numel(fred_group) == height(series_meta)
        series_meta.fred_group = fred_group(:);
    else
        warning('[run_tables_figures] fred_group length (%d) != series_meta rows (%d) — not attached.', ...
            numel(fred_group), height(series_meta));
    end
end

%% -----------------------------------------------------------------------
% 3. SAVE CONSOLIDATED series_meta
%% -----------------------------------------------------------------------
save(fullfile(DATA_INT,'series_meta.mat'), 'series_meta');

fprintf('[run_tables_figures] Saved series_meta.mat\n');

%% -----------------------------------------------------------------------
% 4. DONE — series_meta.mat saved above; figures produced by produce_paper_figures.m
%% -----------------------------------------------------------------------

% tables_figures_sensier_fredmd.m requires series_meta.fred_group.
% fred_group is loaded by run_data_prep.m from series_metadata_original_order.csv.
% Guard: if the column is absent (file was missing during data prep), skip gracefully.
if ~ismember('fred_group', series_meta.Properties.VariableNames)
    warning('[run_tables_figures] series_meta.fred_group not found — Stage 4 presentation skipped.');
    fprintf('[run_tables_figures] Ensure series_metadata_original_order.csv is in the data_raw/ directory\n');
    fprintf('[run_tables_figures] and re-run master_run from Stage 1, or re-run run_data_prep.m alone then this script.\n');
    return;
end

% tables_figures_sensier_fredmd.m is not included in this replication package.
% All paper figures are produced by produce_paper_figures.m (Stage 4b).
% run(fullfile(fileparts(PKG_ROOT), 'tables_figures_sensier_fredmd.m'));

fprintf('\n[run_tables_figures] Done.\n');
