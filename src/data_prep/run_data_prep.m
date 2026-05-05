%% run_data_prep.m
% =========================================================================
% DATA PREPARATION
%
% Loads the FRED-MD vintage, applies stationarity transformations,
% removes outliers, and saves the processed panel for downstream scripts.
%
% INPUTS (from config.m)
%   CSV_IN     - path to FRED-MD vintage CSV
%   HIST_START, HIST_END, WIN_2007_END - sample window endpoints
%   MAIN_THRESH, BALANCED_THRESH, BALANCED_START_CUTOFF - inclusion rules
%   DATA_INT   - folder for saved intermediate files
%
% OUTPUTS (saved to data_intermediate/processed_data.mat)
%   yt           - full transformed panel (T_full x N), starts 1959:03
%   data         - yt with outliers removed (for factor estimation)
%   dates_dt     - datetime vector aligned with yt
%   series       - {1xN} cell of series names
%   tcode        - 1xN transformation codes
%   raw_dates_dt - datetime vector for rawdata rows
%   rawdata      - raw untransformed panel
%
%   yt_rep       - historical window panel (1959:03-1999:12)
%   dates_rep_t  - datetime vector for yt_rep
%
%   yt_ext_2007  - 2007-window panel (1959:03-2007:12)
%   dates_ext_2007
%
%   yt_ext_full  - full-extent panel (1959:03-latest)
%   dates_ext_full
%
%   include_main_hist, include_bal_hist  - historical inclusion flags
%   include_main_2007, include_bal_2007  - 2007-window inclusion flags
%   include_main_full_ext                - full-extent inclusion flag
% =========================================================================

run('config.m');

%% -----------------------------------------------------------------------
% 1. LOAD RAW DATA
%% -----------------------------------------------------------------------
fprintf('\n[run_data_prep] Loading %s\n', CSV_IN);
dum = importdata(CSV_IN, ',');

series  = dum.textdata(1, 2:end);
tcode   = dum.data(1,:);
rawdata = dum.data(2:end,:);

raw_dates_dt = datetime(dum.textdata(3:end,1), 'InputFormat','M/d/yyyy');

final_datevec = datevec(dum.textdata(end,1));
final_month   = final_datevec(2);
final_year    = final_datevec(1);

dates = (1959+1/12 : 1/12 : final_year+final_month/12)';
T     = size(dates,1);
rawdata = rawdata(1:T,:);

%% -----------------------------------------------------------------------
% 2. TRANSFORM TO STATIONARITY
%% -----------------------------------------------------------------------
yt = prepare_missing(rawdata, tcode);

% Remove first two rows (some series first-differenced)
yt          = yt(3:T,:);
dates       = dates(3:T,:);
dates_dt    = raw_dates_dt(3:end);

%% -----------------------------------------------------------------------
% 3. REMOVE OUTLIERS
%% -----------------------------------------------------------------------
[data, ~] = remove_outliers(yt);

%% -----------------------------------------------------------------------
% 4. HISTORICAL WINDOW  (HIST_START : HIST_END)
%% -----------------------------------------------------------------------
rep_idx    = (dates_dt >= HIST_START) & (dates_dt <= HIST_END);
yt_rep     = yt(rep_idx,:);
dates_rep_t = dates_dt(rep_idx);

T_rep        = size(yt_rep,1);
usable_rep   = sum(~isnan(yt_rep),1);
main_min_obs = ceil(MAIN_THRESH    * T_rep);
bal_min_obs  = ceil(BALANCED_THRESH * T_rep);

first_valid_rep = NaT(1,size(yt_rep,2));
for j = 1:size(yt_rep,2)
    idx = find(~isnan(yt_rep(:,j)),1,'first');
    if ~isempty(idx); first_valid_rep(j) = dates_rep_t(idx); end
end

include_main_hist = usable_rep >= main_min_obs;
include_bal_hist  = (usable_rep >= bal_min_obs) & ...
                    (first_valid_rep <= BALANCED_START_CUTOFF);

fprintf('[run_data_prep] Historical window: %s – %s  (%d rows)\n', ...
    string(HIST_START,'yyyy-MM'), string(HIST_END,'yyyy-MM'), T_rep);
fprintf('  Series in main historical sample:     %d\n', sum(include_main_hist));
fprintf('  Series in balanced historical sample: %d\n', sum(include_bal_hist));

%% -----------------------------------------------------------------------
% 5. 2007 WINDOW  (HIST_START : WIN_2007_END)
%% -----------------------------------------------------------------------
idx_2007 = (dates_dt >= HIST_START) & (dates_dt <= WIN_2007_END);
yt_ext_2007    = yt(idx_2007,:);
dates_ext_2007 = dates_dt(idx_2007);

T_2007       = size(yt_ext_2007,1);
usable_2007  = sum(~isnan(yt_ext_2007),1);
main_min_2007 = ceil(MAIN_THRESH    * T_2007);
bal_min_2007  = ceil(BALANCED_THRESH * T_2007);

first_valid_2007 = NaT(1,size(yt_ext_2007,2));
for j = 1:size(yt_ext_2007,2)
    idx = find(~isnan(yt_ext_2007(:,j)),1,'first');
    if ~isempty(idx); first_valid_2007(j) = dates_ext_2007(idx); end
end

include_main_2007 = usable_2007 >= main_min_2007;
include_bal_2007  = (usable_2007 >= bal_min_2007) & ...
                    (first_valid_2007 <= BALANCED_START_CUTOFF);

fprintf('[run_data_prep] 2007 window: %s – %s  (%d rows)\n', ...
    string(HIST_START,'yyyy-MM'), string(WIN_2007_END,'yyyy-MM'), T_2007);
fprintf('  Series in main 2007 sample:     %d\n', sum(include_main_2007));
fprintf('  Series in balanced 2007 sample: %d\n', sum(include_bal_2007));

%% -----------------------------------------------------------------------
% 6. FULL-EXTENT WINDOW  (HIST_START : latest)
%% -----------------------------------------------------------------------
idx_full = (dates_dt >= HIST_START);
yt_ext_full    = yt(idx_full,:);
dates_ext_full = dates_dt(idx_full);

T_full       = size(yt_ext_full,1);
usable_full  = sum(~isnan(yt_ext_full),1);
main_min_full = ceil(MAIN_THRESH * T_full);

first_valid_full = NaT(1,size(yt_ext_full,2));
for j = 1:size(yt_ext_full,2)
    idx = find(~isnan(yt_ext_full(:,j)),1,'first');
    if ~isempty(idx); first_valid_full(j) = dates_ext_full(idx); end
end

include_main_full_ext = usable_full >= main_min_full;

fprintf('[run_data_prep] Full-extent window: %s – %s  (%d rows)\n', ...
    string(HIST_START,'yyyy-MM'), string(dates_ext_full(end),'yyyy-MM'), T_full);
fprintf('  Series in main full-extent sample: %d\n', sum(include_main_full_ext));

%% -----------------------------------------------------------------------
% 6b. MODERN-ERA APPENDIX WINDOW  (MOD_ERA_START : latest)
%
%   Appendix robustness window to study post-1999 volatility dynamics free
%   from the many embedded regimes in the full 1959-present sample.
%   Start date: 1999:01 (first row in yt where dates_dt >= MOD_ERA_START).
%   End date:   last available month in the current FRED-MD vintage.
%
%   Inclusion threshold: same MAIN_THRESH = 0.80 as other windows.
%   No balanced-panel criterion is computed for this window (the exercise
%   does not require factor estimation over the modern-era subsample).
%% -----------------------------------------------------------------------
idx_mod          = (dates_dt >= MOD_ERA_START);
yt_mod           = yt(idx_mod,:);
dates_mod        = dates_dt(idx_mod);

T_mod            = size(yt_mod,1);
usable_mod       = sum(~isnan(yt_mod),1);
main_min_mod     = ceil(MAIN_THRESH * T_mod);
include_main_mod = usable_mod >= main_min_mod;

fprintf('[run_data_prep] Modern-era window: %s – %s  (%d rows)\n', ...
    string(MOD_ERA_START,'yyyy-MM'), string(dates_mod(end),'yyyy-MM'), T_mod);
fprintf('  Series in main modern-era sample: %d\n', sum(include_main_mod));
fprintf('  (min_seg_len_eff for MB = max(ceil(0.15×%d),60) = %d months)\n', ...
    T_mod, max(ceil(0.15*T_mod), 60));

%% -----------------------------------------------------------------------
% 6c. FRED-MD GROUP LABELS
%
%   Loads the FRED-MD native group label for each series
%   (e.g. "Output and Income", "Labor Market", "Housing / Construction")
%   from series_metadata_original_order.csv in the same directory as
%   the input CSV.  Required by Stage 4 (tables_figures_sensier_fredmd.m).
%
%   If the file is absent, fred_group is saved as an array of empty
%   strings and Stage 4 will issue a warning instead of generating
%   grouped figures.
%% -----------------------------------------------------------------------
fred_group = repmat("", 1, length(series));   % default: all empty

META_CSV = fullfile(fileparts(CSV_IN), 'series_metadata_original_order.csv');
if isfile(META_CSV)
    try
        meta_tbl = readtable(META_CSV, 'Delimiter', ',');
        if ismember('series_name', meta_tbl.Properties.VariableNames) && ...
           ismember('fred_group',  meta_tbl.Properties.VariableNames)
            meta_sname  = string(meta_tbl.series_name);
            meta_fgroup = string(meta_tbl.fred_group);
            for jj = 1:length(series)
                hit = find(meta_sname == string(series{jj}), 1);
                if ~isempty(hit)
                    fred_group(jj) = meta_fgroup(hit);
                end
            end
            fprintf('[run_data_prep] fred_group loaded: %d of %d series matched.\n', ...
                sum(fred_group ~= ""), length(series));
        else
            warning('[run_data_prep] series_metadata_original_order.csv is missing series_name or fred_group column.');
        end
    catch ME
        warning('[run_data_prep] Could not read series_metadata_original_order.csv: %s', ME.message);
    end
else
    fprintf('[run_data_prep] NOTE: series_metadata_original_order.csv not found in %s\n', fileparts(CSV_IN));
    fprintf('[run_data_prep]       fred_group will be empty — Stage 4 grouped figures will not be generated.\n');
end

%% -----------------------------------------------------------------------
% 7. SAVE PROCESSED DATA
%% -----------------------------------------------------------------------
if ~exist(DATA_INT,'dir'); mkdir(DATA_INT); end

save(fullfile(DATA_INT,'processed_data.mat'), ...
    'yt', 'data', 'dates_dt', 'series', 'tcode', 'raw_dates_dt', 'rawdata', ...
    'yt_rep',      'dates_rep_t', ...
    'yt_ext_2007', 'dates_ext_2007', ...
    'yt_ext_full', 'dates_ext_full', ...
    'yt_mod',      'dates_mod', ...
    'include_main_hist', 'include_bal_hist', ...
    'include_main_2007', 'include_bal_2007', ...
    'include_main_full_ext', ...
    'include_main_mod', ...
    'fred_group');

fprintf('\n[run_data_prep] Saved processed_data.mat to %s\n', DATA_INT);
fprintf('[run_data_prep] Done.\n\n');
