%% produce_paper_figures.m
% =========================================================================
% PAPER FIGURES — FRED-MD Volatility Break Replication
%
% All figures and Table 4 use the fixed-74-month min-seg, 5-break-cap spec
% (mb_fixed74_results.mat / mb_fixed74_group_summary.csv).
%
%   Figure 1   : Distribution of All MB Break Dates (Historical Window)
%   Figure 2   : Median Recursive Breaks by Sector and Sample Window
%   Figure 3   : Share of Break Transitions Moderating vs. Amplifying
%   Figure 4a  : Mean Break Magnitude by Year — UV, Real vs. Nominal
%   Figure 4b  : Mean Break Magnitude by Year — CV, Real vs. Nominal
%   Figure A1  : Output and Income      — Break Magnitude by Year (Full)
%   Figure A2  : Labor Market           — Break Magnitude by Year (Full)
%   Figure A3  : Orders and Inventories — Break Magnitude by Year (Full)
%   Figure A4  : Housing / Construction — Break Magnitude by Year (Full)
%   Figure A5  : Prices                 — Break Magnitude by Year (Full)
%   Figure A6  : Money and Credit       — Break Magnitude by Year (Full)
%   Figure A7  : Interest & Exch. Rates — Break Magnitude by Year (Full)
%   Table 4    : Sector summary with median breaks column added
%
% Outputs saved to output/figures/ (PNG 300 dpi, PDF, .fig)
%          and output/tables/ (table4_updated.csv)
%
% PREREQUISITES: run_mb_fixed74.m, run_uv_fullsample, run_cm_fullsample
% =========================================================================

run('config.m');
if ~exist(OUT_FIGS,   'dir'); mkdir(OUT_FIGS);   end
if ~exist(OUT_TABLES, 'dir'); mkdir(OUT_TABLES); end

fprintf('\n[produce_paper_figures] Starting (fixed-74 spec)...\n');

% =========================================================================
% 0. LOAD DATA
% =========================================================================

fprintf('[produce_paper_figures] Loading data...\n');

load(fullfile(DATA_INT,'processed_data.mat'), ...
    'dates_rep_t','dates_ext_full','series', ...
    'include_main_hist','include_main_2007','include_main_full_ext', ...
    'fred_group');
load(fullfile(DATA_INT,'uv_full_results.mat'), 'vol_obj_full');
load(fullfile(DATA_INT,'mb_fixed74_results.mat'), 'series_meta_mb_f74');

n_series   = length(series);
fred_group = string(fred_group(:));

% Nominal vs Real classification
real_groups    = ["Housing / Construction","Labor Market", ...
                  "Orders and Inventories","Output and Income"];
nominal_groups = ["Interest and Exchange Rates","Money and Credit", ...
                  "Prices","Stock Market"];

% Groups and labels for sector figures (7 sectors, excluding Stock Market)
sector_groups = ["Output and Income","Labor Market","Orders and Inventories", ...
                 "Housing / Construction","Prices","Money and Credit", ...
                 "Interest and Exchange Rates"];
sector_labels = {'Output & Income','Labor Market','Orders & Inv.', ...
                 'Housing','Prices','Money & Credit','Interest & FX'};
n_sec = numel(sector_groups);

% Figure style constants
fs_title  = 13;
fs_axis   = 11;
fs_tick   = 9;
fs_legend = 9;

clr_uv  = [0.20 0.40 0.75];   % blue  — UV
clr_cv  = [0.85 0.52 0.10];   % amber — CV
clr_gfc = [0.90 0.70 0.70];   % light red  — post-GFC shading
clr_cov = [0.80 0.70 0.90];   % light purple — post-COVID shading

% =========================================================================
% 1. COMPUTE PER-BREAK LOG SIGMA RATIOS (FULL SAMPLE, UV + CV)
%    Used by Figures 4a/4b and A1-A7.
% =========================================================================

fprintf('[produce_paper_figures] Computing per-break sigma ratios (full sample)...\n');

uv_break_years  = [];
uv_break_logsig = [];
uv_break_grp    = strings(0,1);

for j = 1:n_series
    if ~include_main_full_ext(j); continue; end
    zt   = vol_obj_full(:, j);
    bstr = string(series_meta_mb_f74.uvfull_mb_break_dates(j));
    if bstr == ""; continue; end
    [br_yrs, br_lsig] = parse_break_sigma(zt, dates_ext_full, bstr);
    if isempty(br_yrs); continue; end
    uv_break_years  = [uv_break_years;  br_yrs(:)];
    uv_break_logsig = [uv_break_logsig; br_lsig(:)];
    uv_break_grp    = [uv_break_grp;    repmat(fred_group(j), numel(br_yrs), 1)];
end

cv_break_years  = [];
cv_break_logsig = [];
cv_break_grp    = strings(0,1);

cm_full_path = fullfile(DATA_INT,'cm_full_results.mat');
if exist(cm_full_path,'file')
    load(cm_full_path,'cm_supw_resid_break_full','series_meta_cm_full');
    cmfull_p     = series_meta_cm_full.cmfull_p;
    cmfull_found = series_meta_cm_full.cmfull_supw_found;

    for j = 1:n_series
        if ~include_main_full_ext(j); continue; end
        if ~cmfull_found(j) || isempty(cm_supw_resid_break_full{j}); continue; end
        pj      = cmfull_p(j);
        ehat_j  = cm_supw_resid_break_full{j};
        s_j     = sqrt(pi/2) * abs(ehat_j);
        Tj      = length(s_j);
        dates_j = dates_ext_full((pj+1):end);
        if length(dates_j) > Tj; dates_j = dates_j(end-Tj+1:end); end
        bstr = string(series_meta_mb_f74.cvfull_mb_break_dates(j));
        if bstr == ""; continue; end
        [br_yrs, br_lsig] = parse_break_sigma(s_j, dates_j, bstr);
        if isempty(br_yrs); continue; end
        cv_break_years  = [cv_break_years;  br_yrs(:)];
        cv_break_logsig = [cv_break_logsig; br_lsig(:)];
        cv_break_grp    = [cv_break_grp;    repmat(fred_group(j), numel(br_yrs), 1)];
    end
else
    fprintf('[produce_paper_figures] cm_full_results.mat not found — CV magnitude unavailable.\n');
end

fprintf('[produce_paper_figures] UV full: %d break transitions; CV full: %d.\n', ...
    numel(uv_break_years), numel(cv_break_years));

% =========================================================================
% FIGURE 1: Distribution of All MB Break Dates, Historical Window
% =========================================================================

fprintf('[produce_paper_figures] Figure 1...\n');

% Collect every break date (all recursive breaks) for included series
uv_hist_all_yrs = [];
cv_hist_all_yrs = [];

for j = 1:n_series
    if ~include_main_hist(j); continue; end
    bstr_uv = string(series_meta_mb_f74.uv_mb_break_dates_hist(j));
    if bstr_uv ~= ""
        uv_hist_all_yrs = [uv_hist_all_yrs; parse_break_years(bstr_uv)]; %#ok<AGROW>
    end
    bstr_cv = string(series_meta_mb_f74.cv_mb_break_dates_hist(j));
    if bstr_cv ~= ""
        cv_hist_all_yrs = [cv_hist_all_yrs; parse_break_years(bstr_cv)]; %#ok<AGROW>
    end
end

bin_edges = 1959:2:2001;   % 2-year bins; last edge at 2001 keeps 1999-2000 full

fig1 = figure('Color','w','Position',[50 50 880 380],'Visible','off');
ax1  = axes('Parent',fig1);
hold(ax1,'on');

% Offset UV left and CV right so bars don't overlap
if ~isempty(uv_hist_all_yrs)
    [uv_cts, ~] = histcounts(uv_hist_all_yrs, bin_edges);
    bin_ctrs    = bin_edges(1:end-1) + 0.5;   % bin centres
    bar(ax1, bin_ctrs - 0.45, uv_cts, 0.85, ...
        'FaceColor', clr_uv, 'FaceAlpha', 0.85, 'EdgeColor','w', ...
        'DisplayName', sprintf('UV (%d breaks)', numel(uv_hist_all_yrs)));
end
if ~isempty(cv_hist_all_yrs)
    [cv_cts, ~] = histcounts(cv_hist_all_yrs, bin_edges);
    bin_ctrs    = bin_edges(1:end-1) + 0.5;
    bar(ax1, bin_ctrs + 0.45, cv_cts, 0.85, ...
        'FaceColor', clr_cv, 'FaceAlpha', 0.85, 'EdgeColor','w', ...
        'DisplayName', sprintf('CV (%d breaks)', numel(cv_hist_all_yrs)));
end

xline(ax1, 1984, '--k', 'LineWidth',1.4, 'Label','1984', ...
    'HandleVisibility','off', 'LabelVerticalAlignment','middle');
hold(ax1,'off');

xlim(ax1, [1958 2001]);
xlabel(ax1, 'Break Year', 'FontSize', fs_axis);
ylabel(ax1, 'Number of Breaks', 'FontSize', fs_axis);
title(ax1, 'Figure 1: Distribution of Estimated Volatility Break Dates, Historical Window (1959–1999)', ...
    'FontSize', fs_title, 'FontWeight','bold');
legend(ax1, 'show', 'Location','northwest', 'FontSize', fs_legend);
set(ax1, 'FontSize', fs_tick);
grid(ax1,'on'); box(ax1,'on');

save_paper_fig(fig1, OUT_FIGS, 'figure1_breakdate_hist');

% =========================================================================
% FIGURE 2: Median Recursive Breaks by Sector and Sample Window
% =========================================================================

fprintf('[produce_paper_figures] Figure 2...\n');

mb_sum = readtable(fullfile(OUT_TABLES,'mb_fixed74_group_summary.csv'), ...
    'TextType','string');

win_ordered = ["1959:03-1999:12","1959:03-2007:12","1959:03-2026:01"];
win_labels  = {'Historical (1959–99)','Pre-GFC (1959–2007)','Full (1959–2026)'};
win_colors  = [0.25 0.45 0.75; 0.85 0.50 0.10; 0.35 0.65 0.35];

uv_med_mat = NaN(n_sec, 3);
cv_med_mat = NaN(n_sec, 3);

for si = 1:n_sec
    for wi = 1:3
        r_uv = mb_sum.window == win_ordered(wi) & ...
               mb_sum.fred_group == sector_groups(si) & ...
               mb_sum.measure == "UV";
        r_cv = mb_sum.window == win_ordered(wi) & ...
               mb_sum.fred_group == sector_groups(si) & ...
               mb_sum.measure == "CV";
        if any(r_uv), uv_med_mat(si,wi) = mb_sum.median_breaks(r_uv); end
        if any(r_cv), cv_med_mat(si,wi) = mb_sum.median_breaks(r_cv); end
    end
end

x_pos   = 1:n_sec;
bar_w   = 0.25;
offsets = [-bar_w, 0, bar_w];
top_val = max(max([uv_med_mat; cv_med_mat],[],'all','omitnan'), 1);

fig2 = figure('Color','w','Position',[50 50 1100 700],'Visible','off');

for panel = 1:2
    subplot(2,1,panel);
    if panel == 1, med_mat = uv_med_mat; ttl = 'UV — Unconditional Volatility';
    else,          med_mat = cv_med_mat; ttl = 'CV — Conditional Volatility';
    end
    hold on;
    for wi = 1:3
        bar(x_pos + offsets(wi), med_mat(:,wi), bar_w*0.92, ...
            'FaceColor', win_colors(wi,:), 'EdgeColor','w', 'LineWidth',0.5, ...
            'DisplayName', win_labels{wi});
    end
    hold off;
    set(gca,'XTick',x_pos,'XTickLabel',sector_labels,'FontSize',fs_tick);
    xtickangle(30);
    ylabel('Median Number of Breaks','FontSize',fs_axis);
    ylim([0 top_val + 0.5]);
    legend('show','Location','best','FontSize',fs_legend);
    title(ttl,'FontSize',fs_title,'FontWeight','bold');
    grid on; box on;
end

sgtitle('Figure 2: Median Recursive Breaks by Sector and Sample Window', ...
    'FontSize',fs_title,'FontWeight','bold');

save_paper_fig(fig2, OUT_FIGS, 'figure2_median_breaks');

% =========================================================================
% FIGURE 3: Share Moderating vs. Amplifying by Sector (Full Sample)
% =========================================================================

fprintf('[produce_paper_figures] Figure 3...\n');

full_win = "1959:03-2026:01";
uv_mod = NaN(n_sec,1); uv_amp = NaN(n_sec,1);
cv_mod = NaN(n_sec,1); cv_amp = NaN(n_sec,1);

for si = 1:n_sec
    r_uv = mb_sum.window == full_win & mb_sum.fred_group == sector_groups(si) & mb_sum.measure == "UV";
    r_cv = mb_sum.window == full_win & mb_sum.fred_group == sector_groups(si) & mb_sum.measure == "CV";
    if any(r_uv)
        uv_mod(si) = mb_sum.share_moderating(r_uv);
        uv_amp(si) = 1 - mb_sum.share_moderating(r_uv);
    end
    if any(r_cv)
        cv_mod(si) = mb_sum.share_moderating(r_cv);
        cv_amp(si) = 1 - mb_sum.share_moderating(r_cv);
    end
end

clr_mod = [0.25 0.60 0.35];   % green — moderating
clr_amp = [0.80 0.25 0.20];   % red   — amplifying

fig3 = figure('Color','w','Position',[50 50 1100 550],'Visible','off');

for panel = 1:2
    subplot(1,2,panel);
    if panel == 1, mod_v = uv_mod; amp_v = uv_amp; ttl = 'UV';
    else,          mod_v = cv_mod; amp_v = cv_amp; ttl = 'CV';
    end
    b = bar(x_pos, [mod_v, amp_v], 0.65, 'stacked');
    b(1).FaceColor = clr_mod; b(1).EdgeColor = 'w';
    b(2).FaceColor = clr_amp; b(2).EdgeColor = 'w';
    b(1).DisplayName = 'Moderating';
    b(2).DisplayName = 'Amplifying';
    yline(0.5,'--k','LineWidth',1.2,'HandleVisibility','off');
    set(gca,'XTick',x_pos,'XTickLabel',sector_labels,'FontSize',fs_tick);
    xtickangle(30);
    ylabel('Share of Break Transitions','FontSize',fs_axis);
    ylim([0 1]);
    legend({'Moderating','Amplifying'},'Location','best','FontSize',fs_legend);
    title([ttl ' — Share Moderating vs. Amplifying'],'FontSize',fs_title,'FontWeight','bold');
    grid on; box on;
end

sgtitle('Figure 3: Share of Break Transitions Moderating vs. Amplifying by Sector (Full Sample)', ...
    'FontSize',fs_title,'FontWeight','bold');

save_paper_fig(fig3, OUT_FIGS, 'figure3_share_moderating');

% =========================================================================
% FIGURES 4a & 4b: Mean Break Magnitude by Year — Real vs. Nominal
%                  Full Sample, UV and CV
% =========================================================================

fprintf('[produce_paper_figures] Figures 4a and 4b...\n');

% Tag each break record as Real / Nominal
uv_nr = repmat("Unassigned", numel(uv_break_years), 1);
for ii = 1:numel(uv_break_years)
    if     any(real_groups    == uv_break_grp(ii)), uv_nr(ii) = "Real";
    elseif any(nominal_groups == uv_break_grp(ii)), uv_nr(ii) = "Nominal";
    end
end
cv_nr = repmat("Unassigned", numel(cv_break_years), 1);
for ii = 1:numel(cv_break_years)
    if     any(real_groups    == cv_break_grp(ii)), cv_nr(ii) = "Real";
    elseif any(nominal_groups == cv_break_grp(ii)), cv_nr(ii) = "Nominal";
    end
end

yr_range = 1959:2026;

fig4_specs = { ...
    uv_break_years, uv_break_logsig, uv_nr, ...
    'figure4a_magnitude_uv', ...
    'Figure 4a: Mean Break Magnitude by Year — Unconditional Volatility (Full Sample)'; ...
    cv_break_years, cv_break_logsig, cv_nr, ...
    'figure4b_magnitude_cv', ...
    'Figure 4b: Mean Break Magnitude by Year — Conditional Volatility (Full Sample)'};

for panel = 1:2
    brk_yrs   = fig4_specs{panel,1};
    brk_lsig  = fig4_specs{panel,2};
    brk_nr    = fig4_specs{panel,3};
    fig_fname = fig4_specs{panel,4};
    fig_title = fig4_specs{panel,5};

    real_yr_mean = NaN(numel(yr_range),1);
    nom_yr_mean  = NaN(numel(yr_range),1);
    for yi = 1:numel(yr_range)
        yr = yr_range(yi);
        rv = brk_lsig(brk_yrs == yr & brk_nr == "Real");
        nv = brk_lsig(brk_yrs == yr & brk_nr == "Nominal");
        if ~isempty(rv), real_yr_mean(yi) = nanmean(rv); end
        if ~isempty(nv), nom_yr_mean(yi)  = nanmean(nv); end
    end

    has_data = ~isnan(real_yr_mean) | ~isnan(nom_yr_mean);
    yr_plt   = yr_range(has_data);
    real_plt = real_yr_mean(has_data);
    nom_plt  = nom_yr_mean(has_data);

    fig_h = figure('Color','w','Position',[50 50 1100 420],'Visible','off');
    ax    = axes('Parent',fig_h);
    hold(ax,'on');

    bw = 0.40;
    bar(ax, yr_plt - bw/2 - 0.05, real_plt, bw, ...
        'FaceColor',[0.25 0.45 0.75],'EdgeColor','none','DisplayName','Real');
    bar(ax, yr_plt + bw/2 + 0.05, nom_plt,  bw, ...
        'FaceColor',[0.85 0.35 0.25],'EdgeColor','none','DisplayName','Nominal');

    fill_between(ax, 2008, 2019,            clr_gfc);
    fill_between(ax, 2020, max(yr_range)+1, clr_cov);

    yline(ax, 0, '-k', 'LineWidth',0.8, 'HandleVisibility','off');
    xline(ax, 1984, '--k','LineWidth',1.0,'Label','1984', ...
        'HandleVisibility','off','LabelHorizontalAlignment','left');
    xline(ax, 2008, '--k','LineWidth',1.0,'Label','2008', ...
        'HandleVisibility','off','LabelHorizontalAlignment','left');
    xline(ax, 2020, '--k','LineWidth',1.0,'Label','2020', ...
        'HandleVisibility','off','LabelHorizontalAlignment','left');

    patch(ax,NaN,NaN,clr_gfc,'FaceAlpha',1,'EdgeColor','none','DisplayName','Post-GFC (2008–2019)');
    patch(ax,NaN,NaN,clr_cov,'FaceAlpha',1,'EdgeColor','none','DisplayName','Post-COVID (2020–)');

    hold(ax,'off');

    all_v = [real_plt(~isnan(real_plt)); nom_plt(~isnan(nom_plt))];
    if ~isempty(all_v)
        ypad = max(0.1, range(all_v)*0.15);
        ylim(ax, [min(all_v)-ypad, max(all_v)+ypad]);
    end

    legend(ax,'show','Location','best','FontSize',fs_legend);
    xlabel(ax,'Year','FontSize',fs_axis);
    ylabel(ax,'Mean Log \sigma-Ratio','FontSize',fs_axis,'Interpreter','tex');
    title(ax, fig_title,'FontSize',fs_title,'FontWeight','bold');
    xlim(ax,[1959 2027]); grid(ax,'on'); box(ax,'on');
    set(ax,'FontSize',fs_tick);

    save_paper_fig(fig_h, OUT_FIGS, fig_fname);
end

% =========================================================================
% FIGURES A1–A7: Sector-Level Break Magnitude by Year (Full Sample)
% =========================================================================

fprintf('[produce_paper_figures] Figures A1–A7...\n');

sector_fnames = { ...
    'figureA1_output_income', ...
    'figureA2_labor_market', ...
    'figureA3_orders_inventories', ...
    'figureA4_housing_construction', ...
    'figureA5_prices', ...
    'figureA6_money_credit', ...
    'figureA7_interest_exchange'};

sector_ftitles = { ...
    'Figure A1: Output and Income — Mean Break Magnitude by Year (Full Sample, 1959–2026)', ...
    'Figure A2: Labor Market — Mean Break Magnitude by Year (Full Sample, 1959–2026)', ...
    'Figure A3: Orders and Inventories — Mean Break Magnitude by Year (Full Sample, 1959–2026)', ...
    'Figure A4: Housing / Construction — Mean Break Magnitude by Year (Full Sample, 1959–2026)', ...
    'Figure A5: Prices — Mean Break Magnitude by Year (Full Sample, 1959–2026)', ...
    'Figure A6: Money and Credit — Mean Break Magnitude by Year (Full Sample, 1959–2026)', ...
    'Figure A7: Interest and Exchange Rates — Mean Break Magnitude by Year (Full Sample, 1959–2026)'};

for si = 1:n_sec
    grp = sector_groups(si);

    uv_yr_mean = NaN(numel(yr_range),1);
    cv_yr_mean = NaN(numel(yr_range),1);

    uv_grp_mask = uv_break_grp == grp;
    cv_grp_mask = cv_break_grp == grp;

    for yi = 1:numel(yr_range)
        yr = yr_range(yi);
        uv_v = uv_break_logsig(uv_grp_mask & uv_break_years == yr);
        cv_v = cv_break_logsig(cv_grp_mask & cv_break_years == yr);
        if ~isempty(uv_v), uv_yr_mean(yi) = nanmean(uv_v); end
        if ~isempty(cv_v), cv_yr_mean(yi) = nanmean(cv_v); end
    end

    has_data = ~isnan(uv_yr_mean) | ~isnan(cv_yr_mean);
    yr_plt   = yr_range(has_data);
    uv_plt   = uv_yr_mean(has_data);
    cv_plt   = cv_yr_mean(has_data);

    fig_h = figure('Color','w','Position',[50 50 1050 420],'Visible','off');
    ax    = axes('Parent',fig_h);
    hold(ax,'on');

    bw = 0.40;
    bar(ax, yr_plt - bw/2 - 0.05, uv_plt, bw, ...
        'FaceColor', clr_uv, 'EdgeColor','none','DisplayName','UV');
    bar(ax, yr_plt + bw/2 + 0.05, cv_plt, bw, ...
        'FaceColor', clr_cv, 'EdgeColor','none','DisplayName','CV');

    fill_between(ax, 2008, 2019,            clr_gfc);
    fill_between(ax, 2020, max(yr_range)+1, clr_cov);

    yline(ax, 0, '-k', 'LineWidth',0.8, 'HandleVisibility','off');
    xline(ax, 1984,'--k','LineWidth',1.0,'Label','1984', ...
        'HandleVisibility','off','LabelHorizontalAlignment','left');
    xline(ax, 2008,'--k','LineWidth',1.0,'Label','2008', ...
        'HandleVisibility','off','LabelHorizontalAlignment','left');
    xline(ax, 2020,'--k','LineWidth',1.0,'Label','2020', ...
        'HandleVisibility','off','LabelHorizontalAlignment','left');

    patch(ax,NaN,NaN,clr_gfc,'FaceAlpha',1,'EdgeColor','none', ...
        'DisplayName','Post-GFC (2008–2019)');
    patch(ax,NaN,NaN,clr_cov,'FaceAlpha',1,'EdgeColor','none', ...
        'DisplayName','Post-COVID (2020–)');

    hold(ax,'off');

    all_v = [uv_plt(~isnan(uv_plt)); cv_plt(~isnan(cv_plt))];
    yl_left = [-0.3 0.3];  % default if no data
    if ~isempty(all_v)
        ypad    = max(0.1, range(all_v)*0.15);
        yl_left = [min(all_v)-ypad, max(all_v)+ypad];
    end

    % Left axis limits
    yyaxis(ax, 'left');
    ylim(ax, yl_left);
    ylabel(ax,'Mean Log \sigma-Ratio','FontSize',fs_axis,'Interpreter','tex');
    set(ax, 'YColor', 'k');

    % Right y-axis: approx % change in sigma via yyaxis (no secondary toolbar)
    yyaxis(ax, 'right');
    ylim(ax, (exp(yl_left)-1)*100);
    ylabel(ax,'Approx. % Change in \sigma','FontSize',fs_axis,'Interpreter','tex');
    set(ax, 'YColor', [0.35 0.35 0.35]);
    yyaxis(ax, 'left');   % restore active side

    legend(ax,'show','Location','best','FontSize',fs_legend);
    xlabel(ax,'Year','FontSize',fs_axis);
    title(ax, sector_ftitles{si},'FontSize',fs_title,'FontWeight','bold');
    xlim(ax,[1959 2027]); grid(ax,'on'); box(ax,'on');
    set(ax,'FontSize',fs_tick);

    save_paper_fig(fig_h, OUT_FIGS, sector_fnames{si});
    fprintf('  Saved %s\n', sector_fnames{si});
end

% =========================================================================
% TABLE 4: Sector-level summary with median breaks added
% =========================================================================

fprintf('[produce_paper_figures] Table 4...\n');

% Build table rows: sector × window × measure
t4_sector   = strings(0,1);
t4_window   = strings(0,1);
t4_measure  = strings(0,1);
t4_n        = [];
t4_p0b      = [];
t4_p1b      = [];
t4_p2p      = [];
t4_med_brk  = [];
t4_med_lsig = [];
t4_pmod     = [];

for si = 1:n_sec
    for wi = 1:3
        for meas = ["UV","CV"]
            r = mb_sum.window == win_ordered(wi) & ...
                mb_sum.fred_group == sector_groups(si) & ...
                mb_sum.measure == meas;
            if ~any(r); continue; end
            t4_sector(end+1,1)   = sector_groups(si);
            t4_window(end+1,1)   = win_ordered(wi);
            t4_measure(end+1,1)  = meas;
            t4_n(end+1,1)        = mb_sum.n_series(r);
            t4_p0b(end+1,1)      = mb_sum.share_0b(r)     * 100;
            t4_p1b(end+1,1)      = mb_sum.share_1b(r)     * 100;
            t4_p2p(end+1,1)      = mb_sum.share_2plus(r)  * 100;
            t4_med_brk(end+1,1)  = mb_sum.median_breaks(r);
            t4_med_lsig(end+1,1) = mb_sum.median_log_sigma_ratio(r);
            t4_pmod(end+1,1)     = mb_sum.share_moderating(r) * 100;
        end
    end
end

t4 = table(t4_sector, t4_window, t4_measure, t4_n, ...
    t4_p0b, t4_p1b, t4_p2p, t4_med_brk, t4_med_lsig, t4_pmod, ...
    'VariableNames',{ ...
        'Sector','Window','Measure','N', ...
        'Share_0Breaks_pct','Share_1Break_pct','Share_2plus_pct', ...
        'Median_Breaks', ...
        'Median_LogSigmaRatio','Share_Moderating_pct'});

writetable(t4, fullfile(OUT_TABLES,'table4_updated.csv'));
fprintf('[produce_paper_figures] Saved table4_updated.csv (%d rows)\n', height(t4));

fprintf('\n[produce_paper_figures] All figures and Table 4 complete.\n');
fprintf('  Figures → %s\n', OUT_FIGS);
fprintf('  Table 4 → %s\n', fullfile(OUT_TABLES,'table4_updated.csv'));

% =========================================================================
% LOCAL FUNCTIONS
% =========================================================================

function yrs = parse_break_years(bstr)
    yrs = [];
    if bstr == "" || strtrim(bstr) == ""; return; end
    parts = strtrim(strsplit(bstr,'|'));
    parts = parts(~cellfun(@(x) isempty(strtrim(x)), cellstr(parts)));
    if isempty(parts); return; end
    try
        dts = datetime(parts,'InputFormat','dd-MMM-yyyy');
        yrs = year(dts(:));
    catch; end
end

function [br_yrs, br_lsig] = parse_break_sigma(zt, dates, bstr)
    br_yrs  = [];
    br_lsig = [];
    if bstr == "" || strtrim(bstr) == ""; return; end
    parts = strtrim(strsplit(bstr,'|'));
    parts = parts(~cellfun(@(x) isempty(strtrim(x)), cellstr(parts)));
    if isempty(parts); return; end
    try
        break_dates = datetime(parts,'InputFormat','dd-MMM-yyyy');
    catch; return; end
    break_dates = sort(break_dates(:));
    sentinel    = dates(end) + caldays(32);
    boundaries  = [dates(1); break_dates; sentinel];
    for k = 1:numel(break_dates)
        pre  = dates >= boundaries(k)   & dates < boundaries(k+1);
        post = dates >= boundaries(k+1) & dates < boundaries(k+2);
        mp   = nanmean(zt(pre));
        mn   = nanmean(zt(post));
        if isnan(mp)||isnan(mn)||mp<=0||mn<=0; continue; end
        br_yrs(end+1,1)  = year(break_dates(k)); %#ok<AGROW>
        br_lsig(end+1,1) = log(mn/mp);            %#ok<AGROW>
    end
end

function fill_between(ax, x1, x2, col)
    large = 1e6;
    fill(ax,[x1 x2 x2 x1],[-large -large large large],col, ...
        'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
end

function save_paper_fig(fig, out_dir, base_name)
    exportgraphics(fig, fullfile(out_dir,[base_name '.png']), 'Resolution',300);
    exportgraphics(fig, fullfile(out_dir,[base_name '.pdf']), 'ContentType','vector');
    savefig(fig,         fullfile(out_dir,[base_name '.fig']));
    close(fig);
end
