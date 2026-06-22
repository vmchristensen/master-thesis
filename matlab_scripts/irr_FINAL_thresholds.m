% INSIDEOUT_threshold_comparison.m
% Computes IrrLeft, IrrRight, IrrLI, IrrWhole at tau=1
% Both plain mean (q=0) and top-5% threshold (q=0.95)
% Then compares between Left/Right/Bilateral groups
% ------EDIT BEFORE RUNNING------

base_dir_left  = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_dswu';
base_dir_right = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';
results_dir    = '/home/marlene/Desktop/victoria/INSIDEOUT_results';
demo_file      = '/home/marlene/Desktop/victoria/demographics_cleaned_filtered_259.csv';

TAU    = 1;
N_ROIS = 100;

exclude_list = {'sub-EZM20001', 'sub-EZM20002', 'sub-EZM20003', 'sub-EZM20004', ...
                'sub-EZM20005', 'sub-EZM20006', 'sub-EZM20007', 'sub-EZM20008', ...
                'sub-EZM20009', 'sub-EZM20024', 'sub-EZM20028', 'sub-EZM20033', ...
                'sub-EZM20045', 'sub-EZM20046', 'sub-EZM20059', 'sub-EZM20085', ...
                'sub-EZM20086', 'sub-EZM20094', 'sub-EZM20099', 'sub-EZM20104', ...
                'sub-EZM20120', 'sub-EZM20136', 'sub-EZM20137', 'sub-EZM20146', ...
                'sub-EZM20148', 'sub-EZM20152', 'sub-EZM20153', 'sub-EZM20155', ...
                'sub-EZM20168', 'sub-EZM20176', 'sub-EZM20052', 'sub-EZM20057', ...
                'sub-EZM20091', 'sub-EZM20112', 'sub-EZM20115', 'sub-EZM20186', ...
                'sub-EZM20014', 'sub-EZM30214'};

% ------STOP EDIT------

if ~exist(results_dir, 'dir'), mkdir(results_dir); end

%% Get subject list
subj_dirs    = dir(fullfile(base_dir_left, 'sub-EZM*'));
all_subjects = {subj_dirs.name};
subjects     = all_subjects(~ismember(all_subjects, exclude_list));
NSUB         = length(subjects);
fprintf('Subjects after exclusion: %d\n', NSUB);

%% Exclude missing-LQ subject
demo     = readtable(demo_file);
demo.ID  = strtrim(string(demo.ID));
missing_lq = strcat('sub-', demo.ID(isnan(demo.Laterality_Quotient)));
subjects   = subjects(~ismember(subjects, missing_lq));
NSUB       = length(subjects);
fprintf('Subjects after LQ exclusion: %d\n', NSUB);

%% Preallocate
SubjectID     = cell(NSUB, 1);
IrrLeft_q0    = nan(NSUB, 1);
IrrRight_q0   = nan(NSUB, 1);
IrrLI_q0      = nan(NSUB, 1);
IrrWhole_q0   = nan(NSUB, 1);
IrrLeft_q95   = nan(NSUB, 1);
IrrRight_q95  = nan(NSUB, 1);
IrrLI_q95     = nan(NSUB, 1);
IrrWhole_q95  = nan(NSUB, 1);

%% Main loop
for sub = 1:NSUB
    subj = subjects{sub};
    SubjectID{sub} = subj;
    fprintf('[%d/%d] %s\n', sub, NSUB, subj);

    left_file  = fullfile(base_dir_left,  subj, 'R_mean_roi_left.mat');
    right_file = fullfile(base_dir_right, subj, 'R_mean_roi_right.mat');

    if ~exist(left_file,'file') || ~exist(right_file,'file')
        warning('Missing file for %s — skipping', subj);
        continue
    end

    dl   = load(left_file);
    dr   = load(right_file);
    ts_L = dl.R_mean_roi_left;    % 100 x T
    ts_R = dr.R_mean_roi_right;   % 100 x T

    % Validate dimensions
    if size(ts_L,1) ~= N_ROIS || size(ts_R,1) ~= N_ROIS
        warning('Unexpected ROI count for %s: L=%d R=%d — skipping', ...
                subj, size(ts_L,1), size(ts_R,1));
        continue
    end

    Tm = size(ts_L,2);
    if size(ts_R,2) ~= Tm
        warning('Timepoint mismatch for %s — skipping', subj);
        continue
    end
    if Tm <= TAU
        warning('Time series too short for %s — skipping', subj);
        continue
    end

    % ----------------------------------------------------------
    % Compute FSdiff vector for left hemisphere (100 x T)
    % ----------------------------------------------------------
    ts_rev = ts_L(:, end:-1:1);
    FC_fwd = corr(ts_L(:, 1:Tm-TAU)',   ts_L(:, 1+TAU:Tm)');
    FC_rev = corr(ts_rev(:, 1:Tm-TAU)', ts_rev(:, 1+TAU:Tm)');
    FS_fwd = -0.5 * log(1 - FC_fwd.^2 + 1e-12);
    FS_rev = -0.5 * log(1 - FC_rev.^2 + 1e-12);
    ref_L  = (FS_fwd(:) - FS_rev(:)).^2;   % 10000 x 1

    % ----------------------------------------------------------
    % Compute FSdiff vector for right hemisphere (100 x T)
    % ----------------------------------------------------------
    ts_rev = ts_R(:, end:-1:1);
    FC_fwd = corr(ts_R(:, 1:Tm-TAU)',   ts_R(:, 1+TAU:Tm)');
    FC_rev = corr(ts_rev(:, 1:Tm-TAU)', ts_rev(:, 1+TAU:Tm)');
    FS_fwd = -0.5 * log(1 - FC_fwd.^2 + 1e-12);
    FS_rev = -0.5 * log(1 - FC_rev.^2 + 1e-12);
    ref_R  = (FS_fwd(:) - FS_rev(:)).^2;   % 10000 x 1

    % ----------------------------------------------------------
    % Compute FSdiff vector for whole atlas (200 x T)
    % ----------------------------------------------------------
    ts_W   = [ts_L; ts_R];
    ts_rev = ts_W(:, end:-1:1);
    FC_fwd = corr(ts_W(:, 1:Tm-TAU)',   ts_W(:, 1+TAU:Tm)');
    FC_rev = corr(ts_rev(:, 1:Tm-TAU)', ts_rev(:, 1+TAU:Tm)');
    FS_fwd = -0.5 * log(1 - FC_fwd.^2 + 1e-12);
    FS_rev = -0.5 * log(1 - FC_rev.^2 + 1e-12);
    ref_W  = (FS_fwd(:) - FS_rev(:)).^2;   % 40000 x 1

    % ----------------------------------------------------------
    % Plain mean (q=0, all pairs)
    % ----------------------------------------------------------
    IrrLeft_q0(sub)  = mean(ref_L);
    IrrRight_q0(sub) = mean(ref_R);
    IrrWhole_q0(sub) = mean(ref_W);
    IrrLI_q0(sub)    = (IrrLeft_q0(sub) - IrrRight_q0(sub)) / ...
                       (IrrLeft_q0(sub) + IrrRight_q0(sub) + 1e-12);

    % ----------------------------------------------------------
    % Thresholded (q=0.95, top 5% most irreversible pairs)
    % ----------------------------------------------------------
    IrrLeft_q95(sub)  = mean(ref_L(ref_L > quantile(ref_L, 0.95)));
    IrrRight_q95(sub) = mean(ref_R(ref_R > quantile(ref_R, 0.95)));
    IrrWhole_q95(sub) = mean(ref_W(ref_W > quantile(ref_W, 0.95)));
    IrrLI_q95(sub)    = (IrrLeft_q95(sub) - IrrRight_q95(sub)) / ...
                        (IrrLeft_q95(sub) + IrrRight_q95(sub) + 1e-12);

end

%% Save full CSV
T_out = table(SubjectID, ...
    IrrLeft_q0, IrrRight_q0, IrrLI_q0, IrrWhole_q0, ...
    IrrLeft_q95, IrrRight_q95, IrrLI_q95, IrrWhole_q95);

valid = ~cellfun(@isempty, SubjectID) & ~isnan(IrrLeft_q0);
T_out = T_out(valid, :);

outfile = fullfile(results_dir, 'INSIDEOUT_q0_vs_q95_n258.csv');
writetable(T_out, outfile);
fprintf('\nSaved %d subjects to: %s\n', height(T_out), outfile);

%% Merge with demographics for group labels
T_out.ID = strtrim(erase(string(T_out.SubjectID), 'sub-'));
demo.ID  = strtrim(string(demo.ID));
T_merged = innerjoin(T_out, demo(:, {'ID','HEMISPHERIC_2','LI_GEN'}), 'Keys', 'ID');
fprintf('Subjects after demo merge: %d\n', height(T_merged));

%% Group comparisons
groups   = {'Left','Right','Bilateral'};
versions = {
    {'IrrLeft_q0','IrrRight_q0','IrrLI_q0','IrrWhole_q0'},   'PLAIN MEAN (q=0)';
    {'IrrLeft_q95','IrrRight_q95','IrrLI_q95','IrrWhole_q95'}, 'THRESHOLDED (q=0.95)'
};

for ver = 1:2
    vars  = versions{ver,1};
    label = versions{ver,2};

    fprintf('\n=======================================================\n');
    fprintf('VERSION: %s\n', label);
    fprintf('=======================================================\n');

    for v = 1:length(vars)
        vname = vars{v};
        fprintf('\n--- %s ---\n', vname);

        for g = 1:length(groups)
            idx  = strcmp(T_merged.HEMISPHERIC_2, groups{g});
            vals = T_merged.(vname)(idx);
            fprintf('  %s (n=%d): mean=%.6f  sd=%.6f  median=%.6f\n', ...
                groups{g}, sum(idx), mean(vals), std(vals), median(vals));
        end

        % Kruskal-Wallis
        [p_kw, tbl_kw] = kruskalwallis(T_merged.(vname), ...
                                        T_merged.HEMISPHERIC_2, 'off');
        fprintf('  Kruskal-Wallis: chi2=%.3f  df=%d  p=%.4f\n', ...
            tbl_kw{2,5}, tbl_kw{2,3}, p_kw);

        % Pairwise Wilcoxon with Bonferroni correction
        pairs  = {{'Left','Right'}, {'Left','Bilateral'}, {'Right','Bilateral'}};
        raw_p  = nan(1,3);
        for pp = 1:3
            g1     = T_merged.(vname)(strcmp(T_merged.HEMISPHERIC_2, pairs{pp}{1}));
            g2     = T_merged.(vname)(strcmp(T_merged.HEMISPHERIC_2, pairs{pp}{2}));
            raw_p(pp) = ranksum(g1, g2);
        end
        bonf_p = min(raw_p * 3, 1);
        for pp = 1:3
            fprintf('  %s vs %s: p_raw=%.4f  p_bonf=%.4f\n', ...
                pairs{pp}{1}, pairs{pp}{2}, raw_p(pp), bonf_p(pp));
        end
    end
end

%% Spearman correlation with continuous LI
if ismember('LI_GEN', demo.Properties.VariableNames)
    fprintf('\n=======================================================\n');
    fprintf('SPEARMAN CORRELATION WITH CONTINUOUS LI_GEN\n');
    fprintf('=======================================================\n');

    all_vars = [versions{1,1}, versions{2,1}];
    for v = 1:length(all_vars)
        vname = all_vars{v};
        x     = T_merged.(vname);
        y     = T_merged.LI_GEN;
        ok    = ~isnan(x) & ~isnan(y);
        [rho, pval] = corr(x(ok), y(ok), 'Type', 'Spearman');
        fprintf('  %s ~ LI_GEN: rho=%.4f  p=%.4f\n', vname, rho, pval);
    end
end

fprintf('\nDone.\n');


% Binary-level IrrLI check
% Run after T_merged exists in workspace

% If T_merged is not in workspace, reload:
% T_merged = ... (rerun the merge section of the main script)

% Add HEMISPHERIC (binary) to T_merged if not already there
if ~ismember('HEMISPHERIC', T_merged.Properties.VariableNames)
    demo      = readtable('/home/marlene/Desktop/victoria/demographics_cleaned_filtered_259.csv');
    demo.ID   = strtrim(string(demo.ID));
    T_merged2 = innerjoin(T_merged, demo(:, {'ID','HEMISPHERIC'}), 'Keys', 'ID');
else
    T_merged2 = T_merged;
end

% Group means
groups_bin = {'Typical', 'Atypical'};
fprintf('\n--- IrrLI_q0 by binary group ---\n');
for g = 1:length(groups_bin)
    idx  = strcmp(T_merged2.HEMISPHERIC, groups_bin{g});
    vals = T_merged2.IrrLI_q0(idx);
    fprintf('  %s (n=%d): mean=%.6f  sd=%.6f\n', ...
        groups_bin{g}, sum(idx), mean(vals), std(vals));
end

% Wilcoxon rank-sum
typ   = T_merged2.IrrLI_q0(strcmp(T_merged2.HEMISPHERIC, 'Typical'));
atyp  = T_merged2.IrrLI_q0(strcmp(T_merged2.HEMISPHERIC, 'Atypical'));
[p, h, stats] = ranksum(typ, atyp);
fprintf('\n  Wilcoxon ranksum: W=%.0f  p=%.4f\n', stats.ranksum, p);


% Group means
groups_bin = {'Typical', 'Atypical'};
fprintf('\n--- IrrLI_q95 by binary group ---\n');
for g = 1:length(groups_bin)
    idx  = strcmp(T_merged2.HEMISPHERIC, groups_bin{g});
    vals = T_merged2.IrrLI_q95(idx);
    fprintf('  %s (n=%d): mean=%.6f  sd=%.6f\n', ...
        groups_bin{g}, sum(idx), mean(vals), std(vals));
end


% Wilcoxon rank-sum
typ   = T_merged2.IrrLI_q95(strcmp(T_merged2.HEMISPHERIC, 'Typical'));
atyp  = T_merged2.IrrLI_q95(strcmp(T_merged2.HEMISPHERIC, 'Atypical'));
[p, h, stats] = ranksum(typ, atyp);
fprintf('\n  Wilcoxon ranksum: W=%.0f  p=%.4f\n', stats.ranksum, p);
