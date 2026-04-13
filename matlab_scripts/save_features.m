%% =========================================================================
%  FEATURE EXTRACTION — SIMPLE VERSION
%  Predicting language LI from rsfMRI functional connectivity
%
%  FEATURES (all designed as lateralization indices):
%    Block 1: Node-level FC asymmetry    — AI  per ROI         [101 features]
%    Block 2: Node-level FC asymmetry    — raw diff per ROI    [101 features]
%    Block 3: Homotopic FC               — diagonal of LR      [101 features]
%    Block 4: Interhemispheric AI        — outgoing L vs R     [101 features]
%    Block 5: Sample entropy asymmetry   — AI  per ROI         [101 features]
%    Block 6: Sample entropy asymmetry   — raw diff per ROI    [101 features]
%    Block 7: Approx entropy asymmetry   — AI  per ROI         [101 features]
%    Block 8: Approx entropy asymmetry   — raw diff per ROI    [101 features]
%
%  Total: 808 features per subject
%
%  NO graph theory metrics (BCT not required for this script)
%
%  OUTPUT:
%    features_simple.csv         — full feature matrix + LI_GEN
%    feature_block_index.csv     — block names, sizes, descriptions
%
%  REQUIRES: no external toolboxes
% =========================================================================

clear; clc;

%% -------------------------------------------------------------------------
%  USER SETTINGS — edit these paths
% -------------------------------------------------------------------------
base_left_dir  = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_dswu';
base_right_dir = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';
fc_dir         = '/home/marlene/Desktop/victoria/FC_matrices';
csv_path       = '/home/marlene/Desktop/data_victoria/demographics_cleaned_filtered.csv';   % <-- edit this
out_dir        = '/home/marlene/Desktop/victoria/FC_features';
mkdir(out_dir);

% Entropy parameters (standard defaults for fMRI)
SampEn_m = 2;      % template length
SampEn_r = 0.2;    % tolerance (fraction of SD, applied per series)
ApEn_m   = 2;
ApEn_r   = 0.2;

eps_val  = 1e-6;
n_rois   = 101;

%% -------------------------------------------------------------------------
%  SUBJECTS TO EXCLUDE (from batch_fc_matrices_all.m)
% -------------------------------------------------------------------------
exclude_list = {'sub-EZM20001','sub-EZM20002','sub-EZM20003','sub-EZM20004', ...
    'sub-EZM20005','sub-EZM20006','sub-EZM20007','sub-EZM20008', ...
    'sub-EZM20009','sub-EZM20024','sub-EZM20028','sub-EZM20033', ...
    'sub-EZM20045','sub-EZM20046','sub-EZM20059','sub-EZM20085', ...
    'sub-EZM20086','sub-EZM20094','sub-EZM20099','sub-EZM20104', ...
    'sub-EZM20120','sub-EZM20136','sub-EZM20137','sub-EZM20146', ...
    'sub-EZM20148','sub-EZM20152','sub-EZM20153','sub-EZM20155', ...
    'sub-EZM20168','sub-EZM20176','sub-EZM20052','sub-EZM20057', ...
    'sub-EZM20091','sub-EZM20112','sub-EZM20115','sub-EZM20186'};

%% -------------------------------------------------------------------------
%  LOAD GROUP FC MATRICES (already computed)
% -------------------------------------------------------------------------
fprintf('Loading group FC matrices...\n');

tmp = load(fullfile(fc_dir, 'full_fc_left.mat'));
FC_left_flat = tmp.full_fc_left;
subject_ids  = tmp.subject_ids;

tmp = load(fullfile(fc_dir, 'full_fc_right.mat'));
FC_right_flat = tmp.full_fc_right;

tmp = load(fullfile(fc_dir, 'full_fc_matrix.mat'));
FC_lr_flat = tmp.full_fc_matrix;

n_subj = size(FC_left_flat, 1);
fprintf('Loaded %d subjects.\n', n_subj);

FC_L  = reshape(FC_left_flat,  [n_subj, n_rois, n_rois]);
FC_R  = reshape(FC_right_flat, [n_subj, n_rois, n_rois]);
FC_LR = reshape(FC_lr_flat,    [n_subj, n_rois, n_rois]);

%% -------------------------------------------------------------------------
%  LOAD LI TARGET
% -------------------------------------------------------------------------
fprintf('Loading LI_GEN from demographics CSV...\n');
demo     = readtable(csv_path);
ids_clean = regexprep(subject_ids, '^sub-', '');

LI_GEN = NaN(n_subj, 1);
for s = 1:n_subj
    idx = strcmp(demo.ID, ids_clean{s});
    if any(idx)
        val = demo.LI_GEN(idx);
        if ~isempty(val) && ~isnan(val)
            LI_GEN(s) = val;
        end
    end
end

valid_li = ~isnan(LI_GEN);
fprintf('%d / %d subjects have valid LI_GEN.\n', sum(valid_li), n_subj);

FC_L      = FC_L(valid_li,:,:);
FC_R      = FC_R(valid_li,:,:);
FC_LR     = FC_LR(valid_li,:,:);
LI_GEN    = LI_GEN(valid_li);
ids_final = ids_clean(valid_li);
n_valid   = sum(valid_li);

%% -------------------------------------------------------------------------
%  PRE-ALLOCATE FEATURE MATRICES
% -------------------------------------------------------------------------
feat_node_AI        = zeros(n_valid, n_rois);
feat_node_diff      = zeros(n_valid, n_rois);
feat_homotopic      = zeros(n_valid, n_rois);
feat_ihem_AI        = zeros(n_valid, n_rois);
feat_sampen_AI      = zeros(n_valid, n_rois);
feat_sampen_diff    = zeros(n_valid, n_rois);
feat_apen_AI        = zeros(n_valid, n_rois);
feat_apen_diff      = zeros(n_valid, n_rois);

%% =========================================================================
%  MAIN LOOP
% =========================================================================
fprintf('\nExtracting features for %d subjects...\n', n_valid);

for s = 1:n_valid
    if mod(s,25) == 0
        fprintf('  %d / %d\n', s, n_valid);
    end

    sub_id = ['sub-', ids_final{s}];

    % -------------------------------------------------------------------
    %  BLOCKS 1-4: FC asymmetry (from precomputed matrices)
    % -------------------------------------------------------------------
    L  = squeeze(FC_L(s,:,:));
    R  = squeeze(FC_R(s,:,:));
    LR = squeeze(FC_LR(s,:,:));

    % Zero diagonal
    L(1:n_rois+1:end) = 0;
    R(1:n_rois+1:end) = 0;

    % Node strength (mean connectivity per ROI)
    str_L = mean(L, 2);   % 101x1
    str_R = mean(R, 2);

    feat_node_diff(s,:) = (str_L - str_R)';
    feat_node_AI(s,:)   = ((str_L - str_R) ./ (abs(str_L) + abs(str_R) + eps_val))';

    % Homotopic FC (diagonal of interhemispheric matrix)
    feat_homotopic(s,:) = diag(LR)';

    % Interhemispheric outgoing connectivity AI
    ihem_L = mean(LR, 2);
    ihem_R = mean(LR, 1)';
    feat_ihem_AI(s,:) = ((ihem_L - ihem_R) ./ (abs(ihem_L) + abs(ihem_R) + eps_val))';

    % -------------------------------------------------------------------
    %  BLOCKS 5-8: Entropy (from raw ROI time series)
    % -------------------------------------------------------------------
    left_file  = fullfile(base_left_dir,  sub_id, 'R_mean_roi_left.mat');
    right_file = fullfile(base_right_dir, sub_id, 'R_mean_roi_right.mat');

    if ~isfile(left_file) || ~isfile(right_file)
        warning('Time series missing for %s — entropy set to NaN.', sub_id);
        feat_sampen_AI(s,:)   = NaN;
        feat_sampen_diff(s,:) = NaN;
        feat_apen_AI(s,:)     = NaN;
        feat_apen_diff(s,:)   = NaN;
        continue;
    end

    ts_left  = load(left_file);
    ts_right = load(right_file);
    TS_L = ts_left.R_mean_roi_left;    % n_rois x n_timepoints
    TS_R = ts_right.R_mean_roi_right;

    % Compute entropy per ROI
    SampEn_L = zeros(n_rois, 1);
    SampEn_R = zeros(n_rois, 1);
    ApEn_L   = zeros(n_rois, 1);
    ApEn_R   = zeros(n_rois, 1);

    for roi = 1:n_rois
        ts_l = TS_L(roi, :)';
        ts_r = TS_R(roi, :)';

        % Tolerance scaled to each series SD (standard practice)
        r_l = SampEn_r * std(ts_l);
        r_r = SampEn_r * std(ts_r);

        SampEn_L(roi) = compute_sampen(ts_l, SampEn_m, r_l);
        SampEn_R(roi) = compute_sampen(ts_r, SampEn_m, r_r);
        ApEn_L(roi)   = compute_apen(ts_l,   ApEn_m,   ApEn_r * std(ts_l));
        ApEn_R(roi)   = compute_apen(ts_r,   ApEn_m,   ApEn_r * std(ts_r));
    end

    % Asymmetry
    feat_sampen_diff(s,:) = (SampEn_L - SampEn_R)';
    feat_sampen_AI(s,:)   = ((SampEn_L - SampEn_R) ./ (abs(SampEn_L) + abs(SampEn_R) + eps_val))';
    feat_apen_diff(s,:)   = (ApEn_L - ApEn_R)';
    feat_apen_AI(s,:)     = ((ApEn_L - ApEn_R) ./ (abs(ApEn_L) + abs(ApEn_R) + eps_val))';

end
fprintf('Feature extraction complete.\n');

%% =========================================================================
%  ASSEMBLE & SAVE
% =========================================================================
fprintf('\nAssembling feature table...\n');

roi_idx = 1:n_rois;
mk = @(prefix) arrayfun(@(i) sprintf('%s_R%d', prefix, i), roi_idx, 'UniformOutput', false);

names = [mk('nodeAI'), mk('nodeDiff'), mk('homotopic'), mk('ihemAI'), ...
         mk('sampEnAI'), mk('sampEnDiff'), mk('apEnAI'), mk('apEnDiff')];

all_feats = [feat_node_AI, feat_node_diff, feat_homotopic, feat_ihem_AI, ...
             feat_sampen_AI, feat_sampen_diff, feat_apen_AI, feat_apen_diff];

T = array2table(all_feats, 'VariableNames', names);
T.subject_id = ids_final(:);
T.LI_GEN     = LI_GEN;
T = [T(:,end-1), T(:,1:end-2), T(:,end)];   % subject_id first, LI_GEN last

% Save features
out_csv = fullfile(out_dir, 'features_simple.csv');
writetable(T, out_csv);
fprintf('Saved: %s\n', out_csv);

% Save block index
block_info = {
    'nodeAI',       n_rois, 'Node FC asymmetry index (L-R)/(|L|+|R|+eps)';
    'nodeDiff',     n_rois, 'Node FC raw difference L-R';
    'homotopic',    n_rois, 'Homotopic interhemispheric FC (diagonal of LR)';
    'ihemAI',       n_rois, 'Interhemispheric outgoing connectivity AI';
    'sampEnAI',     n_rois, 'Sample entropy asymmetry index';
    'sampEnDiff',   n_rois, 'Sample entropy raw difference L-R';
    'apEnAI',       n_rois, 'Approximate entropy asymmetry index';
    'apEnDiff',     n_rois, 'Approximate entropy raw difference L-R';
};
block_table = cell2table(block_info, 'VariableNames', {'block','n_features','description'});
writetable(block_table, fullfile(out_dir, 'feature_block_index_simple.csv'));

%% =========================================================================
%  SUMMARY
% =========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Subjects:        %d\n', n_valid);
fprintf('ROIs:            %d\n', n_rois);
fprintf('Total features:  %d\n', width(T)-2);
fprintf('Output:          %s\n', out_dir);
fprintf('==============================\n');
disp(block_table);

%% =========================================================================
%  LOCAL FUNCTIONS
%  (placed at end so script can be run as a whole without needing a toolbox)
% =========================================================================

function se = compute_sampen(x, m, r)
% Sample entropy of time series x with template length m and tolerance r.
% Returns NaN if series is too short or has zero variance.
    N = length(x);
    if N < (m+2) || r == 0
        se = NaN; return;
    end
    % Count template matches at length m and m+1
    B = 0; A = 0;
    for i = 1:N-m
        xi = x(i:i+m-1);
        for j = i+1:N-m
            if max(abs(xi - x(j:j+m-1))) < r
                B = B + 1;
                % Check one step further
                if abs(x(i+m) - x(j+m)) < r
                    A = A + 1;
                end
            end
        end
    end
    if B == 0
        se = NaN;
    else
        se = -log(A / B);
    end
end

function ae = compute_apen(x, m, r)
% Approximate entropy of time series x with template length m and tolerance r.
    N = length(x);
    if N < (m+1) || r == 0
        ae = NaN; return;
    end
    ae = phi_apen(x, m, r, N) - phi_apen(x, m+1, r, N);
end

function p = phi_apen(x, m, r, N)
% Helper: compute Phi^m for approximate entropy.
    count = zeros(N-m+1, 1);
    for i = 1:N-m+1
        xi = x(i:i+m-1);
        for j = 1:N-m+1
            if max(abs(xi - x(j:j+m-1))) <= r
                count(i) = count(i) + 1;
            end
        end
    end
    p = mean(log(count / (N-m+1)));
end
