%% =========================================================================
%  FEATURE EXTRACTION — ASYMMETRIC ATLAS VERSION
%  Predicting language LI from rsfMRI functional connectivity
%
%  Atlas: SENS-HCP-Tian (asymmetric)
%    Left  hemisphere: 100 ROIs (99 named + O3_1 unknown)
%    Right hemisphere: 104 ROIs (103 named + O3_1 unknown)
%    Paired homologs:   99 matched L/R pairs (used for all asymmetry features)
%    Unpaired left:      1 ROI — IFSa (local left index 21, language-relevant)
%    Unpaired right:     5 ROIs — IFSp x2, AVI, PUT-VP, pGP (not used in features)
%
%  FEATURES (computed over 99 paired ROIs):
%    Block 1: Node-level FC asymmetry    — AI  per ROI pair    [99 features]
%    Block 2: Node-level FC asymmetry    — raw diff per ROI    [99 features]
%    Block 3: Homotopic FC               — diagonal of LR      [99 features]
%    Block 4: Interhemispheric AI        — outgoing L vs R     [99 features]
%    Block 5: Sample entropy asymmetry   — AI  per ROI pair    [99 features]
%    Block 6: Sample entropy asymmetry   — raw diff per ROI    [99 features]
%    Block 7: Approx entropy asymmetry   — AI  per ROI pair    [99 features]
%    Block 8: Approx entropy asymmetry   — raw diff per ROI    [99 features]
%
%  Total: 792 features per subject (8 x 99)
%
%  OUTPUT:
%    features_asymmetrical.csv        — full feature matrix + LI_GEN
%    feature_block_index.csv          — block names, sizes, descriptions
%
%  REQUIRES: no external toolboxes
% =========================================================================

clear; clc;

%% -------------------------------------------------------------------------
%  USER SETTINGS
% -------------------------------------------------------------------------
base_dir = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_asymmetrical';
fc_dir   = '/home/marlene/Desktop/victoria/FC_matrices_asymmetrical';
csv_path = '/home/marlene/Desktop/data_victoria/demographics_cleaned_filtered.csv';
out_dir  = '/home/marlene/Desktop/victoria/FC_features_asymmetrical';
mkdir(out_dir);

% Entropy parameters (standard defaults for fMRI)
SampEn_m = 2;    % template length
SampEn_r = 0.2;  % tolerance as fraction of SD
ApEn_m   = 2;
ApEn_r   = 0.2;

eps_val  = 1e-6;

%% -------------------------------------------------------------------------
%  ATLAS DIMENSIONS
%  Row ordering in R_mean_roi_full follows lookup table index (1 -> 204):
%    Rows   1-100 = left hemisphere  (99 named + O3_1 at row 55)
%    Rows 101-204 = right hemisphere (103 named + O3_1 at row 157 = 57+100)
% -------------------------------------------------------------------------
n_L     = 100;   % left ROIs
n_R     = 104;   % right ROIs
n_pairs = 99;    % matched homolog pairs

% Local indices within each hemisphere block (1-based)
% Derived from region_lookup.xlsx — 98 named pairs + O3_1 (atlas labels 2141/2142)
paired_left_idx  = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, ...
    18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, ...
    38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, ...
    57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, ...
    76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, ...
    95, 96, 97, 98, 99, 100 ];

paired_right_idx = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, ...
    18, 19, 20, 23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37, 38, 39, ...
    40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, ...
    59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, ...
    78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 95, 96, 97, ...
    98, 99, 100, 101, 102, 103 ];

% Unpaired left ROI: IFSa (local left index 21) — kept as extra feature
ifsa_idx = 21;

%% -------------------------------------------------------------------------
%  SUBJECTS TO EXCLUDE
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
%  LOAD GROUP FC MATRICES
% -------------------------------------------------------------------------
fprintf('Loading group FC matrices...\n');

tmp          = load(fullfile(fc_dir, 'full_fc_left.mat'));
FC_left_flat = tmp.full_fc_left;
subject_ids  = tmp.subject_ids;

tmp           = load(fullfile(fc_dir, 'full_fc_right.mat'));
FC_right_flat = tmp.full_fc_right;

tmp        = load(fullfile(fc_dir, 'full_fc_matrix.mat'));
FC_lr_flat = tmp.full_fc_matrix;

n_subj = size(FC_left_flat, 1);
fprintf('Loaded %d subjects.\n', n_subj);

% Reshape to 3D — note FC_LR is rectangular (n_L x n_R)
FC_L  = reshape(FC_left_flat,  [n_subj, n_L, n_L]);
FC_R  = reshape(FC_right_flat, [n_subj, n_R, n_R]);
FC_LR = reshape(FC_lr_flat,    [n_subj, n_L, n_R]);

%% -------------------------------------------------------------------------
%  LOAD LI TARGET
% -------------------------------------------------------------------------
fprintf('Loading LI_GEN from demographics CSV...\n');
demo      = readtable(csv_path);
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
%  PRE-ALLOCATE FEATURE MATRICES (all sized n_valid x n_pairs = n_valid x 99)
% -------------------------------------------------------------------------
feat_node_AI     = zeros(n_valid, n_pairs);
feat_node_diff   = zeros(n_valid, n_pairs);
feat_homotopic   = zeros(n_valid, n_pairs);
feat_ihem_AI     = zeros(n_valid, n_pairs);
feat_sampen_AI   = zeros(n_valid, n_pairs);
feat_sampen_diff = zeros(n_valid, n_pairs);
feat_apen_AI     = zeros(n_valid, n_pairs);
feat_apen_diff   = zeros(n_valid, n_pairs);

% Extra: IFSa unilateral interhemispheric connectivity (1 feature per subject)
feat_ifsa_ihem   = zeros(n_valid, 1);

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
    L  = squeeze(FC_L(s,:,:));    % 100 x 100
    R  = squeeze(FC_R(s,:,:));    % 104 x 104
    LR = squeeze(FC_LR(s,:,:));   % 100 x 104

    % Zero diagonal of within-hemisphere matrices
    L(1:n_L+1:end) = 0;
    R(1:n_R+1:end) = 0;

    % Slice to paired ROIs only
    L_p  = L(paired_left_idx, paired_left_idx);    % 99 x 99
    R_p  = R(paired_right_idx, paired_right_idx);  % 99 x 99
    LR_p = LR(paired_left_idx, paired_right_idx);  % 99 x 99

    % Node strength (mean connectivity per ROI, over paired ROIs only)
    str_L = mean(L_p, 2);   % 99 x 1
    str_R = mean(R_p, 2);   % 99 x 1

    feat_node_diff(s,:) = (str_L - str_R)';
    feat_node_AI(s,:)   = ((str_L - str_R) ./ (abs(str_L) + abs(str_R) + eps_val))';

    % Homotopic FC (diagonal of paired interhemispheric matrix)
    feat_homotopic(s,:) = diag(LR_p)';

    % Interhemispheric outgoing connectivity AI
    ihem_L = mean(LR_p, 2);   % mean connectivity from each left ROI to all right
    ihem_R = mean(LR_p, 1)';  % mean connectivity from each right ROI to all left
    feat_ihem_AI(s,:) = ((ihem_L - ihem_R) ./ (abs(ihem_L) + abs(ihem_R) + eps_val))';

    % IFSa: mean interhemispheric connectivity (unilateral, left-only)
    feat_ifsa_ihem(s) = mean(LR(ifsa_idx, :));

    % -------------------------------------------------------------------
    %  BLOCKS 5-8: Entropy (from raw ROI time series)
    % -------------------------------------------------------------------
    roi_file = fullfile(base_dir, sub_id, 'R_mean_roi.mat');

    if ~isfile(roi_file)
        warning('R_mean_roi.mat missing for %s — entropy set to NaN.', sub_id);
        feat_sampen_AI(s,:)   = NaN;
        feat_sampen_diff(s,:) = NaN;
        feat_apen_AI(s,:)     = NaN;
        feat_apen_diff(s,:)   = NaN;
        continue;
    end

    loaded = load(roi_file);
    R_all  = loaded.R_mean_roi;   % 204 x T

    % Split into hemispheres then slice to paired ROIs
    TS_left_all  = R_all(1:n_L, :);          % 100 x T
    TS_right_all = R_all(n_L+1:end, :);      % 104 x T

    TS_L = TS_left_all(paired_left_idx, :);   % 99 x T
    TS_R = TS_right_all(paired_right_idx, :); % 99 x T

    % Compute entropy per paired ROI
    SampEn_L = zeros(n_pairs, 1);
    SampEn_R = zeros(n_pairs, 1);
    ApEn_L   = zeros(n_pairs, 1);
    ApEn_R   = zeros(n_pairs, 1);

    for roi = 1:n_pairs
        ts_l = TS_L(roi, :)';
        ts_r = TS_R(roi, :)';

        r_l = SampEn_r * std(ts_l);
        r_r = SampEn_r * std(ts_r);

        SampEn_L(roi) = compute_sampen(ts_l, SampEn_m, r_l);
        SampEn_R(roi) = compute_sampen(ts_r, SampEn_m, r_r);
        ApEn_L(roi)   = compute_apen(ts_l, ApEn_m, ApEn_r * std(ts_l));
        ApEn_R(roi)   = compute_apen(ts_r, ApEn_m, ApEn_r * std(ts_r));
    end

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

roi_idx = 1:n_pairs;
mk = @(prefix) arrayfun(@(i) sprintf('%s_R%d', prefix, i), roi_idx, 'UniformOutput', false);

col_names = [mk('nodeAI'), mk('nodeDiff'), mk('homotopic'), mk('ihemAI'), ...
             mk('sampEnAI'), mk('sampEnDiff'), mk('apEnAI'), mk('apEnDiff'), ...
             {'IFSa_ihem'}];

all_feats = [feat_node_AI, feat_node_diff, feat_homotopic, feat_ihem_AI, ...
             feat_sampen_AI, feat_sampen_diff, feat_apen_AI, feat_apen_diff, ...
             feat_ifsa_ihem];

T = array2table(all_feats, 'VariableNames', col_names);
T.subject_id = ids_final(:);
T.LI_GEN     = LI_GEN;
T = [T(:,end-1), T(:,1:end-2), T(:,end)];   % subject_id first, LI_GEN last

out_csv = fullfile(out_dir, 'features_asymmetrical.csv');
writetable(T, out_csv);
fprintf('Saved: %s\n', out_csv);

% Block index
block_info = {
    'nodeAI',       n_pairs, 'Node FC asymmetry index (L-R)/(|L|+|R|+eps), 99 paired ROIs';
    'nodeDiff',     n_pairs, 'Node FC raw difference L-R, 99 paired ROIs';
    'homotopic',    n_pairs, 'Homotopic interhemispheric FC (diagonal of paired LR matrix)';
    'ihemAI',       n_pairs, 'Interhemispheric outgoing connectivity AI, 99 paired ROIs';
    'sampEnAI',     n_pairs, 'Sample entropy asymmetry index, 99 paired ROIs';
    'sampEnDiff',   n_pairs, 'Sample entropy raw difference L-R, 99 paired ROIs';
    'apEnAI',       n_pairs, 'Approximate entropy asymmetry index, 99 paired ROIs';
    'apEnDiff',     n_pairs, 'Approximate entropy raw difference L-R, 99 paired ROIs';
    'IFSa_ihem',    1,       'IFSa mean interhemispheric FC (unpaired left-only language ROI)';
};
block_table = cell2table(block_info, 'VariableNames', {'block','n_features','description'});
writetable(block_table, fullfile(out_dir, 'feature_block_index_asymmetrical.csv'));

%% =========================================================================
%  SUMMARY
% =========================================================================
fprintf('\n========== SUMMARY ==========\n');
fprintf('Subjects:        %d\n', n_valid);
fprintf('Paired ROIs:     %d\n', n_pairs);
fprintf('Total features:  %d\n', width(T)-2);
fprintf('Output:          %s\n', out_dir);
fprintf('==============================\n');
disp(block_table);

%% =========================================================================
%  LOCAL FUNCTIONS
% =========================================================================

function se = compute_sampen(x, m, r)
% Sample entropy of time series x with template length m and tolerance r.
    N = length(x);
    if N < (m+2) || r == 0
        se = NaN; return;
    end
    B = 0; A = 0;
    for i = 1:N-m
        xi = x(i:i+m-1);
        for j = i+1:N-m
            if max(abs(xi - x(j:j+m-1))) < r
                B = B + 1;
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
