
% batch_SLIC.m
clear; clc;

%% ── 0. SET WORKING DIRECTORY TO SLIC REPO ───────────────────────────────────
% This is critical — SLIC.m loads parc_graymatter.mat and parc_searchlight.mat
% from the current working directory, so we must run from the repo folder
slic_repo_dir = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master';  % ← EDIT IF NEEDED
cd(slic_repo_dir);
fprintf('Working directory set to: %s\n', pwd);

%% ── 1. USER SETTINGS ────────────────────────────────────────────────────────
% Parent directory containing all subject folders
data_dir = '/home/marlene/Desktop/data_victoria/EZM2';

% Brain mask (same for all subjects)
file_mask = '/home/marlene/Desktop/victoria/atlas/grey/atlas_SENS_HCP_Tian_gm.nii'
central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_asymmetrical'

% % using the left  mask
% file_mask = '/home/marlene/Desktop/victoria/atlas/grey/atlas_SENS_HCP_Tian_gm_LEFT_orig_fixed.nii';
% central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_dswu';

% % using the left flipped mask
% file_mask = '/home/marlene/Desktop/victoria/atlas/grey/flipped_atlas_SENS_HCP_Tian_gm_LEFT_fixed.nii';
% central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';

% % using the right mask
% file_mask = '/home/marlene/Desktop/victoria/atlas/grey/atlas_SENS_HCP_Tian_gm_RIGHT_orig_fixed.nii';
% central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_right_dswu';
% 
% % using the right flipped mask
% file_mask = '/home/marlene/Desktop/victoria/atlas/grey/flipped_atlas_SENS_HCP_Tian_gm_RIGHT_fixed.nii';
% central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_right_flipped_dswu';
% 

% SLIC settings
sK = [600]; 
% change this
m  = 40;

% ── TEST MODE ────────────────────────────────────────────────────────────────
test_mode = false;
n_test    = 1;

%% ── DISCOVER SUBJECTS AUTOMATICALLY ─────────────────────────────────────────
fprintf('Scanning for subjects in: %s\n', data_dir);
all_entries = dir(fullfile(data_dir, 'sub-EZM*'));
all_entries = all_entries([all_entries.isdir]);
subject_ids = {all_entries.name};

if isempty(subject_ids)
    error('No subject folders found in %s', data_dir);
end
fprintf('Found %d subjects total.\n', length(subject_ids));

if test_mode
    subject_ids = subject_ids(1:min(n_test, length(subject_ids)));
    fprintf('TEST MODE: running first %d subjects only.\n', length(subject_ids));
else
    fprintf('FULL MODE: running all %d subjects.\n', length(subject_ids));
end

%% ── LOAD MASK ONCE ───────────────────────────────────────────────────────────
fprintf('\nLoading brain mask...\n');
nii_mask = load_untouch_nii(file_mask);
img_mask = double(nii_mask.img);

msk_gray = img_mask ~= 0;
ind_gray = find(msk_gray);
num_gray = length(ind_gray);
siz      = size(msk_gray);

fprintf('Mask size: %d x %d x %d  |  Gray matter voxels: %d\n', ...
        siz(1), siz(2), siz(3), num_gray);

% Save to working directory (where SLIC.m will look for them)
save('parc_graymatter.mat', 'msk_gray', 'ind_gray', 'num_gray', 'siz');
fprintf('Saved parc_graymatter.mat to: %s\n', pwd);

% Build searchlight
r = 1; shape = 'cube';
[riSL, nvSL] = parc_searchlight(r, shape);
save('parc_searchlight.mat', 'riSL', 'nvSL');
fprintf('Saved parc_searchlight.mat to: %s\n', pwd);

%% ── TRACKING ─────────────────────────────────────────────────────────────────
failed_subjects  = {};
skipped_subjects = {};
n_total = length(subject_ids);

%% ── MAIN LOOP ────────────────────────────────────────────────────────────────
for s = 1:n_total
    sub_id = subject_ids{s};
    fprintf('\n========================================\n');
    fprintf('Subject %d / %d : %s\n', s, n_total, sub_id);
    fprintf('========================================\n');

    % Build paths — same structure as working single-subject script
    fmri_path  = fullfile(data_dir, sub_id, 'RS', sprintf('dswuRS_%s.nii', sub_id)); 
    output_dir = fullfile(central_output, sub_id);

    % Check fMRI file exists
    if ~isfile(fmri_path)
        fprintf('WARNING: fMRI file not found, skipping.\n');
        fprintf('  Expected: %s\n', fmri_path);
        skipped_subjects{end+1} = sub_id;
        continue
    end

    try
        %% LOAD fMRI
        fprintf('Loading fMRI...\n');
        nii_fmri   = load_untouch_nii(fmri_path);
        img_rs_4d  = double(nii_fmri.img);
        num_trial  = size(img_rs_4d, 4);
        N          = prod(siz);
        fprintf('fMRI header dims: %d x %d x %d\n', size(img_rs_4d,1), size(img_rs_4d,2), size(img_rs_4d,3));
        fprintf('Mask header dims: %d x %d x %d\n', siz(1), siz(2), siz(3));
        fprintf('fMRI pixdim: %.4f %.4f %.4f\n', nii_fmri.hdr.dime.pixdim(2), nii_fmri.hdr.dime.pixdim(3), nii_fmri.hdr.dime.pixdim(4));
        fprintf('Mask pixdim: %.4f %.4f %.4f\n', nii_mask.hdr.dime.pixdim(2), nii_mask.hdr.dime.pixdim(3), nii_mask.hdr.dime.pixdim(4));
        fprintf('fMRI datatype: %d\n', nii_fmri.hdr.dime.datatype);
        fprintf('Mask datatype: %d\n', nii_mask.hdr.dime.datatype);

        %% VERIFY DIMENSIONS
        siz_fmri = size(img_rs_4d(:,:,:,1));
        if ~isequal(siz_fmri, siz)
            error('DIMENSION MISMATCH: Mask %dx%dx%d vs fMRI %dx%dx%d', ...
                  siz(1),siz(2),siz(3), siz_fmri(1),siz_fmri(2),siz_fmri(3));
        end
        fprintf('✓ Dimensions match: %d x %d x %d  |  Timepoints: %d\n', ...
                siz(1), siz(2), siz(3), num_trial);

        %% PREPROCESS — identical to working script
        img_rs = reshape(img_rs_4d, N, num_trial);

        flag = zeros(N, 1);
        for i = 1:N
            flag(i) = std(img_rs(i, :));
        end
        flag      = flag ~= 0;
        msk_rs    = reshape(flag, siz);
        msk_merge = msk_rs .* msk_gray;
        ind_merge = find(msk_merge == 1);

        assert(~isempty(ind_merge), ...
            'No valid voxels in mask/fMRI overlap.');

        tmp = img_rs(ind_merge, :);
        tmp = tmp - repmat(mean(tmp, 2), [1, num_trial]);
        tmp = tmp ./ repmat(sqrt(sum(tmp .^ 2, 2)), [1, num_trial]);

        [~, index] = ismember(ind_merge, ind_gray);
        img_gray_prep = zeros(num_gray, num_trial);
        img_gray_prep(index, :) = tmp;

        fprintf('Valid gray matter voxels: %d\n', length(ind_merge));

        %% RUN SLIC — identical call to working script
        mkdir(output_dir);
        nK = length(sK);

        for iK = 1:nK
            cK = sK(iK);
            fprintf('\n=== Running SLIC with K = %d ===\n', cK);

            label = SLIC(img_gray_prep, m, cK);

            % Map back to 3D — using fresh variable name to avoid conflicts
            % img_parc = zeros(siz);
            % [~, idx_map] = ismember(ind_merge, ind_gray);
            % valid_labels = label(idx_map);
            % img_parc(ind_merge) = valid_labels;
            % Map back to 3D safely
            % Map back to 3D
            img_parc = zeros(siz);
            img_parc(ind_gray) = label;

            actual_K = length(unique(label));
            fprintf('Requested K = %d  |  Actual clusters = %d\n', cK, actual_K);
            fprintf('Non-zero voxels in parcellation: %d\n', sum(img_parc(:) > 0));

            % Save .mat
            mat_out = fullfile(output_dir, sprintf('SLIC_K%d.mat', cK));
            save(mat_out, 'img_parc', 'actual_K', 'cK');

            % Save .nii — identical to working script
            nii_out     = nii_mask;
            nii_out.img = img_parc;
            nii_out.hdr.dime.dim(1) = 3;
            nii_out.hdr.dime.dim(5) = 1;
            nii_filename = fullfile(output_dir, sprintf('SLIC_K%d.nii', cK));
            save_untouch_nii(nii_out, nii_filename);
            fprintf('Saved: %s\n', nii_filename);
        end

        fprintf('✓ Subject %s done.\n', sub_id);

    catch ME
        fprintf('ERROR on subject %s:\n  %s\n', sub_id, ME.message);
        failed_subjects{end+1} = sub_id;
    end

    % Clear subject variables to free memory
    clear img_rs_4d img_rs img_gray_prep tmp

end

%% ── SUMMARY ──────────────────────────────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('BATCH COMPLETE\n');
fprintf('  Processed : %d subjects\n', n_total - length(failed_subjects) - length(skipped_subjects));
fprintf('  Skipped   : %d  (fMRI file not found)\n', length(skipped_subjects));
fprintf('  Failed    : %d  (error during processing)\n', length(failed_subjects));

if ~isempty(skipped_subjects)
    fprintf('\nSkipped subjects:\n');
    for i = 1:length(skipped_subjects)
        fprintf('  %s\n', skipped_subjects{i});
    end
end

if ~isempty(failed_subjects)
    fprintf('\nFailed subjects:\n');
    for i = 1:length(failed_subjects)
        fprintf('  %s\n', failed_subjects{i});
    end
end
fprintf('========================================\n');