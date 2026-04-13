%% PARAMETERS
%base_slic_dir = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';
base_slic_dir = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_asymmetrical'
base_rs_dir   = '/home/marlene/Desktop/data_victoria/EZM2';

% mask_left          = '/home/marlene/Desktop/victoria/atlas/grey/atlas_SENS_HCP_Tian_gm_LEFT_orig_fixed.nii';
% mask_left_mirrored = '/home/marlene/Desktop/victoria/atlas/grey/flipped_atlas_SENS_HCP_Tian_gm_LEFT_fixed.nii';

K = 600;

subject_dirs = dir(fullfile(base_slic_dir, 'sub-EZM2*'));
subject_dirs = subject_dirs([subject_dirs.isdir]);
fprintf('Found %d EZM2 subjects\n', length(subject_dirs));

%% LOOP OVER SUBJECTS
for s = 1:length(subject_dirs)
    sub_id = subject_dirs(s).name;
    fprintf('\nProcessing %s (%d/%d)...\n', sub_id, s, length(subject_dirs));

    %% DEFINE PATHS FOR THIS SUBJECT
    slic_dir    = fullfile(base_slic_dir, sub_id);
    slic        = fullfile(slic_dir, sprintf('SLIC_K%d.nii', K));
    RS_d        = fullfile(base_rs_dir, sub_id, 'RS', sprintf('dswuRS_%s.nii', sub_id));
    overlap_mat = fullfile(slic_dir, sprintf('left_overlap_results_%d.mat', K));
    out_file    = fullfile(slic_dir, 'R_mean_roi.mat');

    %% CHECK FILES EXIST
    if ~isfile(slic)
        warning('SLIC file missing for %s, skipping.', sub_id);
        continue;
    end
    if ~isfile(RS_d)
        warning('RS file missing for %s, skipping.', sub_id);
        continue;
    end
    if ~isfile(overlap_mat)
        warning('Overlap .mat missing for %s, skipping.', sub_id);
        continue;
    end

    %% LOAD OVERLAP RESULTS
    loaded  = load(overlap_mat);
    results = loaded.results;

    %% CALCULATE SUPERVOXEL SIGNAL
    V = spm_vol(slic);
    S = spm_read_vols(V);

    V1 = spm_vol(RS_d);
    D  = spm_read_vols(V1);

    fprintf('RS dimensions : %d x %d x %d x %d\n', size(D,1), size(D,2), size(D,3), size(D,4));
    fprintf('SLIC dimensions: %d x %d x %d\n',     size(S,1), size(S,2), size(S,3));


    %% COMPUTE R_mean PER CLUSTER
    R_mean = zeros(max(S(:)), size(D,4));

    for i = 1:max(S(:))
        vox_idx = find(S == i);
        if isempty(vox_idx)
            continue;
        end
        [x,y,z] = ind2sub(size(S), vox_idx);
        R = zeros(length(x), size(D,4));
        for j = 1:length(x)
            R(j,:) = squeeze(D(x(j), y(j), z(j), :));
        end
        R_mean(i,:) = mean(R, 1);
    end

    %% COMPUTE R_mean_roi PER ATLAS ROI
    ID_atlas        = results(:,2);
    ID_atlas_unique = unique(ID_atlas);

    R_mean_roi = zeros(length(ID_atlas_unique), size(R_mean,2));

    for i = 1:length(ID_atlas_unique)
        roi          = find(ID_atlas == ID_atlas_unique(i));
        R_mean_roi_j = zeros(length(roi), size(R_mean,2));
        for j = 1:length(roi)
            cluster_id = results(roi(j), 3);
            if cluster_id > 0
                weight = results(roi(j), 5);
                R_mean_roi_j(j,:) = R_mean(cluster_id,:) * weight;
            end
        end
        R_mean_roi(i,:) = sum(R_mean_roi_j, 1);
    end

    %% SAVE
    R_mean_roi = R_mean_roi;
    save(out_file, 'R_mean_roi');
    fprintf('Saved: %s\n', out_file);
end

fprintf('\nDone! Processed %d subjects.\n', length(subject_dirs));