%% paths
% file_mask = '/home/marlene/Desktop/victoria/atlas/grey/atlas_SENS_HCP_Tian_gm_LEFT_orig_fixed.nii';
% central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_dswu';

% file_mask = '/home/marlene/Desktop/victoria/atlas/grey/flipped_atlas_SENS_HCP_Tian_gm_LEFT_fixed.nii';
% central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';

file_mask = '/home/marlene/Desktop/victoria/atlas/grey/atlas_SENS_HCP_Tian_gm.nii'
central_output = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_asymmetrical'

%% load the mask once
mask_nii = load_untouch_nii(file_mask);
mask_nii.img = double(mask_nii.img);
mask_vol = mask_nii.img;
roi_ids  = unique(mask_vol);
roi_ids  = roi_ids(roi_ids ~= 0);

%% loop over subjects
cd(central_output);
data = dir('sub*');
data = data([data.isdir]);

for i = 1:length(data)
    sub_id = data(i).name;
    fprintf('\n=== Subject %d / %d : %s ===\n', i, length(data), sub_id);

    % Skip if already done
    save_path_mat = fullfile(central_output, sub_id, 'left_overlap_results_600.mat');
    save_path_csv = fullfile(central_output, sub_id, 'left_overlap_results_600.csv');
    if isfile(save_path_mat)
        fprintf('Already exists, skipping.\n');
        continue;
    end

    % Load clustering
    clust_path = fullfile(central_output, sub_id, 'SLIC_K600.mat');
    if ~isfile(clust_path)
        fprintf('WARNING: SLIC_K600.mat not found for %s, skipping.\n', sub_id);
        continue;
    end
    clust_data = load(clust_path);
    clust_vol  = clust_data.img_parc;

    % Diagnostic info
    brain_mask = (mask_vol ~= 0);
    clust_mask = (clust_vol ~= 0);
    fprintf('Non-zero voxels in atlas mask : %d\n', sum(brain_mask(:)));
    fprintf('Non-zero voxels in clustering : %d\n', sum(clust_mask(:)));
    fprintf('Voxels in both               : %d\n', sum(brain_mask(:) & clust_mask(:)));

    cluster_ids = unique(clust_vol);
    cluster_ids = cluster_ids(cluster_ids ~= 0);

    % Reset results for this subject
    results = [];

    for r = 1:length(roi_ids)
        roi_id       = roi_ids(r);
        roi_voxels   = (mask_vol == roi_id);
        n_roi_voxels = sum(roi_voxels(:));

        for c = 1:length(cluster_ids)
            cluster_id       = cluster_ids(c);
            cluster_voxels   = (clust_vol == cluster_id);
            overlap_voxels   = roi_voxels & cluster_voxels;
            n_overlap        = sum(overlap_voxels(:));

            if n_overlap > 0
                pct_of_roi       = n_overlap / n_roi_voxels;
                n_cluster_voxels = sum(cluster_voxels(:));
                pct_of_cluster   = n_overlap / n_cluster_voxels;
                results(end+1, :) = [i, roi_id, cluster_id, n_overlap, pct_of_roi, pct_of_cluster];
            end
        end
    end

    % Save to this subject's folder
    save(save_path_mat, 'results');
    header = {'subject_index','roi_id','cluster_id','n_voxels_overlap','pct_of_roi','pct_of_cluster'};
    results_table = array2table(results, 'VariableNames', header);
    writetable(results_table, save_path_csv);
    fprintf('Saved to: %s\n', save_path_mat);
end

disp('All done.');