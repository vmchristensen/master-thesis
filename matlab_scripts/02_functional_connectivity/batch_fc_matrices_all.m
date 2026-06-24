% create FC matrix

%% SUBJECTS TO EXCLUDE FROM GROUP MATRICES (individual files still saved)
exclude_list = {'sub-EZM20001', 'sub-EZM20002', 'sub-EZM20003', 'sub-EZM20004', 'sub-EZM20005', 'sub-EZM20006', 'sub-EZM20007', 'sub-EZM20008', 'sub-EZM20009','sub-EZM20024', 'sub-EZM20028', 'sub-EZM20033', 'sub-EZM20045', 'sub-EZM20046', 'sub-EZM20059', 'sub-EZM20085', 'sub-EZM20086', 'sub-EZM20094', 'sub-EZM20099', 'sub-EZM20104', 'sub-EZM20120', 'sub-EZM20136', 'sub-EZM20137', 'sub-EZM20146', 'sub-EZM20148', 'sub-EZM20152', 'sub-EZM20153', 'sub-EZM20155', 'sub-EZM20168', 'sub-EZM20176', 'sub-EZM20052', 'sub-EZM20057', 'sub-EZM20091', 'sub-EZM20112', 'sub-EZM20115', 'sub-EZM20186'};  % add your subjects here

%% PARAMETERS
base_left_dir  = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_dswu';
base_right_dir = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';
output_dir     = '/home/marlene/Desktop/victoria/FC_matrices';
mkdir(output_dir);

%% DISCOVER SUBJECTS
subject_dirs = dir(fullfile(base_left_dir, 'sub-EZM*'));
subject_dirs = subject_dirs([subject_dirs.isdir]);
fprintf('Found %d subjects\n', length(subject_dirs));

%% PREALLOCATE
n_subjects = length(subject_dirs);
n_rois     = 101;
n_elements = n_rois * n_rois;

full_fc_matrix = zeros(n_subjects, n_elements);
full_fc_left   = zeros(n_subjects, n_elements);
full_fc_right  = zeros(n_subjects, n_elements);
subject_ids    = cell(n_subjects, 1);
valid          = false(n_subjects, 1);

%% LOOP OVER SUBJECTS
for s = 1:n_subjects
    sub_id = subject_dirs(s).name;
    fprintf('\nProcessing %s (%d/%d)...\n', sub_id, s, n_subjects);

    left_file  = fullfile(base_left_dir,  sub_id, 'R_mean_roi_left.mat');
    right_file = fullfile(base_right_dir, sub_id, 'R_mean_roi_right.mat');

    % Check files exist
    if ~isfile(left_file)
        warning('Left file missing for %s, skipping.', sub_id);
        continue;
    end
    if ~isfile(right_file)
        warning('Right file missing for %s, skipping.', sub_id);
        continue;
    end

    % Load
    left_data  = load(left_file);
    right_data = load(right_file);
    R_left     = left_data.R_mean_roi_left;
    R_right    = right_data.R_mean_roi_right;

    % Check dimensions
    if size(R_left,1) ~= n_rois || size(R_right,1) ~= n_rois
        warning('Unexpected dimensions for %s: left=%dx%d, right=%dx%d, skipping.', ...
            sub_id, size(R_left,1), size(R_left,2), size(R_right,1), size(R_right,2));
        continue;
    end

    % Compute FC matrices
    FC_matrix = corr(R_left', R_right');   % 101x101 left vs right
    FC_left   = corr(R_left');             % 101x101 left vs left
    FC_right  = corr(R_right');            % 101x101 right vs right

    % Fisher z-transform
    FC_matrix_z = atanh(FC_matrix);
    FC_left_z   = atanh(FC_left);
    FC_right_z  = atanh(FC_right);

    % Save per-subject FC matrices (ALL subjects, even excluded ones)
    sub_out = fullfile(output_dir, sub_id);
    mkdir(sub_out);
    save(fullfile(sub_out, 'FC_matrix.mat'), 'FC_matrix', 'FC_matrix_z');
    save(fullfile(sub_out, 'FC_left.mat'),   'FC_left',   'FC_left_z');
    save(fullfile(sub_out, 'FC_right.mat'),  'FC_right',  'FC_right_z');

    % Only add to group matrices if NOT in exclude list
    if ~ismember(sub_id, exclude_list)
        full_fc_matrix(s, :) = FC_matrix_z(:)';
        full_fc_left(s, :)   = FC_left_z(:)';
        full_fc_right(s, :)  = FC_right_z(:)';
        subject_ids{s}       = sub_id;
        valid(s)             = true;
        fprintf('Done and added to group matrices: %s\n', sub_id);
    else
        fprintf('Done - individual files saved but excluded from group matrices: %s\n', sub_id);
    end
end

%% REMOVE EMPTY ROWS (skipped or excluded subjects)
full_fc_matrix = full_fc_matrix(valid, :);
full_fc_left   = full_fc_left(valid, :);
full_fc_right  = full_fc_right(valid, :);
subject_ids    = subject_ids(valid);

fprintf('\n%d / %d subjects included in group matrices.\n', sum(valid), n_subjects);

%% SAVE GROUP MATRICES AS .mat
save(fullfile(output_dir, 'full_fc_matrix.mat'), 'full_fc_matrix', 'subject_ids');
save(fullfile(output_dir, 'full_fc_left.mat'),   'full_fc_left',   'subject_ids');
save(fullfile(output_dir, 'full_fc_right.mat'),  'full_fc_right',  'subject_ids');

%% SAVE GROUP MATRICES AS .csv
fc_matrix_table = array2table(full_fc_matrix, 'RowNames', subject_ids);
fc_left_table   = array2table(full_fc_left,   'RowNames', subject_ids);
fc_right_table  = array2table(full_fc_right,  'RowNames', subject_ids);

writetable(fc_matrix_table, fullfile(output_dir, 'full_fc_matrix.csv'), 'WriteRowNames', true);
writetable(fc_left_table,   fullfile(output_dir, 'full_fc_left.csv'),   'WriteRowNames', true);
writetable(fc_right_table,  fullfile(output_dir, 'full_fc_right.csv'),  'WriteRowNames', true);

fprintf('\nAll saved to: %s\n', output_dir);
