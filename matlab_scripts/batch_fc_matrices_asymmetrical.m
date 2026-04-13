%% batch_fc_matrices_asymmetrical.m
% Computes FC matrices from the asymmetric SENS-HCP-Tian atlas.
% Atlas has 99 left ROIs and 103 right ROIs (204 total in-mask).
% All three FC matrices are saved per subject:
%   FC_matrix : 99x103  interhemispheric (left vs right)
%   FC_left   : 99x99   left within-hemisphere
%   FC_right  : 103x103 right within-hemisphere
%
% Paired indices (98 matched L/R homologs) are saved alongside
% for use at the feature extraction stage — no need to rerun
% this script if you want to switch between full and paired features.

%% SUBJECTS TO EXCLUDE FROM GROUP MATRICES (individual files still saved)
exclude_list = {'sub-EZM20001', 'sub-EZM20002', 'sub-EZM20003', 'sub-EZM20004', ...
    'sub-EZM20005', 'sub-EZM20006', 'sub-EZM20007', 'sub-EZM20008', 'sub-EZM20009', ...
    'sub-EZM20024', 'sub-EZM20028', 'sub-EZM20033', 'sub-EZM20045', 'sub-EZM20046', ...
    'sub-EZM20059', 'sub-EZM20085', 'sub-EZM20086', 'sub-EZM20094', 'sub-EZM20099', ...
    'sub-EZM20104', 'sub-EZM20120', 'sub-EZM20136', 'sub-EZM20137', 'sub-EZM20146', ...
    'sub-EZM20148', 'sub-EZM20152', 'sub-EZM20153', 'sub-EZM20155', 'sub-EZM20168', ...
    'sub-EZM20176', 'sub-EZM20052', 'sub-EZM20057', 'sub-EZM20091', 'sub-EZM20112', ...
    'sub-EZM20115', 'sub-EZM20186'};

%% PARAMETERS
base_dir   = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_asymmetrical';
output_dir = '/home/marlene/Desktop/victoria/FC_matrices_asymmetrical';
mkdir(output_dir);

% ROI counts from asymmetric atlas (region_lookup.xlsx)
n_L = 100;   % left ROIs  (atlas indices 1-99)
n_R = 104;  % right ROIs (atlas indices 100-202 in R_mean_roi ordering)

% Paired homolog indices — derived from region_lookup.xlsx
% These index into R_mean_roi rows (1-based, ordered as in the lookup table)
% 98 matched pairs; left ROI 21 (IFSa) has no right homolog
paired_left_idx  = [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, ...
    18, 19, 20, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, ...
    38, 39, 40, 41, 42, 43, 44, 45, 46, 95, 97, 99, 101, 103, 105, 107, 109, 111, 113, ...
    115, 117, 119, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, ...
    147, 149, 151, 153, 155, 157, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, ...
    193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204 ];

paired_right_idx = [ 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, ...
    63, 64, 65, 66, 69, 70, 71, 72, 73, 74, 75, 77, 78, 79, 80, 81, 82, 83, 84, 85, ...
    86, 87, 88, 89, 90, 91, 92, 93, 94, 96, 98, 100, 102, 104, 106, 108, 110, 112, 114, ...
    116, 118, 120, 122, 124, 126, 128, 130, 132, 134, 136, 138, 140, 142, 144, 146, ...
    148, 150, 152, 154, 156, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, ...
    169, 170, 171, 173, 174, 175, 176, 177, 178, 179, 180, 181 ];

% Index of left-only unpaired ROI (IFSa — language relevant, kept separately)
ifsa_idx = 21;

%% DISCOVER SUBJECTS
subject_dirs = dir(fullfile(base_dir, 'sub-EZM*'));
subject_dirs = subject_dirs([subject_dirs.isdir]);
fprintf('Found %d subjects\n', length(subject_dirs));

%% PREALLOCATE GROUP MATRICES
n_subjects       = length(subject_dirs);
n_inter          = n_L * n_R;        % 99x103 = 10197
n_left           = n_L * n_L;        % 99x99  = 9801
n_right          = n_R * n_R;        % 103x103 = 10609

full_fc_matrix   = zeros(n_subjects, n_inter);
full_fc_left     = zeros(n_subjects, n_left);
full_fc_right    = zeros(n_subjects, n_right);
subject_ids      = cell(n_subjects, 1);
valid            = false(n_subjects, 1);

%% MAIN LOOP
for s = 1:n_subjects
    sub_id = subject_dirs(s).name;
    fprintf('\nProcessing %s (%d/%d)...\n', sub_id, s, n_subjects);

    roi_file = fullfile(base_dir, sub_id, 'R_mean_roi.mat');

    if ~isfile(roi_file)
        warning('R_mean_roi.mat missing for %s, skipping.', sub_id);
        continue;
    end

    % Load ROI time series — shape: (n_L + n_R) x timepoints = 202 x T
    loaded = load(roi_file);
    R_all  = loaded.R_mean_roi;

    % Check total ROI count
    if size(R_all, 1) ~= (n_L + n_R)
        warning('Unexpected ROI count for %s: got %d rows, expected %d. Skipping.', ...
            sub_id, size(R_all, 1), n_L + n_R);
        continue;
    end

    % Split into hemispheres
    % Assumes R_mean_roi_full rows are ordered: left ROIs (1:99) then right ROIs (100:202)
    R_left  = R_all(1:n_L, :);           % 99  x T
    R_right = R_all(n_L+1:end, :);       % 103 x T

    %% COMPUTE FC MATRICES
    FC_matrix = corr(R_left', R_right');  % 99x103 interhemispheric
    FC_left   = corr(R_left');            % 99x99  left within-hemisphere
    FC_right  = corr(R_right');           % 103x103 right within-hemisphere

    % Fisher z-transform (atanh; diagonal of within-hemisphere matrices = Inf, expected)
    FC_matrix_z = atanh(FC_matrix);
    FC_left_z   = atanh(FC_left);
    FC_right_z  = atanh(FC_right);

    %% SAVE PER-SUBJECT (all subjects, even excluded ones)
    sub_out = fullfile(output_dir, sub_id);
    mkdir(sub_out);
    save(fullfile(sub_out, 'FC_matrix.mat'), 'FC_matrix', 'FC_matrix_z');
    save(fullfile(sub_out, 'FC_left.mat'),   'FC_left',   'FC_left_z');
    save(fullfile(sub_out, 'FC_right.mat'),  'FC_right',  'FC_right_z');

    %% ADD TO GROUP MATRICES (excluded subjects skipped)
    if ~ismember(sub_id, exclude_list)
        full_fc_matrix(s, :) = FC_matrix_z(:)';
        full_fc_left(s, :)   = FC_left_z(:)';
        full_fc_right(s, :)  = FC_right_z(:)';
        subject_ids{s}       = sub_id;
        valid(s)             = true;
        fprintf('Done and added to group matrices: %s\n', sub_id);
    else
        fprintf('Done - individual saved but excluded from group matrices: %s\n', sub_id);
    end
end

%% TRIM EMPTY ROWS
full_fc_matrix = full_fc_matrix(valid, :);
full_fc_left   = full_fc_left(valid, :);
full_fc_right  = full_fc_right(valid, :);
subject_ids    = subject_ids(valid);

fprintf('\n%d / %d subjects included in group matrices.\n', sum(valid), n_subjects);

%% SAVE GROUP MATRICES AS .mat
% Also save the paired indices so feature extraction scripts know which
% rows/columns correspond to matched L/R homologs
save(fullfile(output_dir, 'full_fc_matrix.mat'), ...
    'full_fc_matrix', 'subject_ids', 'paired_left_idx', 'paired_right_idx', 'ifsa_idx');
save(fullfile(output_dir, 'full_fc_left.mat'),   'full_fc_left',  'subject_ids');
save(fullfile(output_dir, 'full_fc_right.mat'),  'full_fc_right', 'subject_ids');

%% SAVE GROUP MATRICES AS .csv
fc_matrix_table = array2table(full_fc_matrix, 'RowNames', subject_ids);
fc_left_table   = array2table(full_fc_left,   'RowNames', subject_ids);
fc_right_table  = array2table(full_fc_right,  'RowNames', subject_ids);

writetable(fc_matrix_table, fullfile(output_dir, 'full_fc_matrix.csv'), 'WriteRowNames', true);
writetable(fc_left_table,   fullfile(output_dir, 'full_fc_left.csv'),   'WriteRowNames', true);
writetable(fc_right_table,  fullfile(output_dir, 'full_fc_right.csv'),  'WriteRowNames', true);

fprintf('\nAll saved to: %s\n', output_dir);
