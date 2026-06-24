% INSIDEOUT_final_results.m
% Outputs: IrrLeft, IrrRight, IrrLI, IrrWhole (atlas-wide) at tau=1
% Explicitly handles 101 vs 100 ROI ambiguity
% ------EDIT BEFORE RUNNING------

base_dir_left  = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_dswu';
base_dir_right = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';
results_dir    = '/home/marlene/Desktop/victoria/INSIDEOUT_results';

TAU    = 1;
N_ROIS = 100;  % expected ROIs per hemisphere

exclude_list = {'sub-EZM20001', 'sub-EZM20002', 'sub-EZM20003', 'sub-EZM20004', ...
                'sub-EZM20005', 'sub-EZM20006', 'sub-EZM20007', 'sub-EZM20008', ...
                'sub-EZM20009', 'sub-EZM20024', 'sub-EZM20028', 'sub-EZM20033', ...
                'sub-EZM20045', 'sub-EZM20046', 'sub-EZM20059', 'sub-EZM20085', ...
                'sub-EZM20086', 'sub-EZM20094', 'sub-EZM20099', 'sub-EZM20104', ...
                'sub-EZM20120', 'sub-EZM20136', 'sub-EZM20137', 'sub-EZM20146', ...
                'sub-EZM20148', 'sub-EZM20152', 'sub-EZM20153', 'sub-EZM20155', ...
                'sub-EZM20168', 'sub-EZM20176', 'sub-EZM20052', 'sub-EZM20057', ...
                'sub-EZM20091', 'sub-EZM20112', 'sub-EZM20115', 'sub-EZM20186', ...
                'sub-EZM20014'};

% ------STOP EDIT------

if ~exist(results_dir, 'dir'), mkdir(results_dir); end

subj_dirs    = dir(fullfile(base_dir_left, 'sub-EZM*'));
all_subjects = {subj_dirs.name};
subjects     = all_subjects(~ismember(all_subjects, exclude_list));
NSUB         = length(subjects);
fprintf('Subjects after exclusion: %d\n', NSUB);

% Preallocate
SubjectID = cell(NSUB, 1);
IrrLeft   = nan(NSUB, 1);
IrrRight  = nan(NSUB, 1);
IrrLI     = nan(NSUB, 1);
IrrWhole  = nan(NSUB, 1);
N_rois_used_L = nan(NSUB, 1);
N_rois_used_R = nan(NSUB, 1);

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

    dl = load(left_file);
    dr = load(right_file);
    ts_L = dl.R_mean_roi_left;   % should be N_ROIS x T
    ts_R = dr.R_mean_roi_right;  % should be N_ROIS x T

    % --- Handle 101 vs 100 ROI ambiguity ---
    % Row dimension is ROIs, column dimension is timepoints.
    % If the matrix has 101 rows, drop the last row and warn once.
    if size(ts_L, 1) == 101
        if sub == 1
            warning(['R_mean_roi_left has 101 rows. ' ...
                     'Dropping row 101 (assumed summary/noise). ' ...
                     'Verify this is correct before reporting results.']);
        end
        ts_L = ts_L(1:100, :);
    end
    if size(ts_R, 1) == 101
        ts_R = ts_R(1:100, :);
    end

    % Validate dimensions
    nL = size(ts_L, 1);
    nR = size(ts_R, 1);
    N_rois_used_L(sub) = nL;
    N_rois_used_R(sub) = nR;

    if nL ~= N_ROIS || nR ~= N_ROIS
        warning('Unexpected ROI count for %s: L=%d R=%d — skipping', subj, nL, nR);
        continue
    end

    TmL = size(ts_L, 2);
    TmR = size(ts_R, 2);
    if TmL ~= TmR
        warning('Timepoint mismatch for %s: L=%d R=%d — skipping', subj, TmL, TmR);
        continue
    end
    Tm = TmL;

    if Tm <= TAU
        warning('Time series too short for tau=%d in %s — skipping', TAU, subj);
        continue
    end

    % =========================================================
    % INSIDEOUT scalar: mean(FSdiff^2) for a time series matrix
    % Input: ts is N_rois x T (rows=regions, cols=time)
    % =========================================================
    function val = irr_scalar(ts, tau)
        Tm  = size(ts, 2);
        % Forward: correlate ts[:,1:T-tau] with ts[:,1+tau:T]
        FC_fwd = corr(ts(:, 1:Tm-tau)', ts(:, 1+tau:Tm)');
        % Reversed: flip time axis, same operation
        ts_rev = ts(:, end:-1:1);
        FC_rev = corr(ts_rev(:, 1:Tm-tau)', ts_rev(:, 1+tau:Tm)');
        % FS matrices (mutual information approximation)
        FS_fwd = -0.5 * log(1 - FC_fwd.^2 + 1e-12);
        FS_rev = -0.5 * log(1 - FC_rev.^2 + 1e-12);
        % Squared difference, mean over all pairs
        val = mean((FS_fwd(:) - FS_rev(:)).^2);
    end

    % Compute the four quantities
    IrrLeft(sub)  = irr_scalar(ts_L,                  TAU);
    IrrRight(sub) = irr_scalar(ts_R,                  TAU);
    IrrWhole(sub) = irr_scalar([ts_L; ts_R],           TAU);  % 200 x T
    IrrLI(sub)    = (IrrLeft(sub) - IrrRight(sub)) / ...
                    (IrrLeft(sub) + IrrRight(sub) + 1e-12);

end

%% Diagnostic: confirm ROI counts are consistent
fprintf('\nROI count check:\n');
fprintf('  Left  — unique values: ');
disp(unique(N_rois_used_L(~isnan(N_rois_used_L)))');
fprintf('  Right — unique values: ');
disp(unique(N_rois_used_R(~isnan(N_rois_used_R)))');

%% Save
T_out = table(SubjectID, IrrLeft, IrrRight, IrrLI, IrrWhole, ...
              N_rois_used_L, N_rois_used_R);

% Drop subjects with all NaN (missing files)
valid = ~cellfun(@isempty, SubjectID) & ~isnan(IrrLeft);
T_out = T_out(valid, :);

outfile = fullfile(results_dir, 'INSIDEOUT_final_all_components.csv');
writetable(T_out, outfile);
fprintf('\nSaved %d subjects to: %s\n', height(T_out), outfile);

%% Summary statistics
fprintf('\n=== SUMMARY ===\n');
fprintf('n = %d\n', height(T_out));
fprintf('IrrLeft:  mean=%.6f  std=%.6f  min=%.6f  max=%.6f\n', ...
    mean(T_out.IrrLeft), std(T_out.IrrLeft), ...
    min(T_out.IrrLeft),  max(T_out.IrrLeft));
fprintf('IrrRight: mean=%.6f  std=%.6f  min=%.6f  max=%.6f\n', ...
    mean(T_out.IrrRight), std(T_out.IrrRight), ...
    min(T_out.IrrRight),  max(T_out.IrrRight));
fprintf('IrrLI:    mean=%.6f  std=%.6f  min=%.6f  max=%.6f\n', ...
    mean(T_out.IrrLI), std(T_out.IrrLI), ...
    min(T_out.IrrLI),  max(T_out.IrrLI));
fprintf('IrrWhole: mean=%.6f  std=%.6f  min=%.6f  max=%.6f\n', ...
    mean(T_out.IrrWhole), std(T_out.IrrWhole), ...
    min(T_out.IrrWhole),  max(T_out.IrrWhole));

%% Sanity check against existing IrrLI file
old_file = fullfile(results_dir, 'INSIDEOUT_features_final_tau1_m40_leftleftflipped.csv');
if exist(old_file, 'file')
    T_old = readtable(old_file);
    merged = innerjoin(T_out(:,{'SubjectID','IrrLI'}), ...
                       T_old(:,{'SubjectID','IrrLI'}), ...
                       'Keys', 'SubjectID');
    r = corr(merged.IrrLI_T_out, merged.IrrLI_T_old);
    fprintf('\nSanity check: r(IrrLI_new, IrrLI_old) = %.6f\n', r);
    if r < 0.999
        warning('IrrLI correlation < 0.999 — check for differences in ROI handling or exclusions');
    end
end
