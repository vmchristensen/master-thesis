% INSIDEOUT_compute_alltaus.m
% Computes INSIDEOUT irreversibility features for tau = 1, 2, 3
% Saves one CSV with all features for all subjects
% ------EDIT BEFORE RUNNING!!!!!
clear all; close all;

% PATHS
base_dir_left  = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_dswu';
base_dir_right = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/test_SLIC_results_SHT_left_flipped_dswu';
results_dir    = '/home/marlene/Desktop/victoria/INSIDEOUT_results';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

% PARAMETERS
MAX_LAG = 3;   % tau = 1, 2, 3
N_L     = 100; % left hemisphere ROIs
N_R     = 100; % right hemisphere ROIs (flipped symmetric atlas)

% EXCLUSION LIST
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

% SUBJECT LIST
subj_dirs    = dir(fullfile(base_dir_left, 'sub-EZM*'));
all_subjects = {subj_dirs.name};
subjects     = all_subjects(~ismember(all_subjects, exclude_list));
NSUB         = length(subjects);
fprintf('Found %d subjects, %d excluded, %d remaining\n', ...
        length(all_subjects), length(exclude_list), NSUB);

%---- stop edit

% Preallocate output
% For each tau: IrrLeft, IrrRight, IrrLI, IrrInter, IrrInter_LtoR, IrrInter_RtoL, IrrInter_asym
% Total = NSUB rows x (1 + MAX_LAG*7) columns
ncols = 1 + MAX_LAG*7;
results_table = cell(NSUB, ncols);

for sub = 1:NSUB
    fprintf('Processing subject %d/%d: %s\n', sub, NSUB, subjects{sub});
    subj = subjects{sub};

    left_file  = fullfile(base_dir_left,  subj, 'R_mean_roi_left.mat');
    right_file = fullfile(base_dir_right, subj, 'R_mean_roi_right.mat');

    results_table{sub,1} = subj;

    if ~exist(left_file,'file') || ~exist(right_file,'file')
        warning('Missing file for %s — skipping', subj);
        continue
    end

    dl = load(left_file);
    dr = load(right_file);
    ts_left  = dl.R_mean_roi_left;   % 101 x T
    ts_right = dr.R_mean_roi_right;  % 101 x T

    if size(ts_left,2) ~= size(ts_right,2)
        warning('Timepoint mismatch for %s — skipping', subj);
        continue
    end

    Tm = size(ts_left, 2);

    for Tau = 1:MAX_LAG

        if Tm <= Tau
            warning('Time series too short for tau=%d in %s — skipping', Tau, subj);
            continue
        end

        %% --- WITHIN LEFT HEMISPHERE ---
        ts = ts_left;
        [FC_fwd, ~] = corr(ts(:,1:Tm-Tau)',    ts(:,1+Tau:Tm)');
        [FC_rev, ~] = corr(ts(:,Tm:-1:Tau+1)', ts(:,Tm-Tau:-1:1)');
        Itauf     = -0.5 * log(1 - FC_fwd.*FC_fwd);
        Itaur     = -0.5 * log(1 - FC_rev.*FC_rev);
        FSdiff_LL = (Itauf - Itaur).^2;   % squared for magnitude
        IrrLeft   = mean(FSdiff_LL(:));

        %% --- WITHIN RIGHT HEMISPHERE ---
        ts = ts_right;
        [FC_fwd, ~] = corr(ts(:,1:Tm-Tau)',    ts(:,1+Tau:Tm)');
        [FC_rev, ~] = corr(ts(:,Tm:-1:Tau+1)', ts(:,Tm-Tau:-1:1)');
        Itauf     = -0.5 * log(1 - FC_fwd.*FC_fwd);
        Itaur     = -0.5 * log(1 - FC_rev.*FC_rev);
        FSdiff_RR = (Itauf - Itaur).^2;   % squared for magnitude
        IrrRight  = mean(FSdiff_RR(:));

        %% --- INTERHEMISPHERIC ---
        ts = [ts_left; ts_right];   % 202 x T
        [FC_fwd, ~] = corr(ts(:,1:Tm-Tau)',    ts(:,1+Tau:Tm)');
        [FC_rev, ~] = corr(ts(:,Tm:-1:Tau+1)', ts(:,Tm-Tau:-1:1)');
        Itauf = -0.5 * log(1 - FC_fwd.*FC_fwd);
        Itaur = -0.5 * log(1 - FC_rev.*FC_rev);

        % SIGNED difference for directionality
        FSdiff_signed = Itauf - Itaur;

        % Extract off-diagonal blocks from signed matrix
        FSdiff_LR_signed = FSdiff_signed(1:N_L,     N_L+1:end);  % 101x101
        FSdiff_RL_signed = FSdiff_signed(N_L+1:end, 1:N_L);      % 101x101

        % Directional measures from signed blocks
        IrrInter_LtoR = mean(FSdiff_LR_signed(:));  % positive = L drives R
        IrrInter_RtoL = mean(FSdiff_RL_signed(:));  % positive = R drives L
        IrrInter_asym = IrrInter_LtoR - IrrInter_RtoL;  % net directionality

        % Overall interhemispheric magnitude from squared blocks
        FSdiff_LR_sq = FSdiff_LR_signed.^2;
        FSdiff_RL_sq = FSdiff_RL_signed.^2;
        IrrInter = mean([FSdiff_LR_sq(:); FSdiff_RL_sq(:)]);

        %% --- LATERALIZATION INDEX ---
        IrrLI = (IrrLeft - IrrRight) / (IrrLeft + IrrRight);

        %% --- STORE ---
        col = 1 + (Tau-1)*7;
        results_table{sub, col+1} = IrrLeft;
        results_table{sub, col+2} = IrrRight;
        results_table{sub, col+3} = IrrLI;
        results_table{sub, col+4} = IrrInter;
        results_table{sub, col+5} = IrrInter_LtoR;
        results_table{sub, col+6} = IrrInter_RtoL;
        results_table{sub, col+7} = IrrInter_asym;

    end % tau loop
end % subject loop

%% SAVE
header = {'SubjectID'};
for Tau = 1:MAX_LAG
    header{end+1} = sprintf('IrrLeft_tau%d',       Tau);
    header{end+1} = sprintf('IrrRight_tau%d',      Tau);
    header{end+1} = sprintf('IrrLI_tau%d',         Tau);
    header{end+1} = sprintf('IrrInter_tau%d',      Tau);
    header{end+1} = sprintf('IrrInter_LtoR_tau%d', Tau);
    header{end+1} = sprintf('IrrInter_RtoL_tau%d', Tau);
    header{end+1} = sprintf('IrrInter_asym_tau%d', Tau);
end

T = cell2table(results_table, 'VariableNames', header);
outfile = fullfile(results_dir, 'INSIDEOUT_features_alltaus_m5_leftleftflipped.csv');
writetable(T, outfile);
fprintf('Saved: %s\n', outfile);



% other fmri papers use tau =1 with tr = 2
%% EXTRACT TAU=1 COLUMNS AND SAVE AS FINAL CSV

% Keep only tau=1 columns
tau1_cols = {'SubjectID', ...
             'IrrLeft_tau1', ...
             'IrrRight_tau1', ...
             'IrrLI_tau1', ...
             'IrrInter_tau1', ...
             'IrrInter_LtoR_tau1', ...
             'IrrInter_RtoL_tau1', ...
             'IrrInter_asym_tau1'};

T_final = T(:, tau1_cols);

% Rename to drop tau suffix
T_final.Properties.VariableNames = {'SubjectID', ...
                                     'IrrLeft', ...
                                     'IrrRight', ...
                                     'IrrLI', ...
                                     'IrrInter', ...
                                     'IrrInter_LtoR', ...
                                     'IrrInter_RtoL', ...
                                     'IrrInter_asym'};

% Save
outfile_final = fullfile(results_dir, 'INSIDEOUT_features_final_tau1_m40_leftleftflipped.csv');
writetable(T_final, outfile_final);
fprintf('Saved final features (tau=1): %s\n', outfile_final);
