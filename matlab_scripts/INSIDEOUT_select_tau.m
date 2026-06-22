% INSIDEOUT_select_tau.m
% Loads all-tau CSV, selects optimal tau from autocorrelation decay,
% saves final CSV with one set of features per subject
% ------EDIT BEFORE RUNNING!!!!!

base_dir_left = '/home/marlene/Desktop/victoria/slic_individual_light_IMPORTANT/SLIC-individual-light-master/m5_test_SLIC_results_SHT_left_dswu';
results_dir   = '/home/marlene/Desktop/victoria/INSIDEOUT_results';

MAX_LAG   = 3;
AC_THRESH = 0.1;  % autocorrelation decay threshold

% Same exclusion list as Script 1
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

%---- stop edit

%% STEP 1: COMPUTE AUTOCORRELATION ACROSS SUBJECTS
subj_dirs    = dir(fullfile(base_dir_left, 'sub-EZM*'));
all_subjects = {subj_dirs.name};
subjects     = all_subjects(~ismember(all_subjects, exclude_list));
NSUB         = length(subjects);

fprintf('Computing autocorrelation for tau selection...\n');
all_autocorr  = nan(NSUB, MAX_LAG);

for sub = 1:NSUB
    left_file = fullfile(base_dir_left, subjects{sub}, 'R_mean_roi_left.mat');
    if ~exist(left_file,'file'), continue; end

    dl = load(left_file);
    ts = dl.R_mean_roi_left;  % 101 x T
    Tm = size(ts, 2);

    for lag = 1:MAX_LAG
        if Tm > lag
            ac_per_roi = diag(corr(ts(:,1:Tm-lag)', ts(:,1+lag:Tm)'));
            all_autocorr(sub, lag) = mean(ac_per_roi);
        end
    end
end

% Average across subjects (ignore NaNs from missing files)
mean_ac = nanmean(all_autocorr, 1);
fprintf('Mean autocorrelation at tau=1: %.3f  tau=2: %.3f  tau=3: %.3f\n', ...
        mean_ac(1), mean_ac(2), mean_ac(3));

% Select tau
Tauwinner = find(mean_ac < AC_THRESH, 1, 'first');
if isempty(Tauwinner)
    Tauwinner = 1;
    fprintf('Autocorrelation did not drop below %.2f — defaulting to tau=1\n', AC_THRESH);
else
    fprintf('Selected tau = %d\n', Tauwinner);
end

%% STEP 2: EXTRACT TAUWINNER COLUMNS FROM ALL-TAU CSV
infile  = fullfile(results_dir, 'INSIDEOUT_features_alltaus_m5_leftleftflipped.csv');
T_all   = readtable(infile);

% Column names for the winning tau
cols_to_keep = {
    'SubjectID',
    sprintf('IrrLeft_tau%d',       Tauwinner),
    sprintf('IrrRight_tau%d',      Tauwinner),
    sprintf('IrrLI_tau%d',         Tauwinner),
    sprintf('IrrInter_tau%d',      Tauwinner),
    sprintf('IrrInter_LtoR_tau%d', Tauwinner),
    sprintf('IrrInter_RtoL_tau%d', Tauwinner),
    sprintf('IrrInter_asym_tau%d', Tauwinner)
};

T_final = T_all(:, cols_to_keep);

% Rename columns to drop the tau suffix for cleaner R import
T_final.Properties.VariableNames = {
    'SubjectID',
    'IrrLeft',
    'IrrRight',
    'IrrLI',
    'IrrInter',
    'IrrInter_LtoR',
    'IrrInter_RtoL',
    'IrrInter_asym'
};

% Add tau used as reference column
T_final.Tau_used = repmat(Tauwinner, height(T_final), 1);

%% SAVE
outfile = fullfile(results_dir, 'INSIDEOUT_features_final_m5_leftleftflipped.csv');
writetable(T_final, outfile);
fprintf('Saved final features (tau=%d): %s\n', Tauwinner, outfile);