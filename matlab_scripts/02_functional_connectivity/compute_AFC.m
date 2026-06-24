% =============================================================================
% Asymmetric Functional Connectivity (AFC) — Raemaekers et al. (2018)
% 258 subjects (EZM2 + EZM3)
% =============================================================================

clear; clc;

% ── EDIT THESE TWO PATHS ─────────────────────────────────────────────────────
data_dir = '/home/marlene/Desktop/victoria/FC_matrices'; 
output_dir = '/home/marlene/Desktop/victoria/AFC_output';
% ─────────────────────────────────────────────────────────────────────────────

% ── EXACT 258 SUBJECT LIST ───────────────────────────────────────────────────
subject_ids = { ...
    'EZM30009','EZM20010','EZM20011','EZM20015','EZM20016','EZM20017', ...
    'EZM20018','EZM20019','EZM20020','EZM20021','EZM20022','EZM20023', ...
    'EZM20025','EZM20026','EZM20027','EZM20029','EZM20030','EZM20031', ...
    'EZM20032','EZM20034','EZM20035','EZM20036','EZM20037','EZM20038', ...
    'EZM20039','EZM20040','EZM20041','EZM20042','EZM20043','EZM20044', ...
    'EZM20047','EZM20048','EZM20049','EZM20051','EZM20053', ...
    'EZM20054','EZM20055','EZM20056','EZM20058','EZM20060','EZM20061', ...
    'EZM20062','EZM20063','EZM20064','EZM20065','EZM20066','EZM20067', ...
    'EZM20068','EZM20069','EZM20070','EZM20071','EZM20072','EZM20073', ...
    'EZM20074','EZM20075','EZM20076','EZM20077','EZM20078','EZM20079', ...
    'EZM20080','EZM20081','EZM20082','EZM20083','EZM20084','EZM20087', ...
    'EZM20089','EZM20090','EZM20092','EZM20093','EZM20095','EZM20096', ...
    'EZM20097','EZM20098','EZM20100','EZM20101','EZM20102','EZM20103', ...
    'EZM20105','EZM20106','EZM20107','EZM20108','EZM20109','EZM20110', ...
    'EZM20111','EZM20113','EZM20116','EZM20117','EZM20118','EZM20119', ...
    'EZM20121','EZM20122','EZM20124','EZM20125','EZM20126', ...
    'EZM20127','EZM20128','EZM20129','EZM20130','EZM20131','EZM20132', ...
    'EZM20133','EZM20134','EZM20139','EZM20143','EZM20144','EZM20145', ...
    'EZM20151','EZM30152','EZM20154','EZM20156','EZM20157','EZM20158', ...
    'EZM20159','EZM20160','EZM20161','EZM20162','EZM20163','EZM20164', ...
    'EZM20165','EZM20169','EZM20170','EZM20171','EZM20172','EZM20173', ...
    'EZM20174','EZM20175','EZM30176','EZM20177','EZM20180','EZM20181', ...
    'EZM20182','EZM20183','EZM20185','EZM20187','EZM20188','EZM20189', ...
    'EZM20191','EZM20192','EZM20194','EZM20198','EZM20199', ...
    'EZM20200','EZM20201','EZM20202','EZM20204','EZM30002', ...
    'EZM30024','EZM30028','EZM30033','EZM30045','EZM30046','EZM30050', ...
    'EZM30059','EZM30085','EZM30094','EZM30099','EZM30104','EZM30120', ...
    'EZM30123','EZM30136','EZM30137','EZM30141','EZM30146','EZM30148', ...
    'EZM30149','EZM30153','EZM30155','EZM30168','EZM30190','EZM30196', ...
    'EZM30203','EZM30206','EZM30207','EZM30208','EZM30209','EZM30211', ...
    'EZM30212','EZM30213','EZM30214','EZM30215','EZM30216','EZM30217', ...
    'EZM30218','EZM30219','EZM30220','EZM30221','EZM30222','EZM30223', ...
    'EZM30224','EZM30225','EZM30226','EZM30227','EZM30228','EZM30229', ...
    'EZM30230','EZM30231','EZM30232','EZM30233','EZM30234','EZM30235', ...
    'EZM30236','EZM30237','EZM30238','EZM30239','EZM30241','EZM30242', ...
    'EZM30243','EZM30244','EZM30245','EZM30246','EZM30247','EZM30248', ...
    'EZM30249','EZM30251','EZM30252','EZM30253','EZM30254','EZM30255', ...
    'EZM30257','EZM30258','EZM30259','EZM30260','EZM30261','EZM30266', ...
    'EZM30267','EZM30268','EZM30269','EZM30271','EZM30272','EZM30273', ...
    'EZM30274','EZM30275','EZM30276','EZM30277','EZM30278','EZM30281', ...
    'EZM30283','EZM30284','EZM30285','EZM30286','EZM30287','EZM30288', ...
    'EZM30289','EZM30291','EZM30292','EZM30293','EZM30294','EZM30295', ...
    'EZM30296','EZM30297','EZM30298','EZM30299','EZM30301','EZM30302', ...
    'EZM30304' ...
};

n_subj = numel(subject_ids);
fprintf('Processing %d subjects\n', n_subj);

% ── PREALLOCATE ──────────────────────────────────────────────────────────────
n_roi = 100;
n_tri = n_roi * (n_roi - 1) / 2;          % 4950 upper-triangle pairs

[uidx_r, uidx_c] = find(triu(ones(n_roi), 1));   % [4950x1] index arrays

% afc_vectors: [259 x 9900]
%   cols 1:4950      = within-hemisphere  AFC (upper tri of LL - RR)
%   cols 4951:9900   = between-hemisphere AFC (upper tri of LR - LR')
afc_vectors = NaN(n_subj, 2 * n_tri);

failed = {};

% ── MAIN LOOP ────────────────────────────────────────────────────────────────
for s = 1:n_subj
    subj      = subject_ids{s};
    subj_path = fullfile(data_dir, ['sub-' subj]);   % adds 'sub-' prefix

    if ~exist(subj_path, 'dir')
        warning('Folder not found: sub-%s — skipping', subj);
        failed{end+1} = subj;
        continue
    end

    try
        d_LL = load(fullfile(subj_path, 'FC_left.mat'));
        d_RR = load(fullfile(subj_path, 'FC_right.mat'));
        d_LR = load(fullfile(subj_path, 'FC_matrix.mat'));

        LL = d_LL.FC_left_z;     % [100x100] left  x left
        RR = d_RR.FC_right_z;    % [100x100] right x right
        LR = d_LR.FC_matrix_z;   % [100x100] left  x right

        % Within-hemisphere AFC: how much stronger is L-L vs R-R?
        AFC_within = LL - RR;
        AFC_within(1:n_roi+1:end) = NaN;       % clear diagonal

        % Between-hemisphere AFC: LR[i,j] vs its mirror LR[j,i]
        % LR'[i,j] = LR[j,i] = connectivity of left-j with right-i
        AFC_between = LR - LR';
        AFC_between(1:n_roi+1:end) = NaN;      % clear homotopic diagonal

        % Flatten upper triangles using precomputed indices
        lin_idx = sub2ind([n_roi n_roi], uidx_r, uidx_c);
        w = AFC_within(lin_idx);
        b = AFC_between(lin_idx);

        afc_vectors(s, :) = [w; b]';

    catch err
        warning('Error loading sub-%s: %s', subj, err.message);
        failed{end+1} = subj;
    end
end

fprintf('\nCompleted: %d succeeded, %d failed\n', ...
        n_subj - numel(failed), numel(failed));
if ~isempty(failed)
    disp('Failed subjects:'); disp(failed');
end

% ── GROUP-MEAN AFC ───────────────────────────────────────────────────────────
group_mean_vec = mean(afc_vectors, 1, 'omitnan');   % [1 x 9900]

% Rebuild as [100x100] matrices for plotting and inspection
lin_idx     = sub2ind([n_roi n_roi], uidx_r, uidx_c);
lin_idx_T   = sub2ind([n_roi n_roi], uidx_c, uidx_r);  % transposed positions

group_mean_within  = zeros(n_roi);
group_mean_between = zeros(n_roi);

group_mean_within(lin_idx)   =  group_mean_vec(1:n_tri);
group_mean_within(lin_idx_T) =  group_mean_vec(1:n_tri);         % symmetric

group_mean_between(lin_idx)   =  group_mean_vec(n_tri+1:end);
group_mean_between(lin_idx_T) = -group_mean_vec(n_tri+1:end);    % antisymmetric

group_mean_within(1:n_roi+1:end)  = NaN;
group_mean_between(1:n_roi+1:end) = NaN;

% ── AFC-SCORE PER SUBJECT ────────────────────────────────────────────────────
% Regress each subject's AFC vector on the group mean vector
% beta (AFC-score) = how strongly this subject expresses the group pattern
% High = typical left-lateralization; low/negative = atypical

ok_cols    = ~isnan(group_mean_vec);
AFC_scores = NaN(n_subj, 1);

for s = 1:n_subj
    indiv = afc_vectors(s, :);
    ok    = ok_cols & ~isnan(indiv);
    if sum(ok) < 100, continue; end
    X             = [ones(sum(ok), 1), group_mean_vec(ok)'];
    b             = X \ indiv(ok)';
    AFC_scores(s) = b(2);
end

fprintf('\nAFC-score summary (%d subjects):\n', sum(~isnan(AFC_scores)));
fprintf('  Mean=%.4f  SD=%.4f  Min=%.4f  Max=%.4f\n', ...
        mean(AFC_scores,'omitnan'), std(AFC_scores,'omitnan'), ...
        min(AFC_scores,[],'omitnan'), max(AFC_scores,[],'omitnan'));

% ── SIGNIFICANCE TESTING ─────────────────────────────────────────────────────
% One-sample t-test per connection, Bonferroni corrected over 9900 tests
alpha_bonf = 0.05 / (2 * n_tri);

p_within  = NaN(n_tri, 1);
p_between = NaN(n_tri, 1);

for k = 1:n_tri
    w_ok = afc_vectors(~isnan(afc_vectors(:,k)),          k);
    b_ok = afc_vectors(~isnan(afc_vectors(:,n_tri+k)), n_tri+k);
    if numel(w_ok) > 5,  [~, p_within(k)]  = ttest(w_ok); end
    if numel(b_ok) > 5,  [~, p_between(k)] = ttest(b_ok); end
end

sig_within  = p_within  < alpha_bonf;
sig_between = p_between < alpha_bonf;

fprintf('\nSignificant within-hemisphere AFC:  %d / %d\n', ...
        sum(sig_within,'omitnan'),  n_tri);
fprintf('Significant between-hemisphere AFC: %d / %d\n', ...
        sum(sig_between,'omitnan'), n_tri);

% ── SAVE ─────────────────────────────────────────────────────────────────────
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

save(fullfile(output_dir, 'AFC_results_corrected.mat'), ...
     'afc_vectors', ...          % [259 x 9900] — use as features in R
     'group_mean_vec', ...       % [1 x 9900]
     'group_mean_within', ...    % [100x100] for plotting
     'group_mean_between', ...   % [100x100] for plotting
     'AFC_scores', ...           % [259 x 1]
     'p_within',  'p_between', ...
     'sig_within', 'sig_between', ...
     'subject_ids', 'uidx_r', 'uidx_c', ...
     '-v7.3');

% CSV for loading into R
fid = fopen(fullfile(output_dir, 'AFC_scores_corrected.csv'), 'w');
fprintf(fid, 'subject_id,AFC_score\n');
for s = 1:n_subj
    fprintf(fid, '%s,%.6f\n', subject_ids{s}, AFC_scores(s));
end
fclose(fid);

fprintf('\nSaved to %s\n', output_dir);

% ── PLOT ─────────────────────────────────────────────────────────────────────
figure('Name','Group-mean AFC','Position',[100 100 1000 420]);

subplot(1,2,1);
imagesc(group_mean_within, [-0.15 0.15]);
colormap(bluewhitered(256)); colorbar; axis square;
title('Within-hemisphere AFC  (LL - RR)');
xlabel('ROI index'); ylabel('ROI index');

subplot(1,2,2);
imagesc(group_mean_between, [-0.15 0.15]);
colorbar; axis square;
title('Between-hemisphere AFC  (LR - LR^T)');
xlabel('ROI index'); ylabel('ROI index');

% ── HELPER FUNCTION ──────────────────────────────────────────────────────────
function cmap = bluewhitered(n)
    if nargin < 1, n = 256; end
    half = floor(n/2);
    top  = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
    bot  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
    cmap = [top; bot];
end
