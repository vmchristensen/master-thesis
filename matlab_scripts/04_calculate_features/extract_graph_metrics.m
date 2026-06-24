% extract_graph_metrics.m
% Computes graph theory metrics from FC matrices for each subject.
%
% Metrics computed per region (on thresholded, weighted graph):
%   1. Weighted degree         — sum of connection strengths
%   2. Clustering coefficient  — local interconnectedness
%   3. Betweenness centrality  — hub-ness / gateway role
%   4. Local efficiency        — resilience of local communication
%
% For each metric: raw value per region + L-R asymmetry per homotopic pair
%
% No external toolboxes required — all implemented from scratch.
%
% Input:  FC_matrix_sub-EZM*_K800.mat  (FC_matrix_z)
%         region_reference.csv
% Output: graph_metrics_K800.mat + graph_metrics_K800.csv

clear; clc;

%% ── USER SETTINGS ────────────────────────────────────────────────────────────
fc_dir       = '/home/marlene/Desktop/victoria/FC_matrices';
ref_path     = '/home/marlene/Desktop/victoria/region_reference.csv';
output_dir   = '/home/marlene/Desktop/victoria/features';

% fc_dir       = '/home/marlene/Desktop/victoria/FC_matrices_notinflated';
% output_dir   = '/home/marlene/Desktop/victoria/features_notinflated';

K          = 800;

% Proportional threshold: keep top X% of connections per subject
% 0.20 = top 20% of edges. Raised from 0.15 to reduce isolated nodes
% (regions with degree=0) which caused ±1 clipping in asymmetry features.
prop_threshold = 0.20;

test_mode = true;
n_test    = 10;

%% ── LOAD REGION REFERENCE ────────────────────────────────────────────────────
fprintf('Loading region reference...\n');
ref = readtable(ref_path);

left_ref  = ref(strcmp(ref.hemisphere, 'Left'),  :);
right_ref = ref(strcmp(ref.hemisphere, 'Right'), :);

left_paired  = left_ref(strcmp(left_ref.has_homotopic_pair,  'Yes'), :);
right_paired = right_ref(strcmp(right_ref.has_homotopic_pair, 'Yes'), :);

left_paired  = sortrows(left_paired,  'homotopic_pair_id');
right_paired = sortrows(right_paired, 'homotopic_pair_id');

assert(isequal(left_paired.homotopic_pair_id, right_paired.homotopic_pair_id), ...
       'Pair IDs do not align — check region_reference.csv');

n_pairs    = height(left_paired);
n_left     = height(left_ref);
n_right    = height(right_ref);
pair_names = left_paired.abbreviation;

all_left_idx   = left_ref.fc_index;
all_right_idx  = right_ref.fc_index;
pair_left_idx  = left_paired.fc_index;
pair_right_idx = right_paired.fc_index;

fprintf('Left: %d | Right: %d | Pairs: %d\n\n', n_left, n_right, n_pairs);

%% ── DISCOVER SUBJECTS ────────────────────────────────────────────────────────
fprintf('Scanning for FC matrix files...\n');
all_files = dir(fullfile(fc_dir, sprintf('FC_matrix_sub-EZM*_K%d.mat', K)));

if isempty(all_files)
    error('No FC matrix files found in %s', fc_dir);
end
fprintf('Found %d subjects.\n', length(all_files));

if test_mode
    all_files = all_files(1:min(n_test, length(all_files)));
    fprintf('TEST MODE: running first %d subjects only.\n\n', length(all_files));
end

mkdir(output_dir);
n_total = length(all_files);

%% ── PRE-ALLOCATE ─────────────────────────────────────────────────────────────
% Raw metrics: one value per region per subject
deg_L   = NaN(n_total, n_left);   deg_R   = NaN(n_total, n_right);
clust_L = NaN(n_total, n_left);   clust_R = NaN(n_total, n_right);
btwn_L  = NaN(n_total, n_left);   btwn_R  = NaN(n_total, n_right);
leff_L  = NaN(n_total, n_left);   leff_R  = NaN(n_total, n_right);

% Asymmetry: one value per homotopic pair per subject
asym_deg   = NaN(n_total, n_pairs);
asym_clust = NaN(n_total, n_pairs);
asym_btwn  = NaN(n_total, n_pairs);
asym_leff  = NaN(n_total, n_pairs);

subject_ids = cell(n_total, 1);
failed      = {};

%% ── MAIN LOOP ────────────────────────────────────────────────────────────────
for s = 1:n_total
    filename = all_files(s).name;
    sub_id   = regexp(filename, 'sub-EZM\d+', 'match', 'once');
    subject_ids{s} = sub_id;

    fprintf('Subject %d / %d : %s\n', s, n_total, sub_id);

    try
        d  = load(fullfile(fc_dir, filename));
        FC = d.FC_matrix_z;
        n  = size(FC, 1);

        %% THRESHOLD: keep top prop_threshold% of connections (proportional)
        % Work on absolute values (negative FC still represents a connection)
        FC_abs = abs(FC);
        FC_abs(logical(eye(n))) = 0;   % zero diagonal

        % Find threshold value
        upper_tri_vals = FC_abs(triu(true(n), 1));
        thresh_val     = quantile(upper_tri_vals, 1 - prop_threshold);

        % Apply: zero out weak connections, keep sign and magnitude
        W = FC;
        W(FC_abs < thresh_val) = 0;
        W(logical(eye(n))) = 0;

        % Normalise weights to [0, 1] for clustering/efficiency calculations
        W_pos = max(W, 0);   % use positive weights only for these metrics
        w_max = max(W_pos(:));
        if w_max > 0
            W_norm = W_pos / w_max;
        else
            W_norm = W_pos;
        end

        %% COMPUTE METRICS
        deg_vec   = compute_degree(W_norm, n);
        clust_vec = compute_clustering(W_norm, n);
        btwn_vec  = compute_betweenness(W_norm, n);
        leff_vec  = compute_local_efficiency(W_norm, n);

        %% STORE PER-HEMISPHERE VALUES
        valid_L = all_left_idx(all_left_idx <= n);
        valid_R = all_right_idx(all_right_idx <= n);

        for r = 1:n_left
            i = all_left_idx(r);
            if i > n, continue; end
            deg_L(s,r)   = deg_vec(i);
            clust_L(s,r) = clust_vec(i);
            btwn_L(s,r)  = btwn_vec(i);
            leff_L(s,r)  = leff_vec(i);
        end

        for r = 1:n_right
            i = all_right_idx(r);
            if i > n, continue; end
            deg_R(s,r)   = deg_vec(i);
            clust_R(s,r) = clust_vec(i);
            btwn_R(s,r)  = btwn_vec(i);
            leff_R(s,r)  = leff_vec(i);
        end

        %% ASYMMETRY PER HOMOTOPIC PAIR: (L - R) / (L + R)
        for p = 1:n_pairs
            iL = pair_left_idx(p);
            iR = pair_right_idx(p);
            if iL > n || iR > n, continue; end

            asym_deg(s,p)   = safe_asym(deg_vec(iL),   deg_vec(iR));
            asym_clust(s,p) = safe_asym(clust_vec(iL), clust_vec(iR));
            asym_btwn(s,p)  = safe_asym(btwn_vec(iL),  btwn_vec(iR));
            asym_leff(s,p)  = safe_asym(leff_vec(iL),  leff_vec(iR));
        end

    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        failed{end+1} = sub_id;
    end
end

%% ── SAVE .mat ────────────────────────────────────────────────────────────────
output_mat = fullfile(output_dir, sprintf('graph_metrics_K%d.mat', K));

save(output_mat, ...
     'deg_L',   'deg_R',   'clust_L', 'clust_R', ...
     'btwn_L',  'btwn_R',  'leff_L',  'leff_R',  ...
     'asym_deg', 'asym_clust', 'asym_btwn', 'asym_leff', ...
     'subject_ids', 'pair_names', 'left_ref', 'right_ref', ...
     'left_paired', 'right_paired', 'prop_threshold');

fprintf('\nSaved: %s\n', output_mat);

%% ── SAVE .csv ────────────────────────────────────────────────────────────────
% Build column names — append fc_index for any duplicated abbreviation
col_degL   = make_unique_names('deg_L_',    left_ref.abbreviation,  left_ref.fc_index);
col_degR   = make_unique_names('deg_R_',    right_ref.abbreviation, right_ref.fc_index);
col_clustL = make_unique_names('clust_L_',  left_ref.abbreviation,  left_ref.fc_index);
col_clustR = make_unique_names('clust_R_',  right_ref.abbreviation, right_ref.fc_index);
col_btwnL  = make_unique_names('btwn_L_',   left_ref.abbreviation,  left_ref.fc_index);
col_btwnR  = make_unique_names('btwn_R_',   right_ref.abbreviation, right_ref.fc_index);
col_leffL  = make_unique_names('leff_L_',   left_ref.abbreviation,  left_ref.fc_index);
col_leffR  = make_unique_names('leff_R_',   right_ref.abbreviation, right_ref.fc_index);
col_adeg   = strcat('asym_deg_',   pair_names);
col_aclust = strcat('asym_clust_', pair_names);
col_abtwn  = strcat('asym_btwn_',  pair_names);
col_aleff  = strcat('asym_leff_',  pair_names);

data_mat = [deg_L, deg_R, clust_L, clust_R, ...
            btwn_L, btwn_R, leff_L, leff_R, ...
            asym_deg, asym_clust, asym_btwn, asym_leff];

all_cols = [col_degL; col_degR; col_clustL; col_clustR; ...
            col_btwnL; col_btwnR; col_leffL; col_leffR; ...
            col_adeg; col_aclust; col_abtwn; col_aleff];

T_out = array2table(data_mat, 'VariableNames', all_cols);
T_out = addvars(T_out, subject_ids, 'Before', 1, 'NewVariableNames', 'subject_id');

output_csv = fullfile(output_dir, sprintf('graph_metrics_K%d.csv', K));
writetable(T_out, output_csv);
fprintf('Saved: %s\n', output_csv);

%% ── SUMMARY ──────────────────────────────────────────────────────────────────
fprintf('\n========================================\n');
fprintf('DONE\n');
fprintf('  Subjects   : %d processed, %d failed\n', n_total - length(failed), length(failed));
fprintf('  Features   : %d per subject\n', size(data_mat, 2));
fprintf('    Raw L/R  : %d  |  Asymmetry: %d\n', ...
        size([deg_L,deg_R,clust_L,clust_R,btwn_L,btwn_R,leff_L,leff_R], 2), ...
        size([asym_deg,asym_clust,asym_btwn,asym_leff], 2));
fprintf('  Threshold  : top %.0f%% of connections\n', prop_threshold*100);
fprintf('========================================\n');

if ~isempty(failed)
    fprintf('\nFailed subjects:\n');
    for i = 1:length(failed), fprintf('  %s\n', failed{i}); end
end

%% ════════════════════════════════════════════════════════════════════════════
%%  LOCAL FUNCTIONS
%% ════════════════════════════════════════════════════════════════════════════

function d = compute_degree(W, n)
% Weighted degree: sum of weights for each node
    d = sum(W, 2);   % column sum = in-degree; for symmetric W same as row sum
    d = d / (n - 1); % normalise by max possible connections
end

function C = compute_clustering(W, n)
% Weighted clustering coefficient (Onnela et al. 2005)
% C_i = (1/k_i*(k_i-1)) * sum_{j,h} (w_ij * w_ih * w_jh)^(1/3)
% where k_i = number of neighbours (binary degree)
    C = zeros(n, 1);
    for i = 1:n
        neighbours = find(W(i,:) > 0);
        k = length(neighbours);
        if k < 2
            C(i) = 0;
            continue
        end
        % Submatrix of weights among node i and its neighbours
        w_ij = W(i, neighbours)';   % k x 1
        w_ih = W(i, neighbours)';   % k x 1 (same)
        W_sub = W(neighbours, neighbours);  % k x k
        % Geometric mean of triplet weights
        triplet_sum = 0;
        for j = 1:k
            for h = 1:k
                if j ~= h
                    triplet_sum = triplet_sum + ...
                        (w_ij(j) * w_ih(h) * W_sub(j,h))^(1/3);
                end
            end
        end
        C(i) = triplet_sum / (k * (k - 1));
    end
end

function B = compute_betweenness(W, n)
% Betweenness centrality (Brandes algorithm on weighted graph)
% Uses inverse of weights as distances (stronger connection = shorter path)
% Normalised by (n-1)*(n-2)/2

    % Convert weights to distances: d = 1/w (avoid division by zero)
    D = zeros(n);
    nz = W > 0;
    D(nz) = 1 ./ W(nz);
    D(~nz) = Inf;
    D(logical(eye(n))) = 0;

    B = zeros(n, 1);

    for s_node = 1:n
        % Dijkstra from s_node
        dist    = Inf(n, 1);
        dist(s_node) = 0;
        sigma   = zeros(n, 1);   % number of shortest paths
        sigma(s_node) = 1;
        visited = false(n, 1);
        delta   = zeros(n, 1);   % dependency
        pred    = cell(n, 1);    % predecessors

        % Priority queue via sorted list (simple for small n)
        queue = s_node;
        dist_q = 0;

        while ~isempty(queue)
            % Pick node with minimum distance
            [~, min_idx] = min(dist_q);
            u = queue(min_idx);
            queue(min_idx)   = [];
            dist_q(min_idx)  = [];

            if visited(u), continue; end
            visited(u) = true;

            % Relax neighbours
            neighbours = find(D(u,:) < Inf & D(u,:) > 0);
            for v = neighbours
                alt = dist(u) + D(u, v);
                if alt < dist(v)
                    dist(v)  = alt;
                    sigma(v) = sigma(u);
                    pred{v}  = u;
                    queue   = [queue,   v];
                    dist_q  = [dist_q, alt];
                elseif abs(alt - dist(v)) < 1e-10
                    sigma(v) = sigma(v) + sigma(u);
                    pred{v}  = [pred{v}, u];
                end
            end
        end

        % Accumulate dependencies (back-propagation)
        % Process nodes in reverse order of distance
        [~, order] = sort(dist, 'descend');
        for w_node = order'
            if isinf(dist(w_node)), continue; end
            for u = pred{w_node}
                delta(u) = delta(u) + (sigma(u)/sigma(w_node)) * (1 + delta(w_node));
            end
            if w_node ~= s_node
                B(w_node) = B(w_node) + delta(w_node);
            end
        end
    end

    % Normalise
    B = B / ((n-1)*(n-2));
end

function E = compute_local_efficiency(W, n)
% Local efficiency: average efficiency of the subgraph of node i's neighbours
% E_loc_i = E_glob(G_i) where G_i is the subgraph of neighbours of i
    E = zeros(n, 1);
    for i = 1:n
        neighbours = find(W(i,:) > 0);
        k = length(neighbours);
        if k < 2
            E(i) = 0;
            continue
        end
        % Subgraph of neighbours
        W_sub = W(neighbours, neighbours);
        n_sub = length(neighbours);
        % Global efficiency of subgraph = mean of inverse shortest paths
        D_sub = zeros(n_sub);
        nz = W_sub > 0;
        D_sub(nz)  = 1 ./ W_sub(nz);
        D_sub(~nz) = Inf;
        D_sub(logical(eye(n_sub))) = 0;
        % Floyd-Warshall for shortest paths (fast for small subgraphs)
        for via = 1:n_sub
            for u = 1:n_sub
                for v = 1:n_sub
                    if D_sub(u,via) + D_sub(via,v) < D_sub(u,v)
                        D_sub(u,v) = D_sub(u,via) + D_sub(via,v);
                    end
                end
            end
        end
        % Efficiency = mean of 1/d for all pairs (exclude diagonal)
        inv_d = 1 ./ D_sub;
        inv_d(logical(eye(n_sub))) = 0;
        inv_d(isinf(inv_d)) = 0;
        E(i) = sum(inv_d(:)) / (n_sub * (n_sub - 1));
    end
end

function names = make_unique_names(prefix, abbreviations, fc_indices)
% Appends _<fc_index> only for abbreviations that appear more than once.
    [~, ~, ic] = unique(abbreviations);
    counts     = accumarray(ic, 1);
    dup_mask   = counts(ic) > 1;
    names      = strcat(prefix, abbreviations);
    for i = 1:length(abbreviations)
        if dup_mask(i)
            names{i} = sprintf('%s%s_%d', prefix, abbreviations{i}, fc_indices(i));
        end
    end
end

function a = safe_asym(L, R)
% Asymmetry index: (L-R)/(|L|+|R|), bounded to [-1, +1].
% Using |L|+|R| prevents sign cancellation making denominator near zero.
    denom = abs(L) + abs(R);
    if denom > 1e-10
        a = (L - R) / denom;
    else
        a = 0;
    end
end
