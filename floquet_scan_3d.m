function lambda_3d = floquet_scan_3d(p, k_vec, m0_vec, l_vec, n_steps, save_path)
% floquet_scan_3d  Compute Floquet growth rate λ over a 3-D (k, m0, l) grid.
%
%   lambda_3d = floquet_scan_3d(p, k_vec, m0_vec, l_vec, n_steps, save_path)
%   computes λ(k, m0, l) for all combinations using flattened parfor over
%   (k, l) pairs.  Results are reshaped to (n_m0, n_k, n_l) and negative
%   values are clamped to 0.
%
%   p           Parameter struct from floquet_params(f).
%   k_vec       Vector of k wavenumbers.
%   m0_vec      Vector of m0 values.
%   l_vec       Vector of l wavenumbers.
%   n_steps     Time steps per period.
%   save_path   Path to save .mat file (required — data always saved before plotting).
%
%   lambda_3d   (n_m0, n_k, n_l) array of growth rates (≥ 0).

n_k  = length(k_vec);
n_m0 = length(m0_vec);
n_l  = length(l_vec);

% --- flattened (k, l) grid: ndgrid so k varies fastest → matches reshape ---
[K_flat, L_flat] = ndgrid(k_vec, l_vec);   % n_k × n_l
n_pairs = numel(K_flat);

fprintf('  Grid: %d k × %d m0 × %d l = %d (k,l)-pairs, %d steps.\n', ...
    n_k, n_m0, n_l, n_pairs, n_steps);

% --- parallel pool ---
use_parallel = ~isempty(gcp('nocreate'));
if ~use_parallel
    try
        parpool('local');
        use_parallel = true;
    catch
        use_parallel = false;
    end
end

tic;
lambda_by_pair = zeros(n_m0, n_pairs);

if use_parallel
    n_workers = gcp('nocreate').NumWorkers;
    fprintf('  Parallel: %d iterations on %d workers.\n', n_pairs, n_workers);
    q = parallel.pool.DataQueue;
    report_interval = max(1, round(n_pairs / 40));
    afterEach(q, @(~) fprintf('.'));

    parfor pi = 1:n_pairs
        k_val = K_flat(pi);
        l_val = L_flat(pi);
        col = zeros(n_m0, 1);
        for mj = 1:n_m0
            col(mj) = floquet_core(k_val, m0_vec(mj), l_val, p, n_steps);
        end
        lambda_by_pair(:, pi) = col;
        if mod(pi, report_interval) == 0
            send(q, pi);
        end
    end
    fprintf('\n');
else
    for pi = 1:n_pairs
        k_val = K_flat(pi);
        l_val = L_flat(pi);
        for mj = 1:n_m0
            lambda_by_pair(mj, pi) = floquet_core(k_val, m0_vec(mj), l_val, p, n_steps);
        end
        if mod(pi, max(1, round(n_pairs / 20))) == 0
            fprintf('  %d/%d done\n', pi, n_pairs);
        end
    end
end

fprintf('  Computation: %.1f s\n', toc);

% --- reshape: (m0, k, l) ---
lambda_3d = reshape(lambda_by_pair, n_m0, n_k, n_l);
lambda_3d(lambda_3d < 0) = 0;

% --- save data (always) ---
save_dir = fileparts(save_path);
if ~exist(save_dir, 'dir'), mkdir(save_dir); end
save(save_path, 'lambda_3d', 'k_vec', 'm0_vec', 'l_vec');
fprintf('  Saved to %s\n', save_path);
end
