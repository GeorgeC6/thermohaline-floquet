function [lambda, is_cached] = cache_manager(data_path, Ri_vec, omega_vec, acc)
% cache_manager  Load, validate, and merge cached 2-D λmax data.
%
%   [lambda, is_cached] = cache_manager(data_path, Ri_vec, omega_vec, acc)
%
%   acc is a struct of accuracy parameters.  Every field must be present
%   in the cache's |cache_acc| struct and match exactly, or the cache is
%   discarded.  Typical fields: f_val, n_k, n_m0, n_l, n_steps, Pr, Rp,
%   tau, k_vec, m0_vec, l_vec.
%
%   On cache hit, lambda is (n_omega × n_Ri) with 0 at uncached positions;
%   is_cached is a logical mask of the same size.  On miss, both are [].
%
%   The cache file must contain:
%       lambda_max   (n_omega × n_Ri)
%       Ri_vec       or  Ri_vals  (1 × n_Ri)
%       omega_vec    (1 × n_omega)
%       cache_acc    (struct — accuracy metadata)
%
%   Save pattern:
%       cache_acc = acc;  %#ok<NASGU>
%       save(data_path, 'lambda_max', 'Ri_vec', 'omega_vec', 'cache_acc');

lambda    = [];
is_cached = [];

if ~exist(data_path, 'file')
    return;
end

fprintf('  Loading cached data: %s\n', data_path);
ld = load(data_path);

% --- validate cache_acc -----------------------------------------------
if ~isfield(ld, 'cache_acc')
    fprintf('  Cache missing accuracy metadata — recomputing.\n');
    return;
end

acc_fields = fieldnames(acc);
for i = 1:length(acc_fields)
    fn = acc_fields{i};
    if ~isfield(ld.cache_acc, fn)
        fprintf('  Cache missing field ''%s'' — recomputing.\n', fn);
        return;
    end
    if ~isequaln(ld.cache_acc.(fn), acc.(fn))
        fprintf('  Cache ''%s'' mismatch — recomputing.\n', fn);
        return;
    end
end
fprintf('  Accuracy params match — merging grids.\n');

% --- merge: map cached → current (ω, Ri) pairs -----------------------
cached_Ri = [];
if isfield(ld, 'Ri_vec'),  cached_Ri = ld.Ri_vec; end
if isfield(ld, 'Ri_vals'), cached_Ri = ld.Ri_vals; end
if isempty(cached_Ri)
    fprintf('  Cache missing Ri vector — recomputing.\n');
    return;
end

n_omega = length(omega_vec);
n_Ri    = length(Ri_vec);

lambda    = zeros(n_omega, n_Ri);
is_cached = false(n_omega, n_Ri);

[~, oi_c2c] = ismember(ld.omega_vec, omega_vec);
[~, ri_c2c] = ismember(cached_Ri,    Ri_vec);

n_cached = 0;
for oi_c = 1:length(ld.omega_vec)
    oi = oi_c2c(oi_c);
    if oi == 0, continue; end
    for ri_c = 1:length(cached_Ri)
        ri = ri_c2c(ri_c);
        if ri == 0, continue; end
        lambda(oi, ri)    = ld.lambda_max(oi_c, ri_c);
        is_cached(oi, ri) = true;
        n_cached = n_cached + 1;
    end
end

n_total = n_omega * n_Ri;
if n_cached == n_total
    fprintf('  All %d points cached — skipping computation.\n', n_total);
else
    fprintf('  %d/%d points cached, %d new to compute.\n', n_cached, n_total, n_total - n_cached);
end
end
