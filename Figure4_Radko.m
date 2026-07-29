% Figure4_Radko.m
% Reproduce Radko (2019) Figure 4:
%   log10 of the maximal Floquet growth rate λmax as a function of (Ri, ω).
%
% For each (Ri, ω) pair, the 2-D (f=0, l=0) Floquet system is solved over a
% (k, m0) wavenumber grid and the maximum growth rate is recorded.
%
% Reference: Radko (2019), J. Phys. Oceanogr., 49, 2379-2392, Fig. 4.

close all; clear; clc;

fprintf('===== Radko (2019) Fig. 4: log10(λmax) vs (Ri, ω) =====\n');

%% -------------------- user settings --------------------
f_val    = 0;         % Coriolis parameter
n_Ri     = 51;        % Richardson number grid points
n_omega  = 51;        % frequency grid points
n_k      = 51;        % k wavenumber grid for max search
n_m0     = 51;        % m0 grid for max search
n_l      = 21;        % l wavenumber grid
n_steps  = 1000;      % time steps per period

Ri_vec    = linspace(0.1, 10, n_Ri);
omega_vec = linspace(0.1, 10, n_omega);
k_vec     = linspace(-0.5, 0.5, n_k);
m0_vec    = linspace(0, 1.5, n_m0);
l_vec     = linspace(-0.03, 0.03, n_l);

if f_val == 0
    l_vec = 0;  % f=0 → λmax always at l=0
    n_l   = 1;
end

% Fixed parameters
Pr  = 10;
Rp  = 2;
tau = 0.01;

%% -------------------- output directory --------------------
out_dir = fullfile('Figures', 'Figure4', sprintf('f_%.3f', f_val));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
data_path = fullfile(out_dir, 'Figure4_data.mat');

%% -------------------- compute --------------------
% Accuracy struct for cache validation
acc = struct('f_val', f_val, 'n_k', n_k, 'n_m0', n_m0, 'n_l', n_l, ...
    'n_steps', n_steps, 'Pr', Pr, 'Rp', Rp, 'tau', tau, ...
    'k_vec', k_vec, 'm0_vec', m0_vec, 'l_vec', l_vec);

[lambda_max, is_cached] = cache_manager(data_path, Ri_vec, omega_vec, acc);

if isempty(lambda_max)
    lambda_max = zeros(n_omega, n_Ri);
    is_cached = false(n_omega, n_Ri);
end

needs_idx = find(~is_cached);
if ~isempty(needs_idx)
    n_needed = length(needs_idx);
    fprintf('Computing %d (Ri, ω) pairs.  (k,m0,l): %d × %d × %d.\n', ...
        n_needed, n_k, n_m0, n_l);

    new_vals = zeros(n_needed, 1);

    % Parallel pool
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

    if use_parallel
        n_workers = gcp('nocreate').NumWorkers;
        fprintf('  Parallel: %d workers.\n', n_workers);
        interval = max(1, round(n_needed / 20));
        state    = containers.Map({'done', 'n', 't0'}, {0, n_needed, tic});
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) report_progress(state, interval));

        parfor pi = 1:n_needed
            idx = needs_idx(pi);
            [oi, ri] = ind2sub([n_omega, n_Ri], idx);
            p = floquet_params(f_val, 'Ri', Ri_vec(ri), 'omega', omega_vec(oi));

            max_lam = 0;
            for ki = 1:n_k
                k_val = k_vec(ki);
                for mi = 1:n_m0
                    for li = 1:n_l
                        lam = floquet_core(k_val, m0_vec(mi), l_vec(li), p, n_steps);
                        if lam > max_lam, max_lam = lam; end
                    end
                end
            end
            new_vals(pi) = max_lam;

            if mod(pi, interval) == 0
                send(q, pi);
            end
        end
        fprintf('\n');
    else
        fprintf('  Serial mode.\n');
        for pi = 1:n_needed
            idx = needs_idx(pi);
            [oi, ri] = ind2sub([n_omega, n_Ri], idx);
            p = floquet_params(f_val, 'Ri', Ri_vec(ri), 'omega', omega_vec(oi));

            max_lam = 0;
            for ki = 1:n_k
                k_val = k_vec(ki);
                for mi = 1:n_m0
                    for li = 1:n_l
                        lam = floquet_core(k_val, m0_vec(mi), l_vec(li), p, n_steps);
                        if lam > max_lam, max_lam = lam; end
                    end
                end
            end
            new_vals(pi) = max_lam;

            if mod(pi, max(1, round(n_needed / 20))) == 0
                fprintf('  %d/%d done (%.0f s)\n', pi, n_needed, toc);
            end
        end
    end

    fprintf('  Computation: %.1f s\n', toc);

    lambda_max(needs_idx) = new_vals;

    % Save merged result with accuracy metadata
    cache_acc = acc; %#ok<NASGU>
    save(data_path, 'lambda_max', 'Ri_vec', 'omega_vec', 'cache_acc');
    fprintf('  Data saved to %s\n', data_path);
end

%% -------------------- summary statistics --------------------
fprintf('\n--- Summary ---\n');
fprintf('λmax range:  %.6f  to  %.6f\n', min(lambda_max(:)), max(lambda_max(:)));
[~, idx] = max(lambda_max(:));
[oi, ri] = ind2sub(size(lambda_max), idx);
fprintf('Global max λ = %.6f  at Ri = %.3f, ω = %.3f\n', ...
    lambda_max(oi, ri), Ri_vec(ri), omega_vec(oi));

%% -------------------- plot: log10 heatmap --------------------
data_log = log10(lambda_max);
data_log(lambda_max <= 0) = -6;  % stable → floor
data_log(data_log < -6) = -6;
clims = [-6, max(-3, max(data_log(:)))];

figure('Units', 'inches', 'Position', [1 1 8 6.5]);
imagesc(Ri_vec, omega_vec, data_log);
axis xy;

% --- data cursor: UserData + CreateFcn survives .fig save/load ---
set(gcf, 'UserData', struct('type', 'heatmap', 'Ri_vec', Ri_vec, 'omega_vec', omega_vec, 'data_log', data_log));
set(gcf, 'CreateFcn', 'setup_datatip(gcf)');
setup_datatip(gcf);  % also set up immediately for the live figure

set(gca, 'FontSize', 14);
xlabel('$Ri$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 22);
colormap('jet');
clim(clims);
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10}(\lambda_{\max})$';
cb.Label.FontSize = 18;
title(sprintf('Maximal Floquet growth rate $\\log_{10}(\\lambda_{\\max})$ vs $(Ri, \\omega)$, $f=%.1f$', f_val), ...
    'Interpreter', 'latex', 'FontSize', 16);

saveas(gcf, fullfile(out_dir, 'Figure4_Radko.png'));
savefig(gcf, fullfile(out_dir, 'Figure4_Radko.fig'));
fprintf('Figure saved to %s/\n', out_dir);
fprintf('===== Done =====\n');
