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
n_Ri     = 51;        % Richardson number grid points
n_omega  = 51;        % frequency grid points
n_k      = 51;        % k wavenumber grid for max search
n_m0     = 51;        % m0 grid for max search
n_steps  = 1000;       % time steps per period (fewer for speed)

Ri_vec    = linspace(0.5, 10, n_Ri);
omega_vec = linspace(0.05, 1, n_omega);
k_vec     = linspace(-0.5, 0.5, n_k);
m0_vec    = linspace(0, 1.5, n_m0);

% Fixed parameters
Pr  = 10;
Rp  = 2;
tau = 0.01;

%% -------------------- output directory --------------------
out_dir = fullfile('Figures', 'Figure4');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
data_path = fullfile(out_dir, 'Figure4_data.mat');

%% -------------------- compute --------------------
if exist(data_path, 'file')
    fprintf('Loading cached data from %s ...\n', data_path);
    load(data_path, 'lambda_max', 'Ri_vec', 'omega_vec');
else
    % Flatten (Ri, ω) grid for parfor
    [Ri_grid, omega_grid] = ndgrid(Ri_vec, omega_vec);  % n_Ri × n_omega
    n_pairs = numel(Ri_grid);

    fprintf('Grid: %d Ri × %d ω = %d pairs.  (k,m0): %d × %d = %d points.\n', ...
        n_Ri, n_omega, n_pairs, n_k, n_m0, n_k * n_m0);
    fprintf('Total Floquet computations: %d\n', n_pairs * n_k * n_m0);

    lambda_flat = zeros(n_pairs, 1);

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
        interval = max(1, round(n_pairs / 20));
        state    = containers.Map({'done', 'n', 't0'}, {0, n_pairs, tic});
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) report_progress(state, interval));

        parfor pi = 1:n_pairs
            p = floquet_params(0, 'Ri', Ri_grid(pi), 'omega', omega_grid(pi));

            % --- sweep (k, m0) for max λ ---
            max_lam = 0;
            for ki = 1:n_k
                k_val = k_vec(ki);
                for mi = 1:n_m0
                    lam = floquet_core(k_val, m0_vec(mi), 0, p, n_steps);
                    if lam > max_lam, max_lam = lam; end
                end
            end
            lambda_flat(pi) = max_lam;

            if mod(pi, interval) == 0
                send(q, pi);
            end
        end
        fprintf('\n');
    else
        fprintf('  Serial mode.\n');
        for pi = 1:n_pairs
            p = floquet_params(0, 'Ri', Ri_grid(pi), 'omega', omega_grid(pi));

            max_lam = 0;
            for ki = 1:n_k
                k_val = k_vec(ki);
                for mi = 1:n_m0
                    lam = floquet_core(k_val, m0_vec(mi), 0, p, n_steps);
                    if lam > max_lam, max_lam = lam; end
                end
            end
            lambda_flat(pi) = max_lam;

            if mod(pi, max(1, round(n_pairs / 20))) == 0
                fprintf('  %d/%d done (%.0f s)\n', pi, n_pairs, toc);
            end
        end
    end

    fprintf('  Computation: %.1f s\n', toc);

    % Reshape to (n_Ri, n_omega) → transpose for imagesc (n_omega, n_Ri)
    lambda_max = reshape(lambda_flat, n_Ri, n_omega)';  % (omega, Ri)

    save(data_path, 'lambda_max', 'Ri_vec', 'omega_vec', ...
        'k_vec', 'm0_vec', 'n_steps', 'Pr', 'Rp', 'tau');
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
set(gca, 'FontSize', 14);
xlabel('$Ri$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 22);
colormap('jet');
clim(clims);
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10}(\lambda_{\max})$';
cb.Label.FontSize = 18;
title('Maximal Floquet growth rate $\log_{10}(\lambda_{\max})$ vs $(Ri, \omega)$', ...
    'Interpreter', 'latex', 'FontSize', 16);

saveas(gcf, fullfile(out_dir, 'Figure4_Radko.png'));
fprintf('Figure saved to %s/Figure4_Radko.png\n', out_dir);
fprintf('===== Done =====\n');
