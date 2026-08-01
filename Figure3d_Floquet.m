% Figure3d_Floquet.m
% Standalone reproduction of Honor Thesis Fig. 3(d).
% Floquet growth-rate heatmap over (k, m0) at fixed l, using shared floquet_core.
%
% Supports single eta (scalar) or eta_vec (list) — when a list is given,
% one heatmap is produced per eta value.
%
% Usage: adjust n_samples, n_steps, l, f, eta_vec, use_log10 below, then run.

close all; clc;

fprintf('===== Fig. 3(d) -- Floquet Heatmap =====\n');

%% -------------------- user settings --------------------
n_samples = 126;
n_steps   = 2000;
f_val     = 0;         % Coriolis parameter
eta_vec   = [0, 0.1, 0.5, 1, 1.2, 1.5, 1.8, 2];         % second-harmonic amplitude ratio

% l-scanning: f=0 → l=0 always (Coriolis is the only l-symmetry breaker)
if f_val ~= 0
    n_l  = 21;
    l_vec = linspace(-0.03, 0.03, n_l);
else
    n_l  = 1;
    l_vec = 0;
end

use_log10 = true;      % true → log10 colour scale; false → linear
log10_floor = -4;      % truncate log10 below this (only when use_log10 = true)

%% -------------------- shared grids --------------------
k_vec  = linspace(-0.5, 0.5, n_samples);
m0_vec = linspace(0, 1.5, n_samples);

%% -------------------- setup --------------------
for ei = 1:length(eta_vec)
    eta = eta_vec(ei);

    % --- setup parameters for this eta ---
    p = floquet_params(f_val, 'eta', eta);

    if eta > 0
        fprintf('\n--- eta = %.3f (%d/%d) ---\n', eta, ei, length(eta_vec));
        fprintf('  f = %.1f, n_l = %d\n', f_val, n_l);
        fprintf('  au1 = %.4f, av1 = %.4f, au2 = %.4f, av2 = %.4f\n', p.au, p.av, p.au2, p.av2);
    else
        fprintf('\n--- eta = 0 (single-frequency) ---\n');
        fprintf('  f = %.1f, n_l = %d, au = %.4f, av = %.4f\n', f_val, n_l, p.au, p.av);
    end
    fprintf('  Grid: %d × %d = %d points, %d time steps.\n', ...
        n_samples, n_samples, n_samples^2, n_steps);

    %% --- output directory ---
    if eta > 0
        out_dir = fullfile('Figures', 'Figure3d', sprintf('f_%.3f', f_val), sprintf('eta_%.3f', eta));
    else
        out_dir = fullfile('Figures', 'Figure3d', sprintf('f_%.3f', f_val));
    end
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    %% --- compute & save data ---
    data_path = fullfile(out_dir, 'Figure3d_data.mat');

    if exist(data_path, 'file')
        fprintf('  Loading cached data from %s ...\n', data_path);
        ld = load(data_path, 'lambda_mat');
        lambda_mat = ld.lambda_mat;
    else
        tic;
        lambda_mat = zeros(n_samples, n_samples);  % (m0, k)

        % --- parallel pool detection ---
        use_parallel = ~isempty(gcp('nocreate'));
        if ~use_parallel
            try
                parpool('local');
                use_parallel = true;
            catch
                use_parallel = false;
            end
        end

        if use_parallel
            n_workers = gcp('nocreate').NumWorkers;
            fprintf('  Parallel: %d iterations on %d workers.\n', n_samples, n_workers);
            interval = max(1, round(n_samples / 20));
            state    = containers.Map({'done', 'n', 't0'}, {0, n_samples, tic});
            q = parallel.pool.DataQueue;
            afterEach(q, @(~) report_progress(state, interval));

            parfor ki = 1:n_samples
                k_val = k_vec(ki);
                col = zeros(n_samples, 1);
                for mj = 1:n_samples
                    max_lam = 0;
                    for li = 1:n_l
                        lam = floquet_core(k_val, m0_vec(mj), l_vec(li), p, n_steps);
                        if lam > max_lam, max_lam = lam; end
                    end
                    col(mj) = max_lam;
                end
                lambda_mat(:, ki) = col;
                if mod(ki, interval) == 0
                    send(q, ki);
                end
            end
            fprintf('\n');
        else
            fprintf('  Serial mode (%d iterations).\n', n_samples);
            for ki = 1:n_samples
                k_val = k_vec(ki);
                for mj = 1:n_samples
                    max_lam = 0;
                    for li = 1:n_l
                        lam = floquet_core(k_val, m0_vec(mj), l_vec(li), p, n_steps);
                        if lam > max_lam, max_lam = lam; end
                    end
                    lambda_mat(mj, ki) = max_lam;
                end
                if mod(ki, max(1, round(n_samples / 20))) == 0
                    fprintf('  %d/%d done\n', ki, n_samples);
                end
            end
        end
        fprintf('  Elapsed: %.1f s\n', toc);

        save(data_path, 'lambda_mat', 'k_vec', 'm0_vec', 'l_vec', 'f_val', 'eta', 'n_samples', 'n_steps');
        fprintf('  Data saved to %s\n', data_path);
    end

    %% --- summary ---
    [mx_val, idx] = max(lambda_mat(:));
    [m0i, ki] = ind2sub(size(lambda_mat), idx);
    fprintf('  Max lambda = %.6f (log10 = %.2f) at k = %.4f, m0 = %.4f', ...
        mx_val, log10(mx_val), k_vec(ki), m0_vec(m0i));
    % Find optimal l at best (k, m0)
    if n_l > 1
        best_k = k_vec(ki);  best_m0 = m0_vec(m0i);
        lam_best = 0;  li_best = 1;
        for li = 1:n_l
            lam = floquet_core(best_k, best_m0, l_vec(li), p, n_steps);
            if lam > lam_best, lam_best = lam; li_best = li; end
        end
        fprintf(', l = %.4f\n', l_vec(li_best));
    else
        fprintf('\n');
    end

    %% --- plot ---
    if use_log10
        data = log10(lambda_mat);
        data(lambda_mat <= 0) = log10_floor;   % stable → floor
        data(data < log10_floor) = log10_floor;
        clims = [log10_floor, -1];
        cb_label = '$\log_{10}(\lambda)$';
    else
        data = lambda_mat;
        data(data < 0) = 0;
        mx = max(data(:));
        clims = [0, max(1e-6, mx * 1.05)];
        cb_label = '$\lambda$';
    end

    figure('Units', 'inches', 'Position', [1 1 8 6]);
    imagesc(k_vec, m0_vec, data);
    axis xy;
    set(gca, 'FontSize', 14);
    xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 22);
    ylabel('$m_0$', 'Interpreter', 'latex', 'FontSize', 22);
    xticks(-0.5:0.25:0.5);
    yticks(0:0.2:1.5);
    ylim([0 1.2]);
    colormap('jet');
    clim(clims);
    cb = colorbar;
    cb.Label.Interpreter = 'latex';
    cb.Label.String = cb_label;
    cb.Label.FontSize = 18;
    if n_l > 1
        if eta > 0
            title(sprintf('Floquet instability ($\\eta = %.2f$, $l = %.3f$)', ...
                eta, l_vec(li_best)), 'Interpreter', 'latex', 'FontSize', 18);
        else
            title(sprintf('Floquet instability ($l = %.3f$)', l_vec(li_best)), ...
                'Interpreter', 'latex', 'FontSize', 18);
        end
    else
        if eta > 0
            title(sprintf('Floquet instability ($\\eta = %.2f$, $l = 0$)', eta), ...
                'Interpreter', 'latex', 'FontSize', 18);
        else
            title(sprintf('Floquet instability for $l = 0$'), ...
                'Interpreter', 'latex', 'FontSize', 18);
        end
    end

    saveas(gcf, fullfile(out_dir, 'Figure3d_Floquet.png'));
    fprintf('  Figure saved to %s/Figure3d_Floquet.png\n', out_dir);
end  % for ei (eta loop)

fprintf('\n===== Done (%d eta values) =====\n', length(eta_vec));
