% plot_slices.m
% 2-D slice heatmaps of Floquet growth rate λ over (k, m0) and (l, m0).
% Computes the full 3-D grid once per f via floquet_scan_3d, then extracts
% every slice:  n_l heatmaps of (k, m0) and n_k heatmaps of (l, m0).
%
% Runs f = 0, 0.1 using Radko (2019) Eq. 40.

close all; clear; clc;

fprintf('===== 2-D Slice Heatmaps =====\n');

%% -------------------- user settings --------------------
n_k     = 51;
n_m0    = 51;
n_l     = 51;
n_steps = 1000;

k_vec  = linspace(-0.5, 0.5, n_k);
m0_vec = linspace(0, 1.5, n_m0);
l_vec  = linspace(0, 0.2, n_l);

f_vals = [0, 0.1];

use_log10     = true;     % true → log10 colour scale; false → linear
log10_floor   = -4;       % truncate log10 below this (only when use_log10 = true)

%% ==================== loop over f ====================
for fi = 1:length(f_vals)
    f_val = f_vals(fi);
    p = floquet_params(f_val);

    fprintf('\n===== f = %.1f | au = %.4f | av = %.4f =====\n', f_val, p.au, p.av);

    %% -------- output directories --------
    task_dir = fullfile('Figures', 'Slices', sprintf('f_%.3f', f_val));
    km_dir   = fullfile(task_dir, 'km_slices');
    lm_dir   = fullfile(task_dir, 'lm_slices');
    data_dir = fullfile(task_dir, 'data');
    for d = {km_dir, lm_dir, data_dir}
        if ~exist(d{1}, 'dir'), mkdir(d{1}); end
    end
    % Clean old PNGs to avoid mixing runs
    delete(fullfile(km_dir, '*.png'));
    delete(fullfile(lm_dir, '*.png'));

    %% -------- compute & save 3-D grid --------
    data_path = fullfile(data_dir, 'full_3d.mat');
    lambda_3d = floquet_scan_3d(p, k_vec, m0_vec, l_vec, n_steps, data_path);

    % Global max
    [max_global, idx] = max(lambda_3d(:));
    [m0i_best, ki_best, li_best] = ind2sub(size(lambda_3d), idx);
    fprintf('  Global max λ = %.6f  at k = %.4f, m0 = %.4f, l = %.4f\n', ...
        max_global, k_vec(ki_best), m0_vec(m0i_best), l_vec(li_best));

    %% -------- (k, m0) heatmaps: one per l --------
    fprintf('  Generating %d (k,m0) heatmaps...\n', n_l);
    for li = 1:n_l
        l_val = l_vec(li);
        lambda_km = lambda_3d(:, :, li);
        [mx, idx_km] = max(lambda_km(:));
        [~, ki] = ind2sub(size(lambda_km), idx_km);

        [data, clims, cb_label] = prepare_plot_data(lambda_km, use_log10, log10_floor);

        fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 8 6]);
        plot_heatmap(k_vec, m0_vec, data, clims, cb_label, '$k$', '$m_0$', [0 1.2], ...
            sprintf('$f=%.1f,\\; l=%.4f$  ($\\max\\lambda=%.4f$ at $k=%.2f$)', ...
                f_val, l_val, mx, k_vec(ki)));

        fname = sprintf('km_l%04d.png', round(l_val * 10000));
        saveas(fig, fullfile(km_dir, fname));
        close(fig);
    end

    %% -------- (l, m0) heatmaps: one per k --------
    fprintf('  Generating %d (l,m0) heatmaps...\n', n_k);
    for ki = 1:n_k
        k_val = k_vec(ki);
        lambda_lm = squeeze(lambda_3d(:, ki, :));  % (m0, l)
        [mx, idx_lm] = max(lambda_lm(:));
        [~, li] = ind2sub(size(lambda_lm), idx_lm);

        [data, clims, cb_label] = prepare_plot_data(lambda_lm, use_log10, log10_floor);

        fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1 1 8 6]);
        plot_heatmap(l_vec, m0_vec, data, clims, cb_label, '$l$', '$m_0$', [0 1.2], ...
            sprintf('$f=%.1f,\\; k=%.3f$  ($\\max\\lambda=%.4f$ at $l=%.4f$)', ...
                f_val, k_val, mx, l_vec(li)));

        fname = sprintf('lm_k%+04d.png', round(k_val * 1000));
        saveas(fig, fullfile(lm_dir, fname));
        close(fig);
    end

    fprintf('  f = %.1f: %d + %d = %d heatmaps saved.\n', f_val, n_l, n_k, n_l + n_k);
end

fprintf('\n===== Done =====\n');


%% ==================== local helpers ====================

function [data, clims, cb_label] = prepare_plot_data(lambda_in, use_log10, log10_floor)
    if use_log10
        data = log10(lambda_in);
        data(lambda_in <= 0) = log10_floor;
        data(data < log10_floor) = log10_floor;
        clims = [log10_floor, max(log10_floor + 1, max(data(:)))];
        cb_label = '$\log_{10}(\lambda)$';
    else
        data = lambda_in;
        data(data < 0) = 0;
        mx = max(data(:));
        clims = [0, max(1e-6, mx * 1.05)];
        cb_label = '$\lambda$';
    end
end

function plot_heatmap(x_vec, y_vec, data, clims, cb_label, x_lbl, y_lbl, y_lim, ttl)
    imagesc(x_vec, y_vec, data);
    axis xy;
    set(gca, 'FontSize', 14);
    xlabel(x_lbl, 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(y_lbl, 'Interpreter', 'latex', 'FontSize', 20);
    ylim(y_lim);
    colormap('jet');
    if clims(2) > clims(1)
        clim(clims);
    end
    cb = colorbar;
    cb.Label.Interpreter = 'latex';
    cb.Label.String = cb_label;
    cb.Label.FontSize = 16;
    title(ttl, 'Interpreter', 'latex', 'FontSize', 16);
end
