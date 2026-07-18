% Figure3d_Floquet.m
% Standalone reproduction of Honor Thesis Fig. 3(d).
% Floquet growth-rate heatmap over (k, m0) at fixed l, using shared floquet_core.
%
% Usage: adjust n_samples, n_steps, l, f, use_log10 below, then run.

close all; clear; clc;

fprintf('===== Fig. 3(d) -- Floquet Heatmap =====\n');

%% -------------------- user settings --------------------
n_samples = 126;
n_steps   = 2000;
l_val     = 0;         % spanwise wavenumber
f_val     = 0;         % Coriolis parameter

use_log10 = true;      % true → log10 colour scale; false → linear
log10_floor = -4;      % truncate log10 below this (only when use_log10 = true)

%% -------------------- setup --------------------
p = floquet_params(f_val);

k_vec  = linspace(-0.5, 0.5, n_samples);
m0_vec = linspace(0, 1.5, n_samples);

fprintf('f = %.1f, l = %.3f, au = %.4f, av = %.4f\n', f_val, l_val, p.au, p.av);
fprintf('Grid: %d × %d = %d points, %d time steps.\n', ...
    n_samples, n_samples, n_samples^2, n_steps);

%% -------------------- output directory --------------------
out_dir = fullfile('Figures', 'Figure3d');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end

%% -------------------- compute & save data --------------------
data_path = fullfile(out_dir, 'Figure3d_data.mat');

if exist(data_path, 'file')
    fprintf('Loading cached data from %s ...\n', data_path);
    load(data_path, 'lambda_mat', 'k_vec', 'm0_vec');
else
    tic;
    lambda_mat = zeros(n_samples, n_samples);  % (m0, k)

    parfor ki = 1:n_samples
        k_val = k_vec(ki);
        col = zeros(n_samples, 1);
        for mj = 1:n_samples
            col(mj) = floquet_core(k_val, m0_vec(mj), l_val, p, n_steps);
        end
        lambda_mat(:, ki) = col;
    end
    fprintf('Elapsed: %.1f s\n', toc);

    save(data_path, 'lambda_mat', 'k_vec', 'm0_vec', 'l_val', 'f_val', 'n_samples', 'n_steps');
    fprintf('Data saved to %s\n', data_path);
end

%% -------------------- plot --------------------
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

[mx_val, idx] = max(lambda_mat(:));
[m0i, ki] = ind2sub(size(lambda_mat), idx);
fprintf('Max lambda = %.6f (log10 = %.2f) at k = %.4f, m0 = %.4f\n', ...
    mx_val, log10(mx_val), k_vec(ki), m0_vec(m0i));

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
title(sprintf('Floquet instability for $l = %.3f$', l_val), ...
    'Interpreter', 'latex', 'FontSize', 18);

saveas(gcf, fullfile(out_dir, 'Figure3d_Floquet.png'));
fprintf('Figure saved to %s/Figure3d_Floquet.png\n', out_dir);
fprintf('===== Done =====\n');
