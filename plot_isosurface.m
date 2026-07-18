% plot_isosurface.m
% 3-D isosurface visualization of Floquet growth rate λ over (k, m0, l).
% Computes the full 3-D grid via floquet_scan_3d, then renders isosurfaces
% at 20 %, 50 %, 85 % of max λ.
%
% Runs f = 0, 0.1 using Radko (2019) Eq. 40.

close all; clear; clc;

fprintf('===== 3-D Isosurface Visualization =====\n');

%% -------------------- user settings --------------------
n_k     = 51;
n_m0    = 51;
n_l     = 51;
n_steps = 1000;

k_vec  = linspace(-0.5, 0.5, n_k);
m0_vec = linspace(0, 1.5, n_m0);
l_vec  = linspace(0, 0.2, n_l);

f_vals = [0, 0.1];

use_log10 = false;      % isosurface uses linear λ

%% ==================== loop over f ====================
for fi = 1:length(f_vals)
    f_val = f_vals(fi);
    p = floquet_params(f_val);

    fprintf('\n===== f = %.1f | au = %.4f | av = %.4f =====\n', f_val, p.au, p.av);

    %% -------- compute & save 3-D grid --------
    task_dir  = fullfile('Figures', 'Isosurface', sprintf('f_%.3f', f_val));
    data_path = fullfile(task_dir, 'Floquet_3D.mat');
    lambda_3d = floquet_scan_3d(p, k_vec, m0_vec, l_vec, n_steps, data_path);

    [max_lambda, idx] = max(lambda_3d(:));
    [~, k_idx, l_idx] = ind2sub(size(lambda_3d), idx);
    fprintf('  Max λ = %.6f  at k = %.4f, l = %.4f\n', ...
        max_lambda, k_vec(k_idx), l_vec(l_idx));

    %% -------- render isosurface --------
    fprintf('  Rendering 3-D isosurface...\n');

    [K_mg, L_mg, M0_mg] = ndgrid(k_vec, l_vec, m0_vec);
    V = permute(lambda_3d, [2, 3, 1]);

    iso_vals  = [0.20, 0.50, 0.85] * max_lambda;
    colors    = {[0.2 0.4 1.0], [0.0 0.7 0.9], [1.0 0.7 0.1]};
    alphas    = [0.15, 0.45, 0.85];

    fig = figure('Visible', 'on', 'Units', 'inches', 'Position', [1 1 10 8]);
    hold on;

    for i = 1:3
        if iso_vals(i) > 0
            fv = isosurface(K_mg, L_mg, M0_mg, V, iso_vals(i));
            p = patch(fv);
            p.FaceColor = colors{i};
            p.EdgeColor = 'none';
            p.FaceAlpha = alphas(i);
        end
    end

    camlight('headlight');
    lighting gouraud;
    view([50, 25]);
    grid on; box on;
    set(gca, 'FontSize', 14);
    xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 22);
    ylabel('$l$', 'Interpreter', 'latex', 'FontSize', 22);
    zlabel('$m_0$', 'Interpreter', 'latex', 'FontSize', 22);
    title(sprintf('$f=%.1f$  ($\\max\\lambda=%.4f$ at $l=%.3f$)', ...
        f_val, max_lambda, l_vec(l_idx)), 'Interpreter', 'latex', 'FontSize', 16);

    if ~exist(task_dir, 'dir'), mkdir(task_dir); end
    saveas(fig, fullfile(task_dir, 'Floquet_3D_isosurface.png'));
    savefig(fig, fullfile(task_dir, 'Floquet_3D_isosurface.fig'));
    close(fig);
    fprintf('  Figure saved to %s/\n', task_dir);
end

fprintf('\n===== Done =====\n');
