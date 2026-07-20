% explore_3d_scan.m
% Quick exploration: 3-D (k, l, m0) Floquet isosurfaces for selected (Ri, ω)
% pairs.  Used to determine wavenumber ranges for Figure 4 upgrade to 3-D.
%
% If the max λ always occurs at l=0, then the 2-D (k, m0) sweep is sufficient.
% If it shifts to l ≠ 0 for some (Ri, ω), we need full 3-D in Figure4_Radko.m.

close all; clear; clc;

fprintf('===== 3-D Exploration for Selected (Ri, ω) =====\n');

%% -------------------- user settings --------------------
n_k     = 51;
n_m0    = 51;
n_l     = 51;
n_steps = 1000;

% Wavenumber ranges
k_vec  = linspace(-0.5, 0.5, n_k);
m0_vec = linspace(0, 1.5, n_m0);
l_vec  = linspace(0, 0.3, n_l);   % wider l than plot_slices to check for shifts

% Test points: (Ri, ω) — corners of parameter space
test_points = [
    1.0, 0.1;    % low Ri,  low ω
    1.0, 0.5;    % low Ri,  mid ω
    1.0, 0.8;    % low Ri,  high ω
    5.0, 0.1;    % mid Ri,  low ω
    5.0, 0.8;    % mid Ri,  high ω
    9.0, 0.1;    % high Ri, low ω
    9.0, 0.8;    % high Ri, high ω
];

f_val = 0;  % 2-D → 3-D extension, no Coriolis

%% ==================== loop over test points ====================
for ti = 1:size(test_points, 1)
    Ri_val    = test_points(ti, 1);
    omega_val = test_points(ti, 2);

    fprintf('\n===== Ti=%d: Ri=%.1f, ω=%.2f =====\n', ti, Ri_val, omega_val);

    % --- parameter struct ---
    p = floquet_params(f_val, 'Ri', Ri_val, 'omega', omega_val);
    fprintf('  au = %.4f, av = %.4f, T = %.4f\n', p.au, p.av, p.T);

    % --- output directory ---
    out_dir = fullfile('Figures', 'Explore3D', sprintf('Ri%.1f_omega%.2f', Ri_val, omega_val));
    if ~exist(out_dir, 'dir'), mkdir(out_dir); end

    % --- compute 3-D grid ---
    data_path = fullfile(out_dir, 'lambda_3d.mat');
    lambda_3d = floquet_scan_3d(p, k_vec, m0_vec, l_vec, n_steps, data_path);

    % --- analysis ---
    [max_all, idx] = max(lambda_3d(:));
    [m0i, ki, li] = ind2sub(size(lambda_3d), idx);
    fprintf('  Max λ = %.6f  at k=%.4f, m0=%.4f, l=%.4f\n', ...
        max_all, k_vec(ki), m0_vec(m0i), l_vec(li));

    % Compare with l=0 slice
    lambda_l0 = lambda_3d(:, :, 1);  % first l index = l=0
    max_l0 = max(lambda_l0(:));
    [m0i_l0, ki_l0] = find(lambda_l0 == max_l0, 1);
    fprintf('  Max at l=0:  λ = %.6f  at k=%.4f, m0=%.4f\n', ...
        max_l0, k_vec(ki_l0), m0_vec(m0i_l0));
    fprintf('  3-D / l=0 ratio: %.4f  (diff: %+.6f)\n', max_all / max(max_l0, 1e-12), max_all - max_l0);

    % Check if max is at edge of l-range (may need wider range)
    if li == 1
        fprintf('  ⚠ Max at l=0 (left edge of l-range).\n');
    elseif li == n_l
        fprintf('  ⚠ Max at l=%.4f (RIGHT EDGE — consider widening l range).\n', l_vec(li));
    end

    % Check if max is at edge of k-range
    if ki == 1 || ki == n_k
        fprintf('  ⚠ Max at k-edge (k=%.4f) — consider widening k range.\n', k_vec(ki));
    end

    % --- 3-D isosurface ---
    fprintf('  Rendering isosurface...\n');

    [K_mg, L_mg, M0_mg] = ndgrid(k_vec, l_vec, m0_vec);
    V = permute(lambda_3d, [2, 3, 1]);

    iso_vals  = [0.20, 0.50, 0.85] * max_all;
    colors    = {[0.2 0.4 1.0], [0.0 0.7 0.9], [1.0 0.7 0.1]};
    alphas    = [0.15, 0.45, 0.85];

    fig = figure('Visible', 'on', 'Units', 'inches', 'Position', [1 1 10 8]);
    hold on;

    for i = 1:3
        if iso_vals(i) > 0
            fv = isosurface(K_mg, L_mg, M0_mg, V, iso_vals(i));
            pa = patch(fv);
            pa.FaceColor = colors{i};
            pa.EdgeColor = 'none';
            pa.FaceAlpha = alphas(i);
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
    title(sprintf('Ri=%.1f, \\omega=%.2f  (max \\lambda=%.4f at k=%.2f, l=%.3f, m_0=%.2f)', ...
        Ri_val, omega_val, max_all, k_vec(ki), l_vec(li), m0_vec(m0i)), ...
        'Interpreter', 'tex', 'FontSize', 13);

    saveas(fig, fullfile(out_dir, 'isosurface.png'));
    savefig(fig, fullfile(out_dir, 'isosurface.fig'));
    close(fig);

    % --- 2-D slice at l where max occurs ---
    fig2 = figure('Visible', 'on', 'Units', 'inches', 'Position', [1 1 8 6]);
    lambda_km = lambda_3d(:, :, li);
    data = log10(lambda_km);
    data(lambda_km <= 0) = -4;
    data(data < -4) = -4;
    imagesc(k_vec, m0_vec, data);
    axis xy;
    set(gca, 'FontSize', 14);
    xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 22);
    ylabel('$m_0$', 'Interpreter', 'latex', 'FontSize', 22);
    ylim([0 1.5]);
    colormap('jet');
    clim([-4, max(-1, max(data(:)))]);
    cb = colorbar;
    cb.Label.Interpreter = 'latex';
    cb.Label.String = '$\log_{10}(\lambda)$';
    cb.Label.FontSize = 16;
    title(sprintf('Ri=%.1f, \\omega=%.2f, l=%.4f  (\\lambda_{max}=%.4f)', ...
        Ri_val, omega_val, l_vec(li), max_all), 'Interpreter', 'tex', 'FontSize', 14);

    saveas(fig2, fullfile(out_dir, 'km_slice_at_max.png'));
    savefig(fig2, fullfile(out_dir, 'km_slice_at_max.fig'));
    close(fig2);

    fprintf('  Saved to %s/\n', out_dir);
end

fprintf('\n===== Done =====\n');
