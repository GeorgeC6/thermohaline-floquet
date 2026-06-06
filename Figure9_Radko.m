% Figure9_Radko_3D.m
% Reproduction of Fig. 9 from Radko (2019): 3D visualization of growth rate lambda_r
% Uses Floquet theory to find lambda_r as a function of (k, l, m0)
%
% Domain: k in [-3e-3, 3e-3], l in [-3e-3, 3e-3], m0 in [0, 0.01]

close all; clear; clc;

fprintf('===== Radko (2019) Fig. 9 -- 3D Floquet Heatmap =====\n');

%% -------------------- Parameters (Radko 2019, Section 5) --------------------
Ri    = 1.0;      % Mean Richardson number (Typical for Fig 8 & 9)
Pr    = 10;       % Prandtl number
omega = 0.1;      % Background-shear frequency (Near-inertial, typical is 0.1 or 0.2)
f     = 0.1;      % Coriolis parameter

T     = 2*pi/omega; % One period

% Shear amplitudes with planetary rotation (Eq. 40 in Radko 2019)
amp_factor = sqrt(2 * Pr / (Ri * (f^2 + omega^2)));
au = amp_factor * omega; 
av = amp_factor * f;

%% -------------------- Resolution --------------------
% Note: 3D grids grow very fast! 31x31x16 = 15376 points.
% Increase these numbers (e.g. to 50) for publication-quality resolution.
n_k  = 31;    
n_l  = 31;    
n_m0 = 16;    
n_steps = 1000;  % Time steps per period (1000 is usually enough for stability)

% Setup parallel pool
use_parallel = ~isempty(gcp('nocreate'));
if ~use_parallel
    try
        parpool('local');
        use_parallel = true;
    catch
        use_parallel = false;
    end
end

%% -------------------- 3D Grid Generation --------------------
k_vec  = linspace(-3e-3, 3e-3, n_k);
l_vec  = linspace(-3e-3, 3e-3, n_l);
m0_vec = linspace(0.0001, 0.01, n_m0); % Start slightly above 0 to avoid singularity

% lambda_mat will hold the linear growth rate (NOT log10)
lambda_mat = zeros(n_k, n_l, n_m0);

fprintf('Grid: %d x %d x %d = %d points, %d time steps each.\n', ...
    n_k, n_l, n_m0, n_k*n_l*n_m0, n_steps);

tic;

%% -------------------- Floquet Computation (3D Loop) --------------------
if use_parallel
    % Parallel over k-index
    parfor ki = 1:n_k
        k_val = k_vec(ki);
        slice_lambda = zeros(n_l, n_m0);
        for li = 1:n_l
            l_val = l_vec(li);
            for mi = 1:n_m0
                m0_val = m0_vec(mi);
                slice_lambda(li, mi) = compute_radko_lambda(k_val, l_val, m0_val, ...
                    T, n_steps, au, av, omega, Pr, f);
            end
        end
        lambda_mat(ki, :, :) = slice_lambda;
        fprintf('  k-slice %d/%d done.\n', ki, n_k);
    end
else
    % Serial execution
    for ki = 1:n_k
        k_val = k_vec(ki);
        for li = 1:n_l
            l_val = l_vec(li);
            for mi = 1:n_m0
                m0_val = m0_vec(mi);
                lambda_mat(ki, li, mi) = compute_radko_lambda(k_val, l_val, m0_val, ...
                    T, n_steps, au, av, omega, Pr, f);
            end
        end
        fprintf('  k-slice %d/%d done.\n', ki, n_k);
    end
end

elapsed = toc;
fprintf('Elapsed time: %.1f s\n', elapsed);

%% -------------------- Plot 3D Isosurface --------------------
[L_grid, K_grid, M0_grid] = meshgrid(l_vec, k_vec, m0_vec);

figure('Units', 'inches', 'Position', [1 1 9 7]);
hold on;

% Draw nested isosurfaces to mimic the volumetric "alien spider/butterfly" look
% Iso-value 1: Low growth rate (Outer translucent shell)
fv1 = isosurface(K_grid, L_grid, M0_grid, lambda_mat, 0.01);
p1 = patch(fv1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% Iso-value 2: Medium growth rate
fv2 = isosurface(K_grid, L_grid, M0_grid, lambda_mat, 0.03);
p2 = patch(fv2, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% Iso-value 3: High growth rate (Inner solid core)
fv3 = isosurface(K_grid, L_grid, M0_grid, lambda_mat, 0.06);
p3 = patch(fv3, 'FaceColor', 'yellow', 'EdgeColor', 'none', 'FaceAlpha', 0.9);

% Lighting and Camera configuration
camlight('headlight'); 
lighting gouraud;
view([45, 30]); % Adjust viewing angle to match Fig 9
grid on; box on;

% Axes formatting
set(gca, 'FontSize', 14);
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$l$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('$m_0$', 'Interpreter', 'latex', 'FontSize', 22);
xlim([-3e-3, 3e-3]);
ylim([-3e-3, 3e-3]);
zlim([0, 0.01]);

title('3D Visualization of Growth Rate $\lambda_r$', 'Interpreter', 'latex', 'FontSize', 18);
fprintf('===== Done =====\n');


%% ================ Helper Function ================
function lambda = compute_radko_lambda(k, l, m0, T, n_steps, au, av, omega, Pr, f)
    % Compute the REAL growth rate (not log10) for the 4x4 Radko model.
    % Returns 0 if the system is stable.
    
    dt = T / n_steps;
    t_list = linspace(0, T, n_steps + 1);

    % Precompute time-dependent parameters (Eq. 39 and integrals)
    m_t = zeros(size(t_list));
    c_t = zeros(size(t_list));
    Au_t = zeros(size(t_list));
    Av_t = zeros(size(t_list));

    for i = 1:length(t_list)
        t = t_list(i);
        Au_t(i)  = au * sin(omega * t);
        Av_t(i)  = av * cos(omega * t);
        
        % Integrals B_U and B_V
        Bu_t     = (au/omega) * (1 - cos(omega * t));
        Bv_t     = (av/omega) * sin(omega * t);
        
        m_t(i)   = m0 - Bu_t * k - Bv_t * l;  
        c_t(i)   = k^2 + l^2 + m_t(i)^2;     
    end

    if any(c_t < 1e-12)
        lambda = 0;
        return;
    end

    % Monodromy matrix for the 4x4 system [rho, u, v, w]
    Phi = eye(4);

    for i = 1:n_steps
        c      = c_t(i);
        m      = m_t(i);
        Au_val = Au_t(i);
        Av_val = Av_t(i);

        M = zeros(4, 4);

        % Row 1: Density eq: d(rho)/dt = w - c*rho
        M(1,1) = -c;
        M(1,4) = 1;

        % Row 2: u-momentum eq (Eliminating pressure via divergence-free condition)
        M(2,1) = Pr * k * m / c;
        M(2,2) = -Pr * c + f * k * l / c;
        M(2,3) = f * (1 - k^2 / c);
        M(2,4) = 2 * k * (k * Au_val + l * Av_val) / c - Au_val;

        % Row 3: v-momentum eq
        M(3,1) = Pr * l * m / c;
        M(3,2) = -f * (1 - l^2 / c);
        M(3,3) = -Pr * c - f * k * l / c;
        M(3,4) = 2 * l * (k * Au_val + l * Av_val) / c - Av_val;

        % Row 4: w-momentum eq
        M(4,1) = -Pr * (1 - m^2 / c);
        M(4,2) = f * l * m / c;
        M(4,3) = -f * k * m / c;
        M(4,4) = 2 * m * (k * Au_val + l * Av_val) / c - Pr * c;

        % Phi = expm(M * dt) * Phi
        Phi = expm(M * dt) * Phi;
    end

    % Floquet multiplier and actual growth rate extraction
    floquet_multipliers = eig(Phi);
    spectral_radius = max(abs(floquet_multipliers));

    if spectral_radius <= 1
        lambda = 0;      % stable
    else
        lambda = log(spectral_radius) / T;   % actual growth rate lambda_r
    end
end