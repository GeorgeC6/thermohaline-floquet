% Figure3d_Floquet.m
% Standalone reproduction of Fig. 3(d) from the thermohaline-shear instability paper.
% Uses Floquet theory: computes the monodromy matrix over one period T,
% extracts Floquet multipliers mu = eig(Phi(T,0)), and returns
%   lambda = log(max|mu|) / T   (energy growth rate).
%
% Domain: k in [-0.5, 0.5],  m0 in [0, 1.5]
%
% Usage:
%   1. Tune n_samples and n_steps below for speed vs. resolution.
%   2. Run the script.  A heatmap is saved to Figures/Figure3d_Floquet.png

close all; clear; clc;

fprintf('===== Fig. 3(d) -- Floquet Heatmap (standalone) =====\n');

%% -------------------- parameters (Table / Eq. 19) --------------------
Ri    = 2;        % mean Richardson number
Rp    = 2;        % background density ratio
Pr    = 10;       % Prandtl number, fixed
omega = 0.5;      % background-shear frequency
tau   = 0.01;     % diffusivity ratio, fixed
l     = 0.5;      % spanwise wavenumber (non-zero for 3-D, Av coupling)
                    % l for loop
f     = 0;        % Coriolis parameter (2-D case)

T     = 2*pi/omega;          % one period
au    = sqrt(2*Pr*(Rp-1)/Ri); % shear amplitude (Eq. 19)

% time-dependent quantities (Eq. 14)
Au = @(t) au * sin(omega * t);                      % A_U(t)
Av = @(t) au * cos(omega * t);                      % A_V(t)
Bu = @(t) (au/omega) * (1 - cos(omega * t));         % B_U(t) = int_0^t Au
Bv = @(t) (au/omega) * sin(omega * t);              % B_V(t) = int_0^t Av

%% -------------------- resolution --------------------
n_samples = 50;      % grid points in k and m0  (126 for full resolution)
n_steps   = 2000;    % time steps per period   (2000 in paper)
n_samples_slow = n_samples;  % rename for clarity

% try to use parallel pool if available
use_parallel = ~isempty(gcp('nocreate'));
if use_parallel && isempty(gcp('nocreate'))
    try
        parpool('local');
        use_parallel = true;
    catch
        use_parallel = false;
    end
end

%% -------------------- grid --------------------
k_vec  = linspace(-0.5, 0.5, n_samples);
m0_vec = linspace(0, 1.5, n_samples);

lambda_mat = zeros(length(m0_vec), length(k_vec));  % m0 = rows, k = cols

%% -------------------- Floquet computation --------------------
fprintf('Grid: %d x %d = %d points, %d time steps each.\n', ...
    length(k_vec), length(m0_vec), length(k_vec)*length(m0_vec), n_steps);

tic;

if use_parallel
    % Parallel over k-index (each worker does one column)
    parfor ki = 1:length(k_vec)
        k_val = k_vec(ki);
        col_lambda = zeros(length(m0_vec), 1);
        for mj = 1:length(m0_vec)
            col_lambda(mj) = compute_floquet_log10lambda(k_val, m0_vec(mj), ...
                T, n_steps, au, omega, Pr, Rp, tau, l, f);
        end
        lambda_mat(:, ki) = col_lambda;
        fprintf('  k = %.3f  done (%d/%d)\n', k_val, ki, length(k_vec));
    end
else
    % Serial loop
    for ki = 1:length(k_vec)
        k_val = k_vec(ki);
        for mj = 1:length(m0_vec)
            lambda_mat(mj, ki) = compute_floquet_log10lambda(k_val, m0_vec(mj), ...
                T, n_steps, au, omega, Pr, Rp, tau, l, f);
        end
        fprintf('  k = %.3f  done (%d/%d)\n', k_val, ki, length(k_vec));
    end
end

elapsed = toc;
fprintf('Elapsed: %.1f s\n', elapsed);

%% -------------------- save data --------------------
out_dir = 'Figures';
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

% threshold: truncate < -4
lambda_plot = lambda_mat;
lambda_plot(lambda_plot < -4) = -4;

writematrix(lambda_mat, fullfile(out_dir, 'Figure3d_Floquet.csv'));

%% -------------------- plot --------------------
figure('Units', 'inches', 'Position', [1 1 8 6]);
imagesc(k_vec, m0_vec, lambda_plot);
axis xy;
set(gca, 'FontSize', 14);
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$m_0$', 'Interpreter', 'latex', 'FontSize', 22);
xticks(-0.5:0.25:0.5);
yticks(0:0.2:1.5);
ylim([0 1.2]);

colormap('jet');
clim([-4, -1]);
cb = colorbar;
cb.Label.Interpreter = 'latex';
cb.Label.String = '$\log_{10}(\lambda)$';
cb.Label.FontSize = 18;

saveas(gcf, fullfile(out_dir, 'Figure3d_Floquet.png'));
fprintf('Figure saved to Figures/Figure3d_Floquet.png\n');
fprintf('===== Done =====\n');


%% ================ helper function ================
function log10lambda = compute_floquet_log10lambda(k, m0, T, n_steps, ...
        au, omega, Pr, Rp, tau, l, f)
    % Compute log10 of the Floquet growth rate for a single (k, m0) pair.
    % Returns -Inf if the system is stable (lambda <= 0) or singular.

    dt = T / n_steps;
    t_list = linspace(0, T, n_steps + 1);  % n_steps+1 points

    % --- precompute m(t) and c(t) for all time steps ---
    m_t = zeros(size(t_list));
    c_t = zeros(size(t_list));
    Au_t = zeros(size(t_list));
    Av_t = zeros(size(t_list));

    for i = 1:length(t_list)
        t = t_list(i);
        Au_t(i)  = au * sin(omega * t);
        Av_t(i)  = au * cos(omega * t);
        Bu_t     = (au/omega) * (1 - cos(omega * t));
        Bv_t     = (au/omega) * sin(omega * t);
        m_t(i)   = m0 - Bu_t * k - Bv_t * l;  % Eq. 14: m(t) = m0 - B_U*k - B_V*l
        c_t(i)   = k^2 + l^2 + m_t(i)^2;      % c = k^2 + l^2 + m^2
    end

    % --- check for singularity (c = 0 when k=0, m(t)=0) ---
    if any(c_t < 1e-12)
        log10lambda = -Inf;
        return;
    end

    % --- monodromy matrix: Phi(T, 0) = expm(M(t_N)*dt) ... expm(M(t_1)*dt) ---
    Phi = eye(5);

    for i = 1:(length(t_list) - 1) % ！！(-1)
        c      = c_t(i);
        m      = m_t(i);
        Au_val = Au_t(i);
        Av_val = Av_t(i);

        % assemble the 5x5 system matrix M(t) (full 3-D, Eq. from FluidSystemSolver)
        M = zeros(5, 5);

        M(1,1) = -c;
        M(1,5) = 1;

        M(2,2) = -tau * c;
        M(2,5) = Rp;

        M(3,1) = -Pr * k * m / c;
        M(3,2) =  Pr * k * m / c;
        M(3,3) =  (f * k * l / c) - Pr * c;
        M(3,4) =  f - (f * k^2 / c);
        M(3,5) =  2 * (k^2 * Au_val + k * l * Av_val) / c - Au_val;

        M(4,1) = -Pr * l * m / c;
        M(4,2) =  Pr * l * m / c;
        M(4,3) =  (f * l^2 / c) - f;
        M(4,4) = -Pr * c - (f * k * l / c);
        M(4,5) =  2 * (k * l * Au_val + l^2 * Av_val) / c - Av_val;

        M(5,1) = -Pr * (m^2 / c - 1);
        M(5,2) =  Pr * (m^2 / c - 1);
        M(5,3) =  f * l * m / c;
        M(5,4) = -f * k * m / c;
        M(5,5) =  2 * (k * m * Au_val + l * m * Av_val) / c - Pr * c;

        % Phi = expm(M * dt) * Phi   (product-integral, Eq. 12-13)
        Phi = expm(M * dt) * Phi;
    end

    % --- Floquet multipliers and growth rate ---
    floquet_multipliers = eig(Phi);
    spectral_radius = max(abs(floquet_multipliers));

    if spectral_radius <= 1
        log10lambda = -Inf;      % stable (no growth)
    else
        lambda = log(spectral_radius) / T;   % Eq. 12
        log10lambda = log10(lambda);
    end
end
