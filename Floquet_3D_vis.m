% Figure3D_Floquet_5x5.m
% 3D visualization of the 5x5 Thermohaline-Shear Instability Model
% Retains the original 5x5 matrix dynamics but visualizes in (k, l, m0) 3D space.

close all; clear; clc;

fprintf('===== 3D Floquet Heatmap (5x5 Thermohaline Model) =====\n');

%% -------------------- Parameters --------------------
Ri    = 2.0;      % mean Richardson number
Rp    = 2.0;      % density ratio
Pr    = 10;       % Prandtl number
omega = 0.5;      % oscillation frequency
tau   = 0.01;     % diffusivity ratio
f     = 0.1;      % Coriolis parameter (set >0 to see 3D asymmetric rotation effects)

T     = 2*pi/omega; 
au    = sqrt(2*Pr*(Rp-1)/Ri); 
av    = au;       % Assuming equivalent shear amplitude in v-direction

%% -------------------- Resolution & 3D Grid --------------------
n_k   = 51;    
n_l   = 51;    
n_m0  = 16;    
n_steps = 1000; % Time steps per period

k_vec  = linspace(-0.5, 0.5, n_k);
l_vec  = linspace(-0.5, 0.5, n_l);
m0_vec = linspace(0.0001, 0.2, n_m0);

% 尝试开启并行运算
use_parallel = ~isempty(gcp('nocreate'));
if ~use_parallel
    try
        parpool('local');
        use_parallel = true;
    catch
        use_parallel = false;
    end
end

lambda_mat = zeros(n_k, n_l, n_m0); % 注意：这里存的是真实生长率，不是log10

fprintf('Grid: %d x %d x %d = %d points, %d time steps.\n', ...
    n_k, n_l, n_m0, n_k*n_l*n_m0, n_steps);

tic;

%% -------------------- Floquet Computation (3D Loop) --------------------
if use_parallel
    parfor ki = 1:n_k
        k_val = k_vec(ki);
        slice_lambda = zeros(n_l, n_m0);
        for li = 1:n_l
            l_val = l_vec(li);
            for mi = 1:n_m0
                slice_lambda(li, mi) = compute_floquet_lambda(...
                    k_val, l_val, m0_vec(mi), T, n_steps, au, av, omega, Pr, Rp, tau, f);
            end
        end
        lambda_mat(ki, :, :) = slice_lambda;
        fprintf('  k-slice %d/%d done.\n', ki, n_k);
    end
else
    for ki = 1:n_k
        k_val = k_vec(ki);
        for li = 1:n_l
            l_val = l_vec(li);
            for mi = 1:n_m0
                lambda_mat(ki, li, mi) = compute_floquet_lambda(...
                    k_val, l_val, m0_vec(mi), T, n_steps, au, av, omega, Pr, Rp, tau, f);
            end
        end
        fprintf('  k-slice %d/%d done.\n', ki, n_k);
    end
end

elapsed = toc;
fprintf('Elapsed: %.1f s\n', elapsed);

%% -------------------- Dynamic Thresholding & 3D Plot --------------------
max_lambda = max(lambda_mat(:));
fprintf('Maximum growth rate found: %.4f\n', max_lambda);

if max_lambda <= 0
    error('The system is completely stable in this grid range. Nothing to plot!');
end

% 动态决定等值面的阈值（保证不管增长率多大，都能画出三层嵌套模型）
iso_val1 = 0.20 * max_lambda; % 外围低增长率 (半透明)
iso_val2 = 0.50 * max_lambda; % 中间层
iso_val3 = 0.85 * max_lambda; % 核心高增长率 (不透明)

[L_grid, K_grid, M0_grid] = meshgrid(l_vec, k_vec, m0_vec);

figure('Units', 'inches', 'Position', [1 1 9 7]);
hold on;

% 第一层外壳：浅蓝色半透明
fv1 = isosurface(K_grid, L_grid, M0_grid, lambda_mat, iso_val1);
patch(fv1, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.2);

% 第二层中壳：青色半透明
fv2 = isosurface(K_grid, L_grid, M0_grid, lambda_mat, iso_val2);
patch(fv2, 'FaceColor', 'cyan', 'EdgeColor', 'none', 'FaceAlpha', 0.5);

% 第三层内核：黄色不透明
fv3 = isosurface(K_grid, L_grid, M0_grid, lambda_mat, iso_val3);
patch(fv3, 'FaceColor', 'yellow', 'EdgeColor', 'none', 'FaceAlpha', 0.9);

% 渲染光照和视角
camlight('headlight'); 
lighting gouraud;
view([45, 30]);
grid on; box on;

set(gca, 'FontSize', 14);
xlabel('$k$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$l$', 'Interpreter', 'latex', 'FontSize', 22);
zlabel('$m_0$', 'Interpreter', 'latex', 'FontSize', 22);
xlim([min(k_vec), max(k_vec)]);
ylim([min(l_vec), max(l_vec)]);
zlim([min(m0_vec), max(m0_vec)]);

title('3D Growth Rate Visualization', 'Interpreter', 'latex', 'FontSize', 18);
fprintf('===== Done =====\n');


%% ================ Helper Function ================
function lambda = compute_floquet_lambda(k, l, m0, T, n_steps, au, av, omega, Pr, Rp, tau, f)
    % 返回线性的增长率 lambda（如果不稳定）。稳定则返回 0。
    
    dt = T / n_steps;
    t_list = linspace(0, T, n_steps + 1);

    m_t = zeros(size(t_list));
    c_t = zeros(size(t_list));
    Au_t = zeros(size(t_list));
    Av_t = zeros(size(t_list));

    for i = 1:length(t_list)
        t = t_list(i);
        Au_t(i)  = au * sin(omega * t);
        Av_t(i)  = av * cos(omega * t);
        Bu_t     = (au/omega) * (1 - cos(omega * t));
        Bv_t     = (av/omega) * sin(omega * t);
        
        m_t(i)   = m0 - Bu_t * k - Bv_t * l;  
        c_t(i)   = k^2 + l^2 + m_t(i)^2;      
    end

    if any(c_t < 1e-12)
        lambda = 0;
        return;
    end

    Phi = eye(5);

    for i = 1:(length(t_list) - 1)
        c      = c_t(i);
        m      = m_t(i);
        Au_val = Au_t(i);
        Av_val = Av_t(i);

        %  5x5 矩阵
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

        % Forward Euler product-integral
        Phi = expm(M * dt) * Phi;
    end

    floquet_multipliers = eig(Phi);
    spectral_radius = max(abs(floquet_multipliers));

    if spectral_radius <= 1
        lambda = 0;      % stable (no growth)
    else
        lambda = log(spectral_radius) / T;  % 返回真实的 lambda
    end
end