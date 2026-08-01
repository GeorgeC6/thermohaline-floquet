% omega_scan_curves.m
% λmax vs ω for fixed Ri values.
%
% For each fixed Ri, sweep ω over [ω_min, ω_max] and find max λ over
% (k, m0, l) wavenumber grid.  Plot λmax(ω) curves for different Ri on
% the same axes — x = ω, y = λmax, each Ri a different coloured curve.

close all; clear; clc;

fprintf('===== ω-scan: λmax(ω) for fixed Ri =====\n');

%% -------------------- user settings --------------------
f_val    = 0.1;         % Coriolis parameter
Ri_vals  = 0.25;  % fixed Ri values
n_omega  = 100;       % ω grid points
omega_min = 0.01;
omega_max = 7;
n_k      = 51;        % k wavenumber grid for max search
n_m0     = 51;        % m0 grid for max search
n_l      = 21;        % l wavenumber grid
n_steps  = 1000;      % time steps per period

omega_vec = linspace(omega_min, omega_max, n_omega);
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
out_dir = fullfile('Figures', 'OmegaScan', sprintf('f_%.3f', f_val));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
data_path = fullfile(out_dir, 'OmegaScan_data.mat');

%% -------------------- compute --------------------
n_Ri = length(Ri_vals);

% Accuracy struct for cache validation
acc = struct('f_val', f_val, 'n_k', n_k, 'n_m0', n_m0, 'n_l', n_l, ...
    'n_steps', n_steps, 'Pr', Pr, 'Rp', Rp, 'tau', tau, ...
    'k_vec', k_vec, 'm0_vec', m0_vec, 'l_vec', l_vec);

[lambda_max, is_cached] = cache_manager(data_path, Ri_vals, omega_vec, acc);

if isempty(lambda_max)
    lambda_max = zeros(n_omega, n_Ri);
    is_cached = false(n_omega, n_Ri);
end

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

% Compute per Ri — skip Ri values that are fully cached
n_computed = 0;
for ri = 1:n_Ri
    Ri_val = Ri_vals(ri);
    if all(is_cached(:, ri))
        fprintf('  Ri = %.3f — cached, skipping.\n', Ri_val);
        continue;
    end

    fprintf('\n===== Ri = %.3f (%d/%d) =====\n', Ri_val, ri, n_Ri);
    fprintf('  ω range: %.3f – %.3f (%d points)\n', omega_min, omega_max, n_omega);
    n_computed = n_computed + 1;

    lambda_omega = zeros(n_omega, 1);
    tic;

    if use_parallel
        n_workers = gcp('nocreate').NumWorkers;
        fprintf('  Parallel: %d workers.\n', n_workers);
        interval = max(1, round(n_omega / 20));
        state    = containers.Map({'done', 'n', 't0'}, {0, n_omega, tic});
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) report_progress(state, interval));

        parfor oi = 1:n_omega
            p = floquet_params(f_val, 'Ri', Ri_val, 'omega', omega_vec(oi));

            % --- sweep (k, m0, l) for max λ ---
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
            lambda_omega(oi) = max_lam;

            if mod(oi, interval) == 0
                send(q, oi);
            end
        end
        fprintf('\n');
    else
        fprintf('  Serial mode.\n');
        for oi = 1:n_omega
            p = floquet_params(f_val, 'Ri', Ri_val, 'omega', omega_vec(oi));

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
            lambda_omega(oi) = max_lam;

            if mod(oi, max(1, round(n_omega / 20))) == 0
                fprintf('  %d/%d done (%.0f s)\n', oi, n_omega, toc);
            end
        end
    end

    lambda_max(:, ri) = lambda_omega;
    fprintf('  Ri=%.3f: %.1f s,  max λ = %.6f at ω = %.3f\n', ...
        Ri_val, toc, max(lambda_omega), ...
        omega_vec(lambda_omega == max(lambda_omega)));
end  % for ri loop

if n_computed > 0
    cache_acc = acc; %#ok<NASGU>
    save(data_path, 'lambda_max', 'Ri_vals', 'omega_vec', 'cache_acc');
    fprintf('\nData saved to %s (%d Ri computed, %d cached).\n', ...
        data_path, n_computed, n_Ri - n_computed);
else
    fprintf('All %d Ri values already cached — no computation needed.\n', n_Ri);
end

%% -------------------- summary --------------------
fprintf('\n--- Summary ---\n');
for ri = 1:n_Ri
    lam = lambda_max(:, ri);
    [mx, oi] = max(lam);
    fprintf('Ri = %.2f:  peak λmax = %.6f  at ω = %.3f\n', Ri_vals(ri), mx, omega_vec(oi));
end

%% -------------------- plot: λmax vs ω curves --------------------
figure('Units', 'inches', 'Position', [1 1 10 7]);
hold on;

% Colour map for Ri curves
cmap = lines(n_Ri);
h_lines = gobjects(n_Ri, 1);

for ri = 1:n_Ri
    h_lines(ri) = plot(omega_vec, lambda_max(:, ri), ...
        'LineWidth', 2.5, 'Color', cmap(ri, :), ...
        'DisplayName', sprintf('Ri = %.2f', Ri_vals(ri)));
end
hold off;

% --- data cursor: UserData + CreateFcn survives .fig save/load ---
set(gcf, 'UserData', struct('type', 'lines', 'Ri_vals', Ri_vals, 'omega_vec', omega_vec, ...
    'lambda_max', lambda_max, 'h_lines', h_lines));
set(gcf, 'CreateFcn', 'setup_datatip(gcf)');
setup_datatip(gcf);  % also set up immediately for the live figure

set(gca, 'FontSize', 14);
xlabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 22);
ylabel('$\lambda_{\max}$', 'Interpreter', 'latex', 'FontSize', 22);
xlim([omega_min, omega_max]);
legend('Location', 'best', 'FontSize', 14);
grid on; box on;
title(sprintf('Maximal Floquet growth rate $\\lambda_{\\max}(\\omega)$ for fixed $Ri$, $f=%.1f$', f_val), ...
    'Interpreter', 'latex', 'FontSize', 16);

saveas(gcf, fullfile(out_dir, 'OmegaScan_curves.png'));
savefig(gcf, fullfile(out_dir, 'OmegaScan_curves.fig'));
fprintf('Figure saved to %s/\n', out_dir);
fprintf('===== Done =====\n');
