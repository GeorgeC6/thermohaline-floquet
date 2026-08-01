% eta_scan_curves.m
% λmax vs η — scan the second-harmonic amplitude ratio at fixed (k, m0) grid.
%
% For each η, finds max Floquet growth rate over (k, m0) at l=0, f=0.
% Plot λmax(η) as a line curve.

close all; clear; clc;

fprintf('===== η-scan: λmax(η) for second-harmonic amplitude =====\n');

%% -------------------- user settings --------------------
f_val     = 0;         % Coriolis parameter
eta_vec   = linspace(0, 2, 41);  % η values to scan

n_k      = 51;         % k wavenumber grid
n_m0     = 51;         % m0 grid
n_steps  = 1000;       % time steps per period

% Fixed parameters
Pr  = 10;
Rp  = 2;
tau = 0.01;

%% -------------------- shared grids --------------------
k_vec  = linspace(-0.5, 0.5, n_k);
m0_vec = linspace(0, 1.5, n_m0);

% l-scanning: f=0 → l=0 always (Coriolis is the only l-symmetry breaker)
if f_val ~= 0
    n_l   = 21;
    l_vec = linspace(-0.03, 0.03, n_l);
else
    n_l   = 1;
    l_vec = 0;
end

n_eta = length(eta_vec);

fprintf('f = %.1f (using default Ri, ω from floquet_params)\n', f_val);
fprintf('η range: %.3f – %.3f (%d points)\n', eta_vec(1), eta_vec(end), n_eta);
fprintf('Grid: %d k × %d m0 × %d l = %d points, %d steps.\n', n_k, n_m0, n_l, n_k*n_m0*n_l, n_steps);

%% -------------------- output directory --------------------
out_dir = fullfile('Figures', 'EtaScan', sprintf('f_%.3f', f_val));
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
data_path = fullfile(out_dir, 'EtaScan_data.mat');

%% -------------------- compute --------------------
if exist(data_path, 'file')
    fprintf('Loading cached data from %s ...\n', data_path);
    ld = load(data_path, 'lambda_max', 'eta_vec');
    % Validate eta_vec matches
    if isequal(ld.eta_vec(:), eta_vec(:))
        lambda_max = ld.lambda_max;
        fprintf('  Cache valid — %d η points loaded.\n', n_eta);
    else
        fprintf('  eta_vec mismatch — recomputing.\n');
        clear ld;
    end
end

if ~exist('lambda_max', 'var')
    lambda_max = zeros(n_eta, 1);

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
        interval = max(1, round(n_eta / 20));
        state    = containers.Map({'done', 'n', 't0'}, {0, n_eta, tic});
        q = parallel.pool.DataQueue;
        afterEach(q, @(~) report_progress(state, interval));

        parfor ei = 1:n_eta
            eta = eta_vec(ei);
            p = floquet_params(f_val, 'eta', eta);

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
            lambda_max(ei) = max_lam;

            if mod(ei, interval) == 0
                send(q, ei);
            end
        end
        fprintf('\n');
    else
        fprintf('  Serial mode.\n');
        for ei = 1:n_eta
            eta = eta_vec(ei);
            p = floquet_params(f_val, 'eta', eta);

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
            lambda_max(ei) = max_lam;

            if mod(ei, max(1, round(n_eta / 20))) == 0
                fprintf('  %d/%d done (%.0f s)\n', ei, n_eta, toc);
            end
        end
    end

    fprintf('  Computation: %.1f s\n', toc);
    save(data_path, 'lambda_max', 'eta_vec', 'f_val', ...
        'k_vec', 'm0_vec', 'l_vec', 'n_k', 'n_m0', 'n_l', 'n_steps', 'Pr', 'Rp', 'tau');
    fprintf('  Data saved to %s\n', data_path);
end

%% -------------------- summary --------------------
fprintf('\n--- Summary ---\n');
[mx_val, ei_peak] = max(lambda_max);
fprintf('Peak λmax = %.6f at η = %.3f\n', mx_val, eta_vec(ei_peak));
fprintf('η = 0 (single-freq): λmax = %.6f\n', lambda_max(1));

%% -------------------- plot: λmax vs η --------------------
figure('Units', 'inches', 'Position', [1 1 8 6]);
plot(eta_vec, lambda_max, 'LineWidth', 2.5, 'Color', [0 0.4470 0.7410]);
hold on;
plot(eta_vec(ei_peak), mx_val, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
hold off;

set(gca, 'FontSize', 14);
xlabel('$\eta$ (second-harmonic amplitude ratio)', 'Interpreter', 'latex', 'FontSize', 20);
ylabel('$\lambda_{\max}$', 'Interpreter', 'latex', 'FontSize', 22);
xlim([eta_vec(1), eta_vec(end)]);
grid on; box on;
if n_l > 1
    title(sprintf('$\\lambda_{\\max}(\\eta)$ (max over $k,m_0,l$), $f=%.1f$', f_val), ...
        'Interpreter', 'latex', 'FontSize', 15);
else
    title(sprintf('$\\lambda_{\\max}(\\eta)$ (max over $k,m_0$, $l=0$), $f=%.1f$', f_val), ...
        'Interpreter', 'latex', 'FontSize', 15);
end

% --- data cursor ---
set(gcf, 'UserData', struct('type', 'eta_scan', 'eta_vec', eta_vec, 'lambda_max', lambda_max));
set(gcf, 'CreateFcn', 'setup_datatip(gcf)');
setup_datatip(gcf);

saveas(gcf, fullfile(out_dir, 'EtaScan_curves.png'));
savefig(gcf, fullfile(out_dir, 'EtaScan_curves.fig'));
fprintf('Figure saved to %s/\n', out_dir);
fprintf('===== Done =====\n');
