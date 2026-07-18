function lambda = floquet_core(k, m0, l, p, n_steps)
% floquet_core  Floquet growth rate for a single (k, m0, l) point.
%
%   lambda = floquet_core(k, m0, l, p, n_steps) computes the energy growth
%   rate lambda = ln(max|μ|)/T for the 5×5 thermohaline-shear ODE system.
%   p is a parameter struct from floquet_params().
%
%   Convention: AU = au*cos(ωt), AV = -av*sin(ωt)  (π/2 phase shift from Radko 3-D).
%   At f=0, l=0 this exactly recovers the Honor Thesis 2-D case.
%
%   Returns 0 if the system is stable (spectral radius ≤ 1) or singular (c ≈ 0).

dt = p.T / n_steps;
t_list = linspace(0, p.T, n_steps + 1);

% --- precompute time-dependent coefficients ---
m_t  = zeros(size(t_list));
c_t  = zeros(size(t_list));
Au_t = zeros(size(t_list));
Av_t = zeros(size(t_list));

for i = 1:length(t_list)
    t = t_list(i);
    Au_t(i) = p.au * sin(p.omega * t + pi/2);          % = p.au * cos(ωt)
    Av_t(i) = p.av * cos(p.omega * t + pi/2);          % = -p.av * sin(ωt)
    Bu_t    = (p.au / p.omega) * sin(p.omega * t);     % ∫ cos = sin/ω
    Bv_t    = (p.av / p.omega) * (cos(p.omega * t) - 1); % ∫(-sin) = (cos-1)/ω
    m_t(i)  = m0 - Bu_t * k - Bv_t * l;
    c_t(i)  = k^2 + l^2 + m_t(i)^2;
end

% --- singularity guard ---
if any(c_t < 1e-12)
    lambda = 0;
    return;
end

% --- monodromy matrix Φ(T, 0) = ∏ expm(M(t)·dt) ---
Phi = eye(5);

for i = 1:(length(t_list) - 1)
    c  = c_t(i);
    m  = m_t(i);
    Au = Au_t(i);
    Av = Av_t(i);

    % 5×5 system matrix (full 3-D)
    M = zeros(5, 5);

    M(1,1) = -c;
    M(1,5) = 1;

    M(2,2) = -p.tau * c;
    M(2,5) = p.Rp;

    M(3,1) = -p.Pr * k * m / c;
    M(3,2) =  p.Pr * k * m / c;
    M(3,3) =  (p.f * k * l / c) - p.Pr * c;
    M(3,4) =  p.f - (p.f * k^2 / c);
    M(3,5) =  2 * (k^2 * Au + k * l * Av) / c - Au;

    M(4,1) = -p.Pr * l * m / c;
    M(4,2) =  p.Pr * l * m / c;
    M(4,3) =  (p.f * l^2 / c) - p.f;
    M(4,4) = -p.Pr * c - (p.f * k * l / c);
    M(4,5) =  2 * (k * l * Au + l^2 * Av) / c - Av;

    M(5,1) = -p.Pr * (m^2 / c - 1);
    M(5,2) =  p.Pr * (m^2 / c - 1);
    M(5,3) =  p.f * l * m / c;
    M(5,4) = -p.f * k * m / c;
    M(5,5) =  2 * (k * m * Au + l * m * Av) / c - p.Pr * c;

    Phi = expm(M * dt) * Phi;
end

% --- Floquet multipliers → growth rate ---
floquet_multipliers = eig(Phi);
spectral_radius = max(abs(floquet_multipliers));

if spectral_radius <= 1
    lambda = 0;
else
    lambda = log(spectral_radius) / p.T;
end
end
