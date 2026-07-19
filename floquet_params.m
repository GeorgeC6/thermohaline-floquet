function p = floquet_params(f_val, varargin)
% floquet_params  Shared physical parameters and Radko (2019) Eq. 40 amplitudes.
%
%   p = floquet_params(f) returns a struct with defaults (Ri=2, Rp=2, Pr=10,
%   omega=0.5, tau=0.01) plus derived T, amp, au, av for Coriolis parameter f.
%
%   p = floquet_params(f, 'Name', Value, ...) overrides any default before
%   computing derived quantities.  Example:
%       p = floquet_params(0, 'Ri', 5, 'omega', 0.3);
%
%   Convention: AU = au*cos(ωt), AV = -av*sin(ωt)  (π/2 phase-shifted Radko
%   3-D; recovers Honor Thesis 2-D at f=0, l=0).

%% --- defaults ---
Ri    = 2;
Rp    = 2;
Pr    = 10;
omega = 0.5;
tau   = 0.01;

%% --- parse optional name-value overrides ---
for i = 1:2:(length(varargin) - 1)
    name  = varargin{i};
    value = varargin{i + 1};
    switch lower(name)
        case 'ri',    Ri    = value;
        case 'rp',    Rp    = value;
        case 'pr',    Pr    = value;
        case 'omega', omega = value;
        case 'tau',   tau   = value;
        otherwise
            warning('floquet_params: unknown parameter "%s"', name);
    end
end

%% --- assemble struct ---
T   = 2 * pi / omega;
amp = sqrt(2 * Pr * (Rp - 1) / (Ri * (f_val^2 + omega^2)));

p = struct('Ri', Ri, 'Rp', Rp, 'Pr', Pr, 'omega', omega, 'tau', tau, ...
           'f', f_val, 'T', T, 'amp', amp, ...
           'au', amp * omega, 'av', amp * f_val);
end
