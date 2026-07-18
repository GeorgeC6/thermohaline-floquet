function p = floquet_params(f_val)
% floquet_params  Shared physical parameters and Radko (2019) Eq. 40 amplitudes.
%
%   p = floquet_params(f) returns a struct with all fixed parameters,
%   the period T, and the shear amplitudes au, av for Coriolis parameter f.
%
%   Convention: AU = au*cos(ωt), AV = -av*sin(ωt)  (phase-shifted Radko 3-D,
%   equivalent to time translation of Eq. 39; recovers Honor Thesis 2-D at f=0, l=0).

p.Ri    = 2;
p.Rp    = 2;
p.Pr    = 10;
p.omega = 0.5;
p.tau   = 0.01;
p.f     = f_val;

p.T     = 2*pi / p.omega;

% Radko (2019) Eq. 40
p.amp = sqrt(2 * p.Pr * (p.Rp - 1) / (p.Ri * (p.f^2 + p.omega^2)));
p.au  = p.amp * p.omega;
p.av  = p.amp * p.f;
end
