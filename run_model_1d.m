clear all
close all
clc

% Push parameters - b is the Gaussian width parameter. 
% Use a conversion to set a known full-width at half-max
%fwhm = 0.001; %unit m
fwhm = 0.01; %unit m. This is the 'width' of the focus from the 3D printed lens/transducer
b = fwhm/(2*sqrt(2*log(2)));

% Mesh/domain details
L = 0.08;           % Total domain size 
N = ceil(L/fwhm)*36;   % Number of grid points, approximately 36 per wavelength

% Time stepping
T = 0.0049;         % solve until t = T, unit seconds
dtmin = 0.005E-6;   % minimum time step (set to ''false'' to use cfl)
dtout = 10E-6;      % sampling time step for output (decimate solutions)
cfl = 0.1;          % stability condition

% Material properties:
% Sound speed function - step function starting at c0 and transitioning to 
% c1 half-way through the domain
% Density uniformly equal to 1000 kg/m^3
% c0 = 2.136;
% c1 = 2.136;
c = 6;
c0 = c; % 6m/s was approximately measured in my own skin
c1 = c;

% Push location
x0 = L/4;

% Initial condition amplitude
u0 = 1;

% Run solver and differentiate
[usol, x, t] = sw1d(N, L, b, c0, c1, T, cfl, dtmin, dtout, x0, u0, 'gausspulse');
vsol = diff(usol, 1, 1)/(t(2) - t(1));
v=vsol'; %dimensionalize same as OnScale (x,t)

cs=c0;

save(['1d_sim_cs_',num2str(c0,1)],'v','x','t','cs')

% pause(1);
% 
% disp("Hi")