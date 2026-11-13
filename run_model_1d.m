% A function that sets the parameters for and runs the sw1d.m file to 
% simulate a one-dimensional ultrasound pulse through a material

% start with a clean workspace so past results don't interfere
clear % removes variables from current workspace
close all % closes all open figure and plot windows
clc % clears all previous command window text/output

%              Simulating the "Push"
%{ 
    A Gaussian Function (Bell Curve) represents the ultrasound pulse's
    amplitude over time and the spread of its frequencies. A real lens 
    focus is only approximately Gaussian, so this matches the width at 
    half height while treating the rest of the shape as an approximation.

    Push parameters - b is the Gaussian width parameter. 
                    - fwhm is the width of the focus of the transducer
                        ^- either directly measured or from data sheet
    Use a conversion to set a known full-width at half-max
%}
%fwhm = 0.001; % narrower setting (can swap in)
fwhm = 0.01; % unit m
b = fwhm/(2*sqrt(2*log(2)));

%              Simulating the Space
% Mesh/domain details
% N is the discrete points at which data is plotted. Larger N gives a more
% detailed picture of the wave but takes longer to compute
L = 0.08;           % Total domain size 
N = ceil(L/fwhm)*36;   % Number of grid points, approximately 36 per push width

% Time stepping
T = 0.0049;         % solve until t = T, unit s (i.e. simulate T seconds)
dtmin = 0.005E-6;   % minimum time step (set to ''false'' to use cfl)
dtout = 10E-6;      % sampling time step for output (decimate solutions)
cfl = 0.1;          % stability condition - smaller number is safer but slower

%               Material properties:
%{ 
    Sound speed function - step function starting at c0 and transitioning to 
    c1 half-way through the domain
    Two speeds (c0 and c1) allow us to model a boundary between two materials
    Density uniformly equal to 1000 kg/m^3
%}
% c0 = 2.136;
% c1 = 2.136;
c = 6;
c0 = c; % 6m/s was approximately measured in my own skin
c1 = c; % c0 = c1 represents uniform material with no boundary

% Push location (where the ultrasound push happens)
x0 = L/4;
% Initial condition amplitude (how strong the push is)
u0 = 1;

% Run the simulation sw1d.m, passing in the variables defined above 
% Returns usol, the wave's displacement at every position (x) and
% recorded time (t)
% differentiate this result to convert displacement to velocity
[usol, x, t] = sw1d(N, L, b, c0, c1, T, cfl, dtmin, dtout, x0, u0, 'gausspulse');
vsol = diff(usol, 1, 1)/(t(2) - t(1));
v=vsol'; % transpose to dimensionalize same as OnScale (x,t)

cs=c0; % save wave speed used for simulation

% store the resulting simulation in a .mat file (basically a container for
% the variables)
save(['1d_sim_cs_',num2str(c0,1)],'v','x','t','cs')
