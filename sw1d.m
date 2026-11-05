function [usol, x, tsol] = sw1d(N, L, b, c0, c1, T, cfl, dtmin, dtout, x0, u0, pulseshape)
% SW1D - One dimensional shear wave propagation
%
% Solves the wave equation on a one dimensional string vibrating in two
% dimensional space. The domain is divided in half with separate wave
% speeds on either side.
%
% Parameters
% ----------
% N : int
%   Number of grid points
% L : double
%   Domain length in m
% b : double
%   Push spatial width in m
% c0 : double
%   Wave speed in the first half of domain in m/s
% c1 : double
%   Wave speed in the second half of domain in m/s
% T : double
%   End time of simulation
% cfl : double (optional)
%   CFL number for computing the time step, default = 0.25
% dtmin : double | false (optional)
%   Minimum time step, default based on CFL (dtmin = false)
% dtout : double | false (optional)
%   Approx output time step - only stores solutions at decimated rate
% x0 : double (optional)
%   Push location, default = L/4
% u0 : double (optional)
%   Initial condition amplitude, default = 1.0
%
% Returns
% -------
% usol : double of shape [nt, N]
%   Displacement solutions
% x : double of shape [1, N]
%   Domain mesh points
% t : double of shape [1, nt]
%   Time points of solution
%
% Notes
% -----
% * This function solves the dimensionless wave equation using the 
%   following scaling: u ~ b, x ~ b, t ~ b/c0
% * The equations are discretized in the following form: utt = [B*ux]x
% * The time derivative is second order central difference.
% * The spatial derivatives are discretized using the second-order
%   method of Nilsson2007. See references [1] and [2].
% * The boundary conditions are u(0,t) = 0 and u(L,t) = 0
% * The initial conditions are u(x,0) = g(x-x0)*cos(2*pi*(x-x0)/b),
%   where g is a Gaussian of width b, and ut(x,0) = 0.
%
% References
% ----------
% [1] Nillson et al. Stable difference approximations for the elastic wave
%     equation in second order formulation. SIAM J Numer Anal. 45(5) (2007), 
%     pp 1902-1936. 
% [2] Sjogreen and Petersson. A fourth order accurate finite difference 
%     scheme for the elastic wave equation in second order formulation. 
%     Lawrence Livermore National Laboratory, LLNL-JRNL-483427. (2011).
%     Accessed from: https://e-reports-ext.llnl.gov/pdf/491932.pdf
%
% ---

    if nargin < 7
        cfl = 0.25;
        dtmin = false;
        dtout = false;
        x0 = L/4;
        u0 = 1.0;
        pulseshape = 'gausspulse';
    elseif nargin < 8
        dtmin = false;
        dtout = false;
        x0 = L/4;
        u0 = 1.0;
        pulseshape = 'gausspulse';
    elseif nargin < 9
        dtout = false;
        x0 = L/4;
        u0 = 1.0;
        pulseshape = 'gausspulse';
    elseif nargin < 10
        x0 = L/4;
        u0 = 1.0;
        pulseshape = 'gausspulse';
    elseif nargin < 11
        u0 = 1.0;
        pulseshape = 'gausspulse';
    elseif nargin < 12
        pulseshape = 'gausspulse';
    end

    % Physical Grid
    h = L/(N-1);
    x = 0:h:((N-1)*h);

    % Choose the timestep to satisfy the CFL condition
    dt = cfl*h/max(c0, c1);
    if dtmin
        dt = min(dtmin, dt);
    end
    t = 0:dt:T;
    
    % Calculate decimation factor for storing output solutions
    decim = floor(dtout/dt);
    
    % Nondimensionalize the step sizes
    h = h/b;
    dt = c0*dt/b;
    hinv2 = 1/h/h;

    % Allocate storage
    % um, u, and up are solutions at times n-1, n, and n+1, respectively
    nt = length(t);
    usol = zeros(ceil(nt/decim) + 1, N);
    tsol = zeros(ceil(nt/decim) + 1, 1);
    u = zeros(1, N);
    um = zeros(1, N);
    up = zeros(1, N);

    % Wave speed squared, dimensionless
    B = ones(1, N);
    half = floor(N/2) + 1;
    B(half:end) = (c1/c0)^2;

    % Initial condition of the displacement field u(x,0), dimensionless
    % Assume the initial condition for the velocity is u'(x,0) = 0
    fwhm = (2*sqrt(2*log(2)))*b;
    gc = @(x, xi, ui, bi) (ui/bi)*exp(-(x - xi).^2/(2*bi^2)).*cos(2*pi*(x - xi)/fwhm);
    g = @(x, xi, ui, bi) (ui/bi)*exp(-(x - xi).^2/(2*bi^2));
    uinit = zeros(size(x(2:N-1)));
    for i = 1:length(x0)
        if strcmp(pulseshape, 'gausspulse') %Compares if pulse shape is = gausspulse, returns true if true
            uinit = uinit + gc(x(2:N-1), x0(i), u0(i), b(i));
        else
            uinit = uinit + g(x(2:N-1), x0(i), u0(i), b(i));
        end
    end

    % Set initial conditions
    u(2:N-1) = uinit;
    um(2:N-1) = uinit;
    usol(1,:) = u;

    % Time stepping - second order central difference in time, second order 
    % variable coefficient method in space.
    t = dt;
    counter = 1;
%     fprintf(1, ' ---------------------------------------------------\n');
%     fprintf(1, ' SW1D\n');
%     fprintf(1, ' ---------------------------------------------------\n');
%     fprintf(1, ' N = %i points, Nt = %i steps\n', N, nt);
%     fprintf(1, ' Progress:                                         ]\n');
%     fprintf(1, ' ');
    prevprog = 0;
    start_time = tic;
    for n = 2:nt
        % Progress
        prog = floor(n/nt*100/2);
%         if prevprog ~= prog
%             fprintf('*')
%         end
        prevprog = prog;
        
        % Advance solution in time
        Bm = 0.5*(B(2:N-1) + B(1:N-2));
        Bp = 0.5*(B(2:N-1) + B(3:N));
        rhs = hinv2*(Bp.*(u(3:N) - u(2:N-1)) - Bm.*(u(2:N-1) - u(1:N-2)));
        up(2:N-1) = 2*u(2:N-1) - um(2:N-1) + (dt^2)*rhs;
        
        % Update solution
        um = u;
        u = up;
        
        % Store solutions
        if ~mod(n+1, decim)
            usol(counter, :) = u;
            tsol(counter) = t;
            counter = counter + 1;
        end
        t = t + dt;
    end
%     fprintf(1, ']\n');
    elapsed_time = toc(start_time);
%     fprintf(1, ' Solver finished in %6.2f seconds.\n\n', elapsed_time);

    % Redimensionalize usol for output
    usol = usol(1:counter-1, :)*b;
    tsol = tsol(1:counter-1, :)*b/c0;