function [gpvel, conf95, tpeak, best_fit] = xttools_group_velocity_env(x, t, v)% maxdist, niter)

    % Sizes and spacing
    [nx, nt] = size(v);
    dt = t(2) - t(1);

    % Convert to analytic signal and get envelope
    nfft = 2^nextpow2(nt); % smallest power of 2 that is greater than paranthesis
    vspect = fft(v, nfft, 2);
    vspect(:, (nfft/2 + 1):nfft) = 0;
    vspect(:, 2:(nfft/2)) = 2*vspect(:, 2:(nfft/2));
    analytic_signal = ifft(vspect, nfft, 2); % inverse fourier transform
    analytic_signal = analytic_signal(:, 1:nt);
    env = abs(analytic_signal);

    % Time-to-peak for each x,y position
    dt = t(2) - t(1); % unnecessary redefinition? see line 5
    [~, tidx] = max(env, [], 2, 'omitnan');
    tpeak = (tidx - 1)*dt;
    
    % Use RANSAC to fit a linear slope
    % best_fit = fdml_ransac_linear(x, tpeak, maxdist, niter);
    % gpvel = 1./best_fit(1);
    warning('off','all')
    % Use "fit" with robust option for a linear slope
    [xData, yData] = prepareCurveData( x, tpeak );
    ft = fittype( 'poly1' );
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
    opts.Robust = 'LAR';
    [fitresult, gof] = fit( xData, yData, ft, opts );
    best_fit = [fitresult.p1, fitresult.p2];
    gpvel = 1./fitresult.p1;

    % Get 95% confidence interval on the inverse of the fit
    % x = v(t - b)
    n = length(x);
    xhat = gpvel*(tpeak - best_fit(2));
    sigma2 = sum((x' - xhat).^2)/(n - 2);
    sxx = sum((tpeak - mean(tpeak)).^2);
    conf95 = tinv([0.025, 0.975], n-2)*sqrt(sigma2/sxx); 

