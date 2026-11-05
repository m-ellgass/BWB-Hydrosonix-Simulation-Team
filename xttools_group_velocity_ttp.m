function [gpvel, conf95,tpeak, best_fit,pk] = xttools_group_velocity_ttp(x, t, v)%, maxdist, niter)

    % Time-to-peak for each x,y position
    dt = t(2) - t(1);
    [pk, tidx] = max(v, [], 2, 'omitnan');
    tpeak = (tidx - 1)*dt;
    
    % Use RANSAC to fit a linear slope
    % best_fit = fdml_ransac_linear(x, tpeak, maxdist, niter);
    % gpvel = 1./best_fit(1);

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

%     figure
%     set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
%     tt = [min(tpeak), max(tpeak)];
%     plot(tpeak*1000, x*1000, '.b', 'MarkerSize', 8, 'LineWidth', 1)
%     hold on
%     plot(tt*1000, gpvel*(tt - best_fit(2))*1000, '-r', 'LineWidth', 1)
%     plot(tt*1000, (gpvel + conf95(1))*(tt - best_fit(2))*1000, '--m', 'LineWidth', 1)
%     plot(tt*1000, (gpvel + conf95(2))*(tt - best_fit(2))*1000, '--m', 'LineWidth', 1)
%     xlabel('Time (ms)')
%     ylabel('Lateral Position (mm)')
%     box off