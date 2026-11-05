function gpvel = xttools_group_velocity_ncc(x, t, v, kernel_nx, nfilt)

    % Sizes and spacing
    [nx, nt] = size(v);
    dx = x(2) - x(1);
    dt = t(2) - t(1);

    % Convert to analytic signal
    nfft = 2^nextpow2(nt);
    vspect = fft(v, nfft, 2);
    vspect(:, (nfft/2 + 1):nfft) = 0;
    vspect(:, 2:(nfft/2)) = 2*vspect(:, 2:(nfft/2));
    analytic_signal = ifft(vspect, nfft, 2);
    analytic_signal = analytic_signal(:, 1:nt);

    % Tracking and computing normalized cross-correlation
    ncc = zeros(2*nt - 1, nx-kernel_nx);
    for i = 1:(nx - kernel_nx)
        reference_signal = squeeze(analytic_signal(i, :));
        test_signal = squeeze(analytic_signal(i + kernel_nx, :));
        ncc(:, i) = xcorr_td(test_signal, reference_signal, 'coeff');
    end

    % Correlation filter
    nccf = zeros(size(ncc));
    win = window(@gausswin, nfilt);
    win = win/sum(win);
    for i = 1:(2*nt - 1)
        nccf(i, :) = conv(ncc(i, :), win(:), 'same');
    end

    % Coarse estimation followed by subpixel estimation by phase zero-crossing
    [maxncc, maxidx] = max(abs(nccf), [], 1);
    lags = zeros(size(maxncc));
    for i = 1:(nx - kernel_nx)
        % Get phase in a three-point neighborhood of the ncc peak
        c = nccf(:, i);
        neighborhood = maxidx(i) + [-1, 0, 1];
        if any(neighborhood == 0) || any(neighborhood == length(c))
            lags(i) = maxidx(i) - nt;
        else
            % Subpixel estimation by phase zero-crossing
            p = angle(c(neighborhood));
            
            % Analytic solution for zero crossing of a three-point linear fit
            zero_cross = 2/3*sum(p)/(p(1) - p(3));
            
            % Lags estimate for ncc peak
            lags(i) = maxidx(i) - nt + zero_cross;
        end
    end

    gpvel = (dx*kernel_nx)./(lags*dt);