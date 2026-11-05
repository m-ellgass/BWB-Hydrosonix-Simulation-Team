function v = sd_filter_surface_wavefield(v, dx, dt, bwx, bwt)
% SD_FILTER_SURFACE_WAVEFIELD - Apply bandpass filtering to surface data
%
% v = sd_filter_surface_wavefield(v, dx, dy, dt, bwx, bwy, bwt)
%
% Parameters
% ----------
% v : [nx, nt] Surface wave field data (the xt data)
% dx : [double] x-direction spacing in meters
% dt : [double] time sampling interval in seconds
% bwx : [1x2 double] filter cutoffs for x-direction filtering in 1/m
% bwt : [1x2 double] filter cutoffs for temporal filtering in Hz
%
% Returns
% -------
% v : [nx, nt] Filtered surface wave field data (the xt data)
%
% Notes
% -----
% Pass an empty array to bwx or bwt to skip filtering in that dimension.
%
% ---


    % Get sampling rates for each dimension
    sample_rate_t = 1/dt;
    sample_rate_x = 1/dx;

%     % Get the valid range of v
%     first_valid = find(~isnan(surface_z), 1, 'first');
%     last_valid = find(~isnan(surface_z), 1, 'last');
% 
%     % Crop NaN rows from the data
%     v = v(first_valid:last_valid, :);

    % FFT lengths
    [nx, nt] = size(v);
    nfft_t = 2^nextpow2(nt);
    nfft_x = 2^nextpow2(nx);
    
    % Frequency/wavenumber vectors
    f = (0:(nfft_t-1))/nfft_t*sample_rate_t - sample_rate_t/2;
    k = (0:(nfft_x-1))/nfft_x*sample_rate_x - sample_rate_x/2;

    % Temporal filter
    if ~isempty(bwt)

        % HPF - Gaussian
        fc = bwt(1);
        g = exp(-4*log(2)*(f/(2*fc)).^2);
        filter_spectrum_t_hpf = 1 - g;

        % LPF - Gaussian x super Gaussian
        fc = bwt(2);
        f1 = sqrt((fc^2 + sqrt(fc^4 + 4*log(2)*(fc/1.2)^4))/(2*log(2)));
        f2 = 1.2*f1;
        g = exp(-(f/f1).^2);
        s = exp(-(f/f2).^4);
        filter_spectrum_t_lpf = s.*g;

        filter_spectrum_t = fftshift(filter_spectrum_t_hpf.*filter_spectrum_t_lpf);
        vspectrum = fft(v, nfft_t, 2);
        v = real(ifft(vspectrum.*filter_spectrum_t, nfft_t, 2));
        v = v(:, 1:nt);
    end

    % Spatial filter
    if ~isempty(bwx)

        % HPF - Gaussian
        kc = bwx(1);
        g = exp(-4*log(2)*(k/(2*kc)).^2);
        filter_spectrum_x_hpf = 1 - g;

        % LPF - Gaussian x super Gaussian
        kc = bwx(2);
        k1 = sqrt((kc^2 + sqrt(kc^4 + 4*log(2)*(kc/1.2)^4))/(2*log(2)));
        k2 = 1.2*k1;
        g = exp(-(k/k1).^2);
        s = exp(-(k/k2).^4);
        filter_spectrum_x_lpf = s.*g;

        filter_spectrum_x = fftshift(filter_spectrum_x_hpf.*filter_spectrum_x_lpf);
        filter_spectrum_x = reshape(filter_spectrum_x, nfft_x, 1);
        vspectrum = fft(v, nfft_x, 1);
        v = real(ifft(vspectrum.*filter_spectrum_x, nfft_x, 1));
        v = v(1:nx, :);
    end