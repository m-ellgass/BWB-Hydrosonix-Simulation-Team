function [VdB, f, k] = xttools_power_spectrum(v, dt, dx, nfft)
    nfft = 2*ceil(nfft/2);
    fs = 1/dt;
    ks = 1/dx;
    f = (0:(nfft-1))/nfft*fs - fs/2;
    k = (0:(nfft-1))/nfft*ks - ks/2;
    vspect = fft2(v, nfft, nfft);
    VdB = 20*log10(abs(vspect)/max(abs(vspect(:))));
    clear vspect
    VdB = VdB(1:(nfft/2), (nfft/2+1):nfft);
    VdB = fliplr(VdB);
    f = fftshift(f);
    k = fftshift(k);
    f = f(1:nfft/2);
    k = k(1:nfft/2);