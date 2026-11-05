function [pw,pk,snr] = pulse_width(dx,dt,c0,x,vsol,do_fig,dirname)
%Calculate spatial pulse width (m)

% do snr
[pk, ~] = max(abs(vsol), [], 2, 'omitnan');
%noise=1;

nx=size(vsol,1);
i = ceil((((nx/2)*dx)/c0)/dt);
vi = vsol(:,i);

% Find envelope and -6dB push length
vi_envelope = abs(hilbert(vi));
[~, idxmax] = max(vi_envelope);

idxmax=round(.005/dx);%force location input for 1d model

xmax = x(idxmax);
%[a b]=max(vi);
vi_envelope_dB = 20*log10(vi_envelope./max(vi_envelope));

idx1 = find(vi_envelope_dB > -6, 1, 'first');
slope_first = (vi_envelope_dB(idx1) - vi_envelope_dB(idx1 - 1))/dx;

if slope_first==NaN
    slope_first=3.0557e+04;
end

idx_dB1 = x(idx1 - 1) - (vi_envelope_dB(idx1 - 1) + 6)/slope_first;

idx2  = find(vi_envelope_dB > -6, 1, 'last');
slope_last = (vi_envelope_dB(idx2 + 1) - vi_envelope_dB(idx2))/dx;
if slope_first==NaN
    slope_first= -9.4348e+03;
end
idx_dB2 = x(idx2) - (vi_envelope_dB(idx2) + 6)/slope_last;

pw = idx_dB2 - idx_dB1;
%pwi=round(pw/dx);

% Find frequency spectrum
spectrum_ampl = abs(fftshift(fft(vi)));
spectrum_freq = (0:(nx-1))/nx*(1/dx) - (1/dx/2);

pwi=round(.0015/dx);
tmp=circshift(vi_envelope,idxmax+pwi);
noise=tmp(pwi*2:end);

noise2=sqrt((mean(noise))^2);
%sign=sqrt((mean(pk))^2);
sign=.7610; %known input for 1D solution post filter
snr=20*log10(sign/noise2);

%    figure,
%     plot(pk2,'k'),ylim([0 1])
%     hold on,plot(noise2,'b'),ylim([0 1])
%     legend(['signal rms= ',num2str(sign)],['noise rms= ',num2str(noise)])
%     title(['SNR = ',num2str(snr),' dB'])
%     
%     

if do_fig
    %% Plot pulse width and SNR
    figure;
    sgtitle(['Spatial Pulse Width = ',num2str(pw*1000),' mm'])
    subplot(311)
    plot(x*1000, vi, '-k')
    hold on
    plot(x*1000, vi_envelope, '--b')
    xlabel('x (mm)')
    ylabel('Vibration Speed (Norm)')
    
    subplot(312)
    plot(spectrum_freq/1000, 20*log10(spectrum_ampl./max(spectrum_ampl)), '-k')
    xlim([0, 8])
    ylim([-30, 0])
    xlabel('Frequency (1/mm)')
    ylabel('Spectrum (dB)')
    legend(['SNR = ',num2str(snr),' dB'])
    
    subplot(313)
    plot((x - xmax)*1000, vi_envelope_dB, '-b');
    hold on
    plot((x - xmax)*1000, -6*ones(size(x)), '-r');
    plot((x - xmax)*1000, -snr*ones(size(x)), '-g');
    xlabel('x - x_{peak} (mm)')
    ylabel('Envelope (dB)')
    xlim([-5, 5])
    ylim([-30, 0])
    box off
    print(fullfile(dirname,['pulsewidth',num2str(snr,1)]), '-dpng', '-r300')


% figure, 
% plot(vi_envelope,'k'),hold on,
% plot(tmp,'b'),
% plot(noise,'g')

end

    
end

