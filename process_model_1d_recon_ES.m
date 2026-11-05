%% Setup
clear all; close all; clc; 
set(0,'DefaultFigureWindowStyle','docked') %keeps figures docked in desktop version
aa_output=zeros(1,3); %creates output matrix for runs

% Load 1d simulation data
fname='1d_sim_cs_6.mat';

load(fname)
fprintf(1, ['Loading.....',fname,'\n']);
dirname='C:\Users\E-Mily\Desktop\BWB\Outputs';
mkdir(dirname)

%% Choose display options
do_disp=1;
do_figs=1;
%% add noise
do_noise=1;
n_dB=[10];%dB signal to noise, referred to as noise avg amplitude and peak wave displacement
%n_dB=[.01:.01:5]; 
%% Choose TOF processing scheme
do_filt=1;

do_corr=1;%correlation shift (NCC) TOF
corr=zeros(1,size(n_dB,2));

do_env=1;
envo=zeros(1,size(n_dB,2));

do_ttp=1;
ttp=zeros(1,size(n_dB,2));

snr=zeros(1,size(n_dB,2));

% Related to do_corr method
nfilt=10; %correlation filter of some sort.. revisit this
kernel_nx=20; %kernel size

% Related to RANSAC method
maxdist=.1; %unit m
niter=10000;%number of RANSAC iterations

%%set total scan limits
tend=0.1; %m
xend=.05;    %m
%tend=0.0109; %m
%xend=.01;    %m
nxin=256; % set # of x positions (a-scans), roughly on the order of OCT

%set TOF limits
xstarttof = 0.001;
xendtof = 0.02;
%xstarttof = 0.003;
%xendtof = 0.09;

%v(:,160:end)=0;
%% re-size array for analysis
[x,t,vsol]=sim_crop(x,t,v,xend,tend,nxin);

%% Raw video
% figure(10)
% for ii=1:size(t)
% plot(x*10^3,vsol(:,ii)),hold off
% xlabel('position (mm)')
% title(['time = ', num2str(t(ii)*10^3), ' msec'])
% ylim([-1 1])
% drawnow,
% 
% end
% v=v./(max(abs(v(:))));
% t(end)=[];
% v=circshift(v,-size(v,1)/4,1);
%%
nt = size(t,1);
nx = size(x,2);
dx = x(2) - x(1);
dt= t(2)-t(1);

for idb=1:size(n_dB,2)
    if do_noise
        if n_dB(idb)>=1
            % add white noise
            n_lin=1/(10^(n_dB(idb)/20)); %convert input dB to amplitude multiplier
            amp=[-1 1];
        else
            n_lin=1/(10^(n_dB(idb)/20)); %convert input dB to amplitude multiplier
            n_lin=n_lin*(1/n_dB(idb));
            amp=[-1 1];
        end
        rng('default')
        noiseSignal = (amp(2)-amp(1)).*rand(nx,nt) + amp(1);
        v = vsol + (n_lin .* noiseSignal);
        
        %% noise wave video
%         figure(11)
%         for ii=1:size(t)
%             plot(x*10^3,v(:,ii)),hold off
%             xlabel('position (mm)')
%             title(['time = ', num2str(t(ii)*10^3), ' msec'])
%             ylim([-1 1])
%             drawnow,
%             
%         end
    end
    

%%
    %colormap for plotting
    rdbu = flipud(cbrewer('div', 'RdBu', 255));
    blue = flipud(cbrewer('seq', 'Blues', 256));
    
    
    %% Apply filters
    if do_filt
        fprintf(1, 'Temporal Filter.....\n');
        % Initialize filtered surface data
        bwx = [25, 1000];
        bwt = [200, 5000];
        v = sd_filter_surface_wavefield(v, dx, dt, bwx, bwt);
        v=v./(max(max(abs(v))));
        
        %% filtered wave propogation
%         figure(12)
%         for ii=1:size(t)
%             plot(x*10^3,v(:,ii)),hold off
%             xlabel('position (mm)')
%             title(['time = ', num2str(t(ii)*10^3), ' msec'])
%             ylim([-1 1])
%             drawnow,
%             
%         end
    end
    

    %% Calculate Pulse Width
    
    do_fig=0;
    %[pw,pk2,snrtmp]=pulse_width(dx,dt,cs,x,v,do_fig,dirname);
    
    %snr(idb)=snrtmp;
    
    %% display data
    if do_disp
        
        fprintf(1, 'XT/FFT.....\n');
        % 2D FFT
        nfft = 2*4096;
        fmax = 10000;
        kmax = 5000;
        [VdB, f, k] = xttools_power_spectrum(v, dt, dx, nfft);
        fmask = f <= fmax;
        kmask = k <= kmax;
        f = f(fmask);
        k = k(kmask);
        VdB = VdB(kmask, fmask);
        
        % XT Plot
        figure
        set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
        pcolor(t*1000, x*1000, v)
        shading flat
        cb = colorbar;
        colormap(rdbu)
        xlabel('Time (ms)')
        ylabel('Lateral Position (mm)')
        ylabel(cb, 'Surface Vertical Velocity (a.u.)')
        caxis([-1, 1])
        ylim([0, xend*1000])
        if do_figs, print(fullfile(dirname,['xtplot_',num2str(snr(idb),1)]), '-dpng', '-r300')
        end
        %print(fullfile(dirname, 'xtplot'), '-dpng', '-r300')
        
        %
        figure
        set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
        pcolor(f/1000, k/1000, VdB)
        shading flat
        cb = colorbar;
        colormap(blue)
        xlabel('Frequency (kHz)')
        ylabel('Wavenumber (1/mm)')
        ylabel(cb, 'Power Spectrum (dB)')
        caxis([-40, 0])
        if do_figs, print(fullfile(dirname,['xtplot_',num2str(snr(idb),1)]), '-dpng', '-r300')
        end
        %print(fullfile(dirname, 'fkplot'), '-dpng', '-r300')
        
        drawnow
        
        % Show animated plot of wave propagation
        %     decim=2;
        %     figure
        %     for i = 1:decim:(nt-1)
        %         plot((x)*1000, v(:, i))
        %         %xlim([-L/2, L/2]*1000)
        %         ylim([-1.5, 1.5])
        %         xlabel('x (mm)')
        %         ylabel('Vibration Speed (Norm)')
        %         title(['Time = ', num2str(t(i)*1E3), ' ms'])
        %         pause(0.05)
        %     end
        
        
    end
    
    %% Calculate TOF with various techniques
    
    % NCC kernel
    if do_corr
        
        fprintf(1, ['NCC kernel TOF.....\n']);
        
        vi=v(find(x<xendtof),:);
        xi=x(find(x<xendtof));
        
        gpveli = xttools_group_velocity_ncc(x, t, vi, kernel_nx, nfilt);
        
        tmp=find(xi>xstarttof);
        gpvel=mean(gpveli(tmp(1:end-kernel_nx)),'omitnan');
        aa_output(1) = gpvel;
        
        corr(idb)=gpvel;
        
        if do_disp
        % Plot and save XT plot with group velocity displayed
        rdbu = flipud(cbrewer('div', 'RdBu', 255));
        figure
        set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
        subplot(211)
        pcolor(t*1000, x*1000, v)
        hold on
        %plot(tlin*1000, xlin*1000, '--k', 'LineWidth', 1);
        shading flat
        xlabel('Time (ms)')
        ylabel('Lateral Position (mm)')
        caxis([-1 1])
        cb = colorbar;
        colormap(rdbu)
        colormap(cb, rdbu)
        ylabel(cb, 'Surface Vertical Velocity (a.u.)')
        title(sprintf('NCC Group Velocity = %2.3f m/s', gpvel))
        subplot(212)
        plot(xi(1:end-kernel_nx)*1000,gpveli)
        ylim([cs-1 cs+1])
        xlabel('Lateral Position (mm)')
        ylabel('Group Velocity (m/s)')
        if do_figs,print(fullfile(dirname,['corr_',num2str(snr(idb),1)]), '-dpng', '-r300')
        end
        end
    end
    
    
    % Envelope Fit
    if do_env
        
        fprintf(1, ['TOF envelope fit.....\n']);
        
        vi=v(find(x<xendtof),:);
        xi=x(find(x<xendtof));
        tmp=find(xi<xstarttof);
        vi(1:size(tmp,2),:)=[];
        xi(1:size(tmp,2))=[];
        
        [gpvel, conf95, tpeak, best_fit] = xttools_group_velocity_env(xi, t, vi);
        aa_output(2)=gpvel;
        envo(idb)=gpvel;
        if do_disp
        % Plot and save XT plot with group velocity displayed
        rdbu = flipud(cbrewer('div', 'RdBu', 255));
        figure
        set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
        pcolor(t*1000, x*1000, v)
        hold on
        tt = [min(tpeak), max(tpeak)];
        plot(tpeak*1000, xi*1000, '.b', 'MarkerSize', 8, 'LineWidth', 1)
        plot(tt*1000, gpvel*(tt - best_fit(2))*1000, '-k', 'LineWidth', 1)
        plot(tt*1000, (gpvel + conf95(1))*(tt - best_fit(2))*1000, '--g', 'LineWidth', 1)
        plot(tt*1000, (gpvel + conf95(2))*(tt - best_fit(2))*1000, '--g', 'LineWidth', 1)
        %    plot(tlin*1000, xlin*1000, '--k', 'LineWidth', 1);
        shading flat
        xlabel('Time (ms)')
        ylabel('Lateral Position (mm)')
        caxis([-1 1])
        cb = colorbar;
        colormap(rdbu)
        colormap(cb, rdbu)
        ylabel(cb, 'Surface Vertical Velocity (a.u.)')
        title(sprintf('Envelope Group Velocity = %2.3f m/s', gpvel))
        %imgname = fullfile(datadirs{diridx}, sprintf('gpvel_t%i', fidx));
        if do_figs, print(fullfile(dirname,['env_',num2str(snr(idb),1)]), '-dpng', '-r300')
        end
        end
    end
    
    if do_ttp
        
        fprintf(1, ['TOF ttp.....\n']);
        
        vi=v(find(x<xendtof),:);
        xi=x(find(x<xendtof));
        
        
        tmp=find(xi<xstarttof);
        vi(1:size(tmp,2),:)=[];
        xi(1:size(tmp,2))=[];
        
        [gpvel, conf95,tpeak, best_fit,pk] = xttools_group_velocity_ttp(xi, t, vi);
        aa_output(3) = gpvel;
        ttp(idb)=gpvel;
        if do_disp
        % Plot and save XT plot with group velocity displayed
        rdbu = flipud(cbrewer('div', 'RdBu', 255));
        figure
        set(gcf, 'Color', [1, 1, 1], 'PaperUnits', 'centimeters', 'PaperPosition', [0, 0, 8.6, 8.6/1.5]);
        pcolor(t*1000, x*1000, v)
        hold on
        tt = [min(tpeak), max(tpeak)];
        %plot(tpeak*1000, xi*1000, '.b', 'MarkerSize', 8, 'LineWidth', 1)
        plot(tt*1000, gpvel*(tt - best_fit(2))*1000, '-k', 'LineWidth', 1)
        plot(tt*1000, (gpvel + conf95(1))*(tt - best_fit(2))*1000, '--g', 'LineWidth', 1)
        plot(tt*1000, (gpvel + conf95(2))*(tt - best_fit(2))*1000, '--g', 'LineWidth', 1)
        %    plot(tlin*1000, xlin*1000, '--k', 'LineWidth', 1);
        shading flat
        xlabel('Time (ms)')
        ylabel('Lateral Position (mm)')
        caxis([-1 1])
        cb = colorbar;
        colormap(rdbu)
        colormap(cb, rdbu)
        ylabel(cb, 'Surface Vertical Velocity (a.u.)')
        title(sprintf('TTP Group Velocity = %2.3f m/s', gpvel))
        %imgname = fullfile(datadirs{diridx}, sprintf('gpvel_t%i', fidx));
        %print(imgname, '-dpng', '-r300')
        if do_figs,print(fullfile(dirname,['ttp_',num2str(snr(idb),1)]), '-dpng', '-r300')
        end
        end
    end
     
    
    clear v
end