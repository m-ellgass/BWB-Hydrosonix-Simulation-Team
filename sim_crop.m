function [x2,t,vsol] = sim_crop(x,t,vsol,xend,tend,nxin)
%Take the 1D simulation results, crop according to space and time field of
%interest

t(end)=[];

%concatenate to match simulation/experiment wavefields
vsol=circshift(vsol,-size(vsol,1)/4,1);
vsol(find(x>xend),:)=[];
x(find(x>xend))=[];

ratio = nxin/size(x,2);

x2=resample(x,nxin,size(x,2));
v2=resample(vsol,nxin,size(x,2));

vnormfactor = max(max(v2(1:2:end,:)));
vsol=v2./(max(abs(v2(:))));

end

