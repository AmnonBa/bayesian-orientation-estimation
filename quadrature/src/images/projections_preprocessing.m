function [projections,Filt,ds_idx] = projections_preprocessing(projections,d,sigma)
D=size(projections{1},1);
ds_start=floor(D/2)-floor(d/2);
ds_idx=(ds_start+1):(ds_start+d);

if mod(d,2)==0
    k = 2*pi*((-d/2):(d/2-1))/d;  
else
    k = 2*pi*((-(d-1)/2):((d-1)/2))/d; 
end
[kx,ky] = meshgrid(k);
Filt = exp(-sigma^2*(kx.^2+ky.^2)/2);

N=numel(projections);
scaling=d^2/D^2;
for i=1:N
    I=projections{i};
    I_fft=centered_fft2(I);
    I_fft=I_fft(ds_idx,ds_idx);
    I_fft=I_fft.*Filt;
    projections{i}=scaling*centered_ifft2(I_fft);
end

end