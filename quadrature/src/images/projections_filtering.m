function projections = projections_filtering(projections,sigma)

n=sqrt(numel(projections{1}));

if mod(n,2)==0
    k = 2*pi*((-n/2):(n/2-1))/n;  
else
    k = 2*pi*((-(n-1)/2):((n-1)/2))/n; 
end
[kx,ky] = meshgrid(k);

F = exp(-sigma^2*(kx.^2+ky.^2)/2);
for i=1:numel(projections)
    I=reshape(projections{i},[n,n]);
    I_fft = fftshift(fft2(ifftshift(I)));
    I_fft_filtered = I_fft.*F;
    I_filtered = fftshift(ifft2(ifftshift(I_fft_filtered)));
    projections{i} = real(I_filtered(:));
end

end