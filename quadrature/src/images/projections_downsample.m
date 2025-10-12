function projections = projections_downsample(projections, d)
n=sqrt(numel(projections{1}));
start_idx = floor(n/2)-floor(d/2);
slice_idx = (start_idx+1):(start_idx+d);
for i=1:numel(projections)
    I=reshape(projections{i},[n,n]);
    I_fft = fftshift(fft2(ifftshift(I)));
    I_fft_ds = I_fft(slice_idx,slice_idx);
    I_ds = fftshift(ifft2(ifftshift(I_fft_ds)));
%     norm(imag(I_ds),'fro')/norm(I_ds,'fro')
    I_ds = real(I_ds)*(d^2/n^2);
    projections{i} = reshape(I_ds(:),[d,d]);
end
end 


% function I_ds = img_downsample(I,n,ds_size)
% I = reshape(I,[n,n]);
% start_idx = floor(n/2)-floor(ds_size/2);
% slice_idx = (start_idx+1):(start_idx+ds_size);
% I_fft = fftshift(fft2(ifftshift(I)));
% I_fft_ds = I_fft(slice_idx,slice_idx);
% I_ds = fftshift(ifft2(ifftshift(I_fft_ds)));
% I_ds = real(I_ds)*(ds_size^2/n^2);
% I_ds = I_ds(:);
% end