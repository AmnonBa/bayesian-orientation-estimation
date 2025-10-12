function vol_ds = vol_downsample(vol,ds_size)
n = size(vol,1);
scaling = (ds_size^3)/n^3;
start_idx = floor(n/2)-floor(ds_size/2);
slice_idx = (start_idx+1):(start_idx+ds_size);
vol_fft = fftshift(fftn(ifftshift(vol)));
vol_fft_ds = vol_fft(slice_idx,slice_idx,slice_idx);
vol_ds = fftshift(ifftn(ifftshift(vol_fft_ds)));
vol_ds = scaling*real(vol_ds);
end