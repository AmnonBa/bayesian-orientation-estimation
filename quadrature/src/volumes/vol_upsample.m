function Vol = vol_upsample(vol,D)
d = size(vol,1);
scaling = (D^3)/d^3;
start_idx = floor(D/2)-floor(d/2);
insert_idx = (start_idx+1):(start_idx+d);
Vol_fft = zeros(D,D,D);
Vol_fft(insert_idx,insert_idx,insert_idx) = fftshift(fftn(ifftshift(vol)));
Vol = fftshift(ifftn(ifftshift(Vol_fft)));
Vol = scaling*real(Vol);
end