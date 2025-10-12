function I_fft = centered_fft2(I)
I_fft = fftshift(fft2(ifftshift(I)));
end