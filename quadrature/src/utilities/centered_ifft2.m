function I = centered_ifft2(I_fft)
I = fftshift(ifft2(ifftshift(I_fft)));
end