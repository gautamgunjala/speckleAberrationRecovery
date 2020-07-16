function res = fft2c(f)
    res = ifftshift(fft2(fftshift(f)));
end
