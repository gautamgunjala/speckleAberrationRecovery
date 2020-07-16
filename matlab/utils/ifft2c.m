function res = ifft2c(F)
    res = ifftshift(ifft2(fftshift(F)));
end