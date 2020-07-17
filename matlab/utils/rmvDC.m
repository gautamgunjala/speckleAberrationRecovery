function [ res, DC ] = rmvDC( res )
%
% [ res, DC ] = rmvDC( res )
% Separates Fourier transform 'res' into its DC component 'DC' and
% remaining data. Useful when DC term is too large and obscures information
% contained in the rest of the image.
%
    [a,b] = size(res);
    DC = res(floor(a/2)+1,floor(b/2)+1);
    res(floor(a/2)+1,floor(b/2)+1) = 0;
end