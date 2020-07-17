function [ res ] = rmvAxes( res )
%
% [ res, DC ] = rmvAxes( res )
% Deletes axis values from Fourier transform 'res' 
%
    [a,b]   = size(res);
    r       = floor(a/2)+1;
    c       = floor(b/2)+1;
    res(r,:) = 0;
    res(:,c) = 0;
end