function [ Ef, Hs ] = propagate( Ein, lambda, z, ps, zp )
% Propagation via angular spectrum method 
% (as given by Eq. 4-21 in Fourier Optics, Goodman)
%
% [ Ef, Hs ] = propagate( Ein, lambda, z, ps )
% 
% Ein       - input electric field
% lambda    - wavelength of illumination
% z         - propagation distance(s)
%             If z is a scalar, Ef will be an image
%             If z is a vector, Ef will be a cell array
% ps        - (effective) pixel size of sensor
% zp
%
% Ef        - output electric field
% H         - propagation kernel function

[m,n]   = size(Ein);
fx      = [ -(n/2) : 1 : (n/2 -1) ] ./(n*ps);
fy      = [ -(m/2) : 1 : (m/2 -1) ] ./(m*ps);
[Fx,Fy] = meshgrid( fx, fy );
Fx      = ifftshift(Fx);
Fy      = ifftshift(Fy);

Lz      = length(z);
if( Lz == 1 )
    Hs      = exp( -1i*pi*lambda*z .*(Fx.^2 + Fy.^2) );
    Ef      = ifft2( fft2(Ein) .* Hs );
    
else
    Ef      = cell(length(z),1);
    Hs      = cell(length(z),1);
    for ii = 1 : Lz
        H       = exp( -1i*pi*lambda*z(ii) .*(Fx.^2 + Fy.^2) );
        Hs{ii}  = H;
        Ef{ii}  = ifft2( fft2(Ein) .* H );
    end
end
end