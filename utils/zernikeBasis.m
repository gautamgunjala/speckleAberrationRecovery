function [ B , Z ] = zernikeBasis( deg, X, Y )
% Outputs a 3D array of Zernike basis functions, B, up to order deg, 
% evaluated on X and Y (outputs of meshgrid) and a matrix for conversion
% from Zernike polynomials to standard basis polynomials, Z
%
% If z is a vector of Zernike coefficients, and s is a vector of standard
% basis coefficients, then Zz = s and Z\s = z
% 
% Zernike coefficients are ordered by OSA/ANSI index
% Standard basis coefficients are ordered by degree, priority given to x
%       e.g. { 1, x, y, x^2, x*y, y^2, x^3, x^2*y, x*y^2, y^3, ... }
% 
% Radial polynomial formula:
%       http://mathworld.wolfram.com/ZernikePolynomial.html
% Normalization constant:
%       http://www.opt.indiana.edu/vsg/library/vsia/vsia-2000_taskforce/tops4_2.html
% 

[m,n]   = size(X);

% Create m's and n's from Noll indices
npol    = 1/2 * (deg + 1) * (deg + 2);
nn      = zeros(1,npol);
mm      = zeros(1,npol);
pp      = zeros(1,npol);
for ii = 0 : deg
    nn((ii*(ii+1)/2:(ii*(ii+3))/2)+1) = ii;
    mm((ii*(ii+1)/2:(ii*(ii+3))/2)+1) = -ii:2:ii;
    pp((ii*(ii+1)/2:(ii*(ii+3))/2)+1) = ii:-1:0;
end

R       = sqrt(X.^2 + Y.^2);
T       = atan2d(Y,X);

B       = zeros(m,n,npol);
ZE      = zeros(m*n,npol);
SE      = zeros(m*n,npol);

for ii = 1 : npol
    ni      = nn(ii);
    mi      = abs(mm(ii));
    pi      = pp(ii);
    N       = sqrt((2*(ni+1))/(1+1*(mi==0)));
    
    tmp     = zeros(m,n);
    for k = 0 : (ni-mi)/2
        tmp     = tmp + (-1).^k .*nchoosek(ni-k,k) ...
                  .*nchoosek(ni-2*k,((ni-mi)/2)-k) .*(R.^(ni-2*k));
    end
    
    if(mm(ii) > 0)
        tmp     = tmp .* cosd(mi*T);
    elseif(mm(ii) < 0)
        tmp     = tmp .* sind(mi*T);
    end
    
    B(:,:,ii)   = N*tmp;
    ZE(:,ii)    = N*tmp(:);
    tmp2        = X.^pi .* Y.^(ni-pi);
    SE(:,ii)    = tmp2(:);    
end

Z       = SE\ZE;


    




