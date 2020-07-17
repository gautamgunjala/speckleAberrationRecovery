function [ f ] = diffCharOpt( x, c, fx2fy2, FP, pupil )
% optimization function to recover defocus along with diffuser statistical
% properties
% x is a vector of coefficients
% c is a vector of constants: 
%   [ (nx^2*ps*sqrt(2*pi)), (-2*pi^2*fc^2*ps^2), (pi*lambda*(1e-6)*fc^2) ]
% fx2fy2 = fx.^2 + fy.^2
% FP is the smoothed periodogram

s1 = exp((c(2))*(x(3))^2*fx2fy2);
s2 = sin(c(3)*x(2)*fx2fy2);
s = pupil .* ((c(1)*(x(1)/x(3)).* s1 .* s2 ).^2 - FP).^2;

f = sum( s(:) );


end
