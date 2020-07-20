function [ M ] = shiftOperator( maxDeg, a, b, x0, y0 )
% M = shiftOperator( maxDeg, a, b, x0, y0 )
%
% Creates a matrix operator, M, which maps coefficients of a polynomial of 
% degree maxDeg to coefficients of the transformed polynomial 
% g(x,y) = f(ax + x0, by + y0) with respect to the standard basis

% Compute size of matrix
N   = 1/2 * (maxDeg + 1) * (maxDeg + 2);
M   = zeros(N);
I   = eye(N);

for i = 1 : N
    M(:,i)  = getCoefs2DShiftedPoly( I(:,i), maxDeg, x0, y0, a, b );
end