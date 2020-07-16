function [ P ] = polyCoeffEval2D( x, n, A )
% x is an array of coefficients
% n is the size of the image (one dimension, assumed square)
% A is an (n x n x k) array of basis matrices
k = length(x);
x = reshape(x,[1,1,k]);
x = repmat(x,[n,n,1]);
P = sum( A.*x, 3 );
end