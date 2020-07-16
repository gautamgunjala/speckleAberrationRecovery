function [ s ] = getCoefs2DShiftedPoly( c, maxDeg, dx, dy, sx, sy )
% getCoefs2DShiftedPoly
% coefMat: 3 cols [coefs, indX, indY]
% dx, dy: shifts in x, y
% sx, sy: scale in x, y
s = zeros(size(c));
N = maxDeg;

for i = 0 : N
    for j = 0:(N-i)
        idx = (i+j)*(i+j+1)/2 + j + 1; % Formula: N(N+1)/2 + j + 1, N = i+j (total order)
        coef = c( idx );
        for n = 0 : i
            for m = 0 : j
                idxO = (n+m)*(n+m+1)/2 + m + 1;
                s(idxO) = s(idxO) + ...
                    coef*nchoosek(i,n)*sx^n*dx^(i-n) * ... 
                    nchoosek(j,m)*sy^m*dy^(j-m);
            end
        end
    end
end
end