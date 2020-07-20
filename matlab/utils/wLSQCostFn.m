function [ cost, grad, hess ] = wLSQCostFn( wgts, meas, doms, As, En, coef )
% wgts: cell of weights (median filtered measurements)
% meas: cell of measurements
% doms: cell of domains
% As:   cell of precomputed A matrices
% En:   expectation of fitted Rayleigh
% coef: column vector of coefficients

nimg    = length(meas);                 % number of measurements
ncoef   = length(coef);                 % number of coefficients
d       = size(meas{1},1);              % dimension of image (square)
grad    = zeros(ncoef,1);               % initialize gradient
hess    = zeros(length(coef));          % initialize Hessian
cost    = 0;                            % initialize cost

for k = 1 : nimg
    currD   = doms{k}; 
    currM   = meas{k}; 
    currA   = As{k};
    currW   = wgts{k};
    
    Ac      = currA*coef;
    Ac      = Ac(:);
    
    % update cost ---------------------------------------------------------
    diffs   = currM - En*reshape(abs(sin(Ac)),[d,d]);
    tmp     = currD.*( (diffs.^2) ./ currW );
    tmp     = tmp(:);
    cost    = cost + tmp(~isnan(tmp));

    if( nargout > 1 )    
    % update grad -----------------------------------------------------
        t1      = reshape( -2*En*currD.*diffs./(currW.^2), [d^2,1] );
        t1(isnan(t1)) = 0;
        t2      = currD(:).*sign(sin(Ac)).*cos(Ac);
        jac     = zeros(size(currA));
        for j = 1 : size(currA,2)
            jac(:,j) = t2 .* currA(:,j);
        end
        grad    = grad + jac' * t1;
        if( nargout > 2 )
        % update Hessian --------------------------------------------------
            hess = hess + 2*En*currA'* ...
                (repmat((cos(Ac).^2 + diffs(:).*sin(2*Ac).*sign(sin(Ac)))./(currW(:)).^2,[1,ncoef])...
                .*currA);
            
        end     
    end
end