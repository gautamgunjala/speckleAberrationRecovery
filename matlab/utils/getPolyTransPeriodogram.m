function [A,A_plus,A_minus] = getPolyTransPeriodogram( maxDeg, U, V, T )
    A = zeros((maxDeg+1)*(maxDeg+2)/2);
    if nargout > 1
        A_plus = zeros((maxDeg+1)*(maxDeg+2)/2);
        A_minus = zeros((maxDeg+1)*(maxDeg+2)/2);
    end
    for i = 1:size(A,1)
        x = zeros(size(A,1),1);
        x(i) = 1;
        cPlus = getCoefs2DShiftedPoly( x, maxDeg, U, V, 1, 1 );
        cMinus = getCoefs2DShiftedPoly( x, maxDeg, U, V, -1, -1 );
        cEven = (cPlus + cMinus)/2;
        cEven(1) = 0;
        
        A(:,i) = cEven;
        
        if nargout > 1
            cPlus(1) = 0;
            cMinus(1) = 0;
            A_plus(:,i) = cPlus;
            A_minus(:,i) = cMinus;
        end
    end
    
    A = T\A*T;
    if nargout > 1
        A_plus = T\A_plus*T;
        A_minus = T\A_minus*T;
    end
end