function [ err_rep, cost_rep, w_fit_rep ] = aberrationMeasSim( meas, x_true, phiEst, A_sp, pupDom, raylMean, phpx, N_reps, N_iter, step0, N_backtrack )
% Inputs:
%   dataname:       name of dataset
%   phpx:           number of photons per pixel
%   N_reps:         number of initializations per dataset
%   N_iter:         number of iterations per initialization

n_img           = size(meas,3);

%% Add photon noise to measurements
N_ph    = phpx*size(meas,1)*size(meas,2);
for k = 1 : n_img
    I       = meas(:,:,k);
    tmp = sum(I(:));
    I = I./tmp;
%     I = I./sum(I(:));
    I_ph    = I*N_ph;
    I_n     = (I_ph + sqrt(I_ph).*randn(size(I_ph)));
    meas(:,:,k) = I_n / N_ph * tmp;
%     meas(:,:,k) = I_n ./ max(I_n(:));
end
%% 
measCell = cell(n_img,1);
is_meas = full(any(A_sp,2));
for i = 1 : n_img % Loop through angles of illumination
    Iu = abs(fft2c(meas(:,:,i)))./(2*phiEst);
    measCell{i} = Iu(pupDom(:,:,i)~=0);
end
% Data vector
y = cell2mat(cellfun(@(x) vec(x),measCell,'UniformOutput',false));
%% Gradient descent
const       = norm(vec(y));
cost_rep    = NaN(N_reps,1);
err_rep     = NaN(N_reps,1);

w_fit_rep   = zeros(length(x_true),N_reps);

test = false;
if test
    cc = NaN(N_iter,N_reps);
end
for iR = 1:N_reps
    w0 = randn(size(x_true));
    w0 = w0/norm(w0)*norm(x_true);
    w_fit = w0;
    
    if test
%         w_fit = x_true; % Perfect init
%         y = raylMean*abs(sin(A_sp(is_meas,:)*x_true)); % Perfect data
    end
    
    
    Ac      = A_sp(is_meas,:)*w_fit;
    curMinCost = mean( (y - raylMean*abs(sin(Ac))).^2 ); % Evaluate new cost
    
    for i = 1:N_iter
        step = step0; % Initial step size
        % Evaluate gradient
        Ac      = A_sp(is_meas,:)*w_fit;
        sAc     = sin(Ac);
        cAc     = cos(Ac);
        res     = y - raylMean*abs(sAc);
        t1      = -2*raylMean.*res;
        t1 = t1.*sign(sAc).*cAc;
        grad = A_sp(is_meas,:)'*t1;
        for k = 1:N_backtrack % Backtracking line-search using gradient descent
            
            w_fit   = w_fit - step*grad;

            Ac1 = A_sp(is_meas,:)*w_fit; % Evaluate new solution
            curCost = mean( (y - raylMean*abs(sin(Ac1))).^2 ); % Evaluate new cost
            if or(curCost < curMinCost,k == N_backtrack) % Accept or reject
                curMinCost = curCost;
%                 Ac = Ac1; % Update Ac 
                break
            else
                w_fit   = w_fit + step*grad; % Reject and shrink step size
                step    = step*0.5; % Shrink step
            end
        end
        if test
            cc(i,iR) = curMinCost;
            plot(cc(:,1:iR))
            xlim([1 N_iter])
            title([num2str(iR) ' ' num2str(i) ' ' num2str(k) ', ' num2str(min(norm(w_fit + x_true),norm(w_fit - x_true))/norm(x_true)*100) ])
            drawnow()
        end
    end
    
    thisCost    = norm(vec( y - raylMean*(abs(sin(A_sp(is_meas,:)*w_fit))) ))/const;
    if norm(w_fit + x_true) < norm(w_fit - x_true)
        w_fit = -w_fit;
    end
    thisErr     = norm(w_fit - x_true)/norm(x_true)*100;
    
    
    w_fit_rep(:,iR) = w_fit;
    cost_rep(iR)    = thisCost;
    err_rep(iR)     = thisErr;
        
    fprintf(    'Init: %i \t Rel. err: %.2f percent ...\n', ...
            iR, thisErr )

end