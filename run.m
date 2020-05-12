% Gautam Gunjala, Stuart Sherwin
% 4/14/2020
% Noise analysis of aberration estimation
%
% aberration polynomials defined with magnitude (L2-norm) 
%
%               10^pow * pi radians RMS
%
% output is the minimum converged cost of 100 initializations 
%
%

addpath('utils')

if false % path0: rooth path where \data\ is located, and where \outputs\ will be saved
    path0 = 'C:\Users\stuar\Documents\Data Files\Simulated Aberrations Data_2020\';
else
    path0 = pwd; % Save in current directory
    if ~ismember(path0(end),{'\','/'})
        path0 = [path0 '\'];
    end
end
%% Process
N_ds        = 1;%25;                % Number of datasets
N_reps      = 1;%50;               % Number of initializations per dataset
N_iter      = 100;               % Number of iterations per initialization
phpx        = 10^3.8;           % Number of photons per pixel
pow         = -1;               % Aberration magnitude (10^pow)*pi rad RMS
step0       = 2^-14;             % Initial step size
N_backtrack = 10;                % Maximum number of steps in backtracking line search (=1 for no backtracking)

allErrs     = NaN(N_reps,N_ds);
allCosts    = NaN(N_reps,N_ds);

%% Run
t_start0 = tic;
for ii = 1 : N_ds
    fprintf('Processing dataset %i of %i ... \n', ii, N_ds)

    ds_name         = ['dataset' num2str(ii) '.mat'];
    ds_path         = [path0 'data/1e' num2str(pow) '/' ds_name];
    if ii == 1
        %% Precompute terms that will be re-used
        [meas, x_true, phiEst, raylMean, A_sp, nFx, nFy, pupDom]  = precompute(ds_path);
        A_sp = A_sp(:,4:end);
    else
        [meas, x_true, phiEst, raylMean]  = precompute(ds_path);
    end
    x_true = x_true(4:end);
    t_start = tic;
    [errs,costs,~]  = aberrationMeasSim( meas, x_true, phiEst, A_sp, pupDom, raylMean, phpx, N_reps, N_iter, step0, N_backtrack );
    tt = toc(t_start);
    
    [~,idx]         = min(costs);
    argminErr       = errs(idx);
    allErrs(:,ii)   = errs;
    allCosts(:,ii)  = costs;
   
    fprintf('done.\n')
    fprintf('Time elapsed: %.2f seconds \n', tt)
    fprintf('Minimum relative error was %.2f percent \n', min(errs))
    fprintf('Relative error of output was %.2f percent \n', argminErr)
    
    savename        = [out_path 'mag_1e' num2str(pow) '_out'];
    save(savename,'allErrs','allCosts','N_ds','N_reps','N_iter','phpx')
    
end

fprintf('\nTotal time elapsed: %.2f seconds \n', toc(t_start0))