% Gautam Gunjala, Stuart Sherwin
% 4/14/2020
% Noise analysis of aberration estimation
%
% aberration polynomials defined with magnitude (L2-norm) 
%
%               10^pow * pi radians RMS
%
% output is the minimum converged cost among 100 initializations 
%
%

fs          = filesep;
parentDir   = pwd;
seps        = strfind(parentDir, fs);
rootDir     = parentDir(1:seps(end-1));
dataDir     = [parentDir fs 'simData' fs];
outDir      = [parentDir fs 'results' fs];
addpath([rootDir 'matlab' fs 'utils' fs]);


%% Simulation parameters

N_ds        = 5;                % Number of datasets
N_reps      = 5;                % Number of initializations per dataset
N_iter      = 50;               % Number of iterations per initialization
phpx        = 6600;             % Number of photons per pixel
pow         = -1;               % Aberration magnitude (10^pow)*pi rad RMS
step0       = 2^-14;            % Initial step size
N_backtrack = 10;               % Maximum number of steps in backtracking line search (=1 for no backtracking)

allErrs     = NaN(N_reps,N_ds);
allCosts    = NaN(N_reps,N_ds);


%% Create datasets

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
% Set newData to true to create new datasets and/or replace existing ones
newData     = true; 
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

if(newData)   
    fprintf('Creating datasets in simData ...')
    
    dirName     = ['mag_1e' num2str(pow)];
    
    try
        rmdir([dataDir, dirName], 's');
    catch err
    end
        
    mkdir(dataDir, dirName);

    for ii = 1 : N_ds
        makeDatasetSim( 10^pow, ...
                [dataDir dirName fs 'dataset' num2str(ii) '.mat'], 1);
    end
    
    fprintf(' done.\n')
else
    fprintf('Using pre-existing datasets in simData.\n')
end
   
   
%% Run aberration recovery simulation

t_0    = tic;

for ii = 1 : N_ds
    fprintf('Processing dataset %i of %i ... \n', ii, N_ds)

    ds_name     = ['dataset' num2str(ii) '.mat'];
    ds_path     = [dataDir 'mag_1e' num2str(pow) fs ds_name];
    
    if ii == 1
        % Precompute terms that will be re-used
        [meas, x_true, phiEst, raylMean, A_sp, nFx, nFy, pupDom] ...
                                                    = precompute(ds_path);
        A_sp = A_sp(:,4:end);
    else
        [meas, x_true, phiEst, raylMean]  = precompute(ds_path);
    end
    x_true = x_true(4:end);
    t_start = tic;
    [errs,costs,~]  = aberrationMeasSim( meas, x_true, phiEst, A_sp, ...
                                         pupDom, raylMean, phpx, N_reps, ...
                                         N_iter, step0, N_backtrack );
    tt = toc(t_start);
    
    [~,idx]         = min(costs);
    argminErr       = errs(idx);
    allErrs(:,ii)   = errs;
    allCosts(:,ii)  = costs;
   
    fprintf('done.\n')
    fprintf('Time elapsed: %.2f seconds \n', tt)
    fprintf('Minimum relative error was %.2f percent \n', min(errs))
    fprintf('Relative error of output was %.2f percent \n', argminErr)
    
    savename        = [outDir 'mag_1e' num2str(pow) '_out.mat'];
    save(savename,'allErrs','allCosts','N_ds','N_reps','N_iter','phpx')
    
end

fprintf('\nTotal time elapsed: %.2f seconds \n', toc(t_0))