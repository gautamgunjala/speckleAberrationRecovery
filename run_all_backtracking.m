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
out_path = [path0 'outputs\']; % Directory to save outputs
if ~exist(out_path,'dir')
    mkdir(out_path)
end
path1 = [path0 'data\']; % Directory to load data
files = dir(path1);
files(1:2) = [];

%% Process
N_rms       = length(files);  % Number of RMS aberation power levels
N_ds        = 25;               % Number of datasets per aberation strength
N_reps      = 50;               % Number of initializations per dataset
N_iter      = 100;              % Number of iterations per initialization
phpx        = 10^3.8;           % Number of photons per pixel
% pow         = -1;               % Aberration magnitude (10^pow)*pi rad RMS
step0       = 2^-14;             % Initial step size
N_backtrack = 10;                % Maximum number of steps in backtracking line search (=1 for no backtracking)

allErrs     = NaN(N_reps,N_ds);
allCosts    = NaN(N_reps,N_ds);

%% Run
out_path1 = [out_path 'N_iter'  num2str(N_iter) 'N_reps' num2str(N_reps) 'step' num2str(step0) 'N_backtrack' num2str(N_backtrack) 'phpx' num2str(log10(phpx)) '\'];
if ~exist(out_path1,'dir')
    mkdir(out_path1)
end
t_start0 = tic;
for jj = 1 : N_rms
    t_start1 = tic; 
    fn = files(jj).name;
%     if ~strcmp(fn,'1e-1')
%         continue
%     end
    for ii = 1 : N_ds
        fprintf([fn '\n'])
        fprintf('Processing dataset (%i, %i) of (%i, %i) ... \n', jj, ii, N_rms, N_ds)

        ds_name         = ['dataset' num2str(ii) '.mat'];
        ds_path         = [path1 fn '/' ds_name];
%         ds_path         = [path1 '1e' num2str(pow) '/' ds_name];
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

        savename        = [out_path1 'mag_' fn '_out' '.mat'];
        save(savename,'allErrs','allCosts','N_ds','N_reps','N_iter','phpx','step0','N_backtrack')

    end
    fprintf(['\n' 'Data set ' num2str(jj) ' complete\n'])
    fprintf('time elapsed: %.2f seconds \n', toc(t_start1))
end
fprintf('\nTotal time elapsed: %.2f seconds \n', toc(t_start0))