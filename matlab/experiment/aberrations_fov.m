% Gautam Gunjala
% Updated 7/16/2020
% Script for aberration estimation across FOV in experimental data (SHARP)
%
% Output is a collection of Zernike coefficient vectors corresponding to
% every analyzed sub-region of the full-field image

fs          = filesep;
parentDir   = pwd;
seps        = strfind(parentDir, fs);
rootDir     = parentDir(1:seps(end-1));
dataDir     = [rootDir 'data' fs];
addpath([rootDir 'matlab' fs 'utils' fs]);


%% Initialize imaging parameters and coordinate spaces

param   = struct(   'nx',           256, ...        % Pixels in sub-image
                    'ps_m',         15e-9, ...      % Pixel size [m]
                    'lambda_m',     13.5e-9, ...    % Wavelength [m]
                    'NA',           0.33/4 );       % Numerical aperture 
                    
% real & frequency space axes
dfx_mi          = 1/(param.nx *param.ps_m);
[X_m,Y_m]       = meshgrid((-(param.nx/2):1:(param.nx/2-1)) *param.ps_m);
fx_mi           = (-(param.nx/2):1:(param.nx/2-1)) *dfx_mi;
[Fx_mi,Fy_mi]   = meshgrid(fx_mi);

% NA of imaging system / pupil window with cutoff frequency NA/lambda
fc              = param.NA / param.lambda_m;
pupil           = 1*( Fx_mi.^2 + Fy_mi.^2 < fc^2 );
nanPup          = pupil; 
nanPup(nanPup == 0) = NaN;

% Normalize freq coordinates for numerical stability (Fx = fc * nFx)
nfx             = fx_mi ./fc;
ndfx            = dfx_mi ./fc;
[nFx,nFy]       = meshgrid( nfx, nfx );


%% Speckle parameter calibration

% Load calibration (strongly defocused) image and metadata
calImg      = double(imread([dataDir 'SHARP-00.png']));
calDef      = -2.40859;
xroi        = 906:1161;
yroi        = 679:934;

center      = floor(param.nx/2)+1;
radP        = 1 ./ ndfx;
crop1       = center - ceil(radP) - 10;
crop2       = center + ceil(radP) + 9;
pupil_c     = pupil(crop1:crop2, crop1:crop2);
nFx_c       = nFx(crop1:crop2, crop1:crop2);
nFy_c       = nFy(crop1:crop2, crop1:crop2);

imSeg       = calImg(yroi, xroi);
imSegFT     = rmvAxes(abs(fft2c(imSeg))) .*pupil;
imSegFT_c   = imSegFT(crop1:crop2, crop1:crop2);

[sigmaRayl, window, ~]  ...
            = calibrateSpeckle(imSegFT_c, param,[10000; calDef; 2], 0);

        
%% Zernike Polynomial basis
[B,Z]   = zernikeBasis( 5, nFx_c, nFy_c );

%% Load Angle data (full)
series = 12;
[img_full, meta_fp] = Sharp.read(folder,mask_name,series);
[imgH,imgW]   = size(img_full{1});
[~,xroi,yroi]       = Sharp.ROI(img_full, nx, 0, 0);

%% Initialize data storage

FOV_ABERR   = zeros( size(Xs,2), size(Xs,1), 21 );
[sz1,sz2,~] = size(FOV_ABERR);
spectraAll  = cell(sz1,sz2);
rawSpecAll  = cell(sz1,sz2);
raylMeanAll = zeros(sz1,sz2);
pupDomAll   = cell(sz1,sz2);
fmaAll      = cell(sz1,sz2);

%% Gather useful images from dataset

idx     = 1;
U       = 0;
V       = 0;

for i = 1 : size(img_full,1)
    if( meta_fp{i}.ma_arg0 < 0.3 && meta_fp{i}.ma_arg0 ~= 0 )
        idx = [ idx; i ];  
        U   = [ U; meta_fp{i}.ma_arg0.*cosd(meta_fp{i}.ma_arg1) ];
        V   = [ V; meta_fp{i}.ma_arg0.*sind(meta_fp{i}.ma_arg1) ];
    end
end

n_img   = size(idx,1);
img_full= img_full(idx);


%% Main loop

for p = 1 : size(Xs,1)
    
    %h2 = figure;
    
    for q = 1 : size(Xs,2)
        
        if( segments(p,q,3) )
        fprintf('Row: %d \t Column: %d \n',p,q)
        
        %% Get crop of data, and compute DC-suppressed spectra
        img_crop    = Sharp.ROI(img_full, nx, Xs(p,q), Ys(p,q));
        %img_crop    = (img_crop ./ mean(img_crop(:))) .* meanInten;
        spectra     = cell(size(img_crop));

        for i = 1 : size(img_crop,1)
            spectra{i}  = Sharp.remove_dc(abs(fft2c(img_crop{i}./5)));
        end
        
        %% Estimate illumination angles
        Ue      = zeros(size(idx));
        Ve      = zeros(size(idx));

        for i = 1 : size(spectra,1)
            [Hp,Vp]     = findCirc2( spectra{i}, radP );
            freqXYout   = (estIllumAngle(spectra{i}, radP, ...
                           [center+Hp,center+Vp], 0) - center)*ndfx;
            Ue(i)       = freqXYout(1);
            Ve(i)       = freqXYout(2);
        end

        %% Window measurement spectra, process, and crop to final size
        pupDom      = cell(n_img,1);
        spectra2    = spectra;

        for i = 1 : n_img
            mm1         = ( (nFx + Ue(i)).^2 + (nFy + Ve(i)).^2 ) <= 1;
            mm2         = ( (-nFx + Ue(i)).^2 + (-nFy + Ve(i)).^2 ) <= 1;
%             mm1         = ( (nFx + U(i)).^2 + (nFy + V(i)).^2 ) <= 1;
%             mm2         = ( (-nFx + U(i)).^2 + (-nFy + V(i)).^2 ) <= 1;
            pupDom{i}   = double( and( mm1, mm2 ) );

            tmp         = spectra{i} .* pupDom{i};
            spectra2{i} = tmp(crop1:crop2,crop1:crop2);
            spectra2{i} = abs( spectra2{i} ./ (2 * window) );

            pupDom{i}   = pupDom{i}(crop1:crop2,crop1:crop2);
        end

        %% Define evaluation map
        E = diag([1 1 1 0 0 0 0 1 1 1 1 1 0 0 0 0 0 0]);

        n           = size(nFx_c,1);      
        Psi         = zeros(n^2,18);
        Psi(:,1)    = reshape(nFx_c.^2,            [n^2,1]);
        Psi(:,2)    = reshape(nFx_c.*nFy_c,         [n^2,1]);
        Psi(:,3)    = reshape(nFy_c.^2,            [n^2,1]);
        Psi(:,4)    = reshape(nFx_c.^3,            [n^2,1]);
        Psi(:,5)    = reshape(nFx_c.^2.*nFy_c,      [n^2,1]);
        Psi(:,6)    = reshape(nFx_c.*nFy_c.^2,      [n^2,1]);
        Psi(:,7)    = reshape(nFy_c.^3,            [n^2,1]);
        Psi(:,8)    = reshape(nFx_c.^4,            [n^2,1]);
        Psi(:,9)    = reshape(nFx_c.^3.*nFy_c,      [n^2,1]);
        Psi(:,10)   = reshape(nFx_c.^2.*nFy_c.^2,   [n^2,1]);
        Psi(:,11)   = reshape(nFx_c.*nFy_c.^3,      [n^2,1]);
        Psi(:,12)   = reshape(nFy_c.^4,            [n^2,1]);
        Psi(:,13)   = reshape(nFx_c.^5,            [n^2,1]);
        Psi(:,14)   = reshape(nFx_c.^4.*nFy_c,      [n^2,1]);
        Psi(:,15)   = reshape(nFx_c.^3.*nFy_c.^2,   [n^2,1]);
        Psi(:,16)   = reshape(nFx_c.^2.*nFy_c.^3,   [n^2,1]);
        Psi(:,17)   = reshape(nFx_c.*nFy_c.^4,      [n^2,1]);
        Psi(:,18)   = reshape(nFy_c.^5,            [n^2,1]);

        %% Define weights (if doing weighted least squares)
        Wts     = ones(n,n,n_img);
        Wts     = squeeze(num2cell(Wts,[1,2]));

        %% Define shift operators
        S = cell(n_img,1);
        for k = 1 : n_img
%             tmp     = shiftOperator(5,1,1,U(k),V(k));
            tmp     = shiftOperator(5,1,1,Ue(k),Ve(k));
            S{k}    = tmp(4:end,4:end);
        end
        clear tmp

        %% Premultiply all matrices
        A = cell(n_img,1);
        for k = 1 : n_img
            A{k} = Psi*E*S{k};
        end

        %% Number of images to use
        nImgAlg     = 9;
        spectra     = spectra(1:nImgAlg);
        spectra2    = spectra2(1:nImgAlg);
        pupDom2     = pupDom(1:nImgAlg);
        for t = 1 : size(pupDom2,1)
           pupDom2{t}    = Sharp.remove_dc(pupDom2{t});    
        end

        %% Set up forward model
        raylMean    = sigmaRayl .* sqrt(pi/2);    
        FM          = @(w,k) pupDom2{k} .* (raylMean*reshape(abs(sin(A{k}*w)),[n,n]));
        
        %% Save everything
        rawSpecAll{p,q}     = spectra;
        spectraAll{p,q}     = spectra2;
        raylMeanAll(p,q)    = raylMean;
        pupDomAll{p,q}      = pupDom2;
        fmaAll{p,q}         = A;
        
        %% Run multiple initializations
        minCost = Inf;
        argmin  = 0;
        
        nVals       = nImgAlg*n^2;
        costEval    = @(w) wLSQCostFn(Wts,spectra2,pupDom2,A,raylMean,w);
        
        %% Newton's Method--------------------------------------------------------
        mults       = [0.2 0.45 0.7 0.95 1.2];
        
        for trial_group = 1 : 5
            
            mult = mults(trial_group);
        
            for trials = 1 : 25

                w_init      = mult*randn([18,1]);
                w0          = w_init;
                step        = 0.05;
                nitr        = 200;

                for i = 1 : nitr

                    [cost,gradW,hessW] = costEval(w0);
                    %fprintf('Iteration: %d \t Cost: %d \n',i-1,cost/nVals)
                    t   = 1/2*norm(gradW)^2;
                    k   = 0;

                    while( cost - costEval(w0 - (step^k)*(gradW)) < (step^k)*t && k < 12 )
                        k = k + 1;
                    end

                    w0  = w0 - (step^k)*(gradW);    
                end
                
                fprintf('Iteration: %d \t Optimal Cost: %d \n',nitr,wLSQCostFn(Wts,spectra2,pupDom2,A,raylMean,w0)/nVals)
                thisCost    = wLSQCostFn(Wts,spectra2,pupDom2,A,raylMean,w0)/nVals;

                if( thisCost < minCost)
                    argmin  = w0;
                    minCost = thisCost;
                end
                
            end
        end
        %% Plot
        z = Z\[0;0;0;argmin]; z(1:3) = 0;
        if z(5) < 0
            z   = -z;
            argmin  = -argmin; 
        end
        
        FOV_ABERR(p,q,:)    = z;
%         subplot(1,length(xs),q), ...
%                 imagesc(nfx(scrop1:scrop2),nfx(scrop1:scrop2), ...
%                 nanPupCrop.*reshape(Psi*argmin,[n,n]));
%         axis xy image
%         drawnow
        end   
    end
    
end

save('aberrations_full_FOV_angleCorrected_512','FOV_ABERR','rawSpecAll','spectraAll','raylMeanAll','pupDomAll','fmaAll')

