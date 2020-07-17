function [ ] = makeDatasetSim( magnitude, savename, plots )

% inputs:
%   - magnitude ( __*pi radians of RMS wavefront )
%   - string savename
%   - plots (1 for on, 0 for off)

% Gautam Gunjala
% UC Berkeley Computational Imaging Lab
% 02/2020


%% Set up experiment parameters ===========================================

% wavelength of laser
lambda_m    = 13.5e-9;

% sensor dimensions
nx          = 256;
ps_m        = 15e-9;

% real space axes
x_m         = ( -(nx/2) : 1 : (nx/2 -1) ) * ps_m;
[X,Y]       = meshgrid( x_m, x_m );

% frequency space axes
dfx         = 1/(nx*ps_m);
fx          = ( -(nx/2) : 1 : (nx/2 -1) ) * dfx;
[Fx,Fy]     = meshgrid( fx, fx );

% NA of imaging system / pupil window with cutoff frequency NA/lambda
NA          = 0.15; %0.33/4;
fc          = NA / lambda_m;
Pupil       = ( Fx.^2 + Fy.^2 < fc^2 );

% Normalize freq coordinates for numerical stability (Fx = fc * nFx)
% Pupil occupies a radius of 1 in (fx,fy)-space
nfx         = fx ./fc;
[nFx,nFy]   = meshgrid( nfx, nfx );

%% Create weak phase diffuser object ======================================

% adjust weakness of phase
sigma   = 1*ps_m;
scale   = 5e-17;

% create Gaussian
a       = 0.5*(1/sigma^2);
gwin_real = scale*(1/(2*pi*sigma^2))*exp(-a*(X.^2 + Y.^2));
gwin    = scale*exp(-2*pi^2*sigma^2*(Fx.^2+Fy.^2));

% convolve Gaussian with random field to create smooth phase object
temp    = randn(nx);
phix    = ifft2c(fft2c(temp).*fft2c(gwin_real));

% create object transmission function 
E0x     = exp(1i.*phix);
E0u     = fft2c( E0x );

%% Characterize diffuser (Gaussian envelope) ==============================

% Acquire ONE defocused measurement
z       = 1e-6;
Ei      = ifft2c(E0u .*Pupil);
[E1x,H] = propagate(Ei, lambda_m, z, ps_m);
H       = fftshift(H);
Ix      = abs(E1x).^2;

% Compute and process periodogram of intensity measurement 
P       = rmvDC(abs(fft2c(Ix)));
P       = abs(P).^2;
%FP = imgaussfilt(P);
FP      = P;
% Precompute matrices/constants for speed
fx2fy2  = nFx.^2 + nFy.^2;
c       = [ nx^2*ps_m*sqrt(2*pi), ... 
            -2*pi^2*fc^2*ps_m^2, ...
            pi*lambda_m*(1e-6)*fc^2 ];
fxp     = fx((nx/2+1):end)./fc;

% Set up optimization cost function
fun     = @(x) diffCharOpt( x, c, fx2fy2, FP, Pupil );
GUESS   = [1000;1;1];
    
% Run unconstrained quasi-newton algorithm
options = optimset('MaxFunEvals',20000,'MaxIter',2000);
[sol,~] = fminsearch(fun,GUESS,options);

% Check accuracy of recovery
iavg    = radial_avg(FP,nx/2);

%% Diffuser calibration plots

tmp     = c(1) * (sol(1)/sol(3)).*exp(c(2)*(sol(3))^2*fx2fy2);

if(plots)
    
    figure; 
    plot( fxp, sqrt(iavg),'-k','LineWidth',6 )
    hold on

    plot( fxp, abs(tmp((nx/2+1),(nx/2+1):end) .* ... 
          sqrt(sin(c(3)*sol(2)*(fxp.^2)).^2)) ,'-r','LineWidth',2 )
    plot( fxp, abs(tmp((nx/2+1),(nx/2+1):end)),'-b','LineWidth',2 )
    hold off
    title('Simultaneous estimation of window and defocus (chirp')
    ylabel('Magnitude'); 
    xlabel('Normalized frequency');

    legend( 'Measurement radial average', ...
            'Fitted model (damped defocus kernel)', ...
            'Extracted window function')
    grid on
    
end

%% Check Gaussian window estimation + save

gwinEst     = tmp*ps_m^2/(2*nx);
sigmaEst    = sol(3)*ps_m;
scaleEst    = gwinEst((nx/2+1),(nx/2+1));

Pest        = abs(gwinEst.*nx./ps_m^2).^2 .* ...
              (2*sin(c(3)*sol(2)*(fx2fy2))).^2;
PAest       = abs(gwin.*nx./ps_m^2 .* (H - conj(H))).^2;

if(plots)

    % Compare diffuser characterization with true diffuser properties
    figure;
    subplot(1,2,1), imagesc(nfx, nfx, gwin); colorbar
    title(sprintf('True spectral window, sigma = %d, scale = %d', ...
                                                    sigma, scale));
    xlabel('Normalized frequency');
    
    subplot(1,2,2), imagesc(nfx, nfx, gwinEst); colorbar
    title(sprintf('Est. spectral window, sigma = %d, scale = %d', ...
                                                    sigmaEst, scaleEst));
    xlabel('Normalized frequency');

    % Periodogram estimation and comparisons
    figure;
    
    subplot(1,3,1), imagesc(nfx, nfx, P); colorbar
    title('Measured periodogram');
    xlabel('Normalized frequency');
    
    subplot(1,3,2), imagesc(nfx, nfx, PAest); colorbar
    title('Est. periodogram using analytic components');
    xlabel('Normalized frequency');
    
    subplot(1,3,3), imagesc(nfx, nfx, Pest); colorbar
    title('Est. periodogram using estimated components');
    xlabel('Normalized frequency');
    
end


%% Estimate sigma (Rayleigh parameter)
r           = sqrt(rmvDC(P))./sqrt(Pest);
use         = and(Pest/max(Pest(:)) > 0.1,nFx.^2 + nFy.^2 <= 1);
rr          = r(use(:));
K           = length(rr);
sigmaRayl   = sqrt( sum(rr.^2/(2*K)) )*exp(1)*sqrt(K/(K-1))*((K-1)/K)^K;

[N,E]       = histcounts(rr,200);
X           = 0.5*(E(1:end-1) + E(2:end));
if(~plots)
    
    close
    
end

y           = X ./(sigmaRayl^2) .*exp(-X.^2 ./(2*sigmaRayl^2));

raylMean    = sigmaRayl*sqrt(pi/2);

if(plots)
    
    figure; 
    bar(X,N);
    hold on
    plot(X, y*sum(N)/sum(y), 'r', 'Linewidth', 3);
    hold off

    title('Noise distribution estimation')
    ylabel('Counts'); 
    xlabel('Residual noise magnitude');

    legend('Residual noise','Estimated PDF')
    grid on
    
end

%% Simulate measurement acquisition =======================================

% Set up optical system aberrations up to degree maxDeg

maxDeg      = 5;
[B,Z]       = zernikeBasis( maxDeg, nFx, nFy );

zcoefTrue       = randn(21,1);
zcoefTrue(1:3)  = 0;
zcoefTrue       = magnitude * (zcoefTrue ./ norm(zcoefTrue));
scoefTrue       = Z*zcoefTrue;

WEF         = polyCoeffEval2D( zcoefTrue, nx, B );

if(plots)
    
    figure; imagesc(fx,fx,WEF.*Pupil); colorbar
    title('Wavefront error function')
    xlabel('Residual noise magnitude');
    
end

% Define directions of input plane waves by azimuthal angle measured from
% the +x axis (1st entry in row) and angle of deflection from +z axis (2nd
% entry in row) (ALL ANGLES IN DEGREES)

srcs        = [     0       0       ; ...
                    0       0.945   ; ...
                    45      1.340   ; ...
                    90      0.945   ; ...  
                    135     1.340   ; ...
                    180     0.945   ; ...
                    225     1.340   ; ...
                    270     0.945   ; ...
                    315     1.340   ];

% For each directional source, get u and v shifts in normalized frequency
% coordinate space, and produce measured image

numImgs     = size(srcs,1);
U           = zeros(numImgs,1);
V           = zeros(numImgs,1);

E1xs        = zeros(nx,nx,numImgs);
Ixs         = zeros(nx,nx,numImgs);

for i = 1 : numImgs
    src     = srcs(i,:);
    d       = lambda_m / sind( src(2) );
    ffx     = cosd( src(1) ) / d;
    ffy     = sind( src(1) ) / d;
    
    U(i)    = ffx / fc;
    V(i)    = ffy / fc;
    
    sPup    = double( (nFx+U(i)).^2 + (nFy+V(i)).^2 <= 1 );
    sscoefs = getCoefs2DShiftedPoly( scoefTrue, maxDeg, U(i), V(i), 1, 1 );
    sWEF    = polyCoeffEval2D( Z\sscoefs, nx, B );
    E1x     = ifft2c( E0u .* sPup .* exp(1i .* sWEF) );
    Ix      = abs( E1x ).^2;    
    
    E1xs(:,:,i) = E1x;
    Ixs(:,:,i) = Ix;
end
clear src d ffx ffy illum E1x Ix


%% Create stacks of pupil domain masks

Ms = zeros(nx,nx,numImgs);
for i = 1 : numImgs
    mm1 = ( (nFx + U(i)).^2 + (nFy + V(i)).^2 ) <= 1;
    mm2 = ( (-nFx + U(i)).^2 + (-nFy + V(i)).^2 ) <= 1;
    Ms(:,:,i) = double( and( mm1, mm2 ) );
end


%% Save Dataset

meas    = Ixs;
pupDom  = Ms;
phiEst  = gwinEst.*nx./ps_m^2;
z_coef  = zcoefTrue;
s_coef  = scoefTrue;
save(   savename,'meas','pupDom','WEF','phiEst', 'nFx','nFy', ...
        'NA','ps_m','lambda_m','U','V','z_coef','s_coef','raylMean');

