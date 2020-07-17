function [ sigmaRayl, window, def ] = calibrateSpeckle( calImgFT, imgPar, GUESS, varargin )

if nargin == 3
    plots = 1;
else
    plots = 0;
end

const       = 1e6;

nx          = size(calImgFT,1);
fc          = imgPar.NA / imgPar.lambda_m;
ndfx        = 1/(imgPar.nx *imgPar.ps_m) /fc;
nfx         = (-floor(nx/2)+1 : floor(nx/2)) .*ndfx;
[nFx,nFy]   = meshgrid( nfx, nfx );
pupil       = 1*( nFx.^2 + nFy.^2 < 1 );

% Compute periodogram of intensity measurement 
P           = abs(calImgFT ./ const).^2;

% Precompute matrices/constants for speed
fx2fy2      = nFx.^2 + nFy.^2;
c           = [ nx^2*imgPar.ps_m*sqrt(2*pi),  ...
                -2*pi^2*fc^2*imgPar.ps_m^2, ...
                pi*imgPar.lambda_m*(1e-6)*fc^2 ];
fxp         = nfx((nx/2+1):end);

% Set up optimization cost function
fun         = @(x) diffCharOpt(x, c, fx2fy2, P, pupil);
    
% Run unconstrained quasi-newton algorithm
options     = optimset('MaxFunEvals', 20000, 'MaxIter', 2000);
[sol,~]     = fminsearch(fun, GUESS, options);

def         = sol(2);

%% Plot accuracy of recovery

tmp = c(1) * (sol(1)/sol(3)).*exp(c(2)*(sol(3))^2*fx2fy2);

if( plots )
    iavg = radial_avg(P,nx/2);

    fontname = 'Helvetica';
    set(0,'defaultaxesfontname',fontname);
    set(0,'defaulttextfontname',fontname);

    fontsize = 16;
    set(0,'defaultaxesfontsize',fontsize-2);
    set(0,'defaulttextfontsize',fontsize);

    h = figure;
    plot(fxp.*fc*1e-6,sqrt(iavg),'-k','LineWidth',6)
    hold on
    plot(fxp.*fc*1e-6,abs(tmp((nx/2+1),(nx/2+1):end) .* ... 
        sqrt(sin(c(3)*sol(2)*(fxp.^2)).^2)) ,'-r','LineWidth',2)
    plot(fxp.*fc*1e-6,abs(tmp((nx/2+1),(nx/2+1):end)),'-b','LineWidth',2)
    hold off
    xlabel('$\rho$ ($\mu$m$^{-1}$)','interpreter','latex')
    ylabel('$\langle \hat{I_{\O}}(\rho) \rangle$','interpreter','latex')

    legend('Measurement radial average','Fitted model (damped defocus kernel)','Extracted window function')
    grid on
    set(h,'Position',[100,100,800,275])
end

window  = tmp .* const;
P       = P .* const^2; 
Pest = abs(window).^2 .* (2*sin(c(3)*sol(2)*(fx2fy2))).^2;

if( plots )
    % Compare diffuser characterization with true diffuser properties
    figure;
    imagesc(nfx,nfx,window); colorbar

    % Periodogram estimation and comparisons
    figure;
    subplot(1,2,1), imagesc(nfx,nfx,P); colorbar
    title('Measured periodogram');
    ylabel('f_y'); xlabel('f_x');
    subplot(1,2,2), imagesc(nfx,nfx,Pest .* pupil); colorbar
    title('Est. periodogram using estimated components');
    ylabel('f_y'); xlabel('f_x');
end

%% Estimate sigma (Rayleigh)
r           = sqrt(P)./sqrt(Pest);
use         = and(Pest/max(Pest(:)) > 0.1,nFx.^2 + nFy.^2 <= 1);
rr          = r(use(:));
K           = length(rr);
sigmaRayl   = sqrt( sum(rr.^2/(2*K)) )*exp(1)*sqrt(K/(K-1))*((K-1)/K)^K;

if( plots )
    [N,E]       = histcounts(rr,200);
    X           = 0.5*(E(1:end-1) + E(2:end));
    y           = raylpdf(X,sigmaRayl);
    
    fontname = 'Helvetica';
    set(0,'defaultaxesfontname',fontname);
    set(0,'defaulttextfontname',fontname);

    fontsize = 16;
    set(0,'defaultaxesfontsize',fontsize-2);
    set(0,'defaulttextfontsize',fontsize);

    h = figure; 
    set(h,'Position',[100,100,800,275])

    bar(X,N);
    hold on
    plot(X,y*sum(N)/sum(y),'r','Linewidth',3);
    hold off

    title('\textbf{Noise distribution estimation}','interpreter','latex','fontsize',20)
    ylabel('\textbf{Counts}','interpreter','latex'); 
    xlabel('\boldmath$\eta(u,v)$','interpreter','latex');

    legend('Residual noise','Estimated PDF')
    grid on
end
