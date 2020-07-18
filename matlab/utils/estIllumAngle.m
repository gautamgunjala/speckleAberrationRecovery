function [ freqXYout ] = estIllumAngle( FI, radP, freqXY, plots )

I           = ifft2c(FI);
imSz        = size(I);
imSz        = imSz(1:2);
c           = imSz(1)/2 + 1;

XYmid       = [ floor(imSz(2)./2)+1, floor(imSz(1)./2)+1 ];
[xI,yI]     = meshgrid(1:imSz(2),1:imSz(1));

freqDTh     = cart2Pol(freqXY, XYmid);

FIdivG = abs(FI);

%% Circular edge detection
freqDTh2    = circEdge(FIdivG,I,radP,freqDTh, XYmid, xI, yI);
freqXYout   = pol2Cart(freqDTh2,XYmid);

%% Display results
if(plots)
    opp         = [c,c] - (freqXYout - [c,c]);
    figure; imagesc(abs(FI)); axis image xy; 
    viscircles(freqXYout, radP);
    viscircles(opp, radP); 
end