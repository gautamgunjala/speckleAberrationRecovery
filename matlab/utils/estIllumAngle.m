function [ freqXYout ] = estIllumAngle( FI, radP, freqXY, plots )
% Fine illumination angle estimation code based on the following paper:
% https://doi.org/10.1364/AO.57.005434
% 
% Original code written by Regina Eckert
% Modified for this project by Gautam Gunjala  
%

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

end

function [ freqDTh2 ] = circEdge(FIdivG, Ithresh, rad, freqDTh, XYmid, xI, yI, PDI, PD2ind )

if ~exist('PDI','var')
    estPDI=true;
else
    estPDI=false;
end

%Circular edge detection
imSz=size(FIdivG);
numImg=size(FIdivG,3);
centD2=freqDTh(:,1);
theta2=freqDTh(:,2);

%Initialize
chRad=(3:-0.5:-3)'; %Range to scan radius in pixels
%**************************************************************************
chTh=-5:0.25:5; %Angle of circle center from origin in degrees
% chTh=permute(chTh,[3 1 2]); %Make theta stretch into 3rd dim
chDi=-(10/250)*imSz(1):0.5:(10/250)*imSz(1); %Distance of circle center from origin
%Values based on a 20 pixel range on a 250x250 pixel image with circle
%radius = 60 pixels
% chDi=permute(chDi,[4 3 1 2]); %Make distance stretch into 4th dim
%**************************************************************************
radV=chRad + rad;

numTh=length(chTh);
numD=length(chDi);
numR=length(chRad);

dRad1=radV(1:end-1)+diff(radV)./2; %Radii assoc. with first derivative
dRad2=radV(2:end-1);

% PDind=find((chRad(1:end-1)+diff(chRad)./2)<0); 

%Set up angles defining circle
angleV=wrapAngle(0:180/64:360); %(in degrees), row vector
numA=length(angleV);

%Pre expand some of the matrices that don't depend on image
rMat=repmat(radV,[1 numA numTh numD]);
aMat=repmat(angleV,[numR 1 numTh numD]);
rcos=rMat.*cosd(aMat);
rsin=rMat.*sind(aMat);

%Storage
pixDMat=zeros(numTh,numD,numR-1,numImg);
pixD2Mat=zeros(numTh,numD,numR-2,numImg);
iMat=zeros(numImg,2);

maxPD=false(numTh,numD,numImg);
maxPD2=maxPD;

maxBothxy=cell(numImg,2);
adjD=zeros(numImg,1);
adjTh=adjD;

filtXY=zeros(numImg,2,1); %Best x,y coord for each dataset (+overall best) & each metric

% disp('Circular Edge Filter:')
% fprintf('Evaluating image: ')
n=0;
%tic
for ii=1:numImg

    %Print current image number
%     fprintf(repmat('\b',1,n));
    strV=num2str(ii);
%     fprintf(strV);
    n=numel(strV);
    
    %Select circle around the center point
    thetaV=wrapAngle(chTh+theta2(ii));
    distV=chDi+centD2(ii);
    
    done=false;
    while ~done
    
    flagLess=0; %Flag to indicate if distance crosses origin
    validD=true(size(distV));
    if any(distV<0)
        flagLess=1;
        validD=distV>=0; %That is, want distV>=0
    end

    [arcMat, numA2]=calArc(rad, distV, thetaV, angleV);
    arcMat=repmat(arcMat,[numR 1 1 numD]);
    
    %Rad in first dim, angle in 2nd
    %Put theta in 3rd dim, distance in 4th

    %Calculate xC, yC along axial spokes to try
    xC=repmat(distV,[numTh 1]).*cosd(repmat(thetaV',[1 numD]))+XYmid(1);
    yC=repmat(distV,[numTh 1]).*sind(repmat(thetaV',[1 numD]))+XYmid(2);
    
    %Calculate the x,y coordinates of the circle edges
    xR=reshape(rcos(arcMat),[numR numA2 numTh numD])+repmat(permute(xC,[3 4 1 2]),[numR numA2 1 1]);
    yR=reshape(rsin(arcMat),[numR numA2 numTh numD])+repmat(permute(yC,[3 4 1 2]),[numR numA2 1 1]);
    
    %Interpolate the circle values at these locations
    pixVal=interp2(xI,yI,FIdivG(:,:,ii),xR,yR,'spline');
    %pixMat=reshape(pixVal,[numR numA numTh numD]);
    
    %Derivatives
    pixMean=squeeze(mean(pixVal,2)); %Mean value around rim
 
    %1st derivative
    pixD=diff(pixMean);
    
    %2nd derivative
    pixD2=diff(pixMean,2);
    
    %Select which radii to use
    if estPDI
    PDind=find(dRad1<rad); %Only use radii less than the radius
    [~,mPDI]=max(max(pixD(PDind,:),[],2));
    PDI=PDind(mPDI);
    PD2ind=find(dRad2>=dRad1(PDI)+1.2 & dRad2<=dRad1(PDI)+3);
    end
  %   PDI=7;
     
    %Select which radii to use for the second derivative
    %Choose the one with max value, but only larger than the actual radius
    %Only use radii larger than the radius
    [~,mPD2I]=max(max(pixD2(PD2ind,:),[],2));
    PD2I=PD2ind(mPD2I);
    
    %Save it
    iMat(ii,:)=[PDI,PD2I];
    
    %Store
    pixDMat(:,:,:,ii)=permute(pixD,[2 3 1]);
    pixD2Mat(:,:,:,ii)=permute(pixD2,[2 3 1]);
    
    %Deal with case where distV<0 (so don't get circle on other side)
    if ~flagLess
        %No problems with crossing origin
        %1st Derivative
        pd=squeeze(pixD(PDI,:,:)); %Select radius
        mxL=max(pd(:));
        sL=std2(pd);
        maxPD(:,:,ii)=(pd>=mxL-0.1.*sL); %Select region around max

        %2nd Derivative
        pd2=squeeze(pixD2(PD2I,:,:)); %Select radius
        mxL2=max(pd2(:));
        sL2=std2(pd2);
        maxPD2(:,:,ii)=(pd2>=mxL2-0.25.*sL2); %Select region around max
    else
        %Same as above, but ignoring center distances where cross origin
        pd=squeeze(pixD(PDI,:,:));
        sL=std2(pd);
        pd(:,~validD)=0;
        mxL=max(pd(:));  
        maxPD(:,:,ii)=(pd>=mxL-0.1.*sL);

        pd2=squeeze(pixD2(PD2I,:,:));
        sL2=std2(pd2);
        pd2(:,~validD)=0;
        mxL2=max(pd2(:));
        maxPD2(:,:,ii)=(pd2>=mxL2-0.25.*sL2);
    end
    
    %Combine them both into one map
    maxBoth=maxPD(:,:,ii) & maxPD2(:,:,ii);

    if ~any(any(maxBoth))
        %If there's no overlap between the two, just use the 1st derivative
        %map
        maxBoth=maxPD(:,:,ii);
    end
    maxBothxy(ii,:)=[{xC(maxBoth)},{yC(maxBoth)}]; %Save the possible coordinates
    
    %Check if the only triggered pixels are on the extreme edges of the
    %search area - means need to change search area
    [thI,dI]=find(maxBoth);
    unDi=unique(dI); unTh=unique(thI);
    
    adjBnd=false;
    if length(unDi)==1 && (unDi==1 || unDi==numD)
        %Only have one distance and it's on the edge
        adjD(ii)=distV(unDi);
        distV=distV(unDi)+chDi;
        adjBnd=true; %Adjust bounds is true --> iterate
    end
    if length(unTh)==1 && (unTh==1 || unTh==numTh)
        %Only have one theta and it's on the edge --> iterate
        adjTh(ii)=thetaV(unTh);
        thetaV=wrapAngle(chTh+thetaV(unTh));
        adjBnd=true;
    end
    
    if ~adjBnd
        %If not adjusting bounds, iterate
        done=true;    
        tempPt=[maxBothxy{ii,1},maxBothxy{ii,2}]; %Extract x,y coordinates
      
        tempErr=-imageErr(Ithresh(:,:,ii),tempPt(:,1),tempPt(:,2),rad,xI,yI);
    
        %Select the circle center that gives minimum RMSE error
        [~,mI]=max(tempErr,[],1);
        filtXY(ii,1:2)=tempPt(mI,1:2);
        

    else
%         fprintf('\nAdjusting image %i to distance %f pixels, theta %f degrees.\n',ii,adjD(ii),adjTh(ii))
%         fprintf('Evaluating image: ');
        n=0;
    end
    end
    

end

freqDTh2=cart2Pol(filtXY,XYmid);

% fprintf('\n')
%toc

%Take stock of the number of points in each map that we've identified as
%candidates
% numPts(:,1)=squeeze(sum(sum(maxPD)));
% numPts(:,2)=squeeze(sum(sum(maxPD2)));
% numPts(:,3)=squeeze(sum(sum(maxBoth)));

end

function [ arcMat, numA2 ] = calArc( rad, distV, thetaV, angleV)
    %Select arc so don't use center region
    %Only let vary with theta so have the same number of angles
    if mean(distV)<=rad
        phi=acosd(mean(distV)./rad);
    else
        phi=0;
    end
    thBnd2=(thetaV)+phi-180; %Lower bound
    thBnd3=2.*(thetaV) - thBnd2; %Upper bound

    numTh=length(thetaV);
    numA=length(angleV);
    %Deal with angle wrapping issues
    th2w=wrapAngle(thBnd2);    
    th3w=wrapAngle(thBnd3);
    angleA=repmat(angleV,[numTh 1])>=repmat(th2w',[1 numA]);
    angleB=repmat(angleV,[numTh 1])<=repmat(th3w',[1 numA]);

    %Define logical index matrix to select angles
    arcMat=angleA&angleB;
    wrapLoc=repmat(thBnd2',[1 numA])<-180 |repmat(thBnd3',[1 numA])>180;
    arcMat(wrapLoc)=angleA(wrapLoc)|angleB(wrapLoc);
    
    %Make sure number of angles requested across angles is the same (data
    %size)
    numAcap=sum(arcMat,2);
    numA2=mode(numAcap);
    difA=numAcap-numA2; %Need to add or subtract angles
    
    if any(difA<0)
        %Need to add angles
        ind=find(difA<0);
        for kk=1:length(ind) %Loop across all vectors
            numAdd=abs(difA(ind(kk))); %Number to add
            temp=arcMat(ind(kk),:); 
            dT=diff(temp); %Find edges
            even=floor(numAdd/2); %Find even number to add
            odd=mod(numAdd,2); %Find if need to add additional 1
            
            up=find(dT<0); %Upper end (where go 1-->0)
            bot=find(dT>0);%Lower end (where go 0-->1)
            
            if isempty(up)
                up=length(temp); %If none, first in vector is 0 , last is 1 
            end
            if isempty(bot)
                bot=1;  %If none, first is 1, last is 0
            end
            
            upI=up+1:up+even+odd; %Add odd and half of the even to the upper end
            upI(upI>length(temp))=upI(upI>length(temp))-length(temp); %Wrap to front if go off the end
            
            botI=bot:-1:bot-even+1; %Add half of even to front
            botI(botI<1)=botI(botI<1)+length(temp); %Wrap around
            
            temp([upI botI])=1; %Use indices to change value
            arcMat(ind(kk),:)=temp; %Re-store
        end
    end
    if any(difA>0)
        %Need to delete angles
        ind=find(difA>0);
        for kk=1:length(ind)
            numDel=abs(difA(ind(kk)));%Loop across all vectors
            temp=arcMat(ind(kk),:);%Number to delete
            dT=diff(temp);%Find edges
            even=floor(numDel/2);
            odd=mod(numDel,2);
            
            up=find(dT<0);
            bot=find(dT>0);
                        
            if isempty(up)
                up=length(temp);  
            end
            if isempty(bot)
                bot=1;  
            end
            
            upI=up:-1:up-even-odd+1;
            upI(upI<1)=upI(upI<1)+length(temp);
            
            botI=bot+1:bot+even;
            botI(botI>length(temp))=botI(botI>length(temp))-length(temp);
            
            temp([upI botI])=0;
            arcMat(ind(kk),:)=temp;
        end
    end

    arcMat=permute(arcMat,[3 2 1]);

end

function [ rmse ] = imageErr( img, xC, yC, rad, xI, yI)
%Calculates rmse of image filtered by a pupil at centers xC, yC, with
%radius rad on grid xI, yI

%Check if need to expand matrices
if size(img,3)>numel(xC)
    if numel(xC) == 1
        xC=xC.*ones([size(img,3) 1]);
        yC=yC.*ones([size(img,3) 1]);
    else
        error('Size (xC,yC) ~= size(img,3), but xC, yC are vectors. Unclear how to proceed.')
    end
elseif size(img,3)<numel(xC)
    if size(img,3) == 1
        img = repmat(img,[1 1 numel(xC)]);
    else
        error('Size(img,3) ~= size (xC,yC), but img is 3D. Unclear how to proceed.')
    end
end

pupil=zeros(size(img));
for ii=1:size(img,3)
    pupil(:,:,ii)=sqrt((xC(ii)-xI).^2+(yC(ii)-yI).^2)<=rad;
end

Iest=ifft2c(fft2c(img).*pupil);

rmse=squeeze(sqrt(sum(sum((abs(Iest)-img).^2))./numel(img(:,:,1))));

end

function [ out ] = wrapAngle( x )
    out = mod(x+179,360)-179;
end

function [RTh] = cart2Pol( XY, XYmid )
%Convert polar to cartesian
%Output r_theta=[r, theta]; r in pixels, theta in degrees
%Input x_y x,y coordinates in pixels


%Theta, distance of circle center from center frequency
centD=sqrt(sum((XY-repmat(XYmid,[size(XY,1) 1])).^2,2)); %(pixels)
theta=atan2d(XY(:,2)-XYmid(2),XY(:,1)-XYmid(1)); %(degrees)

RTh=[centD, theta];

end

function [XY] = pol2Cart( RTh, XYmid )
%Convert polar to cartesian
%r_theta=[r, theta]; r in pixels, theta in degrees
%Output x,y coordinates in pixels

xC=RTh(:,1).*cosd(RTh(:,2))+XYmid(1);
yC=RTh(:,1).*sind(RTh(:,2))+XYmid(2);
XY=[xC yC];

end





