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
angleV=wrapTo180(0:180/64:360); %(in degrees), row vector
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
    thetaV=wrapTo180(chTh+theta2(ii));
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
        thetaV=wrapTo180(chTh+thetaV(unTh));
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

