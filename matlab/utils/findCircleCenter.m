function [ xp, yp ] = findCircleCenter( I, radP, est, pxAcc )
% Returns the center of circle in off-axis illuminated intensity spectra
% Input:    I       - intensity spectrum
%           radP    - radius of circle (in pixels)
%           est     - estimated circle center (tuple: [x,y])
%           pxAcc   - accuracy of output (in pixels)
%
% For efficient usage, call this function several times with decreasing
% input to pxAcc for a successively finer approximation.
%
% Assume:   radP is known (NA is known)
%           input image is square (n x n)
%
% -------------------------------------------------------------------------
% Regina Eckert, Gautam Gunjala - 05/2018 
% UC Berkeley Computational Imaging Lab
% -------------------------------------------------------------------------

[n,~]   = size(I);
x       = -floor(n/2):1:floor(n/2)-1;
[X,Y]   = meshgrid(x);
cs      = (-5:1:5)*pxAcc;
[cX,cY] = meshgrid(cs);

% For each circle center to test, store x,y-coordinates, norms
% Restrict test points to the right half-plane, flip if necessary
cX      = cX + est(1); cX = cX(:);
cY      = cY + est(2); cY = cY(:);
cY      = cY .* ( 2*(cX > 0) - 1 );
cX      = cX .* ( 2*(cX > 0) - 1 );
cent    = [cX,cY];
cD      = sqrt(sum(cent.^2,2));

% Get angles from circle center to chord endpoints
orth    = get2dOrthVec(cent);
chPt1   = orth .* sqrt(radP.^2 - cD.^2);
chPt2   = -orth .* sqrt(radP.^2 - cD.^2);
chPh1   = atan2d(chPt1(:,2)-cY,chPt1(:,1)-cX);
chPh2   = atan2d(chPt2(:,2)-cY,chPt2(:,1)-cX);
phiRng  = [chPh1,chPh2];

% #########################################################################
% #########################################################################
% PLOTS FOR VERIFICATION
% #########################################################################
% #########################################################################

% for ii = 1:length(cX)
%     
%     cenX    = cent(ii,1);
%     cenY    = cent(ii,2);
% 
%     m1      = (X-cenX).^2 + (Y-cenY).^2 <= radP.^2;
%     m2      = (X+cenX).^2 + (Y+cenY).^2 <= radP.^2;
%     spec    = 1*or(m1,m2);
%     
%     figure; imagesc(spec); axis square xy
%     hold on
%     phis = linspace(phiRng(ii,1),phiRng(ii,2),10);
%     scatter(cenX+251,cenY+251,'r')
%     scatter(chPt1(ii,1)+251,chPt1(ii,2)+251,'g')
%     scatter(chPt2(ii,1)+251,chPt2(ii,2)+251,'b')
%     
%     for jj = 1:10
%         pt = [200*cosd(phis(jj)),200*sind(phis(jj))]+[cenX,cenY];
%         scatter(pt(1)+251,pt(2)+251,'k')
%         quiver(cenX+251,cenY+251,pt(1)-cenX,pt(2)-cenY,0,'linewidth',2,'color',[0,1-0.1*jj,0.1*jj]);
%     end
%     hold off
% end

% Compute new x,y evaluation points and spline interp function values
% Duplicate center points to include circle centered in left half-plane
npts    = size(cent,1);
cent    = [cent; -cent];
phiRng  = [phiRng; phiRng+180];

% #########################################################################
% #########################################################################
% FIXED ALGORITHM PARAMETERS
% h     : numerical derivative spacing, e.g. (f(x+h)-f(x-h))/2h
% nphi  : number of radial lines to evaluate in each phi range
% sigma : parameter of Gaussian smoothing kernel applied to image, I
% #########################################################################
% #########################################################################

h       = 10;
nphi    = 10;
sigma   = 10;

rad     = [radP-h, radP+h, radP+sigma-h, radP+sigma, radP+sigma+h];
rad     = repmat(rad,[2*npts,1,nphi]);
phis    = zeros(2*npts,1,nphi);
for ii = 1:2*npts
    phis(ii,1,:)    = linspace(phiRng(ii,1),phiRng(ii,2),nphi);
end
phis    = repmat(phis,[1,5,1]);
cenX    = repmat(cent(:,1),[1,5,nphi]);
cenY    = repmat(cent(:,2),[1,5,nphi]);

xnew    = cenX + rad.*cosd(phis);
ynew    = cenY + rad.*sind(phis);
[s1,s2,s3] = size(xnew);
numx    = numel(xnew);
epsv    = eps*(1:numx)';

Ism     = imgaussfilt(I,10);
fnew    = griddata(X(:),Y(:),Ism(:),xnew(:)+epsv,ynew(:)+epsv,'cubic');
fnew    = reshape(fnew,[s1,s2,s3]);

% % #########################################################################
% % #########################################################################
% % PLOTS FOR VERIFICATION
% % #########################################################################
% % #########################################################################
% 
% for ii = 1:length(cX) 
%     figure; imagesc(Ism); axis square xy
%     hold on
%     scatter(cent(ii,1)+251,cent(ii,2)+251)
%     xx      = xnew(ii,:,:); xx = xx(:);
%     yy      = ynew(ii,:,:); yy = yy(:);
%     scatter(xx+251,yy+251)       
%     
% end


% Evaluate first and second numerical derivatives at radius radP
% and radP+sigma, respectively
fpAtR   = sum( (0.5/h)*fnew(:,2,:) - (0.5/h)*fnew(:,1,:), 3 );
fppAtRs = sum( (fnew(:,3,:)+fnew(:,5,:)-2*fnew(:,4,:))/(h^2), 3 );
fpAtR   = fpAtR(1:npts) + fpAtR(npts+1:end);
fppAtRs = fppAtRs(1:npts) + fppAtRs(npts+1:end);

[~,ind] = max(abs(fpAtR + fppAtRs));
xp      = cent(ind,1);
yp      = cent(ind,2);

end

function [ orth ] = get2dOrthVec( coord )
% Returns 2d vectors orthogonal to each x,y pair, or (0,-1) if input (0,0)
% Ensures that the cross product of coord and the returned vector points in
% negative z
% Assume:   coord has column vectors [x,y]
n       = size(coord,1);
coord   = coord ./ sqrt(sum(coord.^2,2)); 
tmp     = randn(size(coord));
tmp     = tmp - ( repmat(sum(coord.*tmp,2),[1,2]) .* coord );
orth    = tmp ./ sqrt(sum( (tmp).^2, 2));
c       = cross([coord, zeros(n,1)],[tmp, zeros(n,1)]);
c       = c(:,3);
orth    = orth .* ( 2*(c <= 0) - 1 );
try
    [ind,~]     = find(isnan(orth));
    ind         = unique(ind);
    for ii = 1 : length(ind)
        orth(ind(ii),:) = [-1,0];
    end
catch
end
end