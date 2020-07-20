function [ Hpix, Vpix, minCost ] = findCircGuess( F, radP )

F       = F ./ max(F(:));
spc     = linspace(-radP,radP,50);
[Xs,Ys] = meshgrid(spc,spc);
[m,n]   = size(Xs);
V       = zeros(size(Xs));

for ii = 1 : m*n
    V(ii)   = LossFn(Xs(ii),Ys(ii),F,radP);
end
[minCost,ind]   = min(V(:));
Hpix            = Xs(ind);
Vpix            = Ys(ind);

end

function [ cost ] = LossFn( Hpix, Vpix, F, radP )

[ sx, sy ]  = size(F);
xcen        = floor(sx/2) + 1;
ycen        = floor(sy/2) + 1;

x           = (1:sx) - xcen;
y           = (1:sy) - ycen;
[ X, Y ]    = meshgrid(x,y);

mask        = double(and( (X - Hpix).^2 + (Y - Vpix).^2 >= radP.^2 , ...
                          (-X - Hpix).^2 + (-Y - Vpix).^2 >= radP.^2 ));

tmp         = mask.*F;
cost        = var(sqrt(tmp(:)));

end