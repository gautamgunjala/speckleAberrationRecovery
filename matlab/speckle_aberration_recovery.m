%% Speckle aberration recovery

%%
% Gautam Gunjala and Antoine Wojdyla

%% physical parameters

% physical paramaters

% wavelength
lambda_m = 13.5e-9;

% pixel size
px_m = 15e-9;

% numerical aperture
NA = 0.33/4;

%% Load experimental data

folder = '/Users/awojdyla/speckleAberrationRecovery/data';

N_img = 10;
img = cell(1,10);
for i_f=1:N_img
    filename = sprintf('SHARP-%02.0f.png',i_f-1);
    filepath = [folder,filesep,filename];
    img{i_f} = double(imread(filepath));
end

% illuminations
sx = [+5.000e-02, +0.000e+00, -2.000e-01, +0.000e+00, +2.000e-01, -2.000e-01, +2.000e-01, -2.000e-01, +0.000e+00, +2.000e-01];
sy = [+0.000e+00, +0.000e+00, -2.000e-01, -2.000e-01, -2.000e-01, +0.000e+00, +0.000e+00, +2.000e-01, +2.000e-01, +2.000e-01];


%% display one image

% select image:
i_f = 1;

N_px = 2048;
x_m = (0:px_m:(N_px*px_m))-(N_px*px_m)/2;
imagesc(x_m*1e6,x_m*1e6,img{i_f})
title(sprintf('image %02.0f',i_f))
xlabel('position [um]')
ylabel('position [um]')
colormap gray
axis image

%% read hdf5

folder = '/Users/awojdyla/speckleAberrationRecovery/data';
filename = 'speckle_data.hdf5';
filepath = [folder,filesep,filename];
N_img = 10;
img = cell(1,10);
sx = zeros(1,10);
sy = zeros(1,10);
for i_f=1:N_img
    hier = sprintf('/image%02.0f/',i_f);
    img{i_f} = h5read(filepath,[hier,'data']);
    sx(i_f)  = h5read(filepath,[hier,'sx']);
    sy(i_f)  = h5read(filepath,[hier,'sy']);
end

%% pre-process data

% helper functions

% frequency scale
fs  = @(t) (0:1/((t(2)-t(1))*length(t)):(1/(t(2)-t(1))-1/((t(2)-t(1))*length(t)))) - (1/(t(2)-t(1))-mod(length(t),2)*1/((t(2)-t(1))*length(t)))/2;
% centered Fourier transform
ft  = @(img)  fftshift( ifft2( ifftshift( img ) ) );
% centered inverse Fourier transform
ift = @(IMG)  fftshift(  fft2( ifftshift( IMG ) ) );

