import numpy                 as np
import numpy.fft             as nf

def fft2c( f ):
    # Computes the DC-centered forward 2D Fourier transform 
    return nf.ifftshift(nf.fft2(nf.fftshift(f)))

def ifft2c( F ):
    # Computes the DC-centered backward 2D Fourier transform
    return nf.ifftshift(nf.ifft2(nf.fftshift(F)))

def propagate( E0, lam, z, ps ):
    # Fresnel propagation of electric field, E0, by distance z
    # Assume coherent illumination with wavelength lam
    # Assume sensor (effective) pixel size ps 
    (nx,ny) = E0.shape
    dfx,dfy = (1/(nx*ps), 1/(ny*ps))
    fx      = np.arange(-(nx//2), (nx//2), 1) *dfx
    fy      = np.arange(-(ny//2), (ny//2), 1) *dfy
    Fx,Fy   = np.meshgrid(fx,fy)
    H       = np.exp(1j*2*np.pi/lam*z*np.sqrt(1 - (lam**2)*(Fx**2 + Fy**2)))
    return ifft2c(fft2c(E0)*H), H

def radialAvg( f, b=0 ):
    # Computes radial average of image f in b bins
    # Assume f is square
    # Assume b <= (f.shape)[0]//2
    if b == 0:
        b   = (f.shape)[0]//2

    (nx,ny) = f.shape
    maxR    = nx//2
    x       = np.arange(-maxR, maxR, 1)
    X,Y     = np.meshgrid(x,x)
    R       = np.sqrt(X**2 + Y**2)
    avg     = np.zeros(b)
    rs      = np.arange(maxR+1) -0.5

    for i in np.arange(maxR):
        loc     = ((R >= rs[i]) & (R < rs[i+1]))
        numel   = np.sum(loc)
        avg[i]  = np.sum(f[loc])/numel

    return avg

def rmvDC( F ):
    # Deletes the DC component of input spectrum F
    (nx,ny) = F.shape
    cx,cy   = (nx//2, ny//2)
    DC      = F[cx,cy]
    F[cx,cy]= 0
    return F, DC