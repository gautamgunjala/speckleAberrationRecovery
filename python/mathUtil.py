import numpy                 as np
import numpy.fft             as nf
import scipy                 as sp

from scipy                   import special

def round2( x, m ):
    # Returns the input x rounded to the nearest multiple of m
    return np.rint(x/m)*m

##########################################################################################################################
#####     T R I G O N O M E T R Y     ####################################################################################
##########################################################################################################################

def sind( x ):
    return np.sin( np.deg2rad(x) )

def cosd( x ):
    return np.cos( np.deg2rad(x) )

def tand( x ):
    return np.tan( np.deg2rad(x) )

def asind( x ):
    return np.rad2deg( np.arcsin(x) )

def acosd( x ):
    return np.rad2deg( np.arccos(x) )

def atan2d( y, x ):
    return np.rad2deg( np.arctan2(y,x) )

##########################################################################################################################
#####     C O O R D I N A T E   C O N V E R S I O N     ##################################################################
##########################################################################################################################

def cart2pol( x, y ):
    # Returns tuple of r,t for radius (in same units as x,y)
    # and Theta (in radians) for each corresponding x,y pair
    # Assume x, y have equal dimension
    r   = np.sqrt(x**2 + y**2)
    t   = np.arctan2(y,x)
    return (r,t)

def cart2polD( x, y ):
    # Returns tuple of r,t for radius (in same units as x,y)
    # and Theta (in degrees) for each corresponding x,y pair
    # Assume x, y have equal dimension
    r   = np.sqrt(x**2 + y**2)
    t   = atan2d(y,x)
    return (r,t)

def pol2cart( r, t ):
    # Returns tuple of cartesian coordinates x,y 
    # (in same units as r) for each corresponding r,t pair
    # (t in radians)
    # Assume r, t have equal dimension
    x   = r*np.cos(t)
    y   = r*np.sin(t)
    return (x,y)
    
def polD2cart( r, t ):
    # Returns tuple of cartesian coordinates x,y 
    # (in same units as r) for each corresponding r,t pair
    # (t in degrees)
    # Assume r, t have equal dimension
    x   = r*cosd(t)
    y   = r*sind(t)
    return (x,y)

##########################################################################################################################
#####     P O L Y N O M I A L   T O O L S     ############################################################################
##########################################################################################################################

def zernikeBasis( deg, X, Y ):
    # Creates Zernike basis functions on grid specified by X and Y
    # Polynomials are ordered according to OSA/ANSI standard
    # Assume X,Y output of meshgrid
    (m,n)   = X.shape
    R       = np.sqrt(X**2 + Y**2)
    T       = np.arctan2(Y,X)
    npol    = int(0.5*(deg+1)*(deg+2))
    nn      = np.zeros((1,npol))
    mm      = np.zeros((1,npol))
    B       = np.zeros((npol,m,n))

    for i in np.arange(deg+1):
        nn[0,int(i*(i+1)/2):int((i*(i+3)+2)/2)] = i
        mm[0,int(i*(i+1)/2):int((i*(i+3)+2)/2)] = np.arange(-i,i+1,2)

    for i in np.arange(npol):
        cm      = np.abs(mm[0,i])
        cn      = nn[0,i]
        for k in np.arange(int((cn-cm)/2+1)):
            B[i]    += ((-1)**k *sp.special.comb(cn-k,k) *sp.special.comb(cn-2*k,(cn-cm)/2-k) *(R**(cn-2*k)))

        if( mm[0,i] < 0 ):
            B[i]    *= np.sin(cm*T)
        elif( mm[0,i] > 0 ):
            B[i]    *= np.cos(cm*T)
            
        B[i]    *= np.sqrt((2*(cn+1))/(1+1*(cm==0)))
            
    return B

def standardBasis2D( deg, X, Y ):
    # Creates standard basis functions on grid specified by X and Y
    # Polynomials are ordered by degree, with powers in X decreasing
    # ( i.e. 1, X, Y, X**2, X*Y, Y**2, ... )
    # Assume X,Y output of meshgrid
    (m,n)   = X.shape
    npol    = int(0.5*(deg+1)*(deg+2))
    nn      = np.zeros((1,npol))
    B       = np.zeros((npol,m,n))
    for i in np.arange(deg+1):
        for k in np.arange(i+1):
            B[int(i*(i+1)/2)+k] = (X**(i-k))*(Y**k)
            
    return B

def legendreBasis2D( deg, X, Y ):
    # Creates 2D Legendre polynomial basis with polynomials
    # up to degree deg on grid specified by X and Y
    # Assume X,Y output of meshgrid
    # Polynomials are ordered in increasing degree, with priority
    # given to powers of polynomials in X
    # ( i.e. L0(x)*L0(y), L1(x)*L0(y), L0(x)*L1(y), L2(x)*L0(y), ... 
    # L1(x)*L1(y), L0(x)*L2(y), ... )
    m,n     = X.shape
    npol    = int(0.5*(deg+1)*(deg+2))
    B       = np.zeros((npol,m,n))
    
    for i in np.arange(deg+1):
        for k in np.arange(i+1):
            cx                  = np.zeros(i+1)
            cy                  = np.zeros(i+1)
            cx[i-k]             = 1
            cy[k]               = 1
            B[int(i*(i+1)/2)+k] = np.polynomial.legendre.legval(X,cx)*np.polynomial.legendre.legval(Y,cy)
    
    return B

def zernikeStandardBasisConversion( deg ):
    # Returns a matrix Z such that if s is a vector of standard basis
    # coefficients and z is a vector of Zernike polynomial coefficients
    # in the OSA/ANSI ordering, then Z*s = z and Z\z = s for s and z 
    # with maximal polynomial degree deg
    n       = np.polynomial.legendre.leggauss(2*deg)[0]
    Xn,Yn   = np.meshgrid(n,n)
    npol    = int(0.5*(deg+1)*(deg+2))

    Zleg    = basis2EvalMap( zernikeBasis( deg, Xn, Yn ) )
    Sleg    = basis2EvalMap( standardBasis2D( deg, Xn, Yn ) )
    
    return np.linalg.lstsq(Zleg,Sleg,rcond=None)[0]

def legendreStandardBasisConversion( deg ):
    # Returns a matrix L such that if s is a vector of standard basis
    # coefficients and a is a vector of Legendre polynomial coefficients
    # ordered by degree (x priority), then L*s = a and L\a = s for s and a 
    # with maximal polynomial degree deg
    n       = np.polynomial.legendre.leggauss(2*deg)[0]
    Xn,Yn   = np.meshgrid(n,n)
    npol    = int(0.5*(deg+1)*(deg+2))

    Lleg    = basis2EvalMap( legendreBasis2D( deg, Xn, Yn ) )
    Sleg    = basis2EvalMap( standardBasis2D( deg, Xn, Yn ) )
    
    return np.linalg.lstsq(Lleg,Sleg,rcond=None)[0]

def basis2EvalMap( basis ):
    # Converts basis from ndarray with shape=(n,k,k) to the 2D evaluation map
    # which has shape=(k**2,n) and operates on column vectors of coefficients
    #nimg    = basis.shape[0]
    sz      = np.prod((basis.shape)[1:])
    basis   = basis.reshape((-1,sz),order='F')
    
    return basis.T

def poly2DTransformCoeffs( c, deg, dx, dy, sx=1, sy=1 ):
    # Returns the coefficients of f(sx*x+dx,sy*y+dy) in the standard basis
    # given the coefficients, c, of f(x,y)
    # Assume c is a column vector (shape=(n,1)) where n=(deg+1)*(deg+2)/2
    s       = np.zeros(c.shape)
    for i in np.arange(deg+1):
        for j in np.arange(deg-i+1):
            idx     = int( (i+j)*(i+j+1)/2 + j )
            coef    = c[idx]
            for n in np.arange(i+1):
                for m in np.arange(j+1):
                    idx0    = int( (n+m)*(n+m+1)/2 + m )
                    s[idx0] += coef *sp.special.comb(i,n) *(sx**n) *(dx**(i-n)) *sp.special.comb(j,m) *(sy**m) *(dy**(j-m))
    
    return s

def polyEval2D( coef, basis ):
    # Evaluates coefficients on specified basis and outputs graph
    # Assume coef is a column vector (shape=(n,1)), basis is a 3D
    # array (shape=(n,k,k)) and that coef[i] corresponds to basis[i]
    coef.shape = (-1,1,1)
    return np.sum(basis*coef,axis=0)

def polyShiftOperator( deg, x0, y0, a=1, b=1 ):
    # Returns a matrix operator M which maps standard basis polynomial
    # coefficients of f(x,y) to those of f(ax + x0, by + y0)
    n       = int(0.5*(deg+1)*(deg+2))
    M       = np.zeros((n,n))
    I       = np.identity(n)
    for i in np.arange(n):
        M[:,i]  = poly2DTransformCoeffs(I[:,i],deg,x0,y0,a,b)
    
    return M

def evenStandardPolys( deg ):
    # Returns binary vector indicating even parity of standard basis 
    # polynomials, i.e. f(x,y) = f(-x,-y)
    d   = ()
    for i in np.arange(deg+2):
        for j in np.arange(i):
            d    += (i%2,)
    return d

def nInt2DLegGauss( graph, deg, X, Y ):
    # Integrates graph of 2D function f(x,y) over a rectangular region 
    # using Legendre-Gauss quadrature with accuracy tuned by deg
    # deg should be greater than or equal to the degree of f(x) if f(x) 
    # is a polynomial. Otherwise, larger deg gives better accuracy, but
    # no guarantees!
    # Integration bounds are the min/max of X and Y
    
    xs,xe   = ( min(X.flatten()), max(X.flatten()) )
    ys,ye   = ( min(Y.flatten()), max(Y.flatten()) )
    
    L       = basis2EvalMap( legendreBasis2D( deg, X, Y ) )
    coef    = np.linalg.lstsq(L,np.transpose([graph.flatten()]),rcond=None)
    
    n,w     = np.polynomial.legendre.leggauss(2*deg)
    Xn,Yn   = np.meshgrid(n,n)
    W       = np.transpose([w])*[w]
    graphL  = polyEval2D( coef, legendreBasis2D( deg, ((xe-xs)/2)*Xn +((xe+xs)/2), ((ye-ys)/2)*Yn + ((ye+ys)/2) ) )

    return (xe-xs)*(ye-ys)/4 *np.sum( W*graphL )