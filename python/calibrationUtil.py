import numpy                 as np
import numpy.fft             as nf
import scipy                 as sp
import imgUtil 				 as iu
import mathUtil 			 as mu

def diffCharSSE( x, c, fx2fy2, P ):
	s1 	= np.exp( c[1]*(x[2]**2)*fx2fy2 )
	s2  = np.sin( c[2]*x[1]*fx2fy2 )
	s 	= ((c[0]*(x[0]/x[2])*s1*s2)**2 - P)**2
	return np.sum(s)

def costFnSSE( IuMag, domain, A, raylMu, coef, g=True ):
    nimg    = IuMag.shape[0]
    ncoef   = coef.shape[0]
    d       = IuMag.shape[1]
    grad    = np.zeros((ncoef,1))
    cost    = 0
    out     = ()
    
    for k in np.arange(nimg):
        cD      = domain[k]
        cI      = IuMag[k]
        cA      = A[k]
        Ac      = np.reshape(np.dot(cA,coef),(-1,1),order='F')
        
        # Update cost
        diffs   = cI - raylMu*np.reshape(np.abs(np.sin(Ac)),(d,d),order='F')
        cost    += np.nansum(cD*(diffs**2))
        
        if( g ):
            t1      = np.reshape(-2*raylMu*cD*diffs,(d**2,1),order='F')
            t1[np.isnan(t1)] = 0
            t2      = np.reshape(cD,(-1,1),order='F')*np.sign(np.sin(Ac))
            jac     = np.zeros(cA.shape)
            for j in np.arange(cA.shape[1]):
                jac[:,j] = np.squeeze(t2)
                jac[:,j] *= cA[:,j]
                
            grad    += np.dot(jac.T,t1)
    
    out     += (cost,)
    if( g ):
        out     += (grad,)
    
    return out

def findCirc(F, radP):
	F 		= F / np.max(np.ravel(F))
	spc 	= np.linspace(-radP, radP, 50)
	Xs,Ys 	= np.meshgrid(spc,spc)
	m,n 	= Xs.shape
	Xs,Ys 	= (np.ravel(Xs),np.ravel(Ys))

	V 		= np.zeros(m*n)
	for ii in np.arange(m*n):
	  	V[ii] 	= __findCircLossFn(Xs[ii],Ys[ii],F,radP)

	ind 	= np.argmin(V)
	return (Xs[ind],Ys[ind])

def __findCircLossFn(Hpix, Vpix, F, radP):
    sy,sx 	= F.shape
    xc,yc 	= (np.ceil(sx/2),np.ceil(sy/2))
    x,y 	= (np.arange(sx) - xc,np.arange(sy) - yc)
    X,Y 	= np.meshgrid(x,y)
    
    mask 	= F*(((X - Hpix)**2 + (Y - Vpix)**2 >= radP**2) & ((-X - Hpix)**2 + (-Y - Vpix)**2 >= radP**2))
    return np.var(np.sqrt(np.ravel(mask)))

def findCircleCenter(F, radP, est, pxAcc):
	sy,sx 	= F.shape
	x 		= np.arange(-np.floor(sx/2),np.floor(sx/2))
	X,Y 	= np.meshgrid(x,x)
	cs 		= np.arange(-5,6,1)*pxAcc
	cX,cY 	= np.meshgrid(cs,cs)

	# For each circle center to test, store x,y-coordinates, norms
	# Restrict test points to the right half-plane, flip if necessary
	cX,cY 	= np.ravel(cX + est[0]), np.ravel(cY + est[1])
	cY 		*= 2*(cX > 0) - 1
	cX 		*= 2*(cX > 0) - 1
	cent 	= np.append(cX[:,np.newaxis],cY[:,np.newaxis],axis=1)
	cD 		= np.sqrt(cent[:,0]**2 + cent[:,1]**2)[:,np.newaxis]

	# Get angles from circle center to chord endpoints
	orth    = __get2dOrthVec(cent)
	chPt1   = orth * np.sqrt(-1*(cD**2 - radP**2))
	chPt2   = -orth * np.sqrt(-1*(cD**2 - radP**2))
	chPh1   = mu.atan2d(chPt1[:,1]-cY, chPt1[:,0]-cX)
	chPh2   = mu.atan2d(chPt2[:,1]-cY, chPt2[:,0]-cX)
	phiRng  = np.append(chPh1[:,np.newaxis],chPh2[:,np.newaxis],axis=1)

	# Compute new x,y evaluation points and spline interp function values
	# Duplicate center points to include circle centered in left half-plane
	npts    = cent.shape[0]
	cent    = np.append(cent,-cent,axis=0)
	phiRng  = np.append(phiRng,phiRng+180,axis=0)

	# #########################################################################
	# #########################################################################
	# FIXED ALGORITHM PARAMETERS
	# h     : numerical derivative spacing, e.g. (f(x+h)-f(x-h))/2h
	# nphi  : number of radial lines to evaluate in each phi range
	# sigma : parameter of Gaussian smoothing kernel applied to image, I
	# #########################################################################
	# #########################################################################
	h 		= 10
	nphi 	= 10
	sigma 	= 10

	rad 	= np.array([radP-h, radP+h, radP+sigma-h, radP+sigma, radP+sigma+h])
	rad 	= np.tile(rad,(nphi,2*npts,1))
	phis 	= np.zeros((nphi,2*npts,1))

	for ii in np.arange(2*npts):
	    phis[:,ii,0]    = np.linspace(phiRng[ii,0],phiRng[ii,1],nphi)
	    
	phis    = np.tile(phis,(1,1,5))
	cenX    = np.tile(cent[:,0][:,np.newaxis],(nphi,1,5))
	cenY    = np.tile(cent[:,1][:,np.newaxis],(nphi,1,5))

	xnew    = cenX + rad*mu.cosd(phis)
	ynew    = cenY + rad*mu.sind(phis)
	outShp  = xnew.shape

	epsv    = np.finfo(float).eps*(np.arange(np.prod(xnew.shape))+1)
	Fsm     = sp.ndimage.filters.gaussian_filter(F,10)

	oldPts  = np.append(np.ravel(X)[:,np.newaxis],np.ravel(Y)[:,np.newaxis],axis=1)
	fnew    = sp.interpolate.griddata(oldPts,np.ravel(F),(np.ravel(xnew),np.ravel(ynew)),method='cubic')
	fnew    = np.reshape(fnew,outShp)

	fpAtR   = np.sum( (0.5/h)*fnew[:,:,1] - (0.5/h)*fnew[:,:,0], axis=0 )
	fppAtRs = np.sum( (fnew[:,:,2] + fnew[:,:,4] - 2*fnew[:,:,3])/(h**2), axis=0 )
	fpAtR   = fpAtR[0:npts] + fpAtR[npts:]
	fppAtRs = fppAtRs[0:npts] + fppAtRs[npts:]

	idx     = np.argmax(np.abs(fpAtR + fppAtRs))
	return (cent[idx,0],cent[idx,1])

def __get2dOrthVec(coord):
# Returns 2d vectors orthogonal to each x,y pair, or (0,-1) if input (0,0)
# Ensures that the cross product of coord and the returned vector points in
# negative z
# Assume:   coord has column vectors [x,y]
	n       = coord.shape[0]
	rSum    = np.sqrt((coord**2).sum(axis=1))
	coord   = coord / rSum[:, np.newaxis]

	tmp     = np.random.randn(n,2)
	tmp     -= np.tile((coord*tmp).sum(axis=1)[:,np.newaxis],(1,2))*coord
	orth    = tmp / np.tile(np.sqrt((tmp**2).sum(axis=1))[:,np.newaxis],(1,2))

	coord2  = np.zeros((n,3))
	coord2[:,0:2] = coord
	orth2   = np.zeros((n,3))
	orth2[:,0:2] = orth
	cross   = np.cross(coord2,orth2)[:,2]
	orth    *= np.tile(-1*cross[:,np.newaxis],(1,2))

	nanloc  = np.isnan(orth[:,0])
	orth[nanloc,0] = -1
	nanloc  = np.isnan(orth[:,1])
	orth[nanloc,1] = 0
	
	return orth 
