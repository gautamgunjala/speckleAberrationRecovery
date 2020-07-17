import matplotlib.pyplot     as plt
import matplotlib.image      as mpimg

def imagesc( f, x=None, y=None, clim=None, axis=None, cbar=0 ):
	# Plot f with axes specified by x,y and colorbar range 
	# specified by a tuple clim
	# Default axis plots f as-is, axis="xy" flips y axis
	# Assume figure has been initialized
	if( axis is "xy" ):
		orig 	= "lower" 
	else:
		orig 	= "upper"

	if( x is not None ):
		nx      = f.shape[0]
		psx		= (x[-1]-x[0])/(nx-1)
		xmin 	= x[0] - 0.5*psx
		xmax 	= x[-1] - 0.5*psx

	if( y is not None ):
		ny 		= f.shape[1]
		psy 	= (y[-1]-y[0])/(ny-1)
		ymin 	= y[0] - 0.5*psy
		ymax 	= y[-1] - 0.5*psy

	try:
		ext 	= (xmin,xmax,ymin,ymax)
	except:
		ext 	= None

	plt.imshow(f,clim=clim,origin=orig,extent=ext)

	if( cbar == 1 ):
		plt.colorbar(fraction=0.046, pad=0.04)

