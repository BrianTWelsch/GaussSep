"""
PURPOSE: To test the Gauss Separation routines in the gauss_decomp module,
     this main-level Python program creates a synthetic magnetogram 
     with fields due to currents above and below the magnetogram surface
     by making two calls to the hdip subroutine (in the gauss_decomp module). 
     
     (Each call to the hdip subroutine computes the vector magnetic field 
     in a plane a distance z0 above or below horizontal a dipole [axis 
     parallel to the plane], on a grid defined by two 1D arrays of input 
     points {x1d,y1d}.
     
     NOTE: To make a dipole above the plane, z0 must be *negative*.)
     
     Gaussian separation is then applied to the synthetic magnetogram, and
     the fields separated by gauss_decomp are then graphically compared to 
     each input dipole's field in several images and scatter plots.
     
     Users could vary the dipoles' amplitudes, locations, and orientations, 
     or even add more dipoles.
 
     History: 25/02/13 - N.W. Jarvey & B.T. Welsch, tests completed   
 
"""
import numpy as np
import gauss_decomp as gdc
import matplotlib.pyplot as plt

nx = 128
ny = nx

xrange = 8

# Set this to 1 to plot 2D images of input field components,
# and scatter plots of input vs. reconstructed field components.
show_plots = 1

x1d = 2*xrange*(np.arange(nx)/(nx-1)-.5)
y1d = x1d

# Used to make the arrays 2d. For Bz,Bx,By to index.
x2d,y2d = np.meshgrid(x1d,y1d)
#-------------------------

amplt = 3.0
chilt = -np.pi/12             # + is CCW from right.
x0lt = -1 #-0.71920892*np.sin(chi)
y0lt = 0.5 # 0.71920892*np.cos(chi)
z0lt = 0.4

Bltx_in,Blty_in,Bltz_in = gdc.hdip(x1d,y1d,x0lt,y0lt,z0lt,amplt,chilt)


if (show_plots == 1):
    plt.imshow(Bltz_in, extent=(-8,8,-8,8),origin='lower')
    plt.imshow(Bltx_in, extent=(-8,8,-8,8),origin='lower')
    plt.imshow(Blty_in, extent=(-8,8,-8,8),origin='lower')
    plt.close()

#------------------------------

ampgt = 2.0
chigt = np.pi/6             # + is CCW from right.
x0gt = 1 #-0.71920892*np.sin(chi)
y0gt = -0.5 # 0.71920892*np.cos(chi)
z0gt = -0.75 #2.12

Bgtx_in,Bgty_in,Bgtz_in =gdc.hdip(x1d,y1d,x0gt,y0gt,z0gt,ampgt,chigt)

if (show_plots == 1):
#   plt.imshow(Bz)
    plt.imshow(Bgtz_in, extent=(-8,8,-8,8),origin='lower')
    plt.imshow(Bgtx_in, extent=(-8,8,-8,8),origin='lower')
    plt.imshow(Bgty_in, extent=(-8,8,-8,8),origin='lower')
    plt.close()

#-------------------------
Bztot = Bltz_in + Bgtz_in   
Bxtot = Bltx_in + Bgtx_in 
Bytot = Blty_in + Bgty_in 

#    return GS["bltx"], GS["blty"], GS["bltz"], GS["psilt"], \
#        GS["bgtx"],  GS["bgty"], GS["bgtz"], GS["psigt"], \
#        GS["bx_mean"], GS["by_mean"], GS["bz_mean"]


# Use Gaussian separation to reconstruct the input "lt" & "gt" components
Bltx_rec, Blty_rec, Bltz_rec, psilt_rec, \
Bgtx_rec, Bgty_rec, Bgtz_rec, psigt_rec, \
Bx_mean_rec, By_mean_rec, Bz_mean_rec = \
gdc.gauss_sep( Bxtot, Bytot, Bztot)
 
# If desired, make scatter plots comparing input & reconstructed field 
# components; all points should fall along y = x line. 
if (show_plots == 1):
    max_Blt_in = np.max(np.sqrt(Bltx_in**2 + Blty_in**2 + Bltz_in**2))
    max_Bgt_in = np.max(np.sqrt(Bgtx_in**2 + Bgty_in**2 + Bgtz_in**2))
    plt.figure()
    plt.title(label='Bltx')
    plt.scatter(Bltx_in, Bltx_rec)
    plt.plot([-max_Blt_in,max_Blt_in],[-max_Blt_in,max_Blt_in], linestyle='dashed', color='red' )
    plt.savefig('Bltx_scatter.png',dpi = 200)

    plt.figure()
    plt.title(label='Blty')
    plt.scatter(Blty_in, Blty_rec)
    plt.plot([-max_Blt_in,max_Blt_in],[-max_Blt_in,max_Blt_in], linestyle='dashed', color='red' )
    plt.savefig('Blty_scatter.png',dpi = 200)
    plt.close()

    plt.figure()
    plt.title(label='Bltz')
    plt.scatter(Bltz_in, Bltz_rec)
    plt.plot([-max_Blt_in,max_Blt_in],[-max_Blt_in,max_Blt_in], linestyle='dashed', color='red' )
    plt.savefig('Bltz_scatter.png',dpi = 200)
    plt.close()

    plt.figure()
    plt.title(label='Bgtx')
    plt.scatter(Bgtx_in, Bgtx_rec)
    plt.plot([-max_Bgt_in,max_Bgt_in],[-max_Bgt_in,max_Bgt_in], linestyle='dashed', color='red' )
    plt.savefig('Bgtx_scatter.png',dpi = 200)
    plt.close()

    plt.figure()
    plt.title(label='Bgty')
    plt.scatter(Bgty_in, Bgty_rec)
    plt.plot([-max_Bgt_in,max_Bgt_in],[-max_Bgt_in,max_Bgt_in], linestyle='dashed', color='red' )
    plt.savefig('Bgty_scatter.png',dpi = 200)
    plt.close()

    plt.figure()
    plt.title(label='Bgtz')
    plt.scatter(Bgtz_in, Bgtz_rec)
    plt.plot([-max_Bgt_in,max_Bgt_in],[-max_Bgt_in,max_Bgt_in], linestyle='dashed', color='red' )
    plt.savefig('Bgtz_scatter.png',dpi = 200)
    plt.close()
