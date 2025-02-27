"""
 This module contains subroutines needed for:

  - 2D FFT differentiation in x & y:  xderiv_fft_2d, yderiv_fft_2d

  - solving the 1D Poisson's equation via FFT poisson_fft_1d

  - solving the 2D Poisson's equation via FFT: poisson_fft_2d

  - Poloidal-Toroidal Decomposition of a 2D, 3-comp. vector field: ptd_fft_2d
  
  - Gaussian Separation of a 2D, 3-comp. vector field: gauss_sep
  
  - generating a test magnetic field from a horizontal current dipole: hdip

 History: 2025/02/27 - B.T. Welsch, N.W. Jarvey - First release version

"""

import numpy as np

# This uses FFTs to compute x derivative of 1D input array, 
# which is expected to be floating-point.
#
def xderiv_fft_1d( input_arr, L):

    # Example usage: dfdx_1d = bnp.xderiv_fft_1d( input_1d_arr, xrange)
    
    # Is input array 1d?  We should check for this.
    
    nx = input_arr.shape

    dx = L/nx[0]
    xfreqs = np.fft.fftfreq(nx[0], d = dx) # wavenumbers contain pixel scale

    inspect = np.fft.fft(input_arr)

    # d/dx --> i*k_x*(Fourier inverse of input_arr), k_x = 2*pi*xfreq 
    outspect = (2*np.pi*1j)*xfreqs*inspect
           
    out_arr = np.real(np.fft.ifft(outspect)) 

    return out_arr


#============================
# This uses FFTs to solve the 1D Poisson equation
# which is expected to be floating-point.
#
def poisson_fft_1d( input_arr, L):

    # Example usage: potential_1d = poisson_fft_1d( input_1d_arr, xrange)
    
    # Is input array 1d?  Check for this.
    
    nx = np.size(input_arr)
    dx = L/nx
    
    xfreqs = np.fft.fftfreq(nx, d = dx)
    #print(L,dx,xfreqs)
  
    # Compute 1/freqs -- BEWARE: element [0,0] is zero
    inv_freqs = np.zeros(xfreqs.shape)
    nz_indices = np.nonzero(xfreqs)
    inv_freqs[nz_indices] =   1./xfreqs[nz_indices]**2

    inspect = np.fft.fft(input_arr)
    
    outspect = -inv_freqs*inspect/(2*np.pi)**2
   
    out_arr = np.real(np.fft.ifft(outspect)) 
    
    return out_arr


#============================
# This uses FFTs to compute x derivative of 2D input array, 
# which is expected to be floating-point.
#
def xderiv_fft_2d( input_arr, L):

    # Example usage: dfdx_2d = bnp.xderiv_fft_2d( input_2d_arr, xrange)
    
    # Is input array 2d?  We should check for this.
    
    ny,nx = input_arr.shape

    dx = L/nx
    xfreqs = np.fft.fftfreq(nx, d = dx) # wavenumbers contain pixel scale

    inspect = np.fft.fft2(input_arr)

    # d/dx --> i*k_x*(Fourier inverse of input_arr), k_x = 2*pi*xfreq 
    outspect = np.full_like(inspect,0. + 0.*1j)
    for i in range(0,ny):
        outspect[i,:] = (2*np.pi*1j)*xfreqs[:]*inspect[i,:]
           
    out_arr = np.real(np.fft.ifft2(outspect)) 

    return out_arr


#============================

# This uses FFTs to compute y derivative of 2D input array, 
# which is expected to be floating-point.
#
def yderiv_fft_2d( input_arr, L):

    # Example usage: dfdy_2d = bnp.yderiv_fft_2d( input_2d_arr, yrange)
    
    # Is input array 2d?  Should check for this.
    
    ny,nx = input_arr.shape

    dy = L/ny
    yfreqs = np.fft.fftfreq(ny, d = dy) # wavenumbers contain pixel scale
    
    inspect = np.fft.fft2(input_arr)


    # d/dy --> i*k_y*(Fourier inverse of input_arr), k_y = 2*pi*yfreq 
    outspect = np.full_like(inspect,0. + 0.*1j)
    for i in range(0,nx):
        outspect[:,i] = (2*np.pi*1j)*yfreqs[:]*inspect[:,i]
    out_arr = np.real(np.fft.ifft2(outspect)) 

    return out_arr



#============================
# This uses FFTs to solve the 2D Poisson equation.
# Assumes uniform, 2D grid. 
# 
# NOTE:
#
# INPUTS: 
#  input_arr = 2D, floating-point arr. of surface curl of 2D horiz. vector field
#  xrange = number of points in x component of array
#  yrange = number of points in y component of array
#  
#
# OUTPUT: 
#  out_arr = 
    
def poisson_fft_2d( input_arr, xrange, yrange):

    # Example usage: potential_2d = poisson_fft_2d( input_2d_arr, dx = dx)
    
    ny,nx = input_arr.shape

    BadDims = 0
    if (np.ndim(input_arr) != 2):
        print('Input array is not 2D.')
        BadDims = 1

    if (BadDims == 1):
        return None

    dx = xrange/nx
    dy = yrange/ny
    
    xfreqs = np.fft.fftfreq(nx, d = dx)
    yfreqs = np.fft.fftfreq(ny, d = dy)
  
    # Construct wave number arrays
    xfreqs_2d,yfreqs_2d = np.meshgrid(xfreqs,yfreqs)

    # Radius squared, in wavenumber space
    hfreqs_sq = xfreqs_2d**2 + yfreqs_2d**2

    # Compute 1/hfreqs -- BEWARE: element [0,0] is zero
    inv_hfreqs_sq = np.zeros(hfreqs_sq.shape)
    nz_indices = np.nonzero(hfreqs_sq)
    inv_hfreqs_sq[nz_indices] =   1./hfreqs_sq[nz_indices]

    inspect = np.fft.fft2(input_arr)
        
    outspect = -inv_hfreqs_sq*inspect/(2*np.pi)**2 # expected to be right
    #outspect = -inv_hfreqs*inspect/(8*np.pi) # need a factor of (pi/2) -- why?
   
    out_arr = np.real(np.fft.ifft2(outspect))
        
    return out_arr


"""
============================
 This uses FFTs to compute the Poloidal-Toroidal Decomposition of an
 input 2D array of three-component, divergence-free vector field.

 Input vector field is expected to be floating-point.
 Pixels are assumed to be square: xrange/nx = yrange/ny

 Usage: Btx,Bty,Bpx,Bpy,T,dPdz,P,bx_mean,by_mean,bz_mean = \
    bnp.ptd_fft_2d( bx0, by0, bz0)


 Output is a Dictionary, ptd, containing the keys:
    T:    Toroidal potential (from 2D curl of B_h)
    dPdz: z-derivative of Poloidal potential (from 2D divergence of B_h)
    P:    Poloidal potential (from B_z)
    btx:  x-component of B_toroidal 
    bty:  y-component of B_toroidal 
    bpx:  x-component of B_poloidal
    bpy:  y-component of B_poloidal 

"""

def ptd_fft_2d( bx, by, bz): #, T=T, P=P, dPdz=dPdz):

    nybx, nxbx = bx.shape
    nyby, nxby = by.shape
    nybz, nxbz = bz.shape

    BadSizes = 0
    if (nxbx != nxby):
        print('X extent of B_x & B_y do not match.')
        BadSizes = 1
    if (nxbx != nxbz):
        print('X extent of B_x & B_z do not match.')
        BadSizes = 1

    if (nybx != nyby):
        print('Y extent of B_x & B_y do not match.')
        BadSizes = 1
    if (nybx != nybz):
        print('Y extent of B_x & B_z do not match.')
        BadSizes = 1

    if (BadSizes == 1):
        return None

    xrange = nxbx
    yrange = nybx 
    
    # PTD vector field components & potentials
    #----------------------------------
    ptd = {}

    # Calculating mean field components
    # (not recovered by FFT methods)
    #----------------------------------    
    ptd["bx_mean"] = np.mean(bx)
    ptd["by_mean"] = np.mean(by)
    ptd["bz_mean"] = np.mean(bz)

    # source for T is 2D curl of B_horiz
    #------------------------------------
    dbxdy = yderiv_fft_2d(bx, yrange)
    dbydx = xderiv_fft_2d(by, xrange)

    ptd["T"] = poisson_fft_2d( -( dbydx - dbxdy), xrange, yrange) 
  
    # source for dPdz is 2D divergence of B_horiz
    #-----------------------------------------------
    dbxdx = xderiv_fft_2d(bx, xrange)
    dbydy = yderiv_fft_2d(by, yrange)
    ptd["dPdz"] = poisson_fft_2d( ( dbxdx + dbydy), xrange, yrange) 
  
    # source for P is B_z
    #-----------------------------------------------
    ptd["P"] = poisson_fft_2d( -bz, xrange, yrange) 
  
    # Use T & P potentials to compute toroidal & poloidal parts of B_horiz
    #-----------------------------------------------
    ptd["btx"] =  yderiv_fft_2d(ptd["T"], yrange) 
    ptd["bty"] = -xderiv_fft_2d(ptd["T"], xrange) 
    ptd["bpx"] =  xderiv_fft_2d(ptd["dPdz"], xrange) 
    ptd["bpy"] =  yderiv_fft_2d(ptd["dPdz"], yrange) 

    return ptd["btx"], ptd["bty"], ptd["bpx"], ptd["bpy"], \
        ptd["T"], ptd["dPdz"], ptd["P"], \
        ptd["bx_mean"], ptd["by_mean"], ptd["bz_mean"]
        
        
"""
============================
 This uses FFTs to perform a Gaussian Separation  of an
 input 2D array of three-component, divergence-free vector field.

 Input vector field is expected to be floating-point.
 Pixels are assumed to be square: xrange/nx = yrange/ny
 
 Usage: 
    bxlt, bylt, bzlt, psilt, \
    bxgt, bygt, bzgt, psigt, \
    bx_mean, by_mean, bz_mean = \
    bnp.gauss_sep( bx, by, bz)


 Output is a Dictionary, GS, containing the keys:
     
bxlt: x-comp. of Mag. Field from currents in z < 0 ("lt" = less than)
bylt: y-comp. of Mag. Field from currents in z < 0 ("lt" = less than)
bzlt: z-comp. of Mag. Field from currents in z < 0 ("lt" = less than) 
psilt: Scalar potential for Blt, with Blt = -grad psilt
bxgt: x-comp. of Mag. Field from currents in z > 0 ("gt" = greater than)
bygt: y-comp. of Mag. Field from currents in z > 0 ("gt" = greater than)
bzgt: z-comp. of Mag. Field from currents in z > 0 ("gt" = greater than)
psigt: Scalar potential for Bgt, with Bgt = -grad psigt
bx_mean: x-comp. of mean Mag. Field at z=0
by_mean: y-comp. of mean Mag. Field at z=0
bz_mean: z-comp. of mean Mag. Field at z=0

"""
def gauss_sep( bx, by, bz):
    
    nybx, nxbx = bx.shape
    nyby, nxby = by.shape
    nybz, nxbz = bz.shape

    BadSizes = 0
    if (nxbx != nxby):
        print('X extent of B_x & B_y do not match.')
        BadSizes = 1
    if (nxbx != nxbz):
        print('X extent of B_x & B_z do not match.')
        BadSizes = 1

    if (nybx != nyby):
        print('Y extent of B_x & B_y do not match.')
        BadSizes = 1
    if (nybx != nybz):
        print('Y extent of B_x & B_z do not match.')
        BadSizes = 1

    if (BadSizes == 1):
        return None
    
    xrange = nxbx
    yrange = nybx 
    
    # GS vector field components & potentials
    #----------------------------------
    GS = {}

    # Calculating mean field components
    # (not recovered by FFT methods)
    #----------------------------------    
    GS["bx_mean"] = np.mean(bx)
    GS["by_mean"] = np.mean(by)
    GS["bz_mean"] = np.mean(bz)
  
    # source for Bhgt & Bhlt is 2D divergence of B_horiz
    #-----------------------------------------------------
    dbxdx = xderiv_fft_2d(bx, xrange)
    dbydy = yderiv_fft_2d(by, yrange)

    phi_dirichlet = poisson_fft_2d( -( dbxdx + dbydy), xrange, yrange) 
    
    dirispect = np.fft.fft2(phi_dirichlet)
  
    # source for P is B_z: P(kx,ky) = B_z(kx,ky)/(2*pi*hfreqs)
    #---------------------------------------------------------

    dx = xrange/nxbx
    dy = yrange/nybx
    
    xfreqs = np.fft.fftfreq(nxbx, d = dx)
    yfreqs = np.fft.fftfreq(nybx, d = dy)
  
    # Construct wave number arrays
    xfreqs_2d,yfreqs_2d = np.meshgrid(xfreqs,yfreqs)

    # Radius in wavenumber space
    hfreqs = np.sqrt(xfreqs_2d**2 + yfreqs_2d**2)

    # Compute 1/hfreqs -- BEWARE: element [0,0] is zero
    inv_hfreqs = np.zeros(hfreqs.shape)
    nz_indices = np.nonzero(hfreqs)
    inv_hfreqs[nz_indices] =   1./hfreqs[nz_indices]


    bzspect = np.fft.fft2(bz)
    neumspect = inv_hfreqs * bzspect/(2*np.pi) 
   

    psigtspect = 0.5*(dirispect-neumspect)
    psiltspect = 0.5*(dirispect+neumspect)
    
    # Could take d/dx in Fourier space or (x,y) space.
    # We checked both are consistent. 
    # d/dx in real space is simpler, so comment out d/dx in Fourier (below)
    #-----------------------------------------------------------------------
    #bxgtspect = np.full_like(psigtspect,0. + 0.*1j) 
    #for i in range(0,yrange):
    #    bxgtspect[i,:] = (2*np.pi*1j)*xfreqs[:]*psigtspect[i,:]
    #bxgt_kderiv = np.real(np.fft.ifft2(bxgtspect)) 

    psigt = np.real(np.fft.ifft2(psigtspect))
    GS["psigt"] = psigt
    GS["bxgt"] = -xderiv_fft_2d(psigt,xrange)
    GS["bygt"] = -yderiv_fft_2d(psigt,yrange)

    psilt = np.real(np.fft.ifft2(psiltspect))
    GS["psilt"] = psilt
    GS["bxlt"] = -xderiv_fft_2d(psilt,xrange)
    GS["bylt"] = -yderiv_fft_2d(psilt,yrange)

    bzgtspect = - hfreqs * (2*np.pi) * psigtspect 
    GS["bzgt"] =  np.real(np.fft.ifft2(bzgtspect)) # note no (-)!

    bzltspect = - hfreqs * (2*np.pi) * psiltspect 
    GS["bzlt"] = -np.real(np.fft.ifft2(bzltspect)) 
  
    return GS["bxlt"], GS["bylt"], GS["bzlt"], GS["psilt"], \
        GS["bxgt"],  GS["bygt"], GS["bzgt"], GS["psigt"], \
        GS["bx_mean"], GS["by_mean"], GS["bz_mean"]
        
        
"""
PURPOSE: Computes vector magnetic field in a plane a distance z0
above a dipole whose axis is parallel to the plane, on a grid 
defined by two, 1D arrays of input points {x1d,y1d}.
NOTE: To make a dipole above the plane, z0 must be *negative*.

The dipole's amplitude, location with respect to origin, and
orientation with respect to the x & y axes can be set by user.

ARGUMENTS:
{x1d y1d} = 1D arrays of x & y points, will be used to form a 2D grid
{x0, y0, z0} = coordinates of dipole; z0 is depth, so z0 > 0 is submerged
amp = amplitude of dipole; dipole moment
chi = angle (in radians) dipole's axis makes with the +x axis, increases CCW

OUTPUT: 
Bx,By,Bz: Each a 2d array of mag. field comp. at z=0 due to input dipole
"""
        
def hdip(x1d, y1d, x0, y0, z0, amp, chi):
        
    import numpy as np
        
    nx = np.size(x1d)
    ny = np.size(y1d)
        

    x2d,y2d = np.meshgrid(x1d,y1d)
      
    Bz = np.full_like(x2d, 0)  
    Bx = np.full_like(x2d, 0)
    By = np.full_like(x2d, 0)
    #pot = np.full_like(x2d, 0)  # compute potential, if wanted

    for i in range(0,nx-1):
        for k in range(0,ny-1):
                
            dx = x1d[i]-x0
            dy = y1d[k]-y0

            
            angle = -(chi - np.arctan2(dy,dx))
            rsurf = np.sqrt(dx**2 + dy**2)
            rtot2 = (np.sqrt(rsurf**2 + z0**2))**2
            #rtot3 = (np.sqrt(rsurf**2 + z0**2))**3  # for potential, if wanted
            rtot5 = (np.sqrt(rsurf**2 + z0**2))**5
                
            #------------------------- 
            # Normal Component of the field
            Bz[k,i] = (3*amp/rtot5)*z0*rsurf*np.cos(angle)
            
            # Horizontal (Left-right, in-plane) component of field.
            Bx[k, i] = (amp/rtot5)*(3*dx*rsurf*np.cos(angle)-rtot2*np.cos(chi))

            # Vertical (up-down, in-plane) component of field.
            By[k, i] = (amp/rtot5)*(3*dy*rsurf*np.cos(angle)-rtot2*np.sin(chi))
             
            # Can compute scalar potential, if wanted.
            #pot[k, i] = amp*rsurf*np.cos(angle)/(4*np.pi*rtot3)
            #-------------------------
    return Bx,By,Bz 