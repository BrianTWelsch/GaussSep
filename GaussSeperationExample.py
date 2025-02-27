
"""
Example gaussian seperation of vector megnetogram components into fields
from currents above z = 0 (Bgt) and fields from current below z = 0 (Blt)

HISTORY: Created 24/07/30, N. Jarvey & B. Welsch
            Added comparison metrics 25/01/23, N. Jarvey & B. Welsch
"""

import gauss_decomp as gdc
import numpy as np
from astropy.io import fits
        
# From https://solartheory.ssl.berkeley.edu/~welsch/public/data/Schrijver_etal_2008_NLFFF_input/
# Bvec_File[0].header displays information about contents.
Bvec_File = fits.open('Hinode_SOT_SP_20061212_2030.fits')

Bvec_data = Bvec_File[0].data
Bx0_obs = Bvec_data[0,:,:]
By0_obs = Bvec_data[1,:,:]
Bz0_obs = Bvec_data[2,:,:]
Bvec_File.close()    

#------------------------------
# Do PTD, to get toroidal component of horizontal field 
Btx_obs,Bty_obs,Bpx_obs,Bpy_obs,T_obs,dPdz_obs,P_obs,bx_mean_obs,by_mean_obs,bz_mean_obs = \
gdc.ptd_fft_2d( Bx0_obs, By0_obs, Bz0_obs)

       
#------------------------------
# Do Gauss Sep
Bltx_obs, Blty_obs, Bltz_obs, psilt_obs, \
Bgtx_obs, Bgty_obs, Bgtz_obs, psigt_obs, \
Bx_mean_obs, By_mean_obs, Bz_mean_obs = \
gdc.gauss_sep( Bx0_obs, By0_obs, Bz0_obs)
#------------------------------
    
    
#------------------------------
    
fitsdata = np.dstack((Btx_obs, Bty_obs, \
                      Bltx_obs, Blty_obs, Bltz_obs, \
                      Bgtx_obs, Bgty_obs, Bgtz_obs))
        
hdu_out = fits.PrimaryHDU(fitsdata)
hdu_out.header['COMMENT'] = ''
hdu_out.header['COMMENT'] = 'Arrays are B_tor x,y; B^< x,y,z; B^> x,y,z'
hdu_out.header['COMMENT'] = ''

hdu_out.header['Bx_m_obs'] = bx_mean_obs
hdu_out.header['By_m_obs'] = by_mean_obs
hdu_out.header['Bz_m_obs'] = bz_mean_obs
    
#outfile = 'gauss_decomp_example_output.fits'  # Our output, for comparison
outfile = 'gauss_decomp_trial_output.fits'     # Your output
hdu_out.writeto(outfile)


do_comparison = 1 # set this to zero to *not* run comparison metrics
if (do_comparison == 1):
    example_output = fits.open('gauss_decomp_example_output.fits')
    ref = example_output[0].data
    example_output.close()
    Btx_ref =  ref[:,:,0] 
    Bty_ref =  ref[:,:,1] 
    Bltx_ref =  ref[:,:,2] 
    Blty_ref =  ref[:,:,3] 
    Bltz_ref =  ref[:,:,4]                  
    Bgtx_ref =  ref[:,:,5] 
    Bgty_ref =  ref[:,:,6] 
    Bgtz_ref =  ref[:,:,7]  
    Btx_maxabsdiff = np.max(np.abs(Btx_obs - Btx_ref))
    Bty_maxabsdiff = np.max(np.abs(Bty_obs - Bty_ref))
    Bltx_maxabsdiff =  np.max(np.abs(Bltx_obs - Bltx_ref))
    Blty_maxabsdiff = np.max(np.abs(Blty_obs - Blty_ref))
    Bltz_maxabsdiff = np.max(np.abs(Bltz_obs - Bltz_ref))
    Bgtx_maxabsdiff = np.max(np.abs(Bgtx_obs - Bgtx_ref))
    Bgty_maxabsdiff = np.max(np.abs(Bgty_obs - Bgty_ref))
    Bgtz_maxabsdiff = np.max(np.abs(Bgtz_obs - Bgtz_ref))
    print("diff. for B_tor,x = ", Btx_maxabsdiff)
    print("diff. for B_tor,y = ", Bty_maxabsdiff)
    print("diff. for Blt,x = ", Bltx_maxabsdiff)
    print("diff. for Blt,y = ", Blty_maxabsdiff)
    print("diff. for Blt,z = ", Bltz_maxabsdiff)
    print("diff. for Bgt,x = ", Bgtx_maxabsdiff)
    print("diff. for Bgt,y = ", Bgty_maxabsdiff)
    print("diff. for Bgt,z = ", Bgtz_maxabsdiff)
    
    '''
    Comparing files created by Windows & MacOS, differences were:
    diff. for B_tor,x =  1.7251977624255233e-11
    diff. for B_tor,y =  2.0834889369325538e-11
    diff. for Blt,x =  4.588684987538727e-11
    diff. for Blt,y =  5.070432962384075e-11
    diff. for Blt,z =  2.6830093702301383e-11
    diff. for Bgt,x =  2.2289725620794343e-11
    diff. for Bgt,y =  1.7678303265711293e-11
    diff. for Bgt,z =  2.6396662633487722e-11
    '''