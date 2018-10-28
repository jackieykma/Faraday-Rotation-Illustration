'''

An updated illustration to Faraday rotation during propagation
Old version feature: Illustrated propagation of emission from 10kpc
Faraday rotation occurs en-route, and the pol planes are shown

In this code, the Stokes Q & U are injected (through emission) and rotated (Faraday rotation), which is integrated (or rather, summed) numerically

New features (compared to the old version):
(1) Implemented an "emission profile", similar to the original "Faraday rotation profile"
(2) New diagnostic plots
    - Faraday rotation & emission profiles
    - Stokes Q & U vs lambda2
    - PI vs lambda2
    - Faraday spectrum
    - Polarisation plane plot, at
      - Central freq only
      - Central freq, and two ends of the band
      - All ranges within the band
    - Q versus U trajectory
    - PA versus lambda2

Update list:
v1.0 (12.10.2018) Built up this script
v1.1 (14.10.2018) Fixed a bug in specifying parameter qu_max (which determines the plotting ranges of some panels)
'''

import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.gridspec as gridspec
import matplotlib.cm as mplcm
from matplotlib.collections import LineCollection

########################################
## Settings are done here

## The step interval (in pc) at which the summation of parameters are done
step = 1 ## Suggested value: 1


## The step interval (in units of <step>) at which plots are generated
plot_step = 500 ## Suggested value: 10


## Setting up the rotation and emission profiles...
## Ignore this part
d = np.arange(0, 10001, step) ## Distance in pc
neb = d * 0. ## Profile of Ne * B as function of d
emi = d * 0. ## Profile of emissivity (PI)
pa0 = d * 0.  ## PA0 of the emission


## Change this part to modify the rotation profile
## neb is a numpy array; Unit: [cm-3 uG]
neb[6000:9500] += -0.1   ## A big block
neb[8000:8500] += 0.8  ## A small sharp peak
neb += 0.5 * np.e**(-(d-4500.)**2/2.*0.004**2) ## A Gaussian
neb += -0.1 * np.e**(-(d-2500.)**2/2.*0.001**2) ## Another Gaussian
#for i in np.arange(500, 2000, step):  ## A strange triangular function
#   neb[i] += 0.1 - (i-500)*0.4/1500.


## Change this part to modify the emission profile
## "Emission" as in polarised intensity
## emi is a numpy array; Unit: [Jy pc-1]
## pa0 is a numpy array; Unit: [deg]
emi[9900:10000] += 0.001 ## A background "Faraday simple" emission
pa0[9900:10000] = 60.
emi[8000:8500] += 0.002 ## Creates a Burn slab, since rotation != 0 there
pa0[8000:8500] = -45.


## Specify the frequency band here
## The pol information is evaluated at the centre of each frequency channel.
## Instead of integrating within each channel
## In other words, this assume no significant in-channel Faraday rotation
## Or, equivalently, no significant bandwidth depolarisation within each channel
freq = np.arange(1.e9, 2.e9+1, 1.e6) ## VLA L-band: 1--2 GHz, 1 MHz channel
#freq = np.arange(2.e9, 4.e9+1, 2.e6) ## VLA S-band: 2--4 GHz, 2 MHz channel
#freq = np.arange(1.e9, 3.e9+1, 1.e6) ## ATCA LS-band: 1--3 GHz, 1 MHz channel


## Control which plots to create
plot_neb = True ## Plot ne*B profile; probably a good idea to set this as True
plot_emi = True ## Plot emission profile (same panel as neb)
plot_qu_l = True ## Plot Q&U vs lambda2?
plot_pi_l = True ## Plot PI vs lambda2?
plot_spec = False ## Plot RM-spectrum?
plot_polp = True ## Plot polarisation planes?
plot_q_u = True ## Plot Q vs U trajectory?
plot_pa_l = True ## Plot PA vs lambda2?


## Fine control of plots
save_file = False ## Save output to files instead of just plotting in a window
plot_column = 3 ## How many columns of panels in the plot?
file_ext = '.png' ## Extension of the output file (e.g., .png, .eps)
out_dir = './faraday_output/' ## Where to save the output files? Even if plots will not be saved, this directory will be created to store temporary files for RM-Synthesis
#polp_mode = 'cent' ## Only show centre frequency in polarisation plane plot
polp_mode = 'three' ## Show low, mid, high of frequency band in polarisation plane plot
#polp_mode = 'cont' ## Show continuously the entire band in polarisation plane plot
plot_cmap = 'gnuplot' ## Colour map used for some of the line plots; choose from the supported matplotlib colour maps from, e.g., https://matplotlib.org/examples/color/colormaps_reference.html
l2_range = None ## lambda2-range in some plots; if None --> Use min-max of the band
#l2_tick = None ## Let matplotlib automatically decide lambda2 label
l2_tick = [[0.15**2, 0.20**2, 0.25**2, 0.30**2], [r'$15^2$', r'$20^2$', r'$25^2$', r'$30^2$']] ## Manually put in lambda2 label. First entry is where the labels are put [m2], while second entry is the corresponding text to be put
l2_label = r'$\lambda^2$ (cm$^2$)' ## Manually control the lambda2 label. Only useful if manual l2_tick is set above! If not, probably should set this parameter to <r'$\lambda^2$ (m$^2$)'>
plot_margin = 1.2 ## For some plots (e.g. Q vs U, QU vs l2... The margin of the plot will be <plot_margin> * max data value. Recommended value: 1.2)
shade_neb = True ## Shade the neb plot up till the plotted distance step?
shade_emi = True ## Shade the emi plot up till the plotted distance step?
neb_emi_dline = True ## Plot a vertical line to show the plotted distance step?


## Adjust the figure size
## This will be run before the figure is generated
def plt_def():
   plt.figure(figsize=(10,10))

## Adjust the spacing within the plot
## This will be run right before showing / saving the plots
def plt_adj():
   plt.subplots_adjust(left=0.09, bottom=0.07, right=0.90, top=0.97, wspace=0.35, hspace=0.43)


## Control for RM-Synthesis
## Must be Michael Bell's pyrmsynth installed!
rmsyn_path = '/usr/local/pyrmsynth/rmsynthesis.py'
do_rmclean = True
rmclean_gain = 0.2 ## Gain in RM-Clean
rmclean_niter = 1000 ## Number of iterations for RM-Clean
rmclean_cutoff = 0.001 ## Cutoff in RM-Clean in Jy
fd_range = 1000 ## FD range to form spectrum up to (+/-) [rad m-2]
fd_step = 2 ## Step in Faraday spectrum [rad m-2]

## End of settings
########################################

## Main part of the script starts here
## If changes are made, record here and put as a new version

## First, convert the frequency band into lambda2 space
c = 299792458. ## Speed of light in ms-1
l2 = (c/freq)**2 ## Convert frequency [Hz] to lambda2 [m2]
q = l2 * 0. ## Set up an array to store Stokes Q per channel
u = l2 * 0. ## Set up an array to store Stokes U per channel

## Also, count the number of panels to be made (except for neb & emi)
## Important to change this when new diagnostic panels are added!
plot_panels = np.sum([plot_qu_l, plot_pi_l, plot_spec, plot_polp, plot_q_u, plot_pa_l])
plot_row = (plot_panels-1)/plot_column+1
if plot_neb or plot_emi:
   plot_row += 1

## If l2_range == None, then change xlim to min-max
if l2_range == None:
   l2_range = [np.min(l2), np.max(l2)]

## Create the output file directory, if needed
## Keep consistent format for the output directory name
if save_file or plot_spec:
   if out_dir[-1] != '/':
      out_dir += '/'
   os.system('rm -rf '+out_dir[:-1])
   os.system('mkdir '+out_dir)

## Create the RM-Synthesis parameter file, if needed
if plot_spec:
   import pyfits
   os.system('mkdir '+out_dir+'rmsyn_tmp')
   os.system('mkdir '+out_dir+'rmsyn_tmp/out')
   os.system('mkdir '+out_dir+'rmsyn_tmp/fits')
   os.system('mkdir '+out_dir+'rmsyn_tmp/fits/stokes_q')
   os.system('mkdir '+out_dir+'rmsyn_tmp/fits/stokes_u')
   f = open(out_dir+'rmsynthesis.par', 'w')
   f.write('dec_min -1\ndec_max -1\nra_min -1\nra_max -1\n')
   ## FD range and step size here
   f.write('phi_min '+str(-1*fd_range)+'\nnphi '+str(2*fd_range/fd_step)+'\ndphi '+str(fd_step)+'\n')
   ## RM-Clean parameters here
   f.write('do_clean '+str(do_rmclean)+'\ngain '+str(rmclean_gain)+'\nniter '+str(rmclean_niter)+'\ncutoff '+str(rmclean_cutoff)+'\n')
   ## Output and input directories
   f.write('outputfn '+out_dir+'rmsyn_tmp/out/rm_out\ninput_dir '+out_dir+'rmsyn_tmp/fits')
   f.close()

## Calculate how many plots will be saved, for file numbering
nplots = 1+(np.max(d)-np.min(d))/(plot_step*step)
out_no = 0 ## The file number currently writing to, as part of the file name

## Each for-loop iteration goes through one step in physical distance
for i in np.arange(len(d)-1, -1, -1):
   ## First, impose Faraday rotation to the emission from last step
   ## Convert Q and U into PI and PA for easier calculations
   pi = np.sqrt(q**2 + u**2)
   pa = 0.5*np.arctan2(u, q) ## In radians
   pa += (0.81 * neb[i] * step) * l2
   q = pi * np.cos(2.*pa)
   u = pi * np.sin(2.*pa)
   ## Then, add the extra emission from this step
   q += emi[i] * np.cos(2.*pa0[i])
   u += emi[i] * np.sin(2.*pa0[i])
   ## Polarisation calculation done for this step!
   ## Now, judge if this step should be plotted, as per the settings
   if i % plot_step == 0:
      plot_no = 0 ## The panel number currently writing to
      pan_x = 0 ## This will record which panel to write to next (column)
      pan_y = 0 ## This will record which panel to write to next (row)
      plt_def()
      axarray = []
      gs = gridspec.GridSpec(plot_row, plot_column)
      if plot_neb or plot_emi:
         axarray.append(plt.subplot(gs[pan_y, :])) ## This will be for neb/emi
         axarray[plot_no].set_xlabel('Distance from Observer (pc)')
         if plot_neb: ## First, plot neb to first panel
            axarray[plot_no].axhline(y=0., color='grey', linestyle='--')
            axarray[plot_no].plot(d, neb, lw=1.5, color='b')
            if shade_neb:
               axarray[plot_no].fill_between(d[d >= d[i]], neb[d >= d[i]], alpha=0.3, color='b')
            if neb_emi_dline:
               axarray[plot_no].axvline(x=d[i], color='k', lw=1)
            axarray[plot_no].set_ylabel(r'$n_e B_{||}$ (cm$^{-3} \mu$G)')
            axarray[plot_no].set_ylim([np.min([np.min(neb)*plot_margin, 0.]), np.max([np.max(neb)*plot_margin, 0.])])
         if plot_emi: ## Then, plot emi also to first panel
            if plot_neb: ## If already doing neb, plot on a twin-axis
               ax2 = axarray[plot_no].twinx()
               ax2.plot(d, emi, lw=1.5, color='r')
               if shade_emi:
                  ax2.fill_between(d[d >= d[i]], emi[d >= d[i]], alpha=0.3, color='r')
               ax2.set_ylabel(r'Emissivity (Jy pc$^{-1}$)')
               ## Align this 2nd y-axis with the 1st at y=0
               ## Will probably fail if emi < 0 at some point
               ax2.set_ylim([np.max(emi)*plot_margin*np.min([np.min(neb)*plot_margin, 0.])/np.max([np.max(neb)*plot_margin, 0.]), np.max(emi)*plot_margin])
               ## Align y=0 for both plots
            else: ## If not doing neb, plot on the main axis
               axarray[plot_no].plot(d, emi, lw=1.5, color='r')
               if shade_emi:
                  axarray[plot_no].fill_between(d[d >= d[i]], emi[d >= d[i]], alpha=0.3, color='r')
               if neb_emi_dline:
                  axarray[plot_no].axvline(x=d[i], color='k', lw=1)
               axarray[plot_no].set_ylabel(r'Polarised Emissivity (Jy pc$^{-1}$)')
               axarray[plot_no].set_ylim([np.min([np.min(emi)*plot_margin, 0.]), np.max([np.max(emi)*plot_margin, 0.])])
         plot_no += 1
         pan_y += 1
      if plot_qu_l:
         axarray.append(plt.subplot(gs[pan_y, pan_x]))
         axarray[plot_no].set_xlabel(l2_label)
         axarray[plot_no].set_ylabel(r'Stokes $QU$ (Jy)')
         axarray[plot_no].axhline(y=0., color='grey', linestyle='--')
         axarray[plot_no].plot(l2, q, lw=1.5, color='b', label=r'$Q$')
         axarray[plot_no].plot(l2, u, lw=1.5, color='r', label=r'$U$')
         axarray[plot_no].set_xlim(l2_range)
         axarray[plot_no].set_xticks(l2_tick[0])
         axarray[plot_no].set_xticklabels(l2_tick[1])
         qu_max = np.max([np.max(np.abs(q)), np.max(np.abs(u))])
         axarray[plot_no].set_ylim([-qu_max*plot_margin, qu_max*plot_margin])
         axarray[plot_no].legend(loc='best')
         plot_no += 1
         pan_x += 1
         if pan_x == plot_column: ## If row is finished --> Move to next row
            pan_x = 0
            pan_y += 1
      if plot_pi_l:
         axarray.append(plt.subplot(gs[pan_y, pan_x]))
         axarray[plot_no].set_xlabel(l2_label)
         axarray[plot_no].set_ylabel('PI (Jy)')
         ## Plotting the line in continuous colour map
         ## https://matplotlib.org/examples/pylab_examples/multicolored_line.html
         points = np.array([l2, np.sqrt(q**2+u**2)]).T.reshape(-1, 1, 2)
         segments = np.concatenate([points[:-1], points[1:]], axis=1)
         lc = LineCollection(segments, cmap=plt.get_cmap(plot_cmap), norm=plt.Normalize(np.min(l2), np.max(l2)))
         lc.set_array(l2)
         lc.set_linewidth(1.5)
         axarray[plot_no].add_collection(lc)
         axarray[plot_no].set_xlim(l2_range)
         axarray[plot_no].set_xticks(l2_tick[0])
         axarray[plot_no].set_xticklabels(l2_tick[1])
         pi_max = np.max(np.sqrt(q**2 + u**2))
         if pi_max == 0: ## If PI = 0, the plot will be strange. Set it to some arbitrary number
            pi_max = 0.001
         axarray[plot_no].set_ylim([0., pi_max*plot_margin])
         plot_no += 1
         pan_x += 1
         if pan_x == plot_column: ## If row is finished --> Move to next row
            pan_x = 0
            pan_y += 1
      if plot_pa_l:
         axarray.append(plt.subplot(gs[pan_y, pan_x]))
         axarray[plot_no].set_xlabel(l2_label)
         axarray[plot_no].set_ylabel(r'PA ($^\circ$)')
         ## Plotting the line in continuous colour map
         ## https://matplotlib.org/examples/pylab_examples/multicolored_line.html
         points = np.array([l2, np.degrees(0.5*np.arctan2(u, q))]).T.reshape(-1, 1, 2)
         segments = np.concatenate([points[:-1], points[1:]], axis=1)
         lc = LineCollection(segments, cmap=plt.get_cmap(plot_cmap), norm=plt.Normalize(np.min(l2), np.max(l2)))
         lc.set_array(l2)
         lc.set_linewidth(1.5)
         axarray[plot_no].add_collection(lc)
         axarray[plot_no].set_xlim(l2_range)
         axarray[plot_no].set_xticks(l2_tick[0])
         axarray[plot_no].set_xticklabels(l2_tick[1])
         axarray[plot_no].set_ylim([-90., 90.])
         plot_no += 1
         pan_x += 1
         if pan_x == plot_column: ## If row is finished --> Move to next row
            pan_x = 0
            pan_y += 1
      if plot_q_u:
         axarray.append(plt.subplot(gs[pan_y, pan_x]))
         axarray[plot_no].set_xlabel(r'Stokes $Q$ (Jy)')
         axarray[plot_no].set_ylabel(r'Stokes $U$ (Jy)')
         axarray[plot_no].axhline(y=0., color='grey', linestyle='--')
         axarray[plot_no].axvline(x=0., color='grey', linestyle='--')
         ## Plotting the line in continuous colour map
         ## https://matplotlib.org/examples/pylab_examples/multicolored_line.html
         points = np.array([q, u]).T.reshape(-1, 1, 2)
         segments = np.concatenate([points[:-1], points[1:]], axis=1)
         lc = LineCollection(segments, cmap=plt.get_cmap(plot_cmap), norm=plt.Normalize(np.min(l2), np.max(l2)))
         lc.set_array(l2)
         lc.set_linewidth(1.5)
         axarray[plot_no].add_collection(lc)
         qu_max = np.max([np.max(np.abs(q)), np.max(np.abs(u))])
         axarray[plot_no].set_xlim([-qu_max*plot_margin, qu_max*plot_margin])
         axarray[plot_no].set_ylim([-qu_max*plot_margin, qu_max*plot_margin])
         plot_no += 1
         pan_x += 1
         if pan_x == plot_column: ## If row is finished --> Move to next row
            pan_x = 0
            pan_y += 1
      if plot_spec:
         axarray.append(plt.subplot(gs[pan_y, pan_x]))
         ## Using Michael Bell's RM-Synthesis algorithm to calculate RM spectra
         ## This is needed, as emission at different physical depths can depolarised each other
         ## First, write out .fits files for RM-Synthesis
         ## Set up the header file first
         header = pyfits.Header()
         header['SIMPLE'] = True
         header['NAXIS'] = 4
         header['NAXIS1'] = 1 ## Right Ascension
         header['NAXIS2'] = 1 ## Declination
         header['NAXIS3'] = len(freq) ## Frequency
         header['NAXIS4'] = 1 ## Stokes
         header['OBJECT'] = 'Template Source for RM-Synthesis'
         header['BUNIT'] = 'Jy/beam'
         header['CTYPE1'] = 'RA---SIN'
         header['CRVAL1'] = 0
         header['CDELT1'] = -1
         header['CRPIX1'] = 1
         header['CUNIT1'] = 'deg'
         header['CTYPE2'] = 'DEC--SIN'
         header['CRVAL2'] = 0
         header['CDELT2'] = 1
         header['CRPIX2'] = 1
         header['CUNIT2'] = 'deg'
         header['CTYPE3'] = 'FREQ'
         header['CRVAL3'] = freq[0] ## In Hz
         header['CDELT3'] = freq[1]-freq[0]
         header['CRPIX3'] = len(freq)
         header['CUNIT3'] = 'Hz'
         header['CTYPE4'] = 'STOKES'
         header['CRVAL4'] = 0 ## 2 for Stokes Q; 3 for Stokes U
         header['CDELT4'] = 1
         header['CRPIX4'] = 1
         header['CUNIT4'] = ''
         ## Writing to Stokes Q & U files
         header.update('CRVAL4', 2, comment='Stokes Q')
         data_q = [] ## The output Stokes Q data
         for q_in in q:
            data_q.append([[q_in]])
         data_q = np.asarray([data_q])
         pyfits.writeto(out_dir+'rmsyn_tmp/fits/stokes_q/Q.fits', data_q, header)
         header.update('CRVAL4', 3, comment='Stokes U')
         data_u = [] ## The output Stokes U data
         for u_in in u:
            data_u.append([[u_in]])
         data_u = np.asarray([data_u])
         pyfits.writeto(out_dir+'rmsyn_tmp/fits/stokes_u/U.fits', data_u, header)
         ## Run RM-Synthesis
         os.system('python '+rmsyn_path+' -s '+out_dir+'rmsynthesis.par >> '+out_dir+'rmsyn_tmp/out/rmsynthesis.log')
         ## Finally!! Now we can plot the RM-Spectrum (Amplitude)
         header = pyfits.getheader(out_dir+'rmsyn_tmp/out/rm_out_clean_p.fits')
         polf = pyfits.getdata(out_dir+'rmsyn_tmp/out/rm_out_clean_p.fits')[:,0,0]
         fd = np.arange(header['CRVAL3'], header['CRVAL3']+header['CDELT3']*(header['NAXIS3']), header['CDELT3'])
         axarray[plot_no].plot(fd, polf, color='k', lw=1.5)
         axarray[plot_no].set_xlabel(r'$\phi$ (rad m$^{-2}$)')
         axarray[plot_no].set_ylabel(r'$F(\phi)$ (Jy RMTF$^{-1}$)')
         ## Clean up the files
         os.system('rm -rf '+out_dir+'rmsyn_tmp/fits/stokes_q/*')
         os.system('rm -rf '+out_dir+'rmsyn_tmp/fits/stokes_u/*')
         os.system('rm -rf '+out_dir+'rmsyn_tmp/out/*')
         plot_no += 1
         pan_x += 1
         if pan_x == plot_column: ## If row is finished --> Move to next row
            pan_x = 0
            pan_y += 1
      if plot_polp:
         axarray.append(plt.subplot(gs[pan_y, pan_x], projection='polar'))
         if polp_mode == 'cent':
            ## The "central colour' will be mplcm.get_cmap(plot_cmap)(0.5)
            ## Will only plot angle, no amplitude
            axarray[plot_no].plot([pa[len(pa)/2]-np.pi, pa[len(pa)/2]], [1., 1.], color=mplcm.get_cmap(plot_cmap)(0.5), lw=1.5)
         elif polp_mode == 'three':
            axarray[plot_no].plot([pa[0]-np.pi, pa[0]], [1., 1.], color=mplcm.get_cmap(plot_cmap)(0.9999), lw=1.5)
            axarray[plot_no].plot([pa[len(pa)/2]-np.pi, pa[len(pa)/2]], [1., 1.], color=mplcm.get_cmap(plot_cmap)(0.5), lw=1.5)
            axarray[plot_no].plot([pa[-1]-np.pi, pa[-1]], [1., 1.], color=mplcm.get_cmap(plot_cmap)(0.0001), lw=1.5)
         elif polp_mode == 'cont':
            for i in range(len(pa)):
               axarray[plot_no].plot([pa[i]-np.pi, pa[i]], [1., 1.], color=mplcm.get_cmap(plot_cmap)(0.9999-i*0.9998/(len(pa)-1.)), lw=1.5)
         else:
            raise Exception("Parameter <polp_mode> input not recognised! Please set it to one out of ['cent', 'three', 'cont'].")
         axarray[plot_no].set_rticks([]) ## Remove ticks
         axarray[plot_no].set_thetagrids([0, 45, 90, 135, 180, 225, 270, 315], frac=1.25)
         axarray[plot_no].set_theta_zero_location('N')
         plot_no += 1
         pan_x += 1
         if pan_x == plot_column: ## If row is finished --> Move to next row
            pan_x = 0
            pan_y += 1
      #plt.tight_layout()
      #plt.subplots_adjust(left=0.11, bottom=0.10, right=0.96, top=0.97, wspace=0.54, hspace=0.41)
      plt_adj()
      out_no += 1
      if save_file:
         ## The output file will have a file number attached, so it can be easily sorted
         ## Add in "0"s in front, if necessary, to facilitate sorting
         out_no_str = '0'*(np.int(np.ceil(np.log10(nplots)))-np.int(np.ceil(np.log10(out_no+1))))+str(out_no)
         plt.savefig(out_dir+'plot_'+out_no_str+file_ext)
      else:
         plt.show()

## Cleaning up...
if plot_spec:
   os.system('rm -rf '+out_dir+'rmsynthesis.par')
   os.system('rm -rf '+out_dir+'rmsyn_tmp')


