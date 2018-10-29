# Faraday-Rotation-Illustration  
A Python script that illustrates the effect of Faraday rotation as it traverses through foreground magneto-ionic medium

Watch a sample output on: https://www.youtube.com/watch?v=0zzmx793CxE

!! NOTE !!  
In order to generate Faraday spectra, the RM-Synthesis code by Michael Bell (https://github.com/mrbell/pyrmsynth) is needed.   
Pyfits is also required to create the dummy .FITS files. These are not necessary if RM-Synthesis is not performed


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
