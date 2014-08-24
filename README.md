Monte Carlo Completeness Test
=============================

Python routine to run Monte Carlo completeness test to determine star detection efficiency in CCD images.

- Author:       Navtej Singh
- Contact:      reachnavtej@gmail.com
- Web Site:     http://astro.nuigalway.ie/staff/navtejs
- Organization: CfA@NUIG <http://astro.nuigalway.ie>

This routine was coded as part of research paper "Parallel astronomical data 
processing with Python: Recipes for multicore machines", published in Astronomy 
and Computing. Astro-ph link: http://arxiv.org/abs/1306.0573.

Thanks for downloading completeness code. It include two different versions -
complete_serial.py for serial mode and complete_multi.py for multicore mode execution.
complete_multi.py can also be run in serial mode by setting number of cores to 1.


- Following requirements should be met to run the code.
    
        + A Python 2.4/2.5/2.6/2.7 distribution.
        
        + IRAF software from National Optical Astronomy Observatory (NOAO).
          Download from iraf.noao.edu.
        
        + pyraf python module. Download from
          http://www.stsci.edu/resources/software_hardware/pyraf
          
        + numpy python module. Download from http://numpy.scipy.org/.
        
        + Python multiprocessing module (already in Python >= 2.6). Download
          backport for python 2.4 and 2.5 from
          http://pypi.python.org/pypi/multiprocessing/  
          
        + Parallel python module if running complete_pp.py routine. It can be
          downloaded from www.parallelpython.com  

 ESO's Scisoft include IRAF, pyraf and numpy. Download ESO's Scisoft from ftp://ftp.eso.org/scisoft/scisoft7.7/linux/fedora11/.


- HST WFPC2 test image of M71 (NGC 6838) galactic globular cluster and spatially varying TINY TIM analytical 
  point spread function is included in the archived file. Routine is controlled by the parameters in 
  the configuration file complete.cfg.

         
- Example routine execution:
        $python complete_serial.py data/in.fits[200:711,200:711] data/psf.fits 16 26 50 100
        $python complete_multi.py data/in.fits[200:711,200:711] data/psf.fits 16 26 50 100
        $python complete_pp.py data/in.fits[200:711,200:711] data/psf.fits 16 26 50 100        
         

- Command line options -
    
        $python complete_serial.py --help
        $python complete_multi.py --help
        $python complete_pp.py --help
        

- Program uses IRAF's addstar, daofind and tmatch tasks. All the three tasks
  require input parameters to correctly process the images. 'complete.cfg' is
  configuration file for the program and IRAF task parameters can be set in this
  file. The most important parameters are 'datapars' and 'daopars' - they will
  depend on the input image and psf.
  
  If matplotlib is available, program will plot the completeness graph. Its
  parameters can also be set under [plotter].