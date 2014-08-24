#!/usr/bin/env python

'''
    ------------------------------------------------------------------------    
                STAR DETECTION COMPLETENESS TEST - SERIAL MODE

    Monte Carlo simulation of star detection completeness test. This routine
    take an image and its corresponding psf as input and run iterative tests
    on different star magnitude intervals. Output is a text with magnitude
    interval and detection efficiency. It also plots the completeness test
    graph, if python module matplotlib is available.
    
    Usage: complete_multi.py [options] image psf minmag maxmag nstar niter


    Options:
      --version             show program's version number and exit
      -h, --help            show this help message and exit
      -v, --verbose         print result messages to stdout
      -q, --quiet           don't print result messages to stdout
      -w MAGWIDTH, --magwidth=MAGWIDTH
                            magnitude interval width [default is 0.5 mag]
      -f FILE, --filename=FILE
                            output file name
      -t TEMP, --temp=TEMP  delete temporary files? [default is y]                      

    
    Outputs:
        *.dat: Output file with magnitude and normalized detection efficiency
        *.png: Output plot [Optional: if matplotlib module is available]


    Author:
        Navtej Singh


    Organization:
        Centre for Astronomy, National University of Ireland, Galway, Ireland


    Version:
         20 February 2012    1.0     Original version
    ------------------------------------------------------------------------    
'''

# Load python modules to be used in the routine
import os, sys, subprocess, random
from StringIO import StringIO
from optparse import OptionParser
from ConfigParser import SafeConfigParser


# Check if python module numpy is present
try:
    import numpy as np
except:
    print >> sys.stderr, 'Error: Python module numpy is required. Please install it and try again. Exiting.'
    sys.exit(-1)


# Read complete.cfg configuration file
# ====================================
def readConfig():
    # Read default configuration file
    parser = SafeConfigParser()

    if not parser.read( 'complete.cfg' ):
          print >> sys.stderr, 'Error: Not able to open deconvolve.cfg configuraton file. Exiting.'
          sys.exit( -1 )

    return parser


# Load relevant IRAF packages
# ===========================
def loadPackages():
    iraf.noao(_doprint = 0)
    iraf.noao.digiphot(_doprint = 0)
    iraf.noao.digiphot.daophot(_doprint = 0)    


# Set datapars and daopars parameters    
# ===================================
def setParams():
    parser = readConfig()
    
    # Set datapars and daopars
    for option in parser.options('datapars'):
        iraf.datapars.setParam(option + '.p_value', parser.get('datapars', option))

    for option in parser.options( 'daopars' ):
        iraf.daopars.setParam(option + '.p_value', parser.get('daopars', option))
        

# Plot completeness results if matplotlib is available and 'plot' config variable is yes
def plotter(resfile):
    # Import python matplotlib module   
    try:
        import matplotlib.pylab as plt
    except:
        print >> sys.stderr, '\n Info: Python matplotlib module not found. Skipping plotting.'
        
        return None
  
    else:
        # Open inpur results file  
        try:
            ifile = open(resfile, 'r')
        except:
            print >> sys.stderr, 'Error: Not able to open result file ', resfile
            sys.exit(-1)

        # Read data from the input file    
        idata = ifile.readlines()

        # Close the input file
        try:
              ifile.close()
        except:
              print >> sys.stderr, 'Warning: Not able to close input file ', resfile

        # Read configuration file      
        parser = readConfig()

        # Create and populate python lists
        x, y = [], []
        for value in idata:
            if value[0] != '#':
                  try:
                        value.split()[2]
                  except IndexError:
                        pass
                  else:
                        x.append(value.split()[0])
                        y.append(value.split()[2])

        # Set graph parameters and plot the completeness graph       
        graph = os.path.splitext(resfile)[0] + '.' + parser.get('plotter', 'save_format')

        params = {'backend': 'ps',
		  'font.size': 10,
                  'axes.labelweight': 'medium',
	          'dpi' : 300,
                  'savefig.dpi': 300}
	plt.rcParams.update(params)

        fig = plt.figure()
        plt.title(parser.get('plotter', 'title'), fontweight = 'bold', fontsize = 12)
        plt.xlabel(parser.get('plotter', 'xlabel'))
        plt.ylabel(parser.get('plotter', 'ylabel'))
        plt.axis([float(min(x)) - 0.5, float(max(x)) + 0.5, 0.0, 110])
        plt.grid(parser.get('plotter', 'grid') , linestyle = '-', color = '0.75')
        plt.plot(x, y, parser.get('plotter', 'style'))
        fig.savefig(graph)
    
        return graph


# Create dictionary of magnitude and count. Dictionary is initiallized with 
# input and output counts for each magnitude interval to zero.
def magcntDict(minmag, maxmag, magwidth):
      # Create magnitude intervals between min and max magnitude
      iter = np.arange(minmag, maxmag , magwidth)

      # Create magnitude-count dictionary to be populated by
      # [magnitude, input count, output count]
      magcnt = {}
      for value in iter:
            magcnt[str(round(value, 2))] = [0, 0, 0]

      # Return the initialized python list      
      return magcnt
      

# Main routine to determine detection efficiency. Artificial stars are added to
# the image and IRAF tasks daofind, phot and tjoin are used to detect and match
# the input star list with output list.
def monteCarlo(image, psf, minmag, maxmag, magwidth, nstar, sid, niter, temp):
    print >> sys.stdout, '\n -----> Magnitude : [', minmag, ' - ', maxmag, ']'
    print >> sys.stdout, '\n Adding artificial stars...'
    
    # Read configuration file
    parser = readConfig()
    
    # Use addstar to add stars to the image randomly between the input magnitudes and create nimage number of images
    addimage_prefix = image.rsplit('[')[0].replace('.fits', '.' + str(sid))
    subprocess.call('rm -fr ' + addimage_prefix + '*', shell = True)

    seed = random.uniform(1, 2 * niter)
    
    iraf.unlearn('addstar')
    iraf.addstar(image, photfile = '', psfimage = psf, addimage = addimage_prefix, minmag = minmag, maxmag = maxmag, nstar = nstar, seed = seed, nimage = 1, update = parser.get('addstar', 'update'), verify = parser.get('addstar', 'verify'), mode = parser.get('addstar', 'mode'))
    iraf.flprc()

    print >> sys.stdout, '\n Detecting added artificial stars...'
    
    # Star detection on the image created by addstar
    addimage = addimage_prefix + '.fits'
    addimage_in_coords = addimage_prefix + '.art'

    # Detect stars using IRAF daofind task        
    out_coord_file = addimage.replace('.fits', '.coo.1')
    subprocess.call('rm -f ' + out_coord_file, shell = True)

    iraf.unlearn('daofind')
    iraf.daofind(addimage, output = out_coord_file, boundary = parser.get('daofind', 'boundary'), threshold = parser.get('daofind', 'threshold'), nsigma = parser.get('daofind', 'nsigma'), cache = parser.get('daofind', 'cache'), verify = parser.get('daofind', 'verify'))

    # Determine magnitude of the detected stars
    out_mag_file = addimage.replace('.fits', '.mag.1')
    subprocess.call('rm -f ' + out_mag_file, shell = True)

    iraf.unlearn('phot')
    iraf.phot(addimage, coords = out_coord_file, output = out_mag_file, calgorithm = parser.get('centerpars', 'calgorithm'), salgorithm = parser.get('fitskypars', 'salgorithm'), annulus = parser.get('fitskypars','annulus'), dannulus = parser.get('fitskypars','dannulus'), apertures = parser.get('photpars','apertures'), zmag = parser.get('photpars','zmag'), cache = parser.get('phot', 'cache'), interactive = parser.get('phot', 'interactive'), verify = parser.get('phot', 'verify'), update = parser.get('phot', 'update'))

    # Calculate number of lines in the output photometry file. If zero, skip the next step
    line_cnt = int(subprocess.Popen('cat ' +  out_mag_file + ' | wc -l', shell = True, stdout = subprocess.PIPE).stdout.readline()) - int(subprocess.Popen('cat ' +  out_coord_file + ' | grep "#" | wc -l', shell = True, stdout = subprocess.PIPE).stdout.readline())

    if line_cnt > 0:
          out_mag_tab = out_mag_file + '.tab'
          subprocess.call('rm -f ' + out_mag_tab, shell = True)
          iraf.pconvert(out_mag_file, out_mag_tab, '*')

           # Sort the input tables based on xcenter, ycenter and magnitude
          iraf.tsort(addimage_in_coords, 'c2, c3, c4')
          iraf.tsort(out_mag_tab, 'XCENTER, YCENTER, MAG')

          # Match the input and output list
          matched_tab = addimage.replace('.fits', '.match.list.1')
          subprocess.call('rm -f ' + matched_tab, shell = True)

          iraf.tjoin(addimage_in_coords, out_mag_tab, matched_tab, 'c2, c3, c4', 'XCENTER, YCENTER, MAG', tolerance = parser.get('tjoin', 'xtol') + ',' + parser.get('tjoin', 'ytol') + ',' + parser.get('tjoin', 'magtol'))

          # Determine input and output star count and populate the python list
          inmaglst = iraf.tdump(addimage_in_coords, col = 'c4', Stdout = 1)
          outmaglst = iraf.tdump(matched_tab, col = 'c4', Stdout = 1)
          
          iter = np.arange(minmag, maxmag , magwidth)

          magcnt = {}
          for i in range(len(iter)):
                magcnt[str(round(iter[i], 2))] = [0, 0]

          for value in inmaglst:
                for key in magcnt.keys():
                      try:
                            mag = float(value)
                      except:
                            pass
                      else:
                            if mag >= float(key) and mag < float(key) + magwidth:
                                   magcnt[key][0] += 1
          
          for value in outmaglst:
                for key in magcnt.keys():
                      try:
                            mag = float(value)
                      except:
                            pass
                      else:
                            if mag >= float(key) and mag < float(key) + magwidth:
                                  magcnt[key][1] += 1


    # Delete all the temporary files
    if temp == 'y':
          subprocess.call('rm -f ' + addimage_prefix + '*', shell = True)

    
    # Returning magnitude and match file prefix
    return magcnt


# Worker function for pool/map multicore processing 
# =================================================
def worker(data):
    sid, image, psf, minmag, maxmag, magwidth, nstar, niter, temp = data
    return monteCarlo(image, psf, minmag, maxmag, magwidth, nstar, sid, niter, temp)


# Completeness function
# =====================
def completeness(image, psf, minmag, maxmag, nstar, niter, magwidth, outfile, temp):
    # Load IRAF packages
    loadPackages()
    
    # Set datapars and daopars parameters
    setParams()
    
    # Create python list to map to worker method
    idata = []
    for value in range(niter):  
          idata.append((value, image, psf, minmag, maxmag, magwidth, nstar, niter, temp))  

    # Receive the results in python list
    results = map(worker, idata)    

    # Define a magnitude count dictionary            
    statslst = magcntDict(minmag, maxmag, magwidth)

    # Calculate total input and output count for each magnitude interval
    for key in statslst.keys():
          statslst[key][0] = float(key) + (magwidth / 2)
          for result in results:
                if key in result.keys():
                      statslst[key][1] += result[key][0]
                      statslst[key][2] += result[key][1]
    
    # Write results to output file
    if not outfile:
        outfile = image.rsplit('[')[0].replace('.fits', '.result.dat')

    # Open the output file
    try:
        ofile = open(outfile, 'w')
    except:
        print >> sys.stderr, 'Error: Not able to open output file ', outfile, ' to write results. Exiting.'
        sys.exit(-1)
        
    ofile.write('#-------------------------------------------------------------------------------\n')
    ofile.write('#                      Star Completeness Test Result                            \n')
    ofile.write('#                    =================================                          \n')
    ofile.write('# Magnitude   Stars Detected/Input Stars    Completeness    Detection Efficiency\n')
    ofile.write('#  (mag)              (--)                   (%)                (--)            \n')
    ofile.write('#-------------------------------------------------------------------------------\n')


    # Sort the dictionary on basis of keys and write results to output file
    keylist = statslst.keys()
    keylist.sort()
    for key in keylist:
          if statslst[key][1] != 0:
                ofile.write('%10s%10s%s%6.2f%s%4.2f%s' %(str(statslst[key][0]), str(statslst[key][2]) + '/' + str(statslst[key][1]),  '\t\t', (float(statslst[key][2])/statslst[key][1]) * 100, '\t\t', float(statslst[key][2])/statslst[key][1], '\n'))
                #ofile.write( str( statslst[key][0] ) + '\t\t' + str( statslst[key][2] ) + '/' + str( statslst[key][1] ) + '\t\t' + str( ( float(statslst[key][2])/statslst[key][1] ) * 100 ) + '\t\t' + str( float(statslst[key][2])/statslst[key][1] ) + '\n' )
          else:
                ofile.write(str(statslst[key][0]) + '\t\t' + str(statslst[key][2]) + '/' + str(statslst[key][1]) + '\n')

    # Close the output file
    try:            
          ofile.close()
    except:
          print >> sys.stderr, 'Warning: Not able to close the output file ', outfile


    print >> sys.stdout, '\n Results written to - ', outfile
    
    # If plot_flag is True, plot the completeness graph using matplotlib module
    parser = readConfig()
    
    if parser.get('plotter', 'plot') == 'yes':
       graph = plotter( outfile )
       print >> sys.stdout, '\n Completeness graph - ', graph
    


# Main function for the routine
# =============================
def main(image, psf, minmag, maxmag, nstar, niter, magwidth = 0.5, outfile = None, temp = 'y'):
    # Check if the image exists
    if not os.path.exists(image.rsplit('[')[0]):
        print >> sys.stderr, 'Error: Image ', image, ' does not exist. Exiting.'
        sys.exit(-1)
    
    # Check if the psf exists
    if not os.path.exists(psf):
        print >> sys.stderr, 'Error: PSF ', psf,'does not exist. Exiting.'
        sys.exit(-1)
        
    # Check if the input values are valid
    try:
        minmag = float(minmag) 
        maxmag = float(maxmag)
        nstar  = int(nstar)
        niter = int(niter)
        magwidth = float(magwidth)
    except:
        print >> sys.stderr, 'Error: minmag, maxmag, nstar or niter not a number. Exiting.'
        sys.exit(-1)
    
    completeness(image, psf, minmag, maxmag, nstar, niter, magwidth, outfile, temp) 


# Entry point for python script
# =============================
if __name__ == '__main__':
    # Usage: python complete.py image psf maxmag minmag nstar niter
    usage = "Usage: python %prog [options] image psf minmag maxmag nstar niter"
    description = "Description. Utility to run completeness test on images in multiprocessing mode."
    parser = OptionParser(usage = usage, version = "%prog 1.0", description = description)
    parser.add_option("-v", "--verbose",
                      action="store_true", dest="verbose", default = False,
                      help = "print result messages to stdout"
                      )
    parser.add_option("-q", "--quiet",
                    action="store_false", dest="verbose", default = True,
                    help = "don't print result messages to stdout"
                    )
    parser.add_option("-w", '--magwidth', dest = "magwidth", metavar="MAGWIDTH", 
                    action="store", help = "magnitude interval width [default is 0.5 mag]",
                    default = 0.5
                    )
    parser.add_option("-f", "--filename", dest = "filename",
                    action='store', metavar="FILE", help = "output file name"
                    )
    parser.add_option("-t", "--temp", dest = "temp", metavar="TEMP",
                    action="store", help = "delete temporary files? [default is y]",
                    choices=['y', 'n'], default = 'y'
                    )
    (options, args) = parser.parse_args()

    # Check for number of input arguments
    if len(args) != 6:
        parser.error("Incorrect number of arguments. Check help for further details.")

    print >> sys.stdout, '\n Starting processing...'
    
    # Check verbosity
    if not options.verbose:
        output = StringIO()
        old_stdout = sys.stdout
        sys.stdout = output

    # Check if pyraf module is installed
    try:
        from pyraf import iraf
    except:
        print >> sys.stderr, 'Error: Python module pyraf not found. Exiting.'
        exit(-1)
    
    main(args[0], args[1], args[2], args[3], args[4], args[5], options.magwidth, options.filename, options.temp)
    
    # Reset verbosity
    if not options.verbose:
        sys.stdout = old_stdout
    
    print >> sys.stdout, '\n Process completed successfully.'