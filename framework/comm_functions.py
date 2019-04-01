from termcolor import colored
### Communication functions
def info(string):
    #t = "%d/%d/%d %d:%d:%d"%(time.localtime()[:6])
    print colored("\n>>> MEQSILHOUETTE INFO <<< : %s\n"%(string),'green')

def warn(string):
    #t = "%d/%d/%d %d:%d:%d"%(time.localtime()[:6])
    print colored("\n>> MEQSILHOUETTE WARNING <<: %s\n"%(string),'yellow')

def abort(string,exception=SystemExit):
    #t = "%d/%d/%d %d:%d:%d"%(time.localtime()[:6])
    raise exception(colored("\n>> MEQSILHOUETTE ABORTING <<: %s\n"%(string),'red'))

def print_simulation_summary(ms_dict,im_dict):
    info('MS relevant parameters:')
    print('RA, DEC \t\t = %s deg, %s deg'\
          %(str(ms_dict['RA']),str(ms_dict['DEC'])))
    print('central freq\t\t = %.1f GHz'%ms_dict['nu'])
    print('bandwidth\t\t = %.1f MHz'%(ms_dict['dnu']*1e3))
    print('num channels\t\t = %i'%ms_dict['nchan'])
    print('channel width\t\t = %.3f MHz'%(ms_dict['dnu']*1e3 / ms_dict['nchan']))
    print('observation length\t = %.1f hours'%ms_dict['obslength'])
    print('correlator dump time\t = %.1f seconds'%(ms_dict['tint']))
    print('obs start date/time\t = %s'%ms_dict['StartTime'])

    info('IMAGING relevant parameters:')
    print('output image npixels\t = %i'%im_dict['npix'])
    print('output image cellsize\t = %s'%im_dict['cellsize'])
    print('Stokes parameter\t = %s'%im_dict['stokes'])
    print('image weighting\t\t = %s'%im_dict['weight'])

