# coding: utf-8

# In[ ]:
from Pyxis.ModSupport import *
import mqt
import numpy as np
import pyrap.tables as pt
from im import lwimager 
import subprocess
import glob
from framework.comm_functions import *

def run_wsclean(input_fitsimage,input_fitspol,startvis,endvis):
    msname = II('$MS')

    if input_fitspol == 0:
        subprocess.check_call(['wsclean','-predict','-name',input_fitsimage,'-interval',str(int(startvis)),str(int(endvis)),msname])
    else:
        subprocess.check_call(["wsclean","-predict","-name",input_fitsimage,"-interval",str(int(startvis)),str(int(endvis)),"-pol","I,Q,U,V","-no-reorder",msname])

def copy_between_cols(dest_col, src_col):
    msname = II('$MS')

    tab=pt.table(msname, readonly=False)
    src_data = tab.getcol(src_col)
    tab.putcol(dest_col, src_data)
    tab.close()
        
def run_turbosim(input_fitsimage,output_column,taql_string):

    options = {}
    options['ms_sel.msname'] = II('$MS')
    options['ms_sel.output_column'] = output_column
    if input_fitsimage.endswith(('.fits','.FITS')):
        options['me.sky.siamese_oms_fitsimage_sky'] = 1
        options['fitsimage_sky.image_filename'] = input_fitsimage
        options['fitsimage_sky.pad_factor'] = 2.4
    elif input_fitsimage.endswith(('.txt','.html')):
        options['me.sky.tiggerskymodel'] = 1
        options['tiggerlsm.filename'] = input_fitsimage
    options['ms_sel.tile_size'] = 1000000
    options['ms_sel.ms_taql_str'] = taql_string



    
    mqt.MULTITHREAD = 32 #max number of meqserver threads
    mqt.run(script=II('$FRAMEWORKDIR')+'/turbo-sim.py',
            config=II('$FRAMEWORKDIR')+'/tdlconf.profiles',
            section='turbo-sim',
            job='_simulate_MS',
            options=options)
# removed:                 

def add_pjones(output_column):
    options = {}
    options['ms_sel.msname'] = II('$MS')
    options['ms_sel.output_column'] = output_column
    options['read_ms_model'] = 1 # read existing visibilities from MS
    options['ms_sel.model_column'] = output_column
    options['ms_sel.tile_size'] = 1000000
    options['sim_mode'] = 'sim only'

    options['me.p_enable'] = 1
    options['feed_angle.enable_pa'] = 1 # enable parallactic angle rotation
    options['feed_angle.read_ms'] = 1 # enable reading feed angle from FEED subtable in MS

    mqt.MULTITHREAD = 32 #max number of meqserver threads
    mqt.run(script=II('$FRAMEWORKDIR')+'/turbo-sim.py',
            config=II('$FRAMEWORKDIR')+'/tdlconf.profiles',
            job='_simulate_MS',
            options=options)

'''def add_uvjones(output_column, gterm, dterm, gainR, gainL, leak_ampl_string, leak_phas_string):
    options = {}
    options['ms_sel.msname'] = II('$MS')
    options['ms_sel.output_column'] = output_column
    options['read_ms_model'] = 1 # read existing visibilities from MS
    options['ms_sel.model_column'] = output_column
    options['ms_sel.tile_size'] = 1000000
    options['sim_mode'] = 'sim only'
    if gterm:
        options['me.g_enable'] = 1
    if dterm:
        options['me.p_enable'] = 1
        options['feed_angle.enable_pa'] = 1 # enable parallactic angle rotation
        options['feed_angle.read_ms'] = 1 # enable reading feed angle from FEED subtable in MS

        options['me.d_enable'] = 1
        options['leakage.err-gain.error_model'] = 'ListOfValues'
        options['leakage.err-gain.values_str'] = leak_ampl_string 
        options['leakage.err-phase.error_model'] = 'ListOfValues'
        options['leakage.err-phase.values_str'] = leak_phas_string 

    mqt.MULTITHREAD = 32 #max number of meqserver threads
    mqt.run(script=II('$FRAMEWORKDIR')+'/turbo-sim.py',
            config=II('$FRAMEWORKDIR')+'/tdlconf.profiles',
            job='_simulate_MS',
            options=options)'''


def make_dirty_image_lwimager(im_dict,ms_dict):
        lwimager.make_image(column=ms_dict["datacolumn"],
            dirty_image=II('${OUTDIR>/}${MS:BASE}')+'-dirty_map.fits',
                            dirty=True,**im_dict)


def make_image_wsclean():
    print('todo')

def make_image_pymoresane():
    print('todo')




