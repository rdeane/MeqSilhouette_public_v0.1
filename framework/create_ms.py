# coding: utf-8
from Pyxis.ModSupport import *
import pyrap.tables as pt
import os
from framework.comm_functions import *



def create_ms(msname, input_fits, ms_dict):
    """ create empty MS """
    x.sh(return_simms_string(msname, input_fits, **ms_dict))

    # set STATION names to same as NAME col
    antt = pt.table(os.path.join(msname, 'ANTENNA'), readonly=False,ack=False)
    antt.putcol('STATION', antt.getcol('NAME'))
    antt.close()

    # set FIELD name to input filename (minus extension)
    fieldtab = pt.table(os.path.join(msname,'FIELD'),readonly=False,ack=False)
    fieldtab.putcol('NAME',input_fits.split('/')[-1].split('.')[0])
    fieldtab.close()

    # set SOURCE name to input filename (minus extension)
    srctab = pt.table(os.path.join(msname,'SOURCE'),readonly=False,ack=False)
    srctab.putcol('NAME',input_fits.split('/')[-1].split('.')[0])
    srctab.close()

    # set SPW name to input filename (minus extension)
    spwtab = pt.table(os.path.join(msname,'SPECTRAL_WINDOW'),readonly=False,ack=False)
    spwtab.putcol('NAME',input_fits.split('/')[-1].split('.')[0])
    spwtab.close()

    # INI: Add WEIGHT_SPECTRUM and SIGMA_SPECTRUM columns to the MS
    tab = pt.table(msname,readonly=False)
    data = tab.getcol('DATA')
    tab.addcols(pt.makearrcoldesc('SIGMA_SPECTRUM',value=1.0,shape=[data.shape[1],data.shape[2]],valuetype='float'))
    tab.addcols(pt.makearrcoldesc('WEIGHT_SPECTRUM',value=1.0,shape=[data.shape[1],data.shape[2]],valuetype='float'))
    tab.close()

    info('Measurement Set %s created '%msname)


def return_simms_string(msname, input_fits, RA, DEC, polproducts, antenna_table, \
                        obslength, dnu, tint, nu, StartTime, nchan, \
                        nscan, scan_lag,datacolumn, makeplots):
    
    s = "simms -T VLBA -t casa -n %s -ra %.9fdeg -dec %.9fdeg \
-pl '%s' -st %f -sl %f -slg %f -dt %f -f0 %fGHz -df %fGHz -nc %i  -date %s %s\
    " % ( msname, RA, DEC, polproducts, obslength, obslength/float(nscan), scan_lag,
         tint, nu - (dnu/2.) + (dnu/(float(nchan))/2.), dnu/float(nchan),
          nchan, StartTime, os.path.join(II('$CODEDIR'),antenna_table)) 

    return s

    # NOTE: CASA fails if simms -T option is changed from "VLBA",
    #so this is hard-coded in.
    # It fails as it does not have a reference location for the EHT
    # referencelocation=obs_pos gives error:
    #Error Illegal Measure record in MeasureHolder::fromRecord
    #In converting Position parameter




