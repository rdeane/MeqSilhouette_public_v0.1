""" adds a new antenna to a CASA antenna table
Usage: python add_ant.py

Currently, the inputs are hardcoded, with the AMT as an example. 

Inputs:
- input antenna table name
- output antenna table name
- new antenna dictionary, includes:
    -- Position: Lat [deg], Long [deg, east +ve], height [metres]
    converts to XYZ assuming WGS84 ellipsoid, formulas from Seeber, Satellite Geodesy (2003 ed, sect.2.1.4)
    
    -- antenna mount type (will become important as MeqSilhouette accounts for this)

    -- DISH_DIAMETER [metres]

    -- etc.

    e.g. newAntDict = {'OFFSET'  : [0,0,0],\
              'POSITION': [5627893.7450699108, 1637770.1989314035, -2512490.219404601],\
              'TYPE'    : 'GROUND-BASED',\
              'DISH_DIAMETER': 12,\
              'FLAG_ROW':  0,\
              'MOUNT'   : 'alt-az',\
              'NAME'    : 'AMT',\
              'STATION' : 'AMT'}

Ouputs: new antenna table

"""

import numpy as np
import pyrap.tables as pt
import os
from framework.comm_functions import *
import sys


def latlonh_2_xyz(lat_deg,lon_deg,h_metres):
    """ convert input Lat,Lon,height to XYZ coordinates, assuming WGS84 geoid
    Input units: deg, deg, metres respectively
    Output units: metres"""
    degrad = np.pi/180.
    #WGS84 
    a = 6378137.
    finv = 298.257223563
    f = 1./finv
    e2 = 2.*f - f**2

    latr = lat_deg*degrad
    slat = np.sin(latr)
    slat2 = slat**2
    clat = np.cos(latr)
 
    lonr = lon_deg*degrad
    clon = np.cos(lonr)
    slon = np.sin(lonr)

    nbar = a / np.sqrt(1. - e2*slat2)

    x = (nbar+h_metres) * clat*clon
    y = (nbar+h_metres) * clat*slon
    z = ((1.-e2)*nbar+h_metres) * slat

    return [x,y,z]



def AddAnt(inputAntTableName,outputAntTableName, newAntDict):
    """ add a single row to existing CASA antenna table.
    Use antenna dict as input, see AMT example below"""
    if (os.path.exists(outputAntTableName)):
        abort('Antenna table %s exists!'%outputAntTableName)
    os.system('cp -r %s %s'%(inputAntTableName,outputAntTableName))
    tab=pt.table(outputAntTableName,readonly=False,ack=True)
    tab.addrows(nrows=1)
    info('new antenna row parameter values:')
    for colname in tab.colnames():
        col=tab.getcol(colname)
        print colname + '='
        print newAntDict[colname]
        print ''
        col[-1] = newAntDict[colname]
        tab.putcol(colname,col)
    tab.close() 
    info("Additional antenna %s added to input antenna table '%s'\n\
          and written out to new antenna table '%s'"\
          %(newAntDict['NAME'],inputAntTableName,outputAntTableName))
    
    
AMT_LatLonh = [ -23. - (20/60.) - (31.90)/3600.,  16 + (13/60.) + (31.78)/3600., 2350]
AMT_xyz = latlonh_2_xyz(AMT_LatLonh[0],AMT_LatLonh[1],AMT_LatLonh[2])
AMTAntDict = {'OFFSET'  : [0,0,0],\
              'POSITION': AMT_xyz,\
              'TYPE'    : 'GROUND-BASED',\
              'DISH_DIAMETER': 15,\
              'FLAG_ROW':  0,\
              'MOUNT'   : 'alt-az',\
              'NAME'    : 'AMT',\
              'STATION' : 'AMT'}
#AddAnt('EHT','EHT_AMT',AMTAntDict)

MeerKAT_LatLonh = [-30 - (42./60) - (47.41/3600),    -21 - (26./60) - (38./3600), 1050]
MeerKAT_xyz = latlonh_2_xyz(MeerKAT_LatLonh[0],MeerKAT_LatLonh[1],MeerKAT_LatLonh[2])
MeerKATAntDict = {'OFFSET'  : [0,0,0],\
              'POSITION': MeerKAT_xyz,\
              'TYPE'    : 'GROUND-BASED',\
              'DISH_DIAMETER': 13.5,\
              'FLAG_ROW':  0,\
              'MOUNT'   : 'alt-az',\
              'NAME'    : 'MK',\
              'STATION' : 'MK'}


APEX_LatLonh = [-23 - (00./60) - (21.0/3600),    -67 - (45./60) - (33./3600), 5100]
APEX_xyz = latlonh_2_xyz(APEX_LatLonh[0],APEX_LatLonh[1],APEX_LatLonh[2])

APEXAntDict = {'OFFSET'  : [0,0,0],\
                  'POSITION': APEX_xyz,\
                  'TYPE'    : 'GROUND-BASED',\
                  'DISH_DIAMETER': 12.0,\
                  'FLAG_ROW':  0,\
                  'MOUNT'   : 'alt-az',\
                  'NAME'    : 'APEX',\
                  'STATION' : 'APEX'}



tabin = sys.argv[1]
tabout = sys.argv[2]

AddAnt(tabin,tabout,AMTAntDict) #MeerKATAntDict)


