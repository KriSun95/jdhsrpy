"""
Functions to create proxy AIA channels.
"""

import os
import sys
import sunpy
import sunpy.map

from aiapy.calibrate import (
    register,
    update_pointing,
)
from aiapy.calibrate.utils import get_correction_table, get_pointing_table
import astropy.units as u
from astropy.time import Time
from aiapy.calibrate import degradation
import numpy as np


def prep(smap, table):
    return register(update_pointing(smap, pointing_table=table))

def get_unprepped_file_time(filename_time):
    return Time.strptime(filename_time, "%Y_%m_%dT%H_%M_%S.%fZ")

def create_iron18(dir_094, dir_171, dir_211, outdir, needing_prepped=False, TABLE=None):
    """Takes the 94, 171, 211 channels from SDO/AIA to create an iron18 emission proxy (Del Zanna 2013).
    
    Parameters
    ----------
    dir_*** : str
            The string of the directory with the 94A (dir_094), 171A (dir_171), and 211A (dir_211) files.
            
    outdir : str
            Directory for the output iron18 files.

    needing_prepped : bool
            Is the AIA data in dir_*** needing prepped?
            
    Returns
    -------
    Filenames of the iron 18 files.
    """

    os.makedirs(outdir, exist_ok=True)

    files_094_list = list(os.listdir(dir_094))
    files_094 = [ f for f in files_094_list if f.endswith('.fits')]
    files_094.sort()

    files_171_list = list(os.listdir(dir_171))
    files_171 = [ f for f in files_171_list if f.endswith('.fits')]
    files_171.sort()

    files_211_list = list(os.listdir(dir_211))
    files_211 = [ f for f in files_211_list if f.endswith('.fits')]
    files_211.sort()
    
    co_094 = []
    co_171 = []
    co_211 = []
    output = []
    # aia.lev1.94A_2024_04_17T22_18_23.12Z.image_lev1.fits
    for fn094 in files_094:
        if needing_prepped:
            time_094 = get_unprepped_file_time(fn094[-39:-16])# 2024_04_17T22_18_23.12Z
        else:
            time_094 = Time.strptime(fn094[-39:-16], "%Y_%m_%dT%H_%M_%S.%fZ") # replace when he structure is known

        for fn171 in files_171:
            if needing_prepped:
                time_171 = get_unprepped_file_time(fn171[-39:-16])# 2024_04_17T22_18_23.12Z
            else:
                time_171 = Time.strptime(fn171[-39:-16], "%Y_%m_%dT%H_%M_%S.%fZ") #datetime.datetime.strptime(fn171[3:18], '%Y%m%d_%H%M%S')
            
            if time_094 <= time_171 < time_094 + (12<<u.s):

                for fn211 in files_211:     
                    if needing_prepped:
                        time_211 = get_unprepped_file_time(fn211[-39:-16])# 2024_04_17T22_18_23.12Z
                    else:
                        time_211 = Time.strptime(fn211[-39:-16], "%Y_%m_%dT%H_%M_%S.%fZ") #datetime.datetime.strptime(fn211[3:18], '%Y%m%d_%H%M%S')
                    
                    if time_094 <= time_211 < time_094 + (12<<u.s):
                        co_094.append(fn094)
                        co_171.append(fn171)
                        co_211.append(fn211)
                        break
                break
    files_094 = co_094
    files_171 = co_171
    files_211 = co_211

    _first_map = sunpy.map.Map(os.path.join(dir_094, files_094[0]))
    degs = [degradation(94*u.angstrom, Time(_first_map.meta["T_OBS"], format="isot", scale="utc"), correction_table=get_correction_table("jsoc")), 
            degradation(171*u.angstrom, Time(_first_map.meta["T_OBS"], format="isot", scale="utc"), correction_table=get_correction_table("jsoc")), 
            degradation(211*u.angstrom, Time(_first_map.meta["T_OBS"], format="isot", scale="utc"), correction_table=get_correction_table("jsoc"))]
    del _first_map
    
    d_total = len(files_094)
    for d, (f094, f171, f211) in enumerate(zip(files_094, files_171, files_211)):
        aia_map_094 = sunpy.map.Map(os.path.join(dir_094, f094))

        if needing_prepped:
            if TABLE is None:
                TABLE = get_pointing_table("JSOC", 
                                time_range=[aia_map_094.date - 1 * u.day, 
                                            aia_map_094.date + 1 * u.day])
            aia_map_094 = prep(aia_map_094, TABLE)

        data_094 = aia_map_094.data / aia_map_094.exposure_time
        data_094[data_094 < 0] = 0
        
        aia_map_171 = sunpy.map.Map(os.path.join(dir_171, f171))

        if needing_prepped:
            aia_map_171 = prep(aia_map_171, TABLE)

        data_171 = aia_map_171.data / aia_map_171.exposure_time
        data_171[data_171 < 0] = 0
        
        aia_map_211 = sunpy.map.Map(os.path.join(dir_211, f211))

        if needing_prepped:
            aia_map_211 = prep(aia_map_211, TABLE)
        
        data_211 = aia_map_211.data / aia_map_211.exposure_time
        data_211[data_211 < 0] = 0

        dim0, dim1 = 4096, 4096
        shape_094 = np.shape(data_094)
        dim0delta, dim1delta = dim0-shape_094[0], dim1-shape_094[1]
        data_094 = np.pad(data_094, [(int(dim0delta/2), int(dim0delta/2)), (int(dim1delta/2), int(dim1delta/2))], mode='constant')
        shape_171 = np.shape(data_171)
        dim0delta, dim1delta = dim0-shape_171[0], dim1-shape_171[1]
        data_171 = np.pad(data_171, [(int(dim0delta/2), int(dim0delta/2)), (int(dim1delta/2), int(dim1delta/2))], mode='constant')
        shape_211 = np.shape(data_211)
        dim0delta, dim1delta = dim0-shape_211[0], dim1-shape_211[1]
        data_211 = np.pad(data_211, [(int(dim0delta/2), int(dim0delta/2)), (int(dim1delta/2), int(dim1delta/2))], mode='constant')

        iron_18 = data_094/degs[0] - data_211/(120*degs[2]) - data_171/(450*degs[1])
        iron_18[iron_18 < 0] = 0
        aia_map_fe18 = sunpy.map.Map(iron_18, aia_map_094.meta)
        aia_map_fe18.save(os.path.join(outdir, "fe18".join(f094.split("94A"))))
        
        del aia_map_094
        del aia_map_171
        del aia_map_211
        del aia_map_fe18

        print(f'\r[function: {sys._getframe().f_code.co_name}] Saved {d} submap(s) of {d_total}.        ', end='') 

    return output
