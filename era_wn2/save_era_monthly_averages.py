# 11/01/2018 Evaluate streamfunction, plot for JJA with centre by month on top

from data_handling_updates import month_dic, model_constants as mc, gradients as gr
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from pylab import rcParams
import sh
import pandas as pd


def monthly_mean(filename, filename_out):
    data = xr.open_dataset(filename, chunks={'latitude': 100, 'longitude': 100})

    data.load()

    data_mm = data.resample(time='M').mean()
    
    data_mm.to_netcdf(filename_out)


filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_200.nc'
filename_out = '/scratch/rg419/obs_and_reanalysis/sep_levs_u/era_u_200_mm.nc'
monthly_mean(filename, filename_out)


filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_200.nc'
filename_out = '/scratch/rg419/obs_and_reanalysis/sep_levs_v/era_v_200_mm.nc'
monthly_mean(filename, filename_out)

filename = '/scratch/rg419/obs_and_reanalysis/sep_levs_vo/era_vo_200.nc'
filename_out = '/scratch/rg419/obs_and_reanalysis/sep_levs_vo/era_vo_200_mm.nc'
monthly_mean(filename, filename_out)