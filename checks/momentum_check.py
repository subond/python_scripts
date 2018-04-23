#check momentum budget closure of model
from netCDF4 import Dataset
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from physics import mombudg_2d_an_fn
from data_handling import cell_area

def check_mom_balance(run, years):
    data = mombudg_2d_an_fn(run,'16',years)
    area = cell_area(42,'/scratch/rg419/GFDL_model/GFDLmoistModel/')
    xarea = xr.DataArray(area, [('lat', data.lat), ('lon', data.lon)])
    mom_budg_horiz_av = (data*xarea).sum(('lat','lon'))/xarea.sum(('lat','lon'))
    #mom_budg_sdev =  np.sqrt( (np.square(data - mom_budg_horiz_av)*xarea).sum(('lat','lon'))/xarea.sum(('lat','lon')) )
    #print mom_budg_sdev
    plt.figure()
    mom_budg_horiz_av.dphidx_av.plot()
    mom_budg_horiz_av.fv_av.plot()
    mom_budg_horiz_av.mom_mean.plot()
    mom_budg_horiz_av.mom_eddy.plot()
    mom_budg_horiz_av.ddamp_av.plot()
    mom_budg_horiz_av.mom_sum.plot()
    plt.legend(['dphidx','fv','mmmn','mmed','damp','sum'])
    plt.title(run)
    plt.ylim(-2e-5,2e-5)
    #plt.savefig('/scratch/rg419/plots/momentum_budget/'+run+'.png')
    return mom_budg_horiz_av

def check_mom_vert(run, years, era=True):
    data = mombudg_2d_an_fn(run,'16',years,era=era)
    
    #NB ONLY ERA LEVEL CONSISTENT
    if era==True:
        phalf = [1000., 962.5, 887.5, 775., 650., 550., 450., 350., 275., 225., 175., 125., 85., 60., 40., 25., 15., 0.]
        pdiffs = xr.DataArray(np.diff(phalf), [('pfull', data.pfull )])
    else:
        examplefile = '/scratch/rg419/Data_moist/aquaplanet_10m_diags/np16/run13/atmos_daily.nc'
        exdata = xr.open_dataset( examplefile,decode_times=False)
        pdiffs = -1*xr.DataArray(exdata.phalf.diff('phalf').values, [('pfull', data.pfull )])

    mom_budg_vint = (data * pdiffs*100).sum(('pfull'))/9.8
    mom_budg_zav = mom_budg_vint.mean(('lon'))

    plt.figure()
    mom_budg_zav.dphidx_av.plot()
    mom_budg_zav.fv_av.plot()
    mom_budg_zav.mom_mean.plot()
    mom_budg_zav.mom_eddy.plot()
    mom_budg_zav.ddamp_av.plot()
    mom_budg_zav.mom_sum.plot()
    plt.legend(['dphidx','fv','mmmn','mmed','damp','sum'])
    plt.ylim(-0.2,0.2)
    plt.title(run)
    if era==True:
        plt.savefig('/scratch/rg419/plots/momentum_budget/'+run+'_vint.png')
    else:
        plt.savefig('/scratch/rg419/plots/momentum_budget/'+run+'_vint_model.png')
    return mom_budg_zav

#mb_topo = check_mom_vert('topo_10m_diags',range(12,17))
#mb_flat = check_mom_vert('flat_10m_diags',range(12,17))
mb_ap10 = check_mom_vert('aquaplanet_10m_diags',range(12,17),era=False)
#mb_ap1 = check_mom_vert('aquaplanet_1m_diags',range(12,17))

#mb_topo = check_mom_balance('topo_10m_diags',range(12,17))
#mb_flat = check_mom_balance('flat_10m_diags',range(12,17))
#mb_ap10 = check_mom_balance('aquaplanet_10m_diags',range(12,17))
#mb_ap1 = check_mom_balance('aquaplanet_1m_diags',range(12,17))

plt.show()