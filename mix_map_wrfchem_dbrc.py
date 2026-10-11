# This is a python file
# Creat WRFchem anthropogenic emission input file based on the output files from the WRFchem_emiss_multiprocess.py
# This file maps CB05 in MIX to CBMZ_mosaic in WRFchem
# This file can be excuted with python mix_map_wrfchem.py
import os
from netCDF4 import Dataset,num2date,date2num
import numpy as np
import netCDF4 as nc
import pandas as pd
import time 
from datetime import datetime
import copy

def eachFile(filepath):
    '''get the paths (which is dir+filename) of all the files under the directory of filepath

       para
       ------
       filepath: for example(/Users/cwang161/project/)

       return
       ------
       file_path_list: filepath+filename in the type of list
    '''
    pathDir = os.listdir(filepath)
    file_path_list = list()
    for allDir in pathDir:
        child = os.path.join('%s%s' % (filepath, allDir))
        file_path_list.append(child)
    return file_path_list

############################# Info that need to be changed ###################################
# the dir for the reprojected MIX inventory which is also the output file path in WRFchem_emiss_multiprocess.py
flist = eachFile('./emiss_project_d01/')
# t is 24 or 12, 12 for the period 00-11, 24 for the period 13-23.
t = 24
# the domain number 
domain = 1 
nt = 12 # # of times(months) in MIX files
nx = 142 # # of columns in WRFchem domain
ny = 105 # # of rows in WRFchem domain
nlev = 8 # # of levels in WRFchem domain
res = 50 # resolution in WRFchem domain, units:km
month = 1 # which month in MIX inventory that you want to use, 1 for Jan
year = 2017
########################### Main ####################################################################
## Creat the output file name
if t==12:
    fname = 'wrfchemi_00z_d0' + str(domain) + '.nc' + '_' + str(year) + '_' + str(month)
elif t==24:
    fname = 'wrfchemi_12z_d0' + str(domain) + '.nc' + '_' + str(year) + '_' + str(month)
df_new = Dataset(fname, 'w', format='NETCDF4')

## variables in the cbmz_mosaic package in the registry.chem
#cbmz_wrf = 'e_iso,e_so2,e_no,e_co,e_eth,e_hc3,e_hc5,e_hc8,e_xyl,e_ol2,e_olt,e_oli,e_tol,e_csl,e_hcho,e_ald,e_ket,e_ora2,e_nh3,e_no2,e_ch3oh,e_c2h5oh,e_pm25i,e_pm25j,e_pm_10,e_eci,e_ecj,e_orgi,e_orgj,e_so4i,e_so4j,e_no3i,e_no3j,e_so4c,e_no3c,e_orgc,e_ecc,e_brc1i,e_brc1j,e_brc2i,e_brc2j,e_brc3i,e_brc3j'
cbmz_wrf = 'e_co,e_no,e_no2,e_bigalk,e_bigene,e_c2h4,e_c2h5oh,e_c2h6,e_c3h6,e_c3h8,e_ch2o,e_ch3cho,e_ch3coch3,e_ch3oh,e_mek,e_so2,e_toluene,e_benzene,e_xylene,e_nh3,e_isop,e_c10h16,e_pm25i,e_pm25j,e_eci,e_ecj,e_orgi,e_orgj,e_so4i,e_so4j,e_no3i,e_no3j,e_nh4i,e_nh4j,e_nai,e_naj,e_cli,e_clj,e_co_a,e_orgi_a,e_orgj_a,e_co_bb,e_orgi_bb,e_orgj_bb,e_pm_10,e_c2h2,e_gly,e_sulf,e_macr,e_mgly,e_mvk,e_hcooh,e_honoi,e_brc1i,e_brc1j,e_brc2i,e_brc2j,e_brc3i,e_brc3j,e_dbrci,e_dbrcj'
cbmz_wrf =[x.upper() for x in sorted(i[2:] for i in cbmz_wrf.split(','))]

# Species match  (物种匹配)
# Need one-to-one correspondence. For example: ALD =1*(1*ALD2+1*ALDX); NO=0.9*(1*NOx) NO2=0.1*(1*NOx); numbers outside the
# brackets are from coef_wrf; numbers inside the brackets are from coef_mix.
#mix = [['ALD2','ALDX'],['ETH'],['ETHA'],['ETOH'],['FORM'],['IOLE'],['ISOP'],['MEOH'],['OLE'],['TOL'],['XYL'],['SO2'],['NH3'],\
#       ['CO'],['NOx'],['PM2.5','BC','OC'],['PM10','PM2.5'],['BC'],['OC'],['PAR'],['PAR','UNR'],['PAR'],['PAR','ALDX'],['OLE','TOL']]
#wrf = [['ALD'],['OL2'],['ETH'],['C2H5OH'],['HCHO'],['OLI'],['ISO'],['CH3OH'],
#       ['OLT'],['TOL'],['XYL'],['SO2','SO4I','SO4J','SO4C'],['NH3'],['CO'],['NO','NO2'],
#       ['PM25I','PM25J'],['PM_10'],['ECI','ECJ','ECC'],['ORGI','ORGJ','ORGC'],['KET'],['HC3'],['HC5'],['HC8'],['CSL']]
#coef_mix = [[1, 1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1,-1,-1],[1,-1],[1],[1],[3],[2, 1],[3],[1, 1],[0.5, 1]]
#coef_wrf = [[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[0.95,0.05*0.2,0.05*0.8,0],[1],[1],[0.9,0.1],[0.2,0.8],[1],[0.2,0.8,0],[0.2,0.8,0],[1],[1],[1],[1],[1]]

wrf = [['BIGALK'], ['BIGENE'], ['C2H4'], ['C2H5OH'], ['C2H6'], \
       ['C3H6'], ['C3H8'], ['CH2O'], ['CH3CHO'], ['CH3COCH3'], \
       ['CH3OH'], ['TOLUENE'], ['XYLENE'], ['ISOP'], ['C10H16'], \
       ['C2H2'], ['MVK'], ['CO'], ['NO','NO2','NO3I','NO3J'], ['SO2','SO4I','SO4J'], \
       ['NH3'], ['CO_A'], ['ORGI_A','ORGJ_A'], ['PM_10'], ['PM25I','PM25J'], \
       ['ECI','ECJ'], ['ORGI','ORGJ'], ['BRC1I'], ['BRC1J'], ['BRC2I'], \
       ['BRC2J'], ['BRC3I'], ['BRC3J'], ['DBRCI','DBRCJ']]
mix = [['PAR'], ['IOLE','PAR'], ['ETH'], ['ETOH'], ['ETHA'], \
       ['OLE','PAR'], ['CO'], ['FORM'], ['ALD2','ALDX'], ['CO'], \
       ['MEOH'], ['TOL'], ['XYL'], ['ISOP'], ['TERP'], \
       ['CO'], ['CO'], ['CO'], ['NOx'], ['SO2'], \
       ['NH3'], ['CO'], ['OC'], ['PM10'], ['PM25'], \
       ['BC'], ['OC'], ['BRC1I'], ['BRC1J'], ['BRC2I'], \
       ['BRC2J'], ['BRC3I'], ['BRC3J'], ['DBRC']]

coef_wrf = [[1],[1],[1],[1],[1],\
            [1],[1],[1],[1],[1],\
            [1],[1],[1],[1],[1],\
            [1],[1],[1],[0.9,0.1,0,0],[1,0,0],\
            [1],[1],[0.2,0.8],[1],[0.2,0.8],\
            [0.2,0.8],[0.2,0.8],[1],[1],[1],\
            [1],[1],[1],[0.2,0.8]]
coef_mix = [[5],[1,2],[1],[1],[1],\
            [1,1],[0.0112],[1],[1,1],[0.018],\
            [1],[1],[1],[1],[1],\
            [0.00587],[0.00024],[1],[1],[1],\
            [1],[1],[1],[1],[1],\
            [1],[1],[1],[1],[1],\
            [1],[1],[1],[1]]

# 1.define dimesions
Time = df_new.createDimension('Time', None)
emissions_zdim = df_new.createDimension('emissions_zdim', size=nlev)
west_east = df_new.createDimension('west_east', size=nx)
south_north = df_new.createDimension('south_north', size=ny)
DateStrLen = df_new.createDimension('DateStrLen', size=19)

# 2.define variables
my_dict = {}
# for pollutants
for item in cbmz_wrf: 
    my_dict['E_'+item] = df_new.createVariable('E_' + item, 'f4', dimensions=('Time', 'emissions_zdim', 'south_north', 'west_east'))
# for Times
Times = df_new.createVariable('Times', '|S1', dimensions=('Time','DateStrLen'))
#lon = df_new.createVariable('west_east', 'f4',dimensions=('west_east'))
#lat = df_new.createVariable('south_north','f4',dimensions=('south_north'))

# 3.Add data to pollutants and do the layers distribution
# distribute the emission from Power sector the first six levels according to the following coefficients.
powlist = [0.1,0.1,0.3,0.2,0.2,0.1,0,0]

for i in range(len(mix)):
    # define a null array to store the sum of the data in mix[i], for example: ALD2+ALDX
    copy_array = np.zeros((nt, nlev, ny, nx))
    for j in range(len(mix[i])):
        # get the pollutant name
        var = mix[i][j]
        # get the filename for the pollutant
        print(var+'.nc')
        filename = flist[[idx for idx, s in enumerate(flist) if var+'.nc' == s.split('/')[-1]][0]]
        cid = Dataset(filename)
        # get the variables in the pollutant file
        mix_Var = list(cid.variables)
        mix_Var.remove('Times')
        mix_Var.sort()
        
        # Power sector is the forth variable.
        copy_array[:,0,:,:] += coef_mix[i][j]*(cid[mix_Var[0]][:,0,:,:] + cid[mix_Var[3]][:,0,:,:]\
                                                         * powlist[0] + cid[mix_Var[1]][:,0,:,:] + cid[mix_Var[2]][:,0,:,:] + cid[mix_Var[4]][:,0,:,:] + cid[mix_Var[5]][:,0,:,:] + cid[mix_Var[6]][:,0,:,:]) 
        # distribute the emission from Power sector to different levels
        for k in range(1,len(powlist)):
            copy_array[:,k,:,:] += coef_mix[i][j]*cid[mix_Var[3]][:,0,:,:] * powlist[k]
    
    # allocate the data from mix[i] to different variables in wrf[i]
    for p in range(len(wrf[i])):    
        my_dict['E_'+ wrf[i][p]][:,:,:,:] = coef_wrf[i][p]*copy_array

# get the variables in cbmz_wrf which has no match species in the MIX inventory
diff = set(cbmz_wrf)-set([item for sublist in wrf for item in sublist])
#print('Species without species match in cbmz_wrf:',diff)
for item in diff:
    my_dict['E_' + item][:,:,:,:] = 0.0

# change the negative values in PM25 to 0
#my_dict['E_PM_10']=my_dict['E_PM_10'][:,:,:,:]-my_dict['E_PM25I'][:,:,:,:]-my_dict['E_PM25J'][:,:,:,:]
#my_dict['E_PM_10'][:,:,:,:][my_dict['E_PM_10'][:,:,:,:]<0] = 0
#my_dict['E_PM25I'][:,:,:,:][my_dict['E_PM25I'][:,:,:,:] < 0] = 0
#my_dict['E_PM25J'][:,:,:,:][my_dict['E_PM25J'][:,:,:,:] < 0] = 0    

# 4.Add data to Times variable
time_str =['2017_01_01_00:00:00','2017_01_01_01:00:00',
           '2017_01_01_02:00:00','2017_01_01_03:00:00',
           '2017_01_01_04:00:00','2017_01_01_05:00:00',
           '2017_01_01_06:00:00','2017_01_01_07:00:00',
           '2017_01_01_08:00:00','2017_01_01_09:00:00',
           '2017_01_01_10:00:00','2017_01_01_11:00:00',
           '2017_01_01_12:00:00','2017_01_01_13:00:00',
           '2017_01_01_14:00:00','2017_01_01_15:00:00',
           '2017_01_01_16:00:00','2017_01_01_17:00:00',
           '2017_01_01_18:00:00','2017_01_01_19:00:00',
           '2017_01_01_20:00:00','2017_01_01_21:00:00',
           '2017_01_01_22:00:00','2017_01_01_23:00:00']
time_str = [s.replace('2017', str(year)) for s in time_str]
for i in range(len(time_str[t-12:t])):
    Times[i] = np.array([char for char in time_str[t-12:t][i]])
#print('Times variables:',Times[:])
#lon[:] = np.arange(-nx*res*1000/2+res*1000/2, nx*res*1000/2+res*1000/2, res*1000)
#lat[:] = np.arange(-ny*res*1000/2+res*1000/2, ny*res*1000/2+res*1000/2, res*1000)
    
# 5.Change units and do the time distribution 
Gas = ['CO','NO','NO2','SO2','NH3','CO_A']
VOC = ['BIGALK','BIGENE','C2H4','C2H5OH','C2H6','C3H6','C3H8','CH2O','CH3CHO','CH3COCH3','CH3OH','TOLUENE','XYLENE','ISOP','C10H16','C2H2','MVK']
Aerosol = ['NO3I','NO3J','SO4I','SO4J','ORGI_A','ORGJ_A','PM_10','PM25I','PM25J','ECI','ECJ','ORGI','ORGJ','BRC1I','BRC1J','BRC2I','BRC2J','BRC3I','BRC3J','DBRCI','DBRCJ']

# days in 12 months
monday = [31,28,31,30,31,30,31,31,30,31,30,31]
# molecular weights for Gas variables: CO SO2 NH3 NO NO2
MW = [28, 30, 46, 64, 17, 28]
# Time distribution coefficient for 24 hours from 00 to 23 in UTC
coef_time = [0.0467,0.05366,0.05506,0.05876,0.0602,0.05936,0.05768,0.056,0.05496,0.05852,0.05446,0.04702,0.04532,0.03946,0.02922,0.0265,0.02284,0.01752,0.0174,0.01722,0.01772,0.02266,0.03504,0.04672]

for item in Gas:
    original = copy.deepcopy(my_dict['E_'+item][month-1,:,:,:])
    for i in range(t-12,t):
        # Mg/month(mix) to mol/km2/hr (wrfchem)
        my_dict['E_'+item][i-t+12,:,:,:] =  coef_time[i] * original * 10**6/ res**2 / (monday[month-1]) / MW[Gas.index(item)]
for item in VOC:
    original = copy.deepcopy(my_dict['E_'+item][month-1,:,:,:])
    for i in range(t-12,t):
        # Mmol/month(mix) to mol/km2/hr (wrfchem)
        my_dict['E_'+item][i-t+12,:,:,:] =  coef_time[i] * original * 10**6/ res**2 / (monday[month-1])

for item in Aerosol:
    original = copy.deepcopy(my_dict['E_'+item][month-1,:,:,:])
    for i in range(t-12,t):
        # Mg/month(mix) to ug∙𝒎−𝟐∙s−𝟏 (wrfchem)
        my_dict['E_'+item][i-t+12,:,:,:] =  coef_time[i] * original * 10**6/ res**2 / (monday[month-1]*3600)
    
# 5.Add attributes 
for item in Gas: 
    my_dict['E_'+item].units = 'mol km^-2 hr^-1'
    my_dict['E_'+item].coordinates = 'XLONG XLAT'
    my_dict['E_'+item].description = 'EMISSIONS'
    my_dict['E_'+item].MemoryOrder = 'XYZ'
    my_dict['E_'+item].FieldType = 104
    my_dict['E_'+item].stagger= ''
for item in VOC:
    my_dict['E_'+item].units = 'mol km^-2 hr^-1'
    my_dict['E_'+item].coordinates = 'XLONG XLAT'
    my_dict['E_'+item].description = 'EMISSIONS'
    my_dict['E_'+item].MemoryOrder = 'XYZ'
    my_dict['E_'+item].FieldType = 104
    my_dict['E_'+item].stagger= ''                                                    
for item in Aerosol:
    my_dict['E_'+item].units = 'ug/m3 m/s'
    my_dict['E_'+item].coordinates = 'XLONG XLAT'
    my_dict['E_'+item].description = 'EMISSIONS'
    my_dict['E_'+item].MemoryOrder = 'XYZ'
    my_dict['E_'+item].FieldType = 104
    my_dict['E_'+item].stagger= ''
df_new.close()
