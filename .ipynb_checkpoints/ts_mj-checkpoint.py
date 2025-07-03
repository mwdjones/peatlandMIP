#!/usr/bin/env python
# coding: utf-8
import subprocess
import argparse
import numpy as np 
import netCDF4 as netcdf4 
import matplotlib 
import matplotlib.pyplot as plt 
import pandas as pd

parser = argparse.ArgumentParser(description='')
parser.add_argument("-s","--sim", help="simulation", type=int, nargs='?', default=1)
args = parser.parse_args()

# Function definitions

def ilevel(numlevel):
    return np.arange(numlevel,dtype=float)/float(numlevel-1)
    
def get_history_filenames(idir,htag='h0',yrs=None):
    ttag = ''
    if yrs != None:
        if yrs[0] == yrs[1]:
            ttag = str(yrs[0])
        else:
            ttag = ','.join(str(y) for y in range(yrs[0],yrs[1]+1))
            ttag = '{'+ttag+'}'
            
    command = 'ls '+idir+'*clm2.'+htag+'.'+ttag+'*.nc'
    
    return subprocess.Popen(command,stdout=subprocess.PIPE,shell='True').communicate()[0].strip().decode().split('\n')
    
def adjust_axes_for_horizontal_colorbar(fig,cb_hgt_in,scale=1.0):
    ''' adjust result of plt.subplots to accomodate a 
    horizontal colorbar at the bottom of the figure'''
    cb_hgt = scale * cb_hgt_in
    ax = fig.get_axes()
    ymin=[]
    ymax=[]
    for a in ax:
        x = a.get_position() # (Bbox defined as [[x0,y0],[x1,y1]])
        ymin = np.append(ymin,np.min([x.y0,x.y1]))
        ymax = np.append(ymax,np.max([x.y0,x.y1]))
        #print ('coords: ', x.x0,x.x1,x.y0,x.y1)
        #print ('width, height: ', x.x1-x.x0,x.y1-x.y0)
    ymin = np.min(ymin)
    ymax = np.max(ymax)
 
    yspan = ymax - ymin ; new_yspan = yspan - cb_hgt
    scale_factor = new_yspan / yspan
    for a in ax:
        x = a.get_position()
        x.y0 = scale_factor*(x.y0 - ymin) + (ymin + cb_hgt)
        x.y1 = scale_factor*(x.y1 - ymin) + (ymin + cb_hgt)
        a.set_position(x)

def adjust_axes_for_column_size(col_ax,col_scalar):
    xpad = []
    for n in range(len(col_ax)):
        a = col_ax[n]
        x = a.get_position()
        if n == 0:
            xs = x.x0
        if n == len(col_ax)-1:
            xe = x.x1
        if n > 0 and n < len(col_ax):
            xpad.append(x.x0 - xr)
        xr = x.x1

    xwidth = (xe-xs)-np.sum(xpad)
    for n in range(len(col_ax)):
        a = col_ax[n]
        x = a.get_position()
        x.x0 = xs + xwidth*np.sum(column_scalar[0:n])+np.sum(xpad[0:n])
        x.x1 = x.x0 + xwidth*column_scalar[n]
        a.set_position(x)
        
def create_merged_colormap(cname1,cname2,full=True,half=False):
    from matplotlib.colors import ListedColormap
    N = int(256)
    if half:
        full = False
    if full:
        # this merges two full colormaps
        bottom = plt.colormaps.get_cmap(cname1, N//2)        
        top = plt.colormaps.get_cmap(cname2, N//2)
        newcolors = np.vstack((bottom(np.linspace(0, 1, N//2)),
            top(np.linspace(0, 1, N//2))))
    if half:
        # this merges two half-colormaps
        bottom = plt.get_cmap(name=cname1, lut=N)
        top = plt.get_cmap(name=cname2, lut=N)
        newcolors = np.vstack((bottom(np.linspace(0, 0.5, int(N//2))),
            top(np.linspace(0.5, 1.0, int(N//2)))))

    newcmap = ListedColormap(newcolors)#, name='SoilMoisture')
    return newcmap

def calc_water_table_top_down(sm_sat,zsoi,zisoi,sat_lev = 0.9):
    km,tm = sm_sat.shape
    zwt = np.zeros((tm))
    for t in range(tm):
        kzwt = km-1
        for k in range(km):
            if sm_sat[k,t] >= sat_lev:
                kzwt = k
                break
            
        # interpolate between k and k+1
        k1 = int(kzwt)
        k2 = int(np.max([kzwt-1,0]))
        s1 = sm_sat[k1,t]
        s2 = sm_sat[k2,t]
        
        if s2 >= sat_lev: 
            zwt[t] = zisoi[k1]
        else:
            m = (zsoi[k2]-zsoi[k1])/(s2-s1)
            b = zsoi[k2]-m*s2
            zwt[t] = m*sat_lev+b

    return zwt

def calc_isotherm_top_down(tsoi,zsoi,zisoi,iso = 0.):
    km,tm = tsoi.shape
    isotherm = np.zeros((tm))
    for t in range(tm):
        kzwt = 0
        for k in range(km-1):
            if tsoi[k,t] <= iso and tsoi[k+1,t] > iso:
                kzwt = k
                break
            
        # interpolate between k-1 and k+1
        k1 = int(np.min([km-1,kzwt+1]))
        k2 = int(np.max([kzwt-1,0]))
        s1 = tsoi[k1,t]
        s2 = tsoi[k2,t]
        
        if s2 >= iso: 
            isotherm[t] = zisoi[k2]
        else:
            m = (zsoi[k2]-zsoi[k1])/(s2-s1)
            b = zsoi[k2]-m*s2
            isotherm[t] = m*iso+b
    return isotherm

# Select which figures to output

plot_dir = './'
SavePlotsPng = False
SavePlotsPng = True

plot_suffix = 'png'

file_prefix = 'ts'
figcnt = 0

if SavePlotsPng:
    pngfiles = []

# Read in CTSM history filenames

#@
base_dir = '/glade/derecho/scratch/marielj/archive/'
inum = args.sim
if inum == 1:
    idir = f'{base_dir}hillslope-lce-calib-medlynslope/lnd/hist/'

sitename = 'MEF'

print(idir.split('/')[-1])
    
yr1, yr2 = 2004, 2008
# datm data begins in 2011, not 2004
yr1, yr2 = 2011, 2017
    
h0files = get_history_filenames(idir,htag='h0',yrs=[yr1,yr2])
h1files = get_history_filenames(idir,htag='h1',yrs=[yr1,yr2])
h2files = get_history_filenames(idir,htag='h2',yrs=[yr1,yr2])
h3files = get_history_filenames(idir,htag='h3',yrs=[yr1,yr2])
zfile = get_history_filenames(idir,htag='h0')[0]

if len(h1files[0]) < 1:
    print('no files found in '+idir)
    exit()

# Select hillslope and column

f1 = netcdf4.Dataset(h3files[0], 'r')
nmaxcols = len(f1.dimensions['max_columns_hillslope']) 
nhillslope = len(f1.dimensions['nhillslope']) 
f1.close()

aspect_names = ['North','East','South','West']
naspect = len(aspect_names)

ncol_per_hillslope = int(nmaxcols//nhillslope)

# Read in variables

hm = len(h3files)
for t in range(hm):
    ifile = h3files[t]
    #print(ifile)
    f1 = netcdf4.Dataset(ifile, 'r')
    if t == 0:
        tm0     = len(f1.dimensions['time'])
        nlevsoi = len(f1.dimensions['levsoi'])
        lon = np.asarray(f1.variables['lon'][:,]) 
        lat = np.asarray(f1.variables['lat'][:,])
        pt = [lon[0],lat[0]]
        lonstr = f'{lon[0]:.2f}'
        latstr = f'{lat[0]:.2f}'
        coordlabel = lonstr+'E_'+latstr+'N'
        print(coordlabel)
        h_index= np.asarray(f1.variables['hillslope_index'][:,]) 
        nactualcols = int(np.sum(np.where(h_index > 0,1,0)))
        
        fv = f1.variables['landfrac'].getncattr('_FillValue')
        ref_year = int(f1.variables['time'].units.split()[2].split('-')[0])

        # read in soil layer structure from first h0 file
        z1 = netcdf4.Dataset(zfile, 'r')
        zsoi   = np.squeeze(z1.variables['ZSOI'][:,])[0:nlevsoi,]
        dzsoi  = np.squeeze(z1.variables['DZSOI'][:,])[0:nlevsoi,]
        watsat = np.squeeze(z1.variables['WATSAT'][:,])[0:nlevsoi,]
        z1.close()

        zisoi = zsoi + 0.5*dzsoi
        
        hslp_index = f1.variables['hillslope_index'][:,]
        hslp_elev = f1.variables['hillslope_elev'][:,]
        hslp_dist = f1.variables['hillslope_distance'][:,]
        hslp_area = f1.variables['hillslope_area'][:,]
        hslp_width = f1.variables['hillslope_width'][:,]
        hslp_aspect = f1.variables['hillslope_aspect'][:,]
        hslp_pftndx = f1.variables['pfts1d_itype_veg'][:,]
        hslp_nbedrock = f1.variables['cols1d_nbedrock'][:,]
        hslp_nbedrock -= 1 # convert to 0- index

        pfts1d_wtcol = f1.variables['pfts1d_wtcol'][:,]
        pfts1d_ci = f1.variables['pfts1d_ci'][:,]
        cids = np.unique(pfts1d_ci)
        npfts = int(pfts1d_ci.size/cids.size)
        
        tm = hm * tm0
        time  = np.zeros((tm))
        snow_depth = np.zeros((nactualcols,tm))
        swe = np.zeros((nactualcols,tm))
        qlatflowout = np.zeros((nactualcols,tm))
        qover = np.zeros((nactualcols,tm))
        qdrai = np.zeros((nactualcols,tm))
        h2osfc = np.zeros((nactualcols,tm))
        h2osoi = np.zeros((nactualcols,nlevsoi,tm))
        qvegt = np.zeros((nactualcols,tm))
        qvege = np.zeros((nactualcols,tm))
        qsoil = np.zeros((nactualcols,tm))
        rain = np.zeros((nactualcols,tm))
        snow = np.zeros((nactualcols,tm))
            
    t1, t2 = t*tm0, (t+1)*tm0
    #print(t,t1,t2)

    swe[:,t1:t2]  = f1.variables['H2OSNO'][:,].T
    snow_depth[:,t1:t2]  = f1.variables['SNOW_DEPTH'][:,].T
    qlatflowout[:,t1:t2]  = f1.variables['QLATFLOWOUT'][:,].T
    qover[:,t1:t2]  = f1.variables['QOVER'][:,].T
    qdrai[:,t1:t2]  = f1.variables['QDRAI'][:,].T
    h2osfc[:,t1:t2]  = f1.variables['H2OSFC'][:,].T
    h2osoi[:,:,t1:t2]  = f1.variables['H2OSOI'][:,].T

    for nc in range(nactualcols):

        x = np.sum((f1.variables['QVEGT'][:,].T)[pfts1d_ci==cids[nc],] * np.outer(pfts1d_wtcol[pfts1d_ci==cids[nc],],np.ones(tm0)),axis=0)
        qvegt[nc,t1:t2]  = x
        x = np.sum((f1.variables['QVEGE'][:,].T)[pfts1d_ci==cids[nc],] * np.outer(pfts1d_wtcol[pfts1d_ci==cids[nc],],np.ones(tm0)),axis=0)
        qvege[nc,t1:t2]  = x
        x = np.sum((f1.variables['QSOIL'][:,].T) [pfts1d_ci==cids[nc],]* np.outer(pfts1d_wtcol[pfts1d_ci==cids[nc],],np.ones(tm0)),axis=0)
        qsoil[nc,t1:t2]  = x
        
    rain[:,t1:t2]  = f1.variables['RAIN'][:,].T
    snow[:,t1:t2]  = f1.variables['SNOW'][:,].T

    time[t1:t2] = np.asarray(f1.variables['time'][:,])

    f1.close()
    
snow_depth[snow_depth >= fv] = 0
swe[swe >= fv] = 0
qlatflowout[qlatflowout >= fv] = 0
qvegt[qvegt >= fv] = 0
qvege[qvege >= fv] = 0
qsoil[qsoil >= fv] = 0
h2osfc[h2osfc >= fv] = 0
h2osoi[h2osoi >= fv] = 0
qover[qover >= fv] = 0
qdrai[qdrai >= fv] = 0
rain[rain >= fv] = 0
print('snow depth min/max ',np.min(snow_depth),np.max(snow_depth))
#print('fsds min/max ',np.min(fsds_col),np.max(fsds_col))
        
# convert to mm/day
flux_units = 'mm/day'
sf = 24*3600
qvegt = qvegt * sf
qvege = qvege * sf
qsoil = qsoil * sf
qover = qover * sf
qdrai = qdrai * sf
qlatflowout = qlatflowout * sf
rain = rain * sf
snow = snow * sf


#
#time -= 0.5/365
ytime = time/365 + ref_year
y1 = np.floor(np.min(ytime)).astype(int)
y2 = np.floor(np.max(ytime)).astype(int)
y1 = np.max([y1,yr1])
y2 = np.min([y2,yr2+1])
ytime_ticks = [y for y in range(y1,y2+1)]
ytime_ticklabels = [str(y) for y in range(y1,y2+1)]

for t in range(hm):
    ifile = h2files[t]
    #print(ifile)
    f1 = netcdf4.Dataset(ifile, 'r')
    if t == 0:
        tm0     = len(f1.dimensions['time'])
        tm = hm * tm0
        qdischarge = np.zeros((tm))
        
    t1, t2 = t*tm0, (t+1)*tm0
    #print(t,t1,t2)

    qdischarge[t1:t2]  = sf*np.squeeze(f1.variables['VOLUMETRIC_STREAMFLOW'][:,])

    f1.close()


showObsZwt = True
# read in observed water table
if showObsZwt:
    dpm =[31,28,31,30,31,30,31,31,30,31,30,31]
    obsfile = 'MEF_daily_peatland_water_table.csv'
    sdata = []
    with open(obsfile,'r') as f1:
        for tmp in f1:
            sdata.append(tmp.strip())

    stime, szwt = [],[]
    oref_year = int(sdata[1].split(',')[1].split('-')[0])
    for tmp in sdata[1:]:
        x = tmp.split(',')
        if x[0].strip('"') == 'S2':
            y,m,d = [int(i) for i in x[1].split('-')]
            t = 365*(y-oref_year)+d
            if m > 1:
                t += np.sum(dpm[:m-1])
            stime.append(t)
            szwt.append(float(x[2]))

    otime = np.asarray(stime)/365 + oref_year
    # estimate of surface elevation
    zsfc = 422.27
    ozwt = (zsfc - np.asarray(szwt))

# read in observed streamflow
if 1==1:
    dpm =[31,28,31,30,31,30,31,31,30,31,30,31]
    obsfile = 'Streamflow_daily.csv'
    sdata = []
    with open(obsfile,'r') as f1:
        for tmp in f1:
            sdata.append(tmp.strip())

    stime, sflow = [],[]
    oref_year = int(sdata[1].split(',')[0].split('-')[0])
    for tmp in sdata[1:]:
        x = tmp.split(',')
        if x[1] == 'S2':
            y,m,d = [int(i) for i in x[0].split('-')]
            t = 365*(y-oref_year)+d
            if m > 1:
                t += np.sum(dpm[:m-1])
            stime.append(t)
            sflow.append(float(x[2]))

    oftime = np.asarray(stime)/365 + oref_year
    sf = 1e-3 # convert from L/s to m3/s
    oflow = np.asarray(sflow)*sf

        
##0

##1
# Read in all columns - 2D variable

# Specify variable to read  
varname  = 'SOILLIQ'    ; sm_units = 'kg/m2'
varname2 = 'SOILICE'
varname3 = 'ZWT'


resetMap=True

# get watsat_col from monthly file
f1 = netcdf4.Dataset(h1files[0], 'r')
if 'watsat' in f1.variables:
    watsat_col = f1.variables['watsat'][0,:nlevsoi,].T
f1.close()    

hm = len(h3files)
for t in range(hm):
    ifile = h3files[t]
    #print(ifile)
    f1 = netcdf4.Dataset(ifile, 'r')
    if t == 0:
        tm0       = len(f1.dimensions['time'])
        nlevsoi = len(f1.dimensions['levsoi'])

        # read in soil layer structure from first h0 file
        z1 = netcdf4.Dataset(zfile, 'r')
        zsoi   = np.squeeze(z1.variables['ZSOI'][:,])[0:nlevsoi,]
        dzsoi  = np.squeeze(z1.variables['DZSOI'][:,])[0:nlevsoi,]
        watsat = np.squeeze(z1.variables['WATSAT'][:,])[0:nlevsoi,]
        z1.close()

        zisoi = zsoi + 0.5*dzsoi
        
        hslp_elev = f1.variables['hillslope_elev'][:,]
        hslp_dist = f1.variables['hillslope_distance'][:,]
        hslp_area = f1.variables['hillslope_area'][:,]
        hslp_width = f1.variables['hillslope_width'][:,]
        hslp_aspect = f1.variables['hillslope_aspect'][:,]
        hslp_pftndx = f1.variables['pfts1d_itype_veg'][:,]
        hslp_nbedrock = f1.variables['cols1d_nbedrock'][:,]
        hslp_nbedrock -= 1 # convert to 0- index
        
        if 'watsat' in f1.variables:
            watsat_col = f1.variables['watsat'][0,:nlevsoi,].T
        
        tm = hm * tm0
        time = np.zeros((tm))
        sm = np.zeros((nactualcols,nlevsoi,tm))
        liq_frac = np.zeros((nactualcols,nlevsoi,tm))
        zwt  = np.zeros((nactualcols,tm))
        zwt_perch  = np.zeros((nactualcols,tm))
        frost_table  = np.zeros((nactualcols,tm))
        tsoi = np.zeros((nactualcols,nlevsoi,tm))

    t1, t2 = t*tm0, (t+1)*tm0
    #print(t,t1,t2)

    sm[:,:,t1:t2]   =  f1.variables[varname][:,].T
    liq_frac[:,:,t1:t2]   =  sm[:,:,t1:t2]
    sm[:,:,t1:t2]  += f1.variables[varname2][:,].T
    liq_frac[:,:,t1:t2]   =  liq_frac[:,:,t1:t2]/sm[:,:,t1:t2]
    zwt[:,t1:t2]    = f1.variables[varname3][:,].T
    zwt_perch[:,t1:t2]    = f1.variables['ZWT_PERCH'][:,].T
    frost_table[:,t1:t2]  = f1.variables['FROST_TABLE'][:,].T
    tsoi[:,:,t1:t2] =  f1.variables['TSOI'][:,:nlevsoi,].T

    time[t1:t2] = np.asarray(f1.variables['time'][:,])
    f1.close()

tfrz = 273.15
tsoi -= tfrz

# total soil moisture
smtot = np.sum(sm,axis=1)

# Convert soil moisture to saturation
if varname == 'SOILLIQ':
    calcSaturation = False
    calcSaturation = True 

    print(sm.shape)
    print(watsat_col.shape)
    
    #sm = np.zeros((ncol_per_hillslope,nlevsoi,tm))
    if calcSaturation:
        pvarname = 'Saturation'
        sm_units = ''
        smsattot = np.copy(smtot)
        for k in range(nlevsoi):      
            for n in range(watsat_col.shape[0]):      
                sm[n,k,:] = sm[n,k,:]/(1.e3 * dzsoi[k] * watsat_col[n,k])
                if n==-2:
                    print('watsat col ',k,watsat_col[n,k])
        # total column saturation
        smsattot = np.copy(smtot)
        for n in range(watsat_col.shape[0]):      
            smsattot[n,:] = smtot[n,:]/(1.e3 * np.sum(dzsoi * watsat_col[n,]))

    else:
        pvarname = varname
        sm_units = 'm3/m3' 
        for k in range(nlevsoi):      
            sm[:,k,:] = sm[:,k,:]/(1.e3 * dzsoi[k])

plotH2osfcSM,plotSnowSM = False,False

plotH2osfcSM = True
#plotSnowSM = True

h2osfc_ymax = -0.4

##2
# Soil moisture profiles
nrow, ncol = 1, 3
wx, wy = 10, 6

fig1, ax1 = plt.subplots(nrow, ncol)#, sharex='col', sharey='row')
fax = fig1.get_axes()

fig1.set_size_inches([wx,wy])
fig1.subplots_adjust(hspace=0.5)

# modify bounding box of axes to accomodate colorbar at bottom of plot 
cb_hgt = 0.02
adjust_axes_for_horizontal_colorbar(fig1,cb_hgt,scale=1.)

# add the colorbar axes
cb_w=0.6
cbaxes  = fig1.add_axes([0.2,0.06,cb_w,cb_hgt])

#calcSaturation = False
if calcSaturation:
    plot_title = 'Soil moisture [Saturation]  Lon: '+lonstr+' \ Lat: '+latstr
else:
    plot_title = 'Soil moisture [vwc]  Lon: '+lonstr+' \ Lat: '+latstr

xmin = np.min(ytime)
xmax = np.max(ytime)
ymin = np.max(zsoi[hslp_nbedrock])
kz = 9
ymin = np.max(zsoi[kz-1])

ymin_arr = [1,2,4]
ymax = 0
ymax = -0.4
fax[-1].set_xlabel('time')
for n in range(len(fax)):
    fax[n].set_ylabel('depth')

    fax[n].set_xlim(xmin,xmax)
    #fax[n].set_ylim(ymin,ymax)
    fax[n].set_ylim(ymin_arr[n],ymax)
    fax[n].set_aspect('auto')
    fax[n].set_xticks(ytime_ticks)

w=fig1.suptitle(plot_title,fontsize=16)

lnum=25
# set contour levels
if calcSaturation:
    minfld = 0.2
    maxfld = 1
else:
    minfld = 0.
    maxfld = 0.45
datarange = maxfld - minfld
level=datarange*ilevel(lnum) + minfld
cblevel=datarange*ilevel(7) + minfld
cnorm = matplotlib.colors.BoundaryNorm(level, 256)

#cmap = plt.cm.BrBG
cmap = create_merged_colormap('BrBG','RdBu',half=True)
cmap.set_under('white')

img = []

for n in range(ncol_per_hillslope):
    #fax[n].set_title(col_names[n],fontsize=16)
    #kz = hslp_nbedrock[n]+1
    kz = np.argmin(np.abs(zsoi - ymin_arr[n]))+1
    img.append(ax1[n].contourf(ytime,zsoi[:kz],sm[n,:kz,:],cmap=cmap,
                               levels=level,norm=cnorm,extend='both'))

    ax1[n].plot(ytime,zwt[n,:],linewidth=1,c='y')
    
    if plotH2osfcSM:
        ax1[n].plot(ytime,-1e-3*h2osfc[n,:],linewidth=1,c='g')
    if plotSnowSM:
        ax1[n].plot(ytime,-1e-3*snow_depth[n,:],linewidth=1,c='b')

# observed water table
if showObsZwt:
    lw = 1
    fax[0].plot(otime,ozwt,linewidth=lw,c='w',alpha=1)
    fax[1].plot(otime,ozwt+hslp_elev[1],linewidth=lw,c='w',alpha=1)
    fax[2].plot(otime+0.2,2*ozwt+1,linewidth=lw,c='w',alpha=1)

      
cbar = fig1.colorbar(img[0],ax=fax[0],ticks=cblevel,label='',
                     orientation="horizontal",cax = cbaxes)

if 1==1:
    ff = 5
    fax[0].plot(oftime,ff*ymax*oflow/np.max(oflow),linewidth=1,c='k',alpha=1)
    fax[1].plot(oftime,ff*ymax*oflow/np.max(oflow),linewidth=1,c='k',alpha=1)


if SavePlotsPng:
    figcnt += 1
    fname = plot_dir + file_prefix+'.fig.'+str(int(figcnt))+'.'+plot_suffix 
    pngfiles.append(fname)
    plt.savefig(fname,format='png')

##3
# Soil moisture profiles

fig1, ax1 = plt.subplots(1,1)
fax = fig1.get_axes()
ax1 = [ax1]
fig1.set_size_inches([wx,wy])
fig1.subplots_adjust(hspace=0.5)

# modify bounding box of axes to accomodate colorbar at bottom of plot 
cb_hgt = 0.02
adjust_axes_for_horizontal_colorbar(fig1,cb_hgt,scale=1.)

# add the colorbar axes
cb_w=0.6
cbaxes  = fig1.add_axes([0.2,0.06,cb_w,cb_hgt])

#calcSaturation = False
if calcSaturation:
    plot_title = 'Soil moisture [Saturation]  Lon: '+lonstr+' \ Lat: '+latstr
else:
    plot_title = 'Soil moisture [vwc]  Lon: '+lonstr+' \ Lat: '+latstr

#xmin = np.min(ytime)
#xmax = np.max(ytime)

ymin = np.max(zsoi[hslp_nbedrock])
kz = 9
ymin = np.max(zsoi[kz-1])
#print(ymin)

# column to plot
nplot = 0 # lagg
nplot = 1 # bog
ymin_arr = [1,2,4]
ymax = 0
ymax = -0.4
fax[-1].set_xlabel('time')
for n in [0]:
    fax[n].set_ylabel('depth')

    fax[n].set_xlim(xmin,xmax)
    #fax[n].set_ylim(ymin,ymax)
    fax[n].set_ylim(ymin_arr[nplot],ymax)
    fax[n].set_aspect('auto')
    #fax[n].set_xticks(ytime_ticks)

w=fig1.suptitle(plot_title,fontsize=16)

lnum=25
# set contour levels
if calcSaturation:
    minfld = 0.2
    maxfld = 1
else:
    minfld = 0.
    maxfld = 0.45
    #maxfld = 0.8
datarange = maxfld - minfld
level=datarange*ilevel(lnum) + minfld
cblevel=datarange*ilevel(7) + minfld
cnorm = matplotlib.colors.BoundaryNorm(level, 256)

#cmap = plt.cm.BrBG
cmap = create_merged_colormap('BrBG','RdBu',half=True)
cmap.set_under('white')

img = []

#fax[0].set_title(col_names[nplot],fontsize=16)
#kz = hslp_nbedrock[nplot]+1
kz = np.argmin(np.abs(zsoi - ymin_arr[nplot]))+1
img.append(ax1[0].contourf(ytime,zsoi[:kz],sm[nplot,:kz,:],cmap=cmap,levels=level,norm=cnorm,extend='both'))

ax1[0].plot(ytime,zwt[nplot,:],linewidth=1,c='y') #MWJ: original water table elevation searching method

new_zwt = calc_water_table_top_down(sm[nplot,:kz,:],zsoi[:kz],zisoi[:kz],sat_lev = 0.9) #MWJ: adjusted water table elevation searching method
#ax1[0].plot(ytime,new_zwt,linewidth=1,c='purple')

# calculate composite water table from two methods (i.e. for ice-free / ice-present conditions)
if 1==1: #ice/liq frac based
    # re-use calc_isotherm function to calculate ice_frac contour
    icef = calc_isotherm_top_down(liq_frac[0,:kz,:],zsoi[:kz],zisoi[:kz],iso = 0.5)     
    l1 = np.logical_and(tsoi[0,0,:] < 0.01,icef > new_zwt)
    composite_zwt2 = np.where(l1,icef,new_zwt)
    #ax1[0].plot(ytime,icef,linewidth=1,c='w')#c=',greenyellow')
    ax1[0].plot(ytime,composite_zwt2,linewidth=1,c='greenyellow') #*********************** MWJ: new plotted water table elevation
    np.savetxt("amended-wte.csv", pd.DataFrame([ytime, composite_zwt2]).T, delimiter=",") #save

if plotH2osfcSM:
    ax1[0].plot(ytime,-1e-3*h2osfc[nplot,:],linewidth=1,c='g')
if plotSnowSM:
    ax1[0].plot(ytime,-1e-3*snow_depth[nplot,:],linewidth=1,c='b')

# plot observed water table
if showObsZwt:
    lw = 1
    #fax[0].plot(otime,ozwt,linewidth=lw,c='w',alpha=1)
    zoff = -0.3 # eyeballed...
    #zoff = -hslp_elev[nplot]

    fax[0].plot(otime,ozwt+zoff,linewidth=lw,c='orange',alpha=1)
        
cbar = fig1.colorbar(img[0],ax=fax[0],ticks=cblevel,label='',
                     orientation="horizontal",cax = cbaxes)

# observed streamflow
if 1==1:
    ff = 5
    fax[0].plot(oftime,ff*ymax*oflow/np.max(oflow),linewidth=1,c='k',alpha=1)


if SavePlotsPng:
    figcnt += 1
    fname = plot_dir + file_prefix+'.fig.'+str(int(figcnt))+'.'+plot_suffix 
    pngfiles.append(fname)
    plt.savefig(fname,format='png')

