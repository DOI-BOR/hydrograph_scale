# -*- coding: utf-8 -*-
"""
Hydrograph Scaling Script, Hydrograph Prep  (v1)
@author: USBR (tclarkin 2021)

This script cycles through a user specified NWIS gage, identifies a specified number of largest 1-day events, then
extracts 1-hour flows for use in hydrgrograph scaling.

"""
### Import Libraries ###
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime as dt
from functions import nwis_import,analyze_voldur,init_hydrograph_plot

### Begin User Input ###
# Set Working Directory
#os.chdir("C://Users//tclarkin//Documents//Resources//Complete_scripts//hydrograph_scale")
site_source = ["1954_edit.csv","1964_hydro.csv","1999_hydro.csv","2008_hydro.csv","2015_edit.csv"]

# Hydrograph Settings
duration = 15  # int, days
count = 5       # int, number of hydrographs to extract
sel_wy = None # list or None

### Begin Script ###
# Check for directory
hydrodir = "hydro_prep"
if not os.path.isdir(hydrodir):
    os.mkdir(hydrodir)

# Check if site_source is a list
if isinstance(site_source,list)==False:
    site_source = [site_source]
hydro_out = pd.DataFrame()
top_evs = pd.DataFrame()

# Loop through site sources
for site in site_source:
    # Remove files extension
    if "." in site:
        site_lab = site.split(".")[0]
    else:
        site_lab = site
    print(site_lab)

    # Check for directory
    sitedir = f"{hydrodir}/{site_lab}"
    if not os.path.isdir(sitedir):
        os.mkdir(sitedir)

    # Import daily data
    if len(site)==8 and "." not in site:
        usgs = True
        site_daily = nwis_import(site_source,"dv")

        # Remove data prior to 1989 (instantaneous data)
        site_daily = site_daily.loc[site_daily.index.year>=1989,:]
        findevents = True

        # ID annual maxima for duration
        site_evs = analyze_voldur(site_daily,duration)

        # Select top hydrographs corresponding to user specified count
        site_evs_sort = site_evs["avg_flow"].sort_values(ascending=False)
        top_evs = site_evs.loc[site_evs_sort[0:count].index]

        # If sel_wy specified, add to list
        if sel_wy!=None:
            # If sel_wy is longer than count, return message
            if len(sel_wy)<count:
                print(f'{len(sel_wy)} events specified, while only {count} events are expected')
            sel_evs = pd.DataFrame(columns=site_evs.columns)
            for wy in sel_wy:
                if wy in top_evs.index:
                    print(f"{wy} already in top events")
                    continue
                else:
                    print(f"{wy} not in top events; adding")
                    sel_evs.loc[wy,:] = site_evs.loc[wy,:]
            if len(sel_evs)>0:
                top_evs = top_evs.drop(index=top_evs.tail(len(sel_evs)).index)
                for wy in sel_evs.index:
                    top_evs.loc[wy,:] = sel_evs.loc[wy,:]

        # Retrieve instantaneous data
        hydro = pd.DataFrame(index=range(0,duration*24))
        init_hydrograph_plot(duration*24)
        for wy in top_evs.index:
            print(wy)
            start = top_evs.loc[wy,"start"].strftime("%Y-%m-%d")
            end = top_evs.loc[wy,'end'].strftime('%Y-%m-%d')
            site_inst = nwis_import(site_source,"iv",start=start,end=end)

            # Convert to hourly
            site_hour = site_inst["flow"].resample("H").mean()
            site_hour = site_hour.interpolate()
            site_hour = site_hour.dropna()
            site_hour = site_hour.reset_index(drop=True)

            if len(site_hour) < len(hydro):
                print("Adding timesteps...")
                n = 1
                while len(site_hour) < len(hydro):
                    site_hour = site_hour.append(site_hour.tail(1))
                    n += 1
                print(f"{n} timesteps added")
            if len(site_hour) > len(hydro):
                print(f"Clipping {len(site_hour)-len(hydro)} timesteps")
                site_hour = site_hour[0:duration*24]

            hydro[wy] = site_hour.values.astype("int")
            plt.plot(hydro[wy],label=wy)

    else:
        usgs = False
        # Import data
        site_data = pd.read_csv(site, index_col=0, parse_dates=True)
        # Convert to hourly
        site_hourly = site_data.resample("H").mean().interpolate()
        site_hourly_avg = site_hourly.rolling(duration*24).mean()
        top_evs.loc[site_lab,'start'] = site_hourly_avg.idxmax().item() - dt.timedelta(days=duration)
        top_evs.loc[site_lab,'end'] = site_hourly_avg.idxmax().item()

        # Identify maximum period for selected duration
        site_hourly.reset_index(drop=True,inplace=True)
        site_hourly_avg = site_hourly.rolling(duration*24).mean()
        site_duration_max = site_hourly_avg.idxmax().item()
        # Extract max period
        hydro = site_hourly.iloc[site_duration_max-(duration*24-1):site_duration_max+1]
        hydro.reset_index(drop=True,inplace=True)
        # Append to hydro_out
        hydro_out.loc[:,site_lab] = hydro.astype("int")

        # Plot
        plt.plot(hydro, label=site_lab)

    # Save figure and files
    plt.legend()
    plt.savefig(f"{sitedir}/{site_lab}_input_hydrographs.jpg", bbox_inches='tight', dpi=300)
    hydro.to_csv(f"{sitedir}/{site_lab}_input_hydrographs.csv",index=False)
    top_evs.to_csv(f"{sitedir}/{site_lab}_input_hydrographs_summary.csv")

# Save hydro out
hydro_out.to_csv(f"{hydrodir}/input_hydrographs.csv",index=False)
