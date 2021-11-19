# -*- coding: utf-8 -*-
"""
Created on Fri Dec 06, 2019
Hydrograph Scaling Script  (v4)
@author: tclarkin (USBR 2020)

This script takes a user supplied flood frequency analysis (FFA) curve and/or volume frequency analysis (VFA) curve and scales
user supplied hydrographs to specified return period values.

The Peaks method simply calculates the difference between the FFA and the hydrograph peak and multiplies all values by that ratio.
The Volumes method calculates the slope (m) of a line (intercept = 0) that is used to scale the input hydrograph to the VFA volume:
    The scaling is conducted in linear space, only.

The Peaks & Volumes method:
    First, scales the hydrograph to the Peak, as described above.
    Second, solves for the slope (m) of a line (q_h = q_max, when q_h = q_max) that is used to scale the hydrograph to the VFA volume, while preserving the peak.
    Both values are checked. If the peaks or volumes are not within 0.0001 (or 0.01%) of the FFA/VFA value, an error is returned in the command line

The Balanced method:
    ...
    
A note on notation:
    q is the individual flows, Qbar is the average flow, V is the total volume (sum(q) or Qbar*timesteps), P is peak
    
    
This script requires:
    ffile = flood frequency curve as .csv, with aep (column 1), return period (column 2), and flow (column 3), should include header.
    vfile = volume frequency curve as .csv, with aep (column 1), return period (column 2), and x-day average flow (column 3), should include header.
    hfile = hydrographs to be scaled as .csv, with timestep (hours; column 1), and flows (cfs; columns 2+), should include header.

The user must:
    - Specify the above listed files
    - Verify column numbers (0 indexed)
    - Verify critical duration
    - Verify timestep (note...the program may not run with timesteps other than 1 hour)
    - Select a scaling method ("peaks&volumes", "peaks", "volumes" or "balanced")
    - Select if scaling should be done in log space or not (True/False)
    - Specify return periods to be output
    - Specify which input hydrographs to scale
    - Create "Output" folder

"""
### Import Libraries ###
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from functions import hydroclip,importfreq
from functions import peak_scale,vol_scale,pav_scale
from functions import get_balanced_indices,balanced_scale
from functions import check_hydrographs,check_pmf,init_freq_plot,init_hydrograph_plot,get_linestyles

### Begin User Input ### 
# Set Working Directory
os.chdir("C://Users//tclarkin//Documents//Projects//Anderson_Ranch_Dam//hydrographs")
site = "ard"

# Load FFA Curve (should have same return periods as vfa)
ffile = 'ffa.csv'
ffa_col = "median"

# Load VFA Curve (should have same return periods as ffa)
vfile = 'vfa.csv'
vfa_col = "2880"     # specify the column(s), columns names should be durations in hours (1, 2, etc.), used for volume scaling

# Hydrograph duration
duration = 120   # days

# Load Hydrograph(s)
hfile = 'hydro_prep/ard_input_hydrographs.csv'
hydro_sel = ["2017","1997","1996","2011","1999","1989"]   # List of which hydrographs to scale (column name)

# Select Methods
scaleby = "volumes"           #"peaks&volumes", "peaks", "volumes" or "balanced"
rp_sel = [2,3.333,5,10,20,50,100,200,500,1000,2000,10000,20000,50000,100000]     # Which RPs to produce hydrographs for, must correspond to RPs in ffile vfile#
baseflow = False               # True to include transition from baseflow to scaled hydrograph
plotmin = 10
plotmax = 30000

# Additional settings for Balanced
plot_construct = True       # True to plot evolution of balanced hydrographs.
peak_type = "spiked"        # "peaked" or "spiked"
balanced_cols = ["720","2880"]  # columns to be used for balance hydrographs

# check against PMF (if None, no calculation will be completed)
pmf_peak = None               # int
pmf_volumes = None   # list

### Begin Script ###
## Setup ##
# Check for directory
if not os.path.isdir("output"):
    os.mkdir("output")

# Define duration:
duration = int(duration * 24)

# Import Data:
print("Importing hydrographs and frequency curves")
# Clip hydrographs to critical duration window
hydro_in,hydro_sel = hydroclip(hfile, hydro_sel, duration)

fig, ax = plt.subplots(figsize=(6.25, 4))
for h in range(0, len(hydro_sel)):
    print(h)
    plt.plot(range(0, duration), hydro_in[:, h] * (1 / np.max(hydro_in[:, h])), label=hydro_sel[h])
plt.ylabel("Index Flow")
plt.xlabel("Time (hours)")
plt.legend()
plt.savefig(f'{hfile.split(".")[0]}.jpg', bbox_inches='tight', dpi=300)

# Extract and check frequency curves
ffa, vfa, rps_in, aep = importfreq(ffile, vfile, rp_sel)
vfa_durs = [int(i) for i in vfa.columns[2:]]

# Check that user has input a valid method
if (scaleby != "peaks" and scaleby != "volumes" and scaleby != "peaks&volumes" and scaleby != "balanced") == True:
    sys.exit("Please select scaling method: \"peaks\", \"volumes\", \"peaks&volumes\" or \"balanced\"")
else:
    print("Scaling hydrographs by " + scaleby)

# Create output matrix
output = np.zeros((duration, len(rps_in), hydro_in.shape[1]))

## Analysis ##
# Begin looping through selected hydrographs
for h in range(0, hydro_in.shape[1]):
    print("Beginning hydrograph " + hydro_sel[h])
    hydro = hydro_in[:, h]

    if scaleby == "balanced":
        vfa_durs = [int(i) for i in balanced_cols]  # get duratinos from selected column names
        balanced_idx = get_balanced_indices(vfa_durs, hydro)

    for r in range(0, len(rps_in)):
        print("Beginning return period " + str(rps_in[r]))

        if scaleby == "peaks":
            peak_rp = ffa.loc[r, ffa_col]
            output[:, r, h] = peak_scale(peak_rp, hydro, baseflow)

        if scaleby == "volumes":
            if duration != int(vfa_col):
                print("VFA duration and hydrograph duration different. Conducting partial scaling...")
                vfa_dur = int(vfa_col)
                Q_rp = vfa.loc[r, vfa_col]
                output[:, r, h] = vol_scale(Q_rp, hydro, baseflow, partial=vfa_dur)
            else:
                Q_rp = vfa.loc[r, vfa_col]
                output[:, r, h] = vol_scale(Q_rp, hydro, baseflow)

        if scaleby == "peaks&volumes":
            if duration != int(vfa_col):
                print("VFA duration and hydrograph duration different. Conducting partial scaling...")
                vfa_dur = int(vfa_col)
                peak_rp = ffa.loc[r, ffa_col]
                Q_rp = vfa.loc[r, vfa_col]
                output[:, r, h] = pav_scale(peak_rp, Q_rp, hydro, baseflow, partial=vfa_dur)
            else:
                peak_rp = ffa.loc[r, ffa_col]
                Q_rp = vfa.loc[r, vfa_col]
                output[:, r, h] = pav_scale(peak_rp, Q_rp, hydro, baseflow)

        if scaleby == "balanced":
            peak_rp = ffa.loc[r, ffa_col]
            vfa_rp = vfa.loc[r, balanced_cols]
            output[:, r, h] = balanced_scale(peak_rp, vfa_rp, hydro, balanced_idx, peak_type, baseflow, plot_construct)

    print(f"Scaling of Hydrograph {hydro_sel[h]} complete.")

# Export Hydrograph Ordinates
print("Exporting hydrographs...")
if scaleby == "volumes":
    scale_str = scaleby
else:
    scale_str = f"{scaleby}_{ffa_col}"

for h in range(0, hydro_in.shape[1]):
    df = pd.DataFrame(output[:, :, h], columns=rps_in)
    ofile = f"output/{site}_{hydro_sel[h]}_{scale_str}_scaled.csv"
    df.to_csv(ofile, index=True)
    print(f"   Hydrographs exported to {ofile}")

## Check ##
if not os.path.isdir("check"):
    os.mkdir("check")

# Check values and stats:
hydro_vals, hydro_stats = check_hydrographs(output, ffa[ffa_col], vfa, rps_in, hydro_in)
if pmf_peak is not None:
    pmf_file = f"check/{site}_{scale_str}_PMF_summary.csv"
    hydro_pmf = pd.DataFrame(check_pmf(hydro_vals,pmf_peak,pmf_volumes))
    hydro_pmf.index = list(rps_in)
    hydro_pmf.columns = hydro_sel[0:len(hydro_pmf.columns)]
    hydro_pmf.to_csv(pmf_file)

# Export slices of data
# Create header:
head = ["peak"]
for dur in vfa_durs:
    head.append(str(dur))

# By hydrograph, checking peaks and all durations in raw vfa
for h in range(0, hydro_in.shape[1]):
    val_file = f"check/{site}_{hydro_sel[h]}_{scale_str}_summary_values.csv"
    output_vals = pd.DataFrame(hydro_vals[:, :, h])
    output_vals.index = head
    output_vals.columns = list(rps_in)
    output_vals.to_csv(val_file)

    stat_file = f"check/{site}_{hydro_sel[h]}_{scale_str}_summary_stats.csv"
    output_stats = pd.DataFrame(hydro_stats[:, :, h])
    output_stats.index = head
    output_stats.columns = list(rps_in)
    output_stats.to_csv(stat_file)

# By peaks an durations, including all hydrographs
for d in range(0, len(head)):
    durval_file = f"check/{site}_{scale_str}_{head[d]}_values.csv"
    durval_dat = pd.DataFrame(hydro_vals[d, :, :])
    durval_dat.index = list(rps_in)
    durval_dat.columns = list(hydro_sel[0:len(durval_dat.columns)])
    durval_dat.to_csv(durval_file)

    durstat_file = f"check/{site}_{scale_str}_{head[d]}_stats.csv"
    durstat_dat = pd.DataFrame(hydro_stats[d, :, :])
    durstat_dat.index = list(rps_in)
    durstat_dat.columns = list(hydro_sel[0:len(durstat_dat.columns)])
    durstat_dat.to_csv(durstat_file)

## Plots ##
print("Preparing plots")
if not os.path.isdir("plots"):
    os.mkdir("plots")
## Plot the new scaled hydrographs

# Get line colors and styles
linecolors,linestyles = get_linestyles(len(rps_in))

for h in range(0, hydro_in.shape[1]):
    
    # Initialize figure
    init_hydrograph_plot(duration)
    plt.title(hydro_sel[h])
    
    # Plot frequency hydrographs
    for r in range(0, len(rps_in)):
        plt.plot(range(0, duration), output[:, r, h], linewidth=0.5, linestyle=linestyles[r],color = linecolors[r],label=rps_in[r])
        
    # Plot input hydrograph
    plt.plot(range(0, duration), hydro_in[:, h], color='0', linewidth=1,label="Input Hydrograph")
    
    # Add legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    plt.legend(title="Return Period (yr).",bbox_to_anchor=(1, 0.5), loc='center left',prop={'size': 10})

    # Save figure
    hfig_file = f"plots/{site}_{hydro_sel[h]}_{scale_str}_scaled.jpg"
    plt.savefig(hfig_file, bbox_inches='tight', dpi=300)

## Plot all of the hydrographs together

# Get line colors and styles
linecolors,linestyles = get_linestyles(hydro_in.shape[1],len(rps_in))

# Initialize figure
init_hydrograph_plot(duration)
plt.title("All")

for h in range(0, hydro_in.shape[1]):
    
    # Plot frequency hydrographs
    for r in range(0, len(rps_in)):
        if r==len(rps_in)-1:
            plt.plot(range(0, duration), output[:, r, h], linewidth=0.5, linestyle=linestyles[r],color=linecolors[h],label=hydro_sel[h])
        else:
            plt.plot(range(0, duration), output[:, r, h], linewidth=0.5, linestyle=linestyles[r],color=linecolors[h],label = "_nolegend_")
        
# Plot input hydrographs
for h in range(0, hydro_in.shape[1]):
    if h==0:
        plt.plot(range(0, duration), hydro_in[:, h], color='0', linewidth=1,label="Input Hydrographs")
    else:
        plt.plot(range(0, duration), hydro_in[:, h], color='0', linewidth=1,label="_nolegend_")

# Add legend
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
plt.legend(title="Hydrographs Shapes",bbox_to_anchor=(1, 0.5), loc='center left',prop={'size': 10})

# Save figure
xfile = f"plots/{site}_all_{scale_str}_scaled.jpg"
plt.savefig(xfile, bbox_inches='tight', dpi=300)

## Plot against curves
# Peaks
init_freq_plot()
plt.title(f"{scale_str}_peaks")
if pmf_peak is not None:
    plt.plot(aep,[pmf_peak]*len(aep),color="grey",linestyle="dashed",label="PMF Peak")

for h in range(0, hydro_in.shape[1]):
    plt.plot(aep, hydro_vals[0, :, h], marker='o', linewidth=0,label=hydro_sel[h])
plt.plot(ffa.iloc[:, 0], ffa.iloc[:, 2], color='0', linewidth=2, label="Peak Frequency")
plt.plot(ffa.iloc[:, 0], ffa.iloc[:, 3], color='0', linewidth=0.5, linestyle="dashdot", alpha=1, label='_nolegend_')
plt.plot(ffa.iloc[:, 0], ffa.iloc[:, 4], color='0', linewidth=0.5, linestyle="dashdot", alpha=1, label='_nolegend_')

plt.legend(loc="lower right", title="Hydrograph", prop={'size': 7})
xfile = f"plots/{site}_check_peaks_{scale_str}.jpg"
plt.savefig(xfile, bbox_inches='tight', dpi=300)

# Volumes
for d in range(0, len(vfa_durs)):
    init_freq_plot()
    plt.title(f"{scale_str}_{vfa_durs[d]/24}-day")
    if pmf_volumes is not None:
        plt.plot(aep, [pmf_volumes[d]] * len(aep), color="grey", linestyle="dashed", label="PMF Volume")

    for h in range(0, hydro_in.shape[1]):
        plt.plot(aep[0:len(rp_sel)], hydro_vals[d + 1, :, h], marker='o', linewidth=0,label=hydro_sel[h])

    plt.plot(vfa.iloc[:, 0], vfa.loc[:, str(vfa_durs[d])], color='0', linewidth=2,label=f"{vfa_durs[d]/24}-day Frequency")
    plt.legend(loc="lower right", title="Hydrograph", prop={'size': 7})
    xfile = f"plots/{site}_check_{vfa_durs[d]}_{scale_str}.jpg"
    plt.savefig(xfile, bbox_inches='tight', dpi=300)
