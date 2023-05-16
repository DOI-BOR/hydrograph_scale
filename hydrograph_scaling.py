# -*- coding: utf-8 -*-
"""
Hydrograph Scaling Script  (v4)
@author: USBR (tclarkin 2021)

This script takes a user supplied flood frequency analysis (FFA) curve and/or volume frequency analysis (VFA) curve and
scales user supplied hydrographs to specified return period values.

The Peaks method simply uses the ratio of the FFA and the hydrograph peak to multiply all values.
The Volumes method simply uses the ratio of the VFA and the hydrograph volume to multiply all values.
The Peaks & Volumes method:
    First, uses Peaks method, as described above.
    Second, solves for a multiplication factor (m) to be applied proportionally as flows are lower than the peak.
    Negative flows are set to 1. If volumes are not within 0.0001 (or 0.01%) of the FFA/VFA value, an error is posted.
The Balanced method:
    Uses same method as Peaks & Volumes, but with incremental windows for each duration.
The Beard method:
    Uses the method for balanced scaled hydrograph developed by Leo Beard.

This script requires:
    ffile = flood frequency curve as .csv, with aep in percent (column 0), return period (column 1), median (column 2), lower (column 3), upper (column 4).
        Should include header: ["aep_pct","rp","median","lower","upper"].
    vfile = volume frequency curve as .csv, with aep in percent (column 0), return period (column 1), and x-day average flow (column 2).
        Should include header: ["aep_pct","rp","{hours of duratoin",etc].
    hfile = hydrographs to be scaled as .csv, flows (cfs; columns 1+). Timestep should be 1 hour.
        Should include header: ["{hydrograph name}", etc.]

The user must:
    - Specify the above listed files
    - Verify column names
    - Verify duration (days)
    - Verify timestep (the program may not run with timesteps other than 1 hour)
    - Select a scaling method ("peaks&volumes", "peaks", "volumes" or "balanced")
        - If "balanced", user will need to specify plot_construct, peak_type, and blanaced_cols
    - Specify return periods to be output (list or "all")
    - Specify which input hydrographs to scale (list or "all")

"""
### Import Libraries ###
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from functions import hydroclip, importfreq
from functions import peak_scale, vol_scale, pav_scale
from functions import get_balanced_indices, balanced_scale, beard_balanced_scale
from functions import calc_hydro, check_hydrographs, check_pmf, init_freq_plot, init_hydrograph_plot, get_linestyles

### Begin User Input ###
# Set Working Directory
os.chdir("C://Users//tclarkin//GitHub//hydrograph_scale")
site = "test"

# Load FFA Curve (should have same return periods as vfa)
ffile = 'example/ffa.csv'
ffa_col = "median"

# Load VFA Curve (should have same return periods as ffa)
vfile = 'example/vfa.csv'
vfa_col = "168"  # specify the column(s), columns names should be durations in hours (1, 2, etc.), used for volume scaling

# Hydrograph duration
duration = 7  # days

# Load Hydrograph(s)
hfile = 'example/test_input_hydrographs.csv'
hydro_sel = "all"  # List of which hydrographs to scale (column name) or "all"

# Select Methods
scaleby = "beard"  # "peaks&volumes", "peaks", "volumes", "balanced", or "beard"
rp_sel = "all"  # Which RPs to produce hydrographs for, must correspond to RPs in ffile vfile or "all"
baseflow = False  # True to include transition from baseflow to scaled hydrograph
plotmin = 10
plotmax = 30000

# Additional settings for Balanced
plot_construct = True  # True to plot evolution of balanced hydrographs.
peak_type = "spiked"  # "peaked" or "spiked" ONLY applicable to "balanced"
balanced_cols = ["48", "168"]  # columns to be used for balance hydrographs

# check against PMF (if None, no calculation will be completed)
pmf_peak = None  # int
pmf_volumes = None  # list

### Begin Script ###
## Setup ##
# Check for directory
if not os.path.isdir("output"):
    os.mkdir("output")

# Define duration:
duration = int(duration * 24)

# Import Data:
print("Importing hydrographs and frequency curves")

# Extract and check frequency curves
ffa, vfa, rps_in, aep = importfreq(ffile, vfile, rp_sel)
vfa_durs = [int(i) for i in vfa.columns[2:]]
all_durs = ["peak"]
for i in vfa_durs:
    all_durs.append(i)

# Clip hydrographs to critical duration window
hydro_in, hydro_sel = hydroclip(hfile, hydro_sel, duration)

# Calculate hydrograph parameters
hydro_in_vals = np.zeros((len(vfa_durs) + 1, hydro_in.shape[1]))
for h in range(hydro_in.shape[1]):
    hydro_check = pd.DataFrame(hydro_in[:, h])
    hydro_in_vals[:, h] = calc_hydro(hydro_check, vfa_durs)

hydro_in_vals = pd.DataFrame(hydro_in_vals)
hydro_in_vals.index = all_durs
hydro_in_vals.columns = hydro_sel
hydro_in_vals.to_csv(f'{hfile.split(".")[0]}_vals.csv', index=True, header=True)

fig, ax = plt.subplots(figsize=(6.25, 4))
for h in range(0, len(hydro_sel)):
    print(h)
    plt.plot(range(0, duration), hydro_in[:, h] * (1 / np.max(hydro_in[:, h])), label=hydro_sel[h])
plt.ylabel("Index Flow")
plt.xlabel("Time (hours)")
plt.legend()
plt.savefig(f'{hfile.split(".")[0]}.jpg', bbox_inches='tight', dpi=300)

# Check that user has input a valid method
if (
        scaleby != "peaks" and scaleby != "volumes" and scaleby != "peaks&volumes" and scaleby != "balanced" and scaleby != "beard") == True:
    sys.exit("Please select scaling method: \"peaks\", \"volumes\", \"peaks&volumes\", \"balanced\" or \"beard\"")
else:
    print("Scaling hydrographs by " + scaleby)

# Create output matrix
output = np.zeros((duration, len(rps_in), hydro_in.shape[1]))

## Analysis ##
# Begin looping through selected hydrographs
for h in range(0, hydro_in.shape[1]):
    print("Beginning hydrograph " + hydro_sel[h])
    hydro = hydro_in[:, h]

    if scaleby == "balanced" or scaleby == "beard":
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

        if scaleby == "beard":
            peak_rp = ffa.loc[r, ffa_col]
            vfa_rp = vfa.loc[r, balanced_cols]
            output[:, r, h] = beard_balanced_scale(peak_rp, vfa_rp, hydro, balanced_idx, plot_construct)

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
    hydro_pmf = pd.DataFrame(check_pmf(hydro_vals, pmf_peak, pmf_volumes))
    hydro_pmf.index = list(rps_in)
    hydro_pmf.columns = hydro_sel[0:len(hydro_pmf.columns)]
    hydro_pmf.to_csv(pmf_file)

# Export slices of data

# By hydrograph, checking peaks and all durations in raw vfa
for h in range(0, hydro_in.shape[1]):
    val_file = f"check/{site}_{hydro_sel[h]}_{scale_str}_summary_values.csv"
    output_vals = pd.DataFrame(hydro_vals[:, :, h])
    output_vals.index = all_durs
    output_vals.columns = list(rps_in)
    output_vals.to_csv(val_file)

    stat_file = f"check/{site}_{hydro_sel[h]}_{scale_str}_summary_stats.csv"
    output_stats = pd.DataFrame(hydro_stats[:, :, h])
    output_stats.index = all_durs
    output_stats.columns = list(rps_in)
    output_stats.to_csv(stat_file)

# By peaks an durations, including all hydrographs
for d in range(0, len(all_durs)):
    durval_file = f"check/{site}_{scale_str}_{all_durs[d]}_values.csv"
    durval_dat = pd.DataFrame(hydro_vals[d, :, :])
    durval_dat.index = list(rps_in)
    durval_dat.columns = list(hydro_sel[0:len(durval_dat.columns)])
    durval_dat.to_csv(durval_file)

    durstat_file = f"check/{site}_{scale_str}_{all_durs[d]}_stats.csv"
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
linecolors, linestyles = get_linestyles(len(rps_in))

for h in range(0, hydro_in.shape[1]):

    # Initialize figure
    init_hydrograph_plot(duration)
    plt.title(hydro_sel[h])

    # Plot frequency hydrographs
    for r in range(0, len(rps_in)):
        plt.plot(range(0, duration), output[:, r, h], linewidth=0.5, linestyle=linestyles[r], color=linecolors[r],
                 label=rps_in[r])

    # Plot input hydrograph
    plt.plot(range(0, duration), hydro_in[:, h], color='0', linewidth=1, label="Input Hydrograph")

    # Add legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
    plt.legend(title="Return Period (yr).", bbox_to_anchor=(1, 0.5), loc='center left', prop={'size': 10})

    # Save figure
    hfig_file = f"plots/{site}_{hydro_sel[h]}_{scale_str}_scaled.jpg"
    plt.savefig(hfig_file, bbox_inches='tight', dpi=300)

## Plot all of the hydrographs together

# Get line colors and styles
linecolors, linestyles = get_linestyles(hydro_in.shape[1], len(rps_in))

# Initialize figure
init_hydrograph_plot(duration)
plt.title("All")

for h in range(0, hydro_in.shape[1]):

    # Plot frequency hydrographs
    for r in range(0, len(rps_in)):
        if r == len(rps_in) - 1:
            plt.plot(range(0, duration), output[:, r, h], linewidth=0.5, linestyle=linestyles[r], color=linecolors[h],
                     label=hydro_sel[h])
        else:
            plt.plot(range(0, duration), output[:, r, h], linewidth=0.5, linestyle=linestyles[r], color=linecolors[h],
                     label="_nolegend_")

# Plot input hydrographs
for h in range(0, hydro_in.shape[1]):
    if h == 0:
        plt.plot(range(0, duration), hydro_in[:, h], color='0', linewidth=1, label="Input Hydrographs")
    else:
        plt.plot(range(0, duration), hydro_in[:, h], color='0', linewidth=1, label="_nolegend_")

# Add legend
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
plt.legend(title="Hydrographs Shapes", bbox_to_anchor=(1, 0.5), loc='center left', prop={'size': 10})

# Save figure
xfile = f"plots/{site}_all_{scale_str}_scaled.jpg"
plt.savefig(xfile, bbox_inches='tight', dpi=300)

## Plot against curves
# Peaks
init_freq_plot()
plt.title(f"{scale_str}_peaks")
if pmf_peak is not None:
    plt.plot(aep, [pmf_peak] * len(aep), color="grey", linestyle="dashed", label="PMF Peak")

for h in range(0, hydro_in.shape[1]):
    plt.plot(aep, hydro_vals[0, :, h], marker='o', linewidth=0, label=hydro_sel[h])

plt.plot(ffa.iloc[:, 0], ffa.iloc[:, 2], color='0', linewidth=2, label="Peak Frequency")
plt.plot(ffa.iloc[:, 0], ffa.iloc[:, 3], color='0', linewidth=0.5, linestyle="dashdot", alpha=1, label='_nolegend_')
plt.plot(ffa.iloc[:, 0], ffa.iloc[:, 4], color='0', linewidth=0.5, linestyle="dashdot", alpha=1, label='_nolegend_')

plt.legend(loc="lower right", title="Hydrograph", prop={'size': 7})
xfile = f"plots/{site}_check_peaks_{scale_str}.jpg"
plt.savefig(xfile, bbox_inches='tight', dpi=300)

# Volumes
for d in range(0, len(vfa_durs)):
    init_freq_plot()
    plt.title(f"{scale_str}_{vfa_durs[d] / 24}-day")
    if pmf_volumes is not None:
        plt.plot(aep, [pmf_volumes[d]] * len(aep), color="grey", linestyle="dashed", label="PMF Volume")

    for h in range(0, hydro_in.shape[1]):
        plt.plot(aep[0:len(rps_in)], hydro_vals[d + 1, :, h], marker='o', linewidth=0, label=hydro_sel[h])

    plt.plot(vfa.iloc[:, 0], vfa.loc[:, str(vfa_durs[d])], color='0', linewidth=2,
             label=f"{int(vfa_durs[d] / 24)}-day Frequency")
    plt.legend(loc="lower right", title="Hydrograph", prop={'size': 7})
    xfile = f"plots/{site}_check_{vfa_durs[d]}_{scale_str}.jpg"
    plt.savefig(xfile, bbox_inches='tight', dpi=300)
