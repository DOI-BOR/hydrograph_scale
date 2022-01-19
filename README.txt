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

The hydrograph_scaling.docx describes in detail the mathematics of the scaling script.

All contributions will be licensed as Creative Commons Zero (CC0).