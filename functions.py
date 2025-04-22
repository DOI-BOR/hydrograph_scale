# -*- coding: utf-8 -*-
"""
Hydrograph Scaling Script, Functions  (v1)
@author: USBR (tclarkin 2021)

The supporting functions used in hydrograph scaling and hydorgraph pre scripts
"""
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import probscale
import datetime as dt
import dataretrieval.nwis as nwis


def nwis_import(site, dtype, start=None, end=None, wy="WY"):
    """
    Imports flows from NWIS site
    :param site: str, USGS site number
    :param dtype: str, "dv" or "iv"
    :param start: str, start date (default is None)
    :param end: str, end date (default is None)
    :param wy: str, "WY
    :return: dataframe with date index, dates, flows, month, year and water year
    """
    if dtype == "dv":
        parameter = "00060_Mean"
    elif dtype == "iv":
        parameter = "00060"

    data = pd.DataFrame()

    if (start != None) & (end != None):
        try:
            data = nwis.get_record(sites=site, start=start, end=end, service=dtype, parameterCd="00060")
        except ValueError:
            data["flow"] = np.nan
    else:
        if (start == None) & (end == None):
            try:
                data = nwis.get_record(sites=site, start="1800-01-01", service=dtype, parameterCd="00060")
            except ValueError:
                data["flow"] = np.nan
        else:
            if end == None:
                try:
                    data = nwis.get_record(sites=site, start=start, end="3000-01-01", service=dtype,
                                           parameterCd="00060")
                except ValueError:
                    data["flow"] = np.nan
            if start == None:
                try:
                    data = nwis.get_record(sites=site, start="1800-01-01", end=end, service=dtype, parameterCd="00060")
                except ValueError:
                    data["flow"] = np.nan

    data.index = pd.to_datetime(data.index, utc=True)
    # data = data.tz_localize("UTC")
    end = data.index.max()
    start = data.index.min()

    if dtype == "dv":
        date_index = pd.to_datetime(pd.date_range(start, end, freq="D"), utc=True)
    elif dtype == "iv":
        date_index = pd.to_datetime(pd.date_range(start, end, freq="15T"), utc=True)

    out = pd.DataFrame(index=date_index)
    # out = out.tz_localize("UTC")
    out["flow"] = out.merge(data[parameter], left_index=True, right_index=True, how="left")

    out.loc[out["flow"] == -999999, "flow"] = np.nan

    # Add year, month and wy
    out["year"] = pd.DatetimeIndex(out.index).year
    out["month"] = pd.DatetimeIndex(out.index).month
    out["wy"] = out["year"]
    if wy == "WY":
        out.loc[out["month"] >= 10, "wy"] = out.loc[out["month"] >= 10, "year"] + 1

    return (out)


def analyze_voldur(data, dur):
    """
    This function calculates a rolling mean and then identifies the ann. max. for each WY
    :param data: df, data including at least date, variable
    :param dur: int or "WY", duration to analyze
    :return: df, list of events with date, avg_flow and peak
    """
    var = data.columns[0]
    WYs = data["wy"].unique().astype(int)
    evs = pd.DataFrame(index=WYs)
    if dur == "WY":
        print('Analyzing by WY')
        for wy in WYs:
            if sum(pd.isna(data.loc[data["wy"] == wy, var])) > 365 * 0.1:
                continue
            else:
                if dur == "WY":
                    evs.loc[wy, "annual_sum"] = round(data.loc[data["wy"] == wy, var].sum(), 0)
                    evs.loc[wy, "annual_acft"] = round(data.loc[data["wy"] == wy, var].sum() * 86400 / 43560, 0)
                    evs.loc[wy, "count"] = len(data.loc[data["wy"] == wy, var])
                    max_idx = data.loc[data["wy"] == wy, var].idxmax()
                    evs.loc[wy, "max"] = max_idx
                    evs.loc[wy, f"max_{var}"] = round(data.loc[max_idx, var], 0)
    else:
        dur_data = data[var].rolling(dur, min_periods=int(np.ceil(dur * 0.90))).mean()
        data["wy_shift"] = data["wy"].shift(+(int(dur / 2)))

        for wy in WYs:
            try:
                max_idx = dur_data.loc[data["wy_shift"] == wy].idxmax()
            except ValueError:
                continue
            if pd.isna(max_idx):
                continue
            evs.loc[wy, "start"] = max_idx - dt.timedelta(days=int(dur) - 1)  # place date as start of window
            evs.loc[wy, f"avg_{var}"] = round(dur_data[max_idx], 0)
            evs.loc[wy, "mid"] = max_idx - dt.timedelta(days=int(dur / 2) - 1)  # place date as middle of window
            evs.loc[wy, "end"] = max_idx  # place date as end of window
            evs.loc[wy, "max"] = data.loc[evs.loc[wy, "start"]:evs.loc[wy, "end"], var].idxmax()
            evs.loc[wy, f"max_{var}"] = data.loc[evs.loc[wy, "max"], var]
            evs.loc[wy, "count"] = len(data.loc[data["wy"] == wy, var])

    return (evs)


def hydroclip(hfile, hydro_sel, duration):
    """
    This function opens the file and extracts the selected hydrographs. If hydrograph durations are longer than the
    specified duration, they will be clipped.
    :param hfile: str, .csv file containing hydrographs
    :param hydro_sel: list, columns names of selected hydrographs from hfile
    :param duration: int, desired hydrograph duration in hours
    :return: ND array, selected hydrographs, clipped to duration
    """
    hydro_raw = pd.read_csv(hfile)
    # create output array
    if hydro_sel == "all":
        hydro_sel = hydro_raw.columns
    hydro_clip = np.zeros((duration, len(hydro_sel)))

    # check if hydrographs are of same length as duration
    # if too short, exit
    if len(hydro_raw) < duration:
        sys.exit("Hydrograph duration shorter than critical duration. Check timestep or provide longer hydrographs")

    # if equal, use
    if len(hydro_raw) == duration:
        for h in range(0, len(hydro_sel)):
            if hydro_sel[h] in hydro_raw.columns:
                hydro_clip[:, h] = hydro_raw.loc[:, hydro_sel[h]]
            else:
                print(f"{h} missing from hydrograph column names")
                hydro_clip[:, h] = np.nan
                continue

    # if too long, clip
    else:
        for h in range(0, len(hydro_sel)):
            if hydro_sel[h] in hydro_raw.columns:
                # Find peak critical duration period
                vol_cd = hydro_raw.loc[:, hydro_sel[h]].rolling(duration).mean()
                vol_p_end = vol_cd.idxmax()
                vol_p_beg = vol_p_end - duration + 1
                hydro_clip[:, h] = hydro_raw.loc[vol_p_beg:vol_p_end, hydro_sel[h]].values
            else:
                print(f"{h} missing from hydrograph column names")
                hydro_clip[:, h] = np.nan
                continue
    return hydro_clip, hydro_sel


def importfreq(ffile, vfile, rp_sel):
    """
    This function opens the files containing frequency information, checks for consistent length, and puts in format
    used by script.

    :param ffile: str, .csv file containing FFA
    :param vfile: str, .csv file containing VFA
    :param rp_sel: list, selected return periods
    :param ffa_col: int or str, column index or column name from ffile
    :param vfa_col: int or str, column index or column name from vfile
    :return: 1D arrays: vfa, ffa, rps_in and aep
    """
    # Import FFA & VFA
    ffa_raw = pd.read_csv(ffile)
    vfa_raw = pd.read_csv(vfile)

    # Define rps, check that rps for vfa and ffa are identical
    if all(ffa_raw.iloc[:, 0] == vfa_raw.iloc[:, 0]):
        print("Return Periods Match: extracting selected return periods")
        rps = ffa_raw.loc[:, "rp"]
        if rp_sel == "all":
            ffa = ffa_raw
            vfa = vfa_raw
            aep = (1 / rps) * 100
        else:
            rps = np.array(rps.loc[rps.isin(rp_sel),])
            ffa = ffa_raw.loc[rps.isin(rp_sel),].reset_index(drop=True)
            vfa = vfa_raw.loc[rps.isin(rp_sel),].reset_index(drop=True)
            aep = (1 / rps) * 100
    else:
        sys.exit("Return Periods for FFA and VFA do not match; please correct and try again.")
    return ffa, vfa, rps, aep


def adj_baseflow(output, hydro, vol_method=0):
    """
    This function adjusts beginning and ending flows to transition smoothly from baseflow (the original ending and
    beginning flows on the input hydrograph).

    :param output: 1D array, the already scaled hydrograph
    :param hydro: 1D array, the original, unscaled hydrograph
    :param vol_method: int, options for handling volume corrections...will only be used if > 0.01% error in volume
        0 - None
        1 - Correct volume outside of baseflow periods
        2 - Correct volume outside of baseflow periods and peak
    :return: 1D array, scaled hydrograph with baseflow
    """
    # First, determine duration:
    duration = len(output)
    v_rp = output.sum()
    start = int(duration / 12)

    # Second, scale starting baseflow
    # Check if scaled starting flow are greater than raw starting flow
    t_ratio = output[start] / hydro[start]  # transition ratio
    # If smaller, downscale
    if output[start] < hydro[start]:
        t_steps = np.linspace(t_ratio, t_ratio, start)  # transition timesteps
    # If equal or larger, linear ratios from 1 to ratio
    else:
        t_steps = np.linspace(1, t_ratio, start)

    # Scale start of hydrograph for baseflow transition
    output[0:start] = np.ceil(t_steps * hydro[0:start])

    # Third, scale endign baseflow
    end = duration - start
    # Check if scaled starting flow are greater than raw starting flow
    t_ratio = output[end] / hydro[end]  # transtion ratio
    # If smaller, downscale
    if output[end] < hydro[end]:
        t_steps = np.linspace(t_ratio, t_ratio, start)  # transition timesteps
    # If equal or larger, linear ratios from 1 to ratio
    else:
        t_steps = np.linspace(t_ratio, 1, start)

    # Scale end of hydrograph for baseflow transition
    output[end:duration] = np.ceil(t_steps * hydro[end:duration])

    # Check if volume correction needed
    if abs(check_error(output.sum(), v_rp)) > 0.0001:  # greater than 0.0001 or 0.01% in CFS
        print(f"Baseflow adjusted volumes do not agree: {round(abs(check_error(output.sum(), v_rp)), 6)} error")
        if vol_method > 0:
            print(f"Correcting volumes using method {vol_method}")
            # Check total volume and baseflow volume
            v_bf = sum([sum(output[0:start]), sum(output[end:duration])])

            if vol_method == 1:
                # Correct volume, only modifying flows outside of the baseflow transition
                v_ratio = (v_rp - v_bf) / sum(output[start:end])
                output[start:end] = np.ceil(v_ratio * output[start:end])

            if vol_method == 2:
                # Correct volume,  modifying flows outside of the baseflow AND peak
                # Add flow of the peak to volume
                peak_idx = output.argmax()
                peak_rp = output[peak_idx]
                v_bfp = v_bf + peak_rp

                # Correct volume, only modifying flows outside of the baseflow AND peak
                v_ratio = (v_rp - v_bfp) / (sum(output[start + 1:end - 1]) - peak_rp)
                output[start + 1:peak_idx - 1] = np.ceil(v_ratio * output[start + 1:peak_idx - 1])
                output[peak_idx + 1:end - 1] = np.ceil(v_ratio * output[peak_idx + 1:end - 1])
                output[peak_idx] = peak_rp

                # Check volume again...
                if abs(check_error(output.sum(), v_rp)) > 0.0001:  # greater than 0.0001 or 0.01% in CFS
                    print(f"Volumes do not agree: {round(abs(check_error(output.sum(), v_rp)), 6)} error")
                    print("Baseflow adjustment not recommended.")
                else:
                    print("Volume correction successful.")
                # Check peak again
                if abs(check_error(output.max(), peak_rp)) > 0.0001:  # greater than 0.0001 or 0.01% in CFS
                    print(f"Peaks do not agree: {round(abs(check_error(output.max(), peak_rp)), 6)} error")
                    print("Baseflow adjustment not recommended.")
                else:
                    print("Peak preservation successful.")
        else:
            print("Volume correction recommended.")
    else:
        print("Baseflow adjustment successful.")

    return output


def check_error(estimate, real):
    """
    Calculates error (decimal)
    :param estimate: numeric, estimated value
    :param real: numeric, real value
    :return: numeric, error (decimal)
    """
    error = (estimate - real) / real
    return error


def peak_scale(peak_rp, hydro, baseflow=False):
    """
    The concept of peak scaling is simple: multiply the entire hydrograph by the ratio of the desired peak over the
    input hydrograph peak.

    :param peak_rp: float, the desired peak
    :param hydro: 1D array, the input hydrograph
    :param baseflow: boolean, option to smooth beginning and end of hydrographs for baseflow transition
    :return: 1D array, the scaled hydrograph
    """
    # First, scale by peak
    peak_hydro = hydro.max()  # peak of hydrograph
    p_ratio = peak_rp / peak_hydro
    output = np.ceil(p_ratio * hydro)

    # Second, fix baseflow transition (if selected)
    if baseflow:
        output = adj_baseflow(output, hydro, vol_method=0)

    # Third, check if peak is correct
    if abs(check_error(output.max(), peak_rp)) > 0.0001:  # greater than 0.0001 or 0.01%
        print(f"Peaks do not agree: {round(abs(check_error(output.max(), peak_rp)), 6)} error")
    else:
        print("Peak scaling successful.")
    return output


def vol_scale(Q_rp, hydro, baseflow=False, partial=False):
    """
    The concept of volume scaling is simple: multiply the entire hydrograph by the ratio of the desired volume over the
    input hydrograph volume.

    :param peak_rp: float, the desired peak
    :param hydro: 1D array, the input hydrograph
    :param baseflow: boolean, option to smooth beginning and end of hydrographs for baseflow transition
    :param partial: boolean, option to scale entire hydrograph based on calculations for portion of hydrograph
    :return: 1D array, the scaled hydrograph
    """
    # First, scale by volume
    # If hydrograph duration = volume duration
    if partial == False:
        Q_hydro = np.mean(hydro)
    # If not, determine mean flow for duration
    else:
        hydro_dur = pd.DataFrame(hydro)
        vol_dur = hydro_dur.rolling(partial).mean()
        vol_dur_end = int(vol_dur.idxmax() + 1)
        vol_dur_beg = int(vol_dur_end - int(partial))
        Q_hydro = np.mean(hydro[vol_dur_beg:vol_dur_end])

    v_ratio = Q_rp / Q_hydro
    output = np.ceil(v_ratio * hydro)

    # Second, fix baseflow transition (if selected) and correct volumes
    if baseflow:
        output = adj_baseflow(output, hydro, vol_method=1)

    # Third, check if volume is correct
    if abs(check_error(output.mean(), Q_rp)) > 0.0001:  # greater than 0.0001 or 0.01% in CFS
        print(f"Volumes do not agree: {round(abs(check_error(output.mean(), Q_rp)), 6)} error")
    else:
        print("Volume scaling successful.")
    return output


def pav_scale(peak_rp, Q_rp, hydro, baseflow=False, partial=False):
    """
    The concept of peak and volume scaling is scaling the hydrograph to the peak, then correct the rest of the
    hydrograph ordinates to have the correct volume.

    :param peak_rp: float, the desired peak
    :param Q_rp: float, the desired volume (average flow)
    :param hydro: 1D array, the input hydrograph
    :param baseflow: boolean, option to smooth beginning and end of hydrographs for baseflow transition
    :param partial: boolean, option to scale entire hydrograph based on calculations for portion of hydrograph
    :return: 1D array, the scaled hydrograph
    """
    # First, scale by peak
    output = peak_scale(peak_rp, hydro, baseflow=False)

    # Second, use equation 5 to scale volume
    N = len(output)
    q_max = output.max()
    # If hydrograph duration = volume duration
    if partial == False:
        m = (Q_rp * N - q_max * N) / sum(q_max - output)

    # If not, determine mean flow for duration
    else:
        hydro_dur = pd.DataFrame(hydro)
        vol_dur = hydro_dur.rolling(partial).mean()
        vol_dur_end = int(vol_dur.idxmax() + 1)
        vol_dur_beg = int(vol_dur_end - int(partial))
        m = (Q_rp * partial - q_max * partial) / sum(q_max - output[vol_dur_beg:vol_dur_end])

    output = np.ceil(m * (q_max - output) + q_max)

    # Third, fix baseflow transition (if selected) and correct volumes
    if baseflow:
        output = adj_baseflow(output, hydro, vol_method=0)

    # Fix negative flows
    if (any(output < 1)) or (pd.isna(output).any()):
        output[output < 1] = 1
        output[output == np.nan] = 1

    # Fourth, check if peak and volume are correct
    if abs(check_error(output.mean(), Q_rp)) > 0.0001:  # greater than 0.0001 or 0.01% in CFS
        print(f"Volumes do not agree: {round(abs(check_error(output.mean(), Q_rp)), 6)} error")
    else:
        print("Volume scaling successful.")
    if abs(check_error(output.max(), peak_rp)) > 0.0001:  # greater than 0.0001 or 0.01%
        print(f"Peaks do not agree: {round(abs(check_error(output.max(), peak_rp)), 6)} error")
    else:
        print("Peak scaling successful.")
    return output


def peakonly_scale(peak_rp,hydro):
    """
    Scales the hydrograph to the peak, slowly transitioning to no adjustment at the start and end (low flows).

    :param peak_rp: float, the desired peak
    :param hydro: 1D array, the input hydrograph
    :return: 1D array, the scaled hydrograph
    """
    # First, select the hydrograph
    output = hydro.copy()

    # Second, invert equation used in pav_scale
    q_max_idx = output.argmax()
    q_max = output[q_max_idx]
    q_min1 = output[0]
    q_min2 = output[-1]

    m1 = (peak_rp-q_min1)/(q_max-q_min1)
    m2 = (peak_rp-q_min2)/(q_max-q_min2)

    output[:q_max_idx] = np.ceil(m1*(output[:q_max_idx]-q_min1)+q_min1)
    output[q_max_idx:] = np.ceil(m2*(output[q_max_idx:]-q_min2)+q_min2)

    # If any values are less than q_min, set to q_min
    output[output<min(q_min1,q_min2)] = min(q_min1,q_min2)

    # Check if peak is correct

    if abs(check_error(output.max(), peak_rp)) > 0.0001:  # greater than 0.0001 or 0.01%
        print(f"Peaks do not agree: {round(abs(check_error(output.max(), peak_rp)), 6)} error")
    else:
        print("Peak scaling successful.")
    return output

def volumeonly_scale(Q_rp,hydro):
    """
    Scales the hydrograph from the start and end (low flows), while leaving the peak unchanged.

    :param Q_rp: float, the desired average flow (volume)
    :param hydro: 1D array, the input hydrograph
    :return: 1D array, the scaled hydrograph
    """
    # First, select the hydrograph
    output = hydro.copy()

    # Second, invert equation used in pav_scale
    q_max_idx = output.argmax()
    q_max = output[q_max_idx]
    q_min = output.min()
    N = len(output)

    b = (N*(Q_rp-q_min)*(q_max-q_min)-sum((output-q_min)**2))/sum((output-q_min)*(q_max-q_min)-(output-q_min)**2)
    a = (1-b)/(q_max-q_min)

    output = np.ceil(a*(output-q_min)**2+b*(output-q_min)+q_min)

    # Check if volume is correct
    if abs(check_error(output.mean(), Q_rp)) > 0.0001:  # greater than 0.0001 or 0.01%
        print(f"Volumes do not agree: {round(abs(check_error(output.max(), Q_rp)), 6)} error")
    else:
        print("Volume scaling successful.")
    return output

def modified_pav_scale(peak_rp, Q_rp, hydro):
    """
    The concept of peak and volume scaling is scaling the hydrograph to the peak, then correct the rest of the
    hydrograph ordinates to have the correct volume.

    :param peak_rp: float, the desired peak
    :param Q_rp: float, the desired volume (average flow)
    :param hydro: 1D array, the input hydrograph
    :return: 1D array, the scaled hydrograph
    """
    # First, scale peak only
    output = peakonly_scale(peak_rp, hydro)

    # Second, scale volume only
    output = volumeonly_scale(Q_rp,output)

    # Fix negative flows
    if (any(output < 1)) or (pd.isna(output).any()):
        output[output < 1] = 1
        output[output == np.nan] = 1

    # Fourth, check if peak and volume are correct
    if abs(check_error(output.mean(), Q_rp)) > 0.0001:  # greater than 0.0001 or 0.01% in CFS
        print(f"Volumes do not agree: {round(abs(check_error(output.mean(), Q_rp)), 6)} error")
    else:
        print("Volume scaling successful.")
    if abs(check_error(output.max(), peak_rp)) > 0.0001:  # greater than 0.0001 or 0.01%
        print(f"Peaks do not agree: {round(abs(check_error(output.max(), peak_rp)), 6)} error")
    else:
        print("Peak scaling successful.")
    return output

def get_balanced_indices(vfa_durs, hydro):
    """
    This function identifies the beginning and end of the volume frequency windows

    :param vfa_durs: list of int, durations in hours
    :param hydro: 1D array, hydrograph
    :return: 2D array, durations, beginning and ending indices
    """
    # Create output table
    vfa_tab = np.zeros((len(vfa_durs), 3))

    # Define periods for balanced scaling
    hydro = pd.DataFrame(hydro)
    for d in range(0, len(vfa_durs)):
        print(f"Defining window for {vfa_durs[d]} hour duration")
        vol_dur = hydro.rolling(vfa_durs[d]).mean()
        vol_dur_end = int(vol_dur.idxmax() + 1)
        vol_dur_beg = int(vol_dur_end - int(vfa_durs[d]))
        print("Begins at timestep", vol_dur_beg)
        print("Ends at timestep", vol_dur_end)
        vfa_tab[d, :] = [int(vfa_durs[d]), int(vol_dur_beg), int(vol_dur_end)]

    return vfa_tab


def balanced_scale(peak_rp, vfa_rp, hydro, balanced_idx, type="peaked", baseflow=True, plot_construct=False):
    """
    The concept of balanced is similar to peak and volume scaling is simple, except instead of reserving only the peak,
    we reserve successively larger windows for the shorter durations

    :param peak_rp: float, the desired peak
    :param vfa_rp: 1D array, volumes for each duration
    :param hydro: 1D array, the input hydrograph
    :param balanced_idx: 2D array, durations, beginning and ending indices
    :param type: str, "peaked" or "spiked"
    :param baseflow: boolean, option to smooth beginning and end of hydrographs for baseflow transition
    :param plot_construct: boolean, show plots of balanced hydrograph construction
    :return: 1D array, the scaled hydrograph
    """
    if plot_construct:
        # Check for plot construct directory
        if not os.path.isdir("pc"):
            os.mkdir("pc")
        itr = 1

    vfa_durs = balanced_idx[:, 0]

    # First, scale the peak
    peak_idx = hydro.argmax()
    N = duration = len(hydro)
    # Scale to peak
    output = peak_scale(peak_rp, hydro)
    # If type is spiked, only adjust peak hour, adjust back remaining ordinates
    if type == "spiked" or type=="beard":
        output[0:peak_idx] = hydro[0:peak_idx]
        output[peak_idx + 1:duration] = hydro[peak_idx + 1:duration]

    # If plot_construct is true, plot construction of balanced hydrograph
    if plot_construct:
        fig, ax = plt.subplots(figsize=(6.25, 4))
        plt.ylabel('Flow ($ft^3$/s)')
        plt.xlabel('Time (hour)')
        plt.xlim(0, duration)
        ax.grid()

        # Plot locations
        plt.fill_between([peak_idx - 0.5, peak_idx + 0.5], [peak_rp, peak_rp], color="grey", alpha=0.3,
                         label="Peak Location")
        for d in range(0, len(vfa_durs)):
            plt.fill_between([balanced_idx[d, 1], balanced_idx[d, 2]], [vfa_rp[d], vfa_rp[d]], alpha=0.3,
                             label=f"{vfa_durs[d]} Window")

        plt.plot(hydro, color="black", label="0. Raw Hydrograph")
        plt.plot(output, label=f"1. Peak Scaled")
        linestyles = ['--', '-.', ':', '-','--', '-.', ':', '-']

    for d in range(0, len(vfa_durs)):
        print(f"Scaling to {vfa_durs[d]} timestep curve")
        Q_rp = vfa_rp[d]

        # Handle the next duration below peak, simple scaling by volumes (excluding peak)
        if d == 0:
            N = balanced_idx[d, 0]
            if type == "spiked" or type=="beard":
                # Scale volume outside of peak
                Q_0 = output[peak_idx]
                q_h1 = output[int(balanced_idx[d, 1]):peak_idx]
                q_h2 = output[peak_idx + 1:int(balanced_idx[d, 2])]
                Q_h = sum(q_h1) + sum(q_h2)
                v_ratio = (Q_rp * N - Q_0) / Q_h

                # Apply v_ratio
                if type =="spiked":
                    output[0:peak_idx] = np.ceil(v_ratio * output[0:peak_idx])
                    output[peak_idx + 1:duration] = np.ceil(v_ratio * output[peak_idx + 1:duration])
                if type =="beard":
                    output[int(balanced_idx[d,1]):peak_idx] = np.ceil(v_ratio * output[int(balanced_idx[d, 1]):peak_idx])
                    output[peak_idx+1:int(balanced_idx[d,2])] = np.ceil(v_ratio * output[peak_idx+1:int(balanced_idx[d, 2])])

            if type == "peaked":
                # Apply peak and volume scaling to entire hydrograph
                q_max = output.max()
                q_h = output[int(balanced_idx[d, 1]):int(balanced_idx[d, 2])]

                m = (Q_rp * N - q_max * N) / sum(q_max - q_h)
                # Apply equation to modify Q values
                output = np.ceil(m * (q_max - output) + q_max)

        # Handle longer durations
        else:
            # Define durations
            N1 = balanced_idx[d - 1, 0]
            N2 = balanced_idx[d, 0] - N1
            # Define portion of hydrograph BEING scaled
            q_h = np.concatenate((output[int(balanced_idx[d,1]):int(balanced_idx[d-1, 1])],
                                  output[int(balanced_idx[d-1, 2]):int(balanced_idx[d, 2])]))

            # Define portion of hydrograph NOT being scaled
            q_0 = output[int(balanced_idx[d-1, 1]):int(balanced_idx[d-1, 2])]
            q_max = np.min(q_0)
            Q_0 = sum(q_0)

            # If all values == q_max or type = "beard", use simple multiplication only
            if all(q_h == q_max) or type=="beard":
                m = (Q_rp * (N1 + N2) - q_0.sum()) / (q_h.sum())

                output[int(balanced_idx[d,1]):int(balanced_idx[d-1,1])] = np.ceil(
                    m * (output[int(balanced_idx[d,1]):int(balanced_idx[d-1, 1])]))
                output[int(balanced_idx[d-1,2]):int(balanced_idx[d,2])] = np.ceil(
                    m * (output[int(balanced_idx[d-1, 2]):int(balanced_idx[d,2])]))

            # Else, do usual linear factor multiplication
            else:
                # Calculate scaling factor (m)
                m = (Q_rp * (N1 + N2) - q_max * N2 - Q_0) / sum(q_max - q_h)

                # Next, apply equation to modify Q values
                output[0:int(balanced_idx[d - 1, 1])] = np.ceil(
                    m * (q_max - output[0:int(balanced_idx[d - 1, 1])]) + q_max)
                output[int(balanced_idx[d - 1, 2]):duration] = np.ceil(
                    m * (q_max - output[int(balanced_idx[d - 1, 2]):duration]) + q_max)

            # Use smooth transition from time 0 to start of next shorter duration
            if baseflow == False:
                print("Baseflow transition must be applied for balanced method.")
            output = adj_baseflow(output, hydro, vol_method=0)

        if plot_construct:
            itr += 1
            plt.plot(output, linestyle=linestyles[d], label=f"{itr}. {int(balanced_idx[d, 0])} Scaled")

        # Fix negative flows
        if (any(output < 1)) or (pd.isna(output).any()):
            output[output < 1] = 1
            output[output == np.nan] = 1

            if plot_construct:
                itr += 1
                plt.plot(output, linestyle=linestyles[d], label=f"{itr}. {int(balanced_idx[d, 0])} Scaled (Correction)")

        # Finally, check volume
        if abs(check_error(np.mean(output[int(balanced_idx[d, 1]):int(balanced_idx[d, 2])]),
                           vfa_rp[d])) > 0.0001:  # greater than 0.0001 or 0.01% in CFS
            print("Volumes do not agree: ", np.round(
                abs(check_error(np.mean(output[int(balanced_idx[d, 1]):int(balanced_idx[d, 2])]), vfa_rp[d])), 6),
                  " error")
        else:
            print("Volume correction successful.")
        # Recheck peaks
        if abs(check_error(output.max(), peak_rp)) > 0.0001:  # greater than 0.0001 or 0.01%
            print(f"Peaks do not agree: {round(abs(check_error(output.max(), peak_rp)), 6)} error")
        else:
            print("Peak scaling successful.")

    if plot_construct:
        plt.plot(output, color="black", linewidth=0.5, label=f"F. Final Hydrograph")
        plt.legend()
        plt.savefig(f"pc/inputpeak_{int(hydro.max())}_scalepeak_{int(peak_rp)}_plot_construct.jpg", bbox_inches='tight',
                    dpi=300)
        plt.close()

    return output


def calc_hydro(hydro_check, vfa_durs):
    """
    This function calculates the average flow for user specified durations
    :param hydro_check: df, the hydrograph being evaluated
    :param vfa_durs: list, the durations used for calcs
    :return: list, average flow for durations
    """
    hydro_vals = np.zeros((len(vfa_durs) + 1))
    # Check peak
    hydro_vals[0] = round(hydro_check.max(), 0)
    # Loop through volumes
    for d in range(0, len(vfa_durs)):
        hydro_vals[d + 1] = round(hydro_check.rolling(int(vfa_durs[d])).mean().max(), 0)

    return hydro_vals

def check_hydrographs(output, ffa, vfa, rps_in, hydro_in):
    """
    This function checks the peaks and volumes of the hydrograph against the ffa and vfa

    :param output: 3D array, set of scaled hydrographs
    :param ffa: 1D array, peaks in order of rps_in
    :param vfa: 2D array, volumes
    :param rps_in: 1D array, return periods
    :param hydro_in: 2D array, input hydrographs
    :return: table of values, table of statistics
    """
    print("Preparing hydrograph stats")

    # Remove extra columns from vfa, identify durations
    vfa = vfa.iloc[:, 2:]
    vfa_durs = [int(i) for i in vfa.columns]

    # Calculate stats
    hydro_vals = np.zeros((len(vfa_durs) + 1, len(rps_in), hydro_in.shape[1]))
    hydro_stats = np.zeros((len(vfa_durs) + 1, len(rps_in), hydro_in.shape[1]))

    for h in range(0, hydro_in.shape[1]):
        for r in range(0, len(rps_in)):
            hydro_check = pd.DataFrame(output[:, r, h])
            hydro_vals[:, r, h] = calc_hydro(hydro_check, vfa_durs)
            hydro_stats[0, r, h] = np.round(check_error(hydro_vals[0, r, h], ffa[r]), 4)

            for d in range(0, len(vfa_durs)):
                hydro_stats[d + 1, r, h] = round(check_error(hydro_vals[d + 1, r, h], vfa.iloc[r, d]), 4)

    return hydro_vals, hydro_stats


def check_pmf(hydro_vals, peak, volumes):
    """
    This function checks the peaks and volumes of the hydrograph against the ffa and vfa

    :param values: 3D array, set of scaled hydrographs
    :param peak: int, PMF peak
    :param volumes: 1D array, PMF volumes (corresponding to durations in hydro_vals)
    :return: 3D array, boolean
    """
    print("Comparing hydrograph values to PMF")

    dd, rr, hh = hydro_vals.shape
    hydro_pmf = np.empty(hydro_vals.shape, dtype=bool)

    for d in range(0, dd):
        if d == 0:
            comparison = peak
        else:
            comparison = volumes[d - 1]
        for r in range(0, rr):
            for h in range(0, hh):
                hydro_pmf[d, r, h] = bool(hydro_vals[d, r, h] >= comparison)

    hydro_pmf_summary = np.empty((rr, hh), dtype='<U16')
    for r in range(0, rr):
        for h in range(0, hh):
            if (hydro_pmf[0, r, h] == True):
                if (hydro_pmf[1, r, h] == True):
                    hydro_pmf_summary[r, h] = "Peak & Volume"
                else:
                    hydro_pmf_summary[r, h] = "Peak"
            elif (hydro_pmf[1, r, h] == True):
                hydro_pmf_summary[r, h] = "Volume"
            else:
                hydro_pmf_summary[r, h] = "N/A"

    return hydro_pmf_summary


def init_hydrograph_plot(duration):
    """
    Initializes plot for hydrographs

    :param duration: int, duration in hours for x-axis
    :return: figure
    """
    fig, ax = plt.subplots(figsize=(6.25, 4))
    plt.ylabel('Flow ($ft^3$/s)')
    plt.xlabel('Time (hour)')
    plt.xlim(0, duration)
    ax.grid()
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))


def init_freq_plot():
    """
    Initializes plot for frequency comparison

    :return: figure
    """
    fig, ax = plt.subplots(figsize=(6.25, 4))
    plt.ylabel('Flow ($ft^3$/s)')
    plt.xlabel('AEP (%)')
    plt.yscale('log')
    ax.set_xscale('prob')
    plt.xlim(0.001, 1)
    plt.gca().invert_xaxis()
    ax.grid()
    ax.grid(which='minor', linestyle=':', linewidth='0.1', color='black')
    ax.yaxis.set_major_formatter(mpl.ticker.StrMethodFormatter('{x:,.0f}'))


def get_linestyles(n, m=0):
    """
    This function provides line styles for n x m lines being plotted

    :param n:
    :param m:
    :return:
    """
    baselines = ['-', '--', '-.', ':'] * 10
    basecolors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3',
                  '#ff7f00'] * 8  # https://colorbrewer2.org/?type=qualitative&scheme=Set1&n=9

    linestyles = list()
    linecolors = list()

    if m == 0:
        for nn in range(0, n):
            linestyles.append(baselines[nn])
            linecolors.append(basecolors[nn])

    else:
        for nn in range(0, n):
            linecolors.append(basecolors[nn])
        for mm in range(0, m):
            linestyles.append(baselines[mm])

    return linecolors, linestyles
