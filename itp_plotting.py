#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
itp_plotting.py

Code to plot temperature or salinity vs depth for
data from Ice Tethered Profilers
Data available from Woods Hole Oceanographic Institute:
http://www.whoi.edu/page.do?pid=20781

Created on Tuesday, April 9 15:50 2019

@author: Mikhail
"""

import numpy as np
import matplotlib.pyplot as plt
# Pandas used for reading in csv files
import pandas as pd
# for creating a dictionary of lists
from collections import defaultdict
# For checking if a file exists
import os

pwd = "/home/mschee/Documents/2509_Arctic/ITP_data/"

# Constructs a specific filepath for an itp data file
def make_filepath(itp_num, instance, pwd):
    itp_num_str = str(itp_num).zfill(3)
    inst_str = str(instance).zfill(4)
    filepath = pwd + 'itp' + itp_num_str + 'data/itp' + str(itp_num) + 'grd' + inst_str + '.dat'
    return filepath

# Puts data into arrays
def in_array(df, x_str):
    x = df[x_str].values
    return x

# Read in csv file and return array of data
#   (year day pressure(dbar) temperature(C) salinity oxygen(umol/kg))
def import_data_file(itp_num, instance, pwd):
    # generate filepath
    filepath = make_filepath(itp_num, instance, pwd)
    # import to a data frame
    df = pd.read_csv(filepath, engine='python', delim_whitespace=True, header=0, skiprows=2, skipfooter=1)
    return df

# Returns pressure, temperature, and salinity arrays from one file
def p_t_s(itp_num, instance, pwd):
    this_df = import_data_file(itp_num, instance, pwd)
    pressures = in_array(this_df, 'pressure(dbar)')
    temps = in_array(this_df, 'temperature(C)')
    salinities = in_array(this_df, 'salinity')
    return pressures, temps, salinities

def force_aspect(ax, ratio=1):
    xleft, xrite = ax.get_xlim()
    ybott, ytopp = ax.get_ylim()
    ax.set_aspect(abs((xrite-xleft)/(ybott-ytopp))*ratio)

# plots data from a list of itps and instances
def plot_itps(itp_num_inst_list, pwd):
    # Find itp number
    itp_num = itp_num_inst_list[0][0]
    plt_title = 'Ice Thetherd Profiler (ITP' + str(itp_num) + ')'
    # Define plot aspects
    fg, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    fg.suptitle(plt_title)
    ax1.set_title('Temperature Profile')
    ax1.set_xlabel('Temperature (C)')
    ax1.set_ylabel('Pressure (dbar)')
    ax2.set_title('Salinity Profile')
    ax2.set_xlabel('Salinity (psu)')
    ax2.set_ylabel('Pressure (dbar)')
    # Plot temp and salintiy profiles for each instance in list
    for ele in itp_num_inst_list:
        num = ele[0]
        inst = ele[1]
        pressures, temps, salinities = p_t_s(num, inst, pwd)
        ax1.plot(temps, pressures, '-', label=str(inst))
        ax2.plot(salinities, pressures, '-', label=str(inst))
    #plt.legend()
    force_aspect(ax1, ratio=2)
    force_aspect(ax2, ratio=2)
    plt.gca().invert_yaxis()
    #plt.tight_layout() # makes suptitle intersect with plots
    plt.show()

def make_itp_list(itp_num, min_inst, max_inst, pwd):
    itp_list = []
    for i in range(min_inst, max_inst):
        # make filepath
        filepath = make_filepath(itp_num, i, pwd)
        if os.path.isfile(filepath):
            itp_list.append((itp_num, i))
    return itp_list

# Adds to a list in a dictionary
def add_to_dict_list(dict1, ind_data, dep_data):
    len_ = len(ind_data)
    for i in range(len_):
        dict1[ind_data[i]].append(dep_data[i])
    return dict1

# finds the mean and standard deviation of a dictionary
def find_mean_profile(dict2):
    ind_data = sorted([*dict2])     # get list of keys
    dep_data = []
    std_data = []
    for i in range(len(ind_data)):
        p_mean = np.mean(dict2.get(ind_data[i]))
        p_stdv = np.std(dict2.get(ind_data[i]))
        dep_data.append(p_mean)
        std_data.append(p_stdv)
    return ind_data, dep_data, std_data

# finds and returns array of mean temp and salinty vs pressure
#   from a list of itps
def find_itp_mean(itp_num_inst_list, pwd):
    t_dict = defaultdict(list)
    s_dict = defaultdict(list)
    for ele in itp_num_inst_list:
        num = ele[0]
        inst = ele[1]
        pressures, temps, salinities = p_t_s(num, inst, pwd)
        t_dict = add_to_dict_list(t_dict, pressures, temps)
        s_dict = add_to_dict_list(s_dict, pressures, salinities)
    mean_p_t, mean_t_profile, std_t = find_mean_profile(t_dict)
    mean_p_s, mean_s_profile, std_s = find_mean_profile(s_dict)
    return mean_p_t, mean_t_profile, std_t, mean_p_s, mean_s_profile, std_s

# plots mean and standard deviation from itp list
def plot_itp_mean(itp_num_inst_list, pwd):
    # Find itp number
    itp_num = itp_num_inst_list[0][0]
    plt_title = 'Ice Thetherd Profiler (ITP' + str(itp_num) + ')'
    # Define plot aspects
    fg, (ax1, ax2) = plt.subplots(1,2, sharey=True)
    fg.suptitle(plt_title)
    ax1.set_title('Mean Temperature Profile')
    ax1.set_xlabel('Temperature (C)')
    ax1.set_ylabel('Pressure (dbar)')
    ax2.set_title('Mean Salinity Profile')
    ax2.set_xlabel('Salinity (psu)')
    ax2.set_ylabel('Pressure (dbar)')
    # Find mean temp and salinity profiles
    p_t, t_prof, t_std, p_s, s_prof, s_std = find_itp_mean(itp_num_inst_list, pwd)
    # Plot 2 standard deviations as shaded area
    t_profile = np.asarray(t_prof)
    t_stdev = np.asarray(t_std)
    ax1.fill_betweenx(p_t, t_profile-2.*t_stdev, t_profile+2.*t_stdev, label='2 stdv', alpha=0.25)
    s_profile = np.asarray(s_prof)
    s_stdev = np.asarray(s_std)
    ax2.fill_betweenx(p_s, s_profile-2.*s_stdev, s_profile+2.*s_stdev, label='2 stdv', alpha=0.25)
    # Plot temp and salintiy mean profiles
    ax1.plot(t_prof, p_t, '-', label='Mean')
    ax1.legend()
    ax2.plot(s_prof, p_s, '-', label='Mean')
    ax2.legend()
    force_aspect(ax1, ratio=2)
    force_aspect(ax2, ratio=2)
    plt.gca().invert_yaxis()
    #plt.tight_layout() # makes suptitle intersect with plots
    plt.show()

itp_list89 = make_itp_list(89, 168, 678, pwd)

#plot_itps(itp_list89, pwd)

plot_itp_mean(itp_list89, pwd)