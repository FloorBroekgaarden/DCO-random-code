#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 23 12:23:45 2017

@author: Fabian Gittins & David Perkins
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import jinja2 as jj2
import pandas

"""""""""
Global variables
"""""""""

webpage_titles = []         # Holds the webpage figure titles
webpage_descriptions = []   # Holds the webpage figure descriptions
webpage_imagepaths = []     # Holds the webpage figure image paths

"""
Choices used instead of calling from a source file
"""

nbins = 50                                  # Number of bins to be used in 
                                            # histogram distributions

initialdata_name = './allInitial.dat'       # File path for initial data
mergingdata_name = './allMergers.dat'       # File path for merging data
programOptions_name = './output0/programOptions.txt'# File path for program options

"""
HTML options
"""

html_choice = True                          # Choice on whether to output all 
                                            # figures to an HTML webpage 
html_path = '/home/fgittins/www_html/'      # This holds the path to where the 
                                            # HTML webpage will be created
webpage_name = 'test.html'                  # Name of the webpage to be created
webpage_title = 'Test'                      # Title to appear at the top of the
                                            # webpage

file_extension = '.png'                     # File extension for the plots
output_path = './graphs/'                   # Output path to the directory where
                                            # the plots will be saved

"""
Exclusions
"""

Hubble = True                               # Choice on whether to exclude
                                            # binaries which merge in over a
                                            # Hubble time
RLOF = False                                 # Choice on whether to exclude
                                            # binaries where the secondary
                                            # has undergone RLOF after a CEE
Pessimistic = False                          # Choice on whether to exclude
                                            # binaries which use the optimistic
                                            # CE prescription


"""""""""
General functions
"""""""""

def m1DCOm2DCO_symmetrisation(data):
    """
    DCO component masses are symmetrised using this function, with the 
    symmetrised data being returned as arrays.
    """
    m1_raw = data['M1'].as_matrix()
    m2_raw = data['M2'].as_matrix()
    
    m1 = []
    m2 = []
    
    # While loop and if statements to sort m1 and m2 of the data into their
    # symmetrised lists
    for j in range(len(m1_raw)):
        if m1_raw[j] < m2_raw[j]:
            m1.append(m2_raw[j])
            m2.append(m1_raw[j])
        else:
            m1.append(m1_raw[j])
            m2.append(m2_raw[j])
    
    m1DCO_sym = np.array(m1)
    m2DCO_sym = np.array(m2)
    
    return m1DCO_sym, m2DCO_sym

def histogram2d_mask(x, y, Nbins):
    """
    Takes the inputs for the x and y axes, x and y respectively, and returns the
    associated masked 2D histograms required values to create a 2D histogram 
    plot, with Nbins being an integer or array defining the number of bins or 
    bin edges.
    """
    Hist, xedges, yedges = np.histogram2d(x, y, bins = Nbins)
    
    # Hist needs to be rotated and flipped
    Hist = np.rot90(Hist)
    Hist = np.flipud(Hist)
    # Adds on an extremely small value to avoid log(0) errors, these are later 
    # masked
    Hist_noZeros = Hist + 1e-11
    
    # Masking zeros
    Hmasked = np.ma.masked_where(Hist_noZeros <= 1e-10, Hist_noZeros)
    # Masks pixels with a value less than 1e-10
    
    return Hmasked, xedges, yedges

def summaryTable(mergingData_raw):
    # Creating truth tables for DCOs which belong to a specific population
    BNSall_truth = np.logical_and(mergingData_raw['stellarType1'] == 13, 
                                  mergingData_raw['stellarType2'] == 13)
    BHNSall_truth = (np.logical_and(mergingData_raw['stellarType1'] == 14, 
                                    mergingData_raw['stellarType2'] == 13) + 
                    np.logical_and(mergingData_raw['stellarType1'] == 13, 
                                   mergingData_raw['stellarType2'] == 14))
    BBHall_truth = np.logical_and(mergingData_raw['stellarType1'] == 14, 
                                  mergingData_raw['stellarType2'] == 14)
    
    # Pulling out the specific DCO population data
    BNSall_data = mergingData_raw[BNSall_truth]
    BHNSall_data = mergingData_raw[BHNSall_truth]
    BBHall_data = mergingData_raw[BBHall_truth]
    
    total_DCO = len(mergingData_raw)
    total_BNS = len(BNSall_data)
    total_BHNS = len(BHNSall_data)
    total_BBH = len(BBHall_data)

    
    # Checks to see if binaries merge in a Hubble time
    HubbleFlag_data = mergingData_raw['mergesInHubbleTimeFlag']
    # Creates a truth table of events merging in less than a Hubble time
    HubbleFlag_truth = np.logical_not(HubbleFlag_data == 0)
    # Pulls out data that has mergers in less than a Hubble time
    merging_data1 = mergingData_raw[HubbleFlag_truth]
    
    # Creating truth tables for DCOs which belong to a specific population
    BNSHubble_truth = np.logical_and(merging_data1['stellarType1'] == 13, 
                                     merging_data1['stellarType2'] == 13)
    BHNSHubble_truth = (np.logical_and(merging_data1['stellarType1'] == 14, 
                                       merging_data1['stellarType2'] == 13) + 
                        np.logical_and(merging_data1['stellarType1'] == 13, 
                                       merging_data1['stellarType2'] == 14))
    BBHHubble_truth = np.logical_and(merging_data1['stellarType1'] == 14, 
                                     merging_data1['stellarType2'] == 14)
    
    # Pulling out the specific DCO population data
    BNSHubble_data = merging_data1[BNSHubble_truth]
    BHNSHubble_data = merging_data1[BHNSHubble_truth]
    BBHHubble_data = merging_data1[BBHHubble_truth]
    
    Hubble_DCO = len(merging_data1)
    Hubble_BNS = len(BNSHubble_data)
    Hubble_BHNS = len(BHNSHubble_data)
    Hubble_BBH = len(BBHHubble_data)
    
    
    RLOFFlag_data = merging_data1['doubleCommonEnvelopeFlag']
    # Creates truth table of events which don't undergo RLOF immediately after 
    # a CEE
    RLOFFlag_truth = np.logical_not(RLOFFlag_data == 1)
    # Pulls out data which don't undergo RLOF immediately after a CEE
    merging_data2 = merging_data1[RLOFFlag_truth]
    
    # Creating truth tables for DCOs which belong to a specific population
    BNSRLOF_truth = np.logical_and(merging_data2['stellarType1'] == 13, 
                                   merging_data2['stellarType2'] == 13)
    BHNSRLOF_truth = (np.logical_and(merging_data2['stellarType1'] == 14, 
                                     merging_data2['stellarType2'] == 13) + 
                        np.logical_and(merging_data2['stellarType1'] == 13, 
                                       merging_data2['stellarType2'] == 14))
    BBHRLOF_truth = np.logical_and(merging_data2['stellarType1'] == 14, 
                                   merging_data2['stellarType2'] == 14)
    
    # Pulling out the specific DCO population data
    BNSRLOF_data = merging_data2[BNSRLOF_truth]
    BHNSRLOF_data = merging_data2[BHNSRLOF_truth]
    BBHRLOF_data = merging_data2[BBHRLOF_truth]
    
    RLOF_DCO = len(merging_data2)
    RLOF_BNS = len(BNSRLOF_data)
    RLOF_BHNS = len(BHNSRLOF_data)
    RLOF_BBH = len(BBHRLOF_data)
    
    
    OptimisticFlag_data = merging_data2['optimisticCEFlag']
    # Creates truth table of events which don't use the optimistic CE 
    # prescription
    PessimisticFlag_truth = np.logical_not(OptimisticFlag_data == 1)
    # Pulls out data which don't use the optimistic CE prescription
    merging_data = merging_data2[PessimisticFlag_truth]
    
    # Creating truth tables for DCOs which belong to a specific population
    BNSPessimistic_truth = np.logical_and(merging_data['stellarType1'] == 13, 
                                          merging_data['stellarType2'] == 13)
    BHNSPessimistic_truth = (np.logical_and(merging_data['stellarType1'] == 14, 
                                            merging_data['stellarType2'] == 13) + 
                            np.logical_and(merging_data['stellarType1'] == 13, 
                                           merging_data['stellarType2'] == 14))
    BBHPessimistic_truth = np.logical_and(merging_data['stellarType1'] == 14, 
                                          merging_data['stellarType2'] == 14)
    
    # Pulling out the specific DCO population data
    BNSPessimistic_data = merging_data[BNSPessimistic_truth]
    BHNSPessimistic_data = merging_data[BHNSPessimistic_truth]
    BBHPessimistic_data = merging_data[BBHPessimistic_truth]
    
    Pessimistic_DCO = len(merging_data)
    Pessimistic_BNS = len(BNSPessimistic_data)
    Pessimistic_BHNS = len(BHNSPessimistic_data)
    Pessimistic_BBH = len(BBHPessimistic_data)
    
    return [[total_DCO, total_BNS, total_BHNS, total_BBH], 
            [Hubble_DCO, Hubble_BNS, Hubble_BHNS, Hubble_BBH], 
            [RLOF_DCO, RLOF_BNS, RLOF_BHNS, RLOF_BBH], 
            [Pessimistic_DCO, Pessimistic_BNS, Pessimistic_BHNS, Pessimistic_BBH]]

def save_and_add_to_html(plot_title, plot_description):
    """
    Saves each plot and, if necessary, adds them to the HTML webpage.
    """
    plt.savefig(output_path + plot_title + file_extension)

    if html_choice == True:
        plt.savefig(html_path + plot_title + '.png')
        webpage_titles.append(plot_title)
        webpage_descriptions.append(plot_description)
        webpage_imagepaths.append(plot_title + '.png')


"""""""""
Plotting functions
"""""""""

"""
Full population plots
"""

def a_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots PDF of semi-major axis for the initial data.
    """
    plt.figure()
    plt.hist(initial_data['separation'], histtype = 'step', normed = True, 
             bins = nbins)
    plt.xlabel(r'Semi-major axis, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'PDF($a$ / AU)', fontsize = 16)

    plot_title = 'Initial_separation_check'
    plot_description = 'PDF of initial separations'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots PDF of primary mass for the initial data.
    """
    plt.figure(figsize = (8, 6))
    plt.hist(initial_data['mass1'], histtype = 'step', normed = True,
             bins = 3. * nbins, range = [0, 150])
    plt.xlabel(r'Primary ZAMS mass, $M_{1, {\rm ZAMS}}$ ($M_{\odot}$)', 
               fontsize = 16)
    plt.ylabel(r'PDF($M_{1, {\rm ZAMS}}$ / $M_\odot$)', fontsize = 16)

    plot_title = 'Primary_mass_check'
    plot_description = 'PDF of primary ZAMS mass'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def q_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots PDF of mass ratio for the initial data.
    """
    q = initial_data['mass2'] / initial_data['mass1']
    
    plt.figure()
    plt.hist(q, histtype = 'step', normed = True, bins = nbins, range = [0, 1])
    plt.xlabel(r'ZAMS Mass ratio, $q_{\rm ZAMS}$', fontsize = 16)
    plt.ylabel(r'PDF($q_{\rm ZAMS}$)', fontsize = 16)

    plot_title = 'ZAMS_mass_ratio_check'
    plot_description = 'PDF of mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def eccentricity_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                      BBH_data):
    """
    Plots PDF of eccentricity for the initial data.
    """
    plt.figure()
    plt.hist(initial_data['eccentricity'], histtype = 'step', normed = True, 
             bins = nbins, range=[0, 1])
    plt.xlabel(r'Eccentricity, $e$', fontsize = 16)
    plt.ylabel(r'PDF($e$)', fontsize = 16)

    plot_title = 'Eccentricity_check'
    plot_description = 'PDF of eccentricity'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def chirp_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots DCO chirp mass distribution for the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    chirp_mass = ((m1 * m2) ** (3.0 / 5.0)) / ((m1 + m2) ** (1.0 / 5.0))

    x_max = round(max(chirp_mass)) + 5.0

    plt.figure()
    plt.hist(chirp_mass.as_matrix(), histtype = 'step', normed = True,
             bins = nbins, range = [0, x_max])
    plt.xlim(0, x_max)
    plt.xlabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'PDF($M_{\rm c}$ / $M_{\odot}$)', fontsize = 16)
    plt.yscale('log', nonposy = 'clip')
    
    plot_title = 'Chirp_mass_plot'
    plot_description = 'Chirp mass distribution of whole population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def qDCO_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots distribution of DCO mass ratio of all DCOs for the merging data.
    """
    q = merging_data['M2'] / merging_data['M1']
    
    plt.figure(figsize = (8, 6))
    plt.hist(q.as_matrix(), histtype = 'step', normed = True, bins = nbins, 
             range = [0, 10])
    plt.xlabel(r'DCO mass ratio, $q_{\rm DCO}$', fontsize = 16)
    plt.ylabel(r'PDF($q_{\rm DCO}$)', fontsize = 16)
    
    plot_title = 'DCO_mass_ratio_plot'
    plot_description = 'Unsymmetrised DCO mass ratio of whole population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def tc_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots distribution of coallescence times for the merging data.
    """
    Myrs = 60 * 60 * 24 * 365 * 1000000
    
    t = (merging_data['tc'] + 1) / Myrs # + 1 to avoid 0 in logaritmic tc scale. This is negligible as tc is recorded in seconds in the merging data
    # Creates a set of bin edges, such that the histogram will appear to have
    # "equal" bin lengths on a logarithmic scale. Bins start at minimum value of
    # t, and finish at the maximum value of t, split into nbins bins
    logBins = np.logspace(np.floor(np.log10(np.amin(t))), 
                          np.ceil(np.log10(np.amax(t))), nbins)
    
    tc_range = [-2, 5] 
        # Default tc range, currently used for when individual DCO populations
        # aren't identified
    
    plt.figure(figsize = (8, 6))
    plt.hist(t.as_matrix(), histtype = 'step', bins = logBins, range = tc_range)
    # pl.xlim(1,1e5) # Can be used to set x axis limits, normally used to ignore DCOs with tc=0
    
    plt.xlabel(r'Coalescence time, $t_{\rm c}$ (Myr)', fontsize = 16)
    plt.ylabel(r'Number of DCOs',fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(0.1, 1e6) # Lower limit set to 0.1 so bins with a value of 1 can be identified with ease
    
    plot_title = 'tc_plot'
    plot_description = 'Coalescence time distribution for whole population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def compare_m1DCOm2DCO_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                            BBH_data):
    """
    Plots comparitive distributions of the component masses for all DCOs in the 
    merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    # Checks to set upper limit of mass axis to either 100 or 50, so as to keep
    # the bins for both distributions equal
    m1_max = max(m1)
    m2_max = max(m2)
    maxDCO_mass = max([m1_max, m2_max])
    if maxDCO_mass > 50: rangeMax = 100
    else: rangeMax = 50
    
    plt.figure(figsize = (8, 6))
    plt.hist(m1.as_matrix(), histtype = 'step', color = 'b', bins = nbins, 
             range = [0, rangeMax], label = r'$M_{1}$')
    plt.hist(m2.as_matrix(), histtype = 'step', color = 'r', bins = nbins, 
             range = [0, rangeMax], label = r'$M_{2}$')
    plt.xlabel(r'Mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of COs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'DCO_primary_mass_and_DCO_secondary_mass_distributions'
    plot_description = 'Comparative distribution of DCO component masses'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def compare_m1ZAMSm2ZAMS_DCOplot(initial_data, merging_data, BNS_data, 
                                 BHNS_data, BBH_data):
    """
    Plots comparitive distributions of the merging DCO component's ZAMS masses 
    for the merging data.
    """
    m1ZAMS = merging_data['M1ZAMS']
    m2ZAMS = merging_data['M2ZAMS']
    
    plt.figure(figsize = (8, 6))
    plt.hist(m1ZAMS.as_matrix(), histtype = 'step', color = 'b',  bins = nbins, 
             range = [0, 150], label = r'$M_{1,{\rm ZAMS}}$')
    plt.hist(m2ZAMS.as_matrix(), histtype = 'step', color = 'r',  bins = nbins, 
            range = [0, 150], label = r'$M_{2,{\rm ZAMS}}$')
    plt.xlabel(r'Mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of COs',fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'DCO_primary_ZAMS_mass_and_DCO_secondary_ZAMS_mass_distributions'
    plot_description = 'Comparative distribution of DCO component ZAMS masses'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def aInitial_DCOplot(initial_data, merging_data, BNS_data, BHNS_data, 
                     BBH_data):
    """
    Plots the distribution of the initial separation of all merging DCOs in the 
    merging data.
    """
    a = merging_data['separationInitial']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(a)))) to 1e3, with a number of bins equal
    # to nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(a))), 3, nbins)
    
    plt.figure()
    plt.hist(a.as_matrix(), histtype = 'step', bins = logBins)
    plt.xlabel(r'Initial separation, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of DCOs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin = 0.5)
    
    plot_title = 'DCO_initial_separation_distribution'
    plot_description = 'Distribution of DCO initial separations'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def aDCO_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots the distribution of the separation at DCO formation of all merging 
    DCOs in the merging data.
    """
    aDCO = merging_data['separationDCOFormation']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(aDCO)))) to 
    # 1e(np.floor(np.log10(np.amax(aDCO)))), with a number of bins equal to 
    # nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(aDCO))), 
                          np.ceil(np.log10(np.amax(aDCO))), nbins)
    
    plt.figure(figsize = (7, 6))
    plt.hist(aDCO.as_matrix(), histtype = 'step', bins = logBins)
    plt.xlabel(r'Separation at DCO formation, $a_{\rm DCO}$ (AU)', 
              fontsize = 16)
    plt.ylabel(r'Number of DCOs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin = 0.5)
    
    plot_title = 'Separation_at_DCO_formation_distribution'
    plot_description = 'Distribution of DCO formation separations'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def periPreSN_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots the distribution of the periapsis just before the second SN 
    for all DCOs in the merging data. 
    """
    a = merging_data['separationPrior2ndSN']
    e = merging_data['eccentricityPrior2ndSN']
    
    r = a * (1 - e) # Calculates periapsis
    
    plt.figure(figsize = (8, 6))
    plt.hist(r.as_matrix(), histtype = 'step', bins = nbins)
    plt.xlabel(r'Periapsis, $r_{\rm p}$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of DCOs', fontsize = 16)
    
    plot_title = 'DCO_periapsis_just_before_the_second_SN_distribution'
    plot_description = 'Distribution of DCO pre second supernova periapsies'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def periPostSN_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots the distribution of the periapsis at DCO formation for all DCOs in 
    the merging data.
    """
    a = merging_data['separationDCOFormation']
    e = merging_data['eccentricityDCOFormation']
    
    r = a * (1 - e) # Calculates periapsis
    
    plt.figure(figsize = (8, 6))
    plt.hist(r.as_matrix(), histtype = 'step', bins = nbins)
    plt.xlabel(r'Periapsis, $r_{\rm p}$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of DCOs', fontsize = 16)
    
    plot_title = 'DCO_periapsis_at_DCO_formation_distribution'
    plot_description = 'Distribution of periapsies at DCO formation'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def periZAMS_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots the distribution of the initial periapsis for all DCOs in the 
    merging data.
    """
    a = merging_data['separationInitial']
    e = merging_data['eccentricityInitial']
    
    r = a * (1 - e) # Calculates periapsis
    
    plt.figure(figsize = (8, 6))
    plt.hist(r.as_matrix(), histtype = 'step', bins = nbins)
    plt.xlabel(r'Periapsis, $r_{\rm p}$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of DCOs', fontsize = 16)
    
    plot_title = 'DCO_initial_periapsis_distribution'
    plot_description = 'Distribution of DCO initial periapsies'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1ZAMSm2ZAMS_densityPlot_initial(initial_data, merging_data, BNS_data, 
                                     BHNS_data, BBH_data):
    """
    Plots the 2D distribution of ZAMS masses for the initial data.
    """
    m1 = initial_data['mass1']
    m2 = initial_data['mass2']
    
    # Obtains masked 2D histogram data, and the appropriate x and y binnings.
    Hmasked, xedges, yedges = histogram2d_mask(m1, m2, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Primary ZAMS mass, $M_{1,{\rm ZAMS}}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'Secondary ZAMS mass, $M_{2,{\rm ZAMS}}$ ($M_\odot$)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10}\left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'Seconary_ZAMS_mass_against_primary_ZAMS_mass_for_initial_binaries_distribution'
    plot_description = '2D density distribution of all initial binaries, of m2ZAMS against m1ZAMS'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def totalZAMS_q_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                            BBH_data):
    """
    Plots the 2D density distribution of total ZAMS mass and ZAMS mass ratio 
    for all initial binaries in the initial data.
    """
    M = initial_data['mass1'] + initial_data['mass2']

    q = initial_data['mass2'] / initial_data['mass1']
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(M, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Total ZAMS mass, $M_{\rm total,ZAMS}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'ZAMS mass ratio, $q_{\rm ZAMS}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'ZAMS_mass_ratio_against_total_ZAMS_mass_for_initial_binaries_distribution'
    plot_description = '2D density distribution of all initial binaries, of total ZAMS mass against mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1ZAMS_q_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                         BBH_data):
    """
    Plots the 2D density distribution of primary ZAMS mass and ZAMS mass ratio 
    for all initial binaries in the initial data.
    """
    m1 = initial_data['mass1']
    m2 = initial_data['mass2']
    q = m2 / m1
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(m1, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Primary ZAMS mass, $M_{1,{\rm ZAMS}}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'ZAMS mass ratio, $q_{\rm ZAMS}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'ZAMS_mass_ratio_against_primary_ZAMS_mass_distribution'
    plot_description = '2D density distribution of all initial binaries, m1ZAMS against mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1ZAMSm2ZAMS_densityPlot_merging(initial_data, merging_data, BNS_data, 
                                     BHNS_data, BBH_data):
    """
    Plots the 2D distribution of ZAMS masses for the merging data.
    """
    m1 = merging_data['M1ZAMS']
    m2 = merging_data['M2ZAMS']
    
    # Obtains masked 2D histogram data, and the appropriate x and y binnings.
    Hmasked, xedges, yedges = histogram2d_mask(m1, m2, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Primary ZAMS mass, $M_{1,{\rm ZAMS}}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'Secondary ZAMS mass, $M_{2,{\rm ZAMS}}$ ($M_\odot$)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10}\left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'Secondary_ZAMS_mass_against_primary_ZAMS_mass_for_merging_binaries_distribution'
    plot_description = '2D density distribution of all merging binaries, of m2ZAMS against m1ZAMS'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1DCOm2DCO_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                           BBH_data):
    """
    Plots the 2D density distribution of component masses for all merging DCOs 
    in the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(m1, m2, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'DCO primary mass, $M_{1}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'DCO secondary mass, $M_{2}$ ($M_\odot$)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10}\left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'DCO_secondary_mass_against_DCO_primary_mass_distribution'
    plot_description = '2D density distribution of all merging binaries, of m2DCO against m1DCO'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def chirp_qDCO_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                           BBH_data):
    """
    Plots the 2D density distribution of chirp mass and DCO mass ratio for all 
    merging DCOs in the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    chirp_mass = ((m1 * m2) ** (3.0 / 5.0)) / ((m1 + m2) ** (1.0 / 5.0))

    q = m2 / m1
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(chirp_mass, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Chirp mass, $M_{\rm c}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'DCO mass ratio, $q_{\rm DCO}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'DCO_mass_ratio_against_chirp_mass_distribution'
    plot_description = '2D density distribution of all merging binaries, of chirp mass against unsymmetrised DCO mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def totalDCO_qDCO_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                              BBH_data):
    """
    Plots the 2D density distribution of total DCO mass and DCO mass ratio for 
    all merging DCOs in the merging data.
    """
    M = merging_data['M1'] + merging_data['M2']

    q = merging_data['M2'] / merging_data['M1']
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(M, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Total DCO mass, $M_{\rm total, DCO}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'DCO mass ratio, $q_{\rm DCO}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'DCO_mass_against_total_DCO_mass_ratio_distribution'
    plot_description = '2D density distribution of all merging binaries, of total DCO mass against unsymmetrised DCO mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def chirp_q_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                        BBH_data):
    """
    Plots the 2D density distribution of chirp mass and ZAMS mass ratio for all 
    merging DCOs in the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    chirp_mass = ((m1 * m2) ** (3.0 / 5.0)) / ((m1 + m2) ** (1.0 / 5.0))

    q = merging_data['M2ZAMS'] / merging_data['M1ZAMS']
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(chirp_mass, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Chirp mass, $M_{\rm c}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'ZAMS mass ratio, $q_{\rm ZAMS}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'ZAMS_mass_ratio_against_chirp_mass_for_merging_binaries_distribution'
    plot_description = '2D density distribution of all merging binaries, of chirp mass against ZAMS mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def totalDCO_q_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                           BBH_data):
    """
    Plots the 2D density distribution of total DCO mass and ZAMS mass ratio for 
    all merging DCOs in the merging data.
    """
    M = merging_data['M1'] + merging_data['M2']

    q = merging_data['M2ZAMS'] / merging_data['M1ZAMS']
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(M, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Total DCO mass, $M_{\rm total, DCO}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'ZAMS mass ratio, $q_{\rm ZAMS}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'ZAMS_mass_ratio_against_total_DCO_mass_distribution'
    plot_description = '2D density distribution of all merging binaries, of total DCO mass against ZAMS mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def totalZAMS_q_DCOdensityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                               BBH_data):
    """
    Plots the 2D density distribution of total ZAMS mass and ZAMS mass ratio 
    for all merging DCOs in the merging data.
    """
    M = merging_data['M1ZAMS'] + merging_data['M2ZAMS']

    q = merging_data['M2ZAMS'] / merging_data['M1ZAMS']
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(M, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Total ZAMS mass, $M_{\rm total, ZAMS}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'ZAMS mass ratio, $q_{\rm ZAMS}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'ZAMS_mass_ratio_against_total_ZAMS_mass_for_merging_binaries_distribution'
    plot_description = '2D density distribution of all merging binaries, of total ZAMS mass against ZAMS mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1ZAMS_qDCO_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                            BBH_data):
    """
    Plots the 2D density distribution of primary ZAMS mass and DCO mass ratio 
    for all merging DCOs in the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    
    q = m2 / m1
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(m1, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Primary ZAMS mass, $M_{1, {\rm ZAMS}}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'DCO mass ratio, $q_{\rm DCO}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'DCO_mass_ratio_against_primary_ZAMS_mass_distribution'
    plot_description = '2D density distribution of all merging binaries, of m1ZAMS against unsymmetrised DCO mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1ZAMS_q_DCOdensityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                            BBH_data):
    """
    Plots the 2D density distribution of primary ZAMS mass and ZAMS mass ratio 
    for all merging DCOs in the merging data.
    """
    m1 = merging_data['M1ZAMS']
    m2 = merging_data['M2ZAMS']
    q = m2 / m1
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(m1, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Primary ZAMS mass, $M_{1, {\rm ZAMS}}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'ZAMS mass ratio, $q_{\rm ZAMS}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'ZAMS_mass_ratio_against_primary_ZAMS_mass_for_merging_binaries_distribution'
    plot_description = '2D density distribution of all merging binaries, of m1ZAMS against ZAMS mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1DCO_qDCO_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                           BBH_data):
    """
    Plots the 2D density distribution of primary DCO mass and DCO mass ratio 
    for all merging DCOs in the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    
    q = m2 / m1
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(m1, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'DCO primary mass, $M_{1}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'DCO mass ratio, $q_{\rm DCO}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'DCO_mass_ratio_against_DCO_primary_mass_distribution'
    plot_description = '2D density distribution of all merging binaries, of m1DCO against unsymmetrised DCO mass ratio'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def etaChirp_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                         BBH_data):
    """
    Plots the 2D density distribution of the symmetric mass ratio against chirp 
    mass for all merging DCOs in the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    
    eta = m1 * m2 * (m1 + m2) ** (-2)
    chirp_mass = ((m1 * m2) ** (3.0 / 5.0)) / ((m1 + m2) ** (1.0 / 5.0))
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(chirp_mass, eta, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Chirp mass, $M_{\rm c}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'Symmetric mass ratio, $\eta$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'Symmetric_mass_ratio_against_chirp_mass_for_merging_binaries_distribution'
    plot_description = '2D density distribution of all merging binaries, of symmetric mass ratio against chirp mass'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def aDCOEccentricityDCO_densityPlot(initial_data, merging_data, BNS_data, 
                                    BHNS_data, BBH_data):
    """
    Plots the 2D density distribution of the semi-major axes and eccentricity 
    of the DCOs at DCO formation for the merging data.
    """
    aDCO = merging_data['separationDCOFormation']
    eccentricityDCO = merging_data['eccentricityDCOFormation'] + 1e-15
    
    logbins_aDCO = np.logspace(np.floor(np.log10(np.amin(aDCO))), 
                               np.ceil(np.log10(np.amax(aDCO))), nbins)
    logbins_eDCO = np.logspace(np.floor(np.log10(np.amin(eccentricityDCO))), 
                               np.ceil(np.log10(np.amax(eccentricityDCO))), 
                                      nbins)
    
    logbins = [logbins_aDCO, logbins_eDCO]
    
    Hmasked, xedges, yedges = histogram2d_mask(aDCO, eccentricityDCO, logbins)

    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Eccentricity, $e$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    
    plot_title = 'Eccentricity_at_DCO_formation_against_separation_at_DCO_formation_distribution'
    plot_description = '2D density distribution of all merging binaries, of separation at DCO formation against eccentricity at DCO formation'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def symmetric_m1DCOm2DCO_densityPlot(initial_data, merging_data, BNS_data, 
                                     BHNS_data, BBH_data):
    """
    Plots the 2D density distribution of the symmetrised component masses 
    (m1 >= m2) for all DCOs, in the merging data.
    """
    # Creating the symmetrised arrays for the DCO component masses
    m1DCO_sym, m2DCO_sym = m1DCOm2DCO_symmetrisation(merging_data)
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(m1DCO_sym, m2DCO_sym, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'DCO symmetrised primary mass, $M_{1}$ ($M_\odot$)', fontsize = 16)
    plt.ylabel(r'DCO symmetrised secondary mass, $M_{2}$ ($M_\odot$)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'DCO_symmetrised_secondary_mass_against_DCO_symmetrised_primary_mass_distribution'
    plot_description = 'Symmetrised 2D density distribution of all merging binaries, of m1DCO against m2DCO'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def m1ZAMSSemi_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                           BBH_data):
    """
    Plots 2D density distribution of the more massive initial component's mass 
    against the initial separation of the binary, taken from the merging data.
    """
    m1ZAMS = merging_data['M1ZAMS']
    a = merging_data['separationInitial']
    
    # Creates a set of bin edges which appear equally spaced on a log scale. This is intended
    # to be used for the y (a) axis
    logBins_a = np.logspace(np.floor(np.log10(np.amin(a))), 
                            np.ceil(np.log10(np.amax(a))), nbins)
    
    # Allows for the x axes bins and y axes bins to be passed and used independetly of one another
    binning = [nbins, logBins_a]
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(m1ZAMS, a, binning)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, Hmasked, cmap = 'summer')
    plt.xlabel(r'Primary ZAMS mass, $M_{1, {\rm ZAMS}}$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'DCO initial separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\mathrm{Number \, of \, DCOs}$', fontsize = 16)
    plt.yscale('log', nonposy = 'clip')
    
    plot_title = 'DCO_initial_separation_against_primary_ZAMS_mass_distribution'
    plot_description = '2D density distribution of all merging binaries, of m1ZAMS against ZAMS separation'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def totalDCOmass_chirp_densityPlot(initial_data, merging_data, BNS_data, 
                                   BHNS_data, BBH_data):
    """
    Plots 2D density distribution of total mass and chirp mass for DCOs in the 
    merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    
    total_mass = m1 + m2
    
    chirp_mass = ((m1 * m2) ** (3.0 / 5.0)) / ((m1 + m2) ** (1.0 / 5.0))
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(chirp_mass, total_mass, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, Hmasked, cmap = 'summer')
    plt.xlabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Total DCO mass, $M_{\rm total}$ ($M_{\odot}$)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\mathrm{Number \, of \, DCOs}$', 
                       fontsize = 16)
    
    plot_title = 'Total_DCO_mass_against_chirp_mass_distribution'
    plot_description = '2D density distribution of all merging binaries, of total DCO mass against chirp mass'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def aPreSN_ePreSN_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                              BBH_data):
    """
    Plots the 2D density distribution of the orbital separation and 
    eccentricity just before the second supernova, from the merging data.
    """
    a = merging_data['separationPrior2ndSN']
    e = merging_data['eccentricityPrior2ndSN']
    
    # Creates a set of bin edges which appear as equal size with a logarithmic 
    # axis. This is intended for the x axis
    xlogBins = np.logspace(-5, np.log10(np.amax(a)), nbins)
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(a, e, [xlogBins, nbins])
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Separation, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'Eccentricity, $e$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    plt.xscale('log', nonposy = 'clip') # Sets the x axis as a logarithmic scale
    
    plot_title = 'Separation_just_before_second_SN_against_eccentricity_just_before_second_SN_for_merging_binaries_distribution'
    plot_description = '2D density distribution of all merging binaries, of binary separation prior to the second SN against binary eccentricity prior to the second SN'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def VrelVkick_densityPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                          BBH_data):
    """
    Plots the 2D density distribution of the kick velocity and component 
    object's relative velocity for the merging.
    """
    Vrel = merging_data['relativeVelocity2ndSN']
    Vkick = merging_data['drawnKick2']
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(Vkick, Vrel, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure(figsize = (8, 6))
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'Kick velocity, $V_{{\rm kick}, 2}$', fontsize = 16)
    plt.ylabel(r'Relative velocity, $V_{\rm rel}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    
    plot_title = 'Relative_velocity_against_kick_velocity_for_merging_binaries_distribution'
    plot_description = "2D density distribution of all merging binaries, of the binary component's relative velocity against the second SN kick velocity"
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def chirp_separationZAMS_densityPlot(initial_data, merging_data, BNS_data, 
                                     BHNS_data, BBH_data):
    """
    Plots the 2D density distribution of chirp mass against initial separation 
    for the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    a = merging_data['separationInitial']
    
    chirp_mass = ((m1 * m2) ** 0.6) / ((m1 + m2) ** 0.2)
    
    # Creates a set of bin edges which appear as equal size with a logarithmic axis. This
    # is intended for the x (a) axis
    logBins_a = np.logspace(np.floor(np.log10(np.amin(a))), 
                            np.ceil(np.log10(np.amax(a))), nbins)
    
    # Allows for the x axes bins and y axes bins to be passed and used independetly of one another
    binning = [logBins_a, nbins]
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(a, chirp_mass, binning)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'DCO separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    plt.xscale('log', nonposy = 'clip') # Log scales the x axis
    
    plot_title = 'Chirp_mass_against_DCO_separation_distribution'
    plot_description = '2D density distribution of all merging binaries, of chirp mass against initial separation'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def chirp_separationDCO_densityPlot(initial_data, merging_data, BNS_data, 
                                    BHNS_data, BBH_data):
    """
    Plots the 2D density distribution of chirp mass against DCO separation for 
    the merging data.
    """
    m1 = merging_data['M1']
    m2 = merging_data['M2']
    a = merging_data['separationDCOFormation']
    
    chirp_mass = ((m1 * m2) ** 0.6) / ((m1 + m2) ** 0.2)
    
    # Creates a set of bin edges which appear as equal size with a logarithmic axis. This
    # is intended for the x (a) axis
    logBins_a = np.logspace(np.floor(np.log10(np.amin(a))), 
                            np.ceil(np.log10(np.amax(a))), nbins)
    
    # Allows for the x axes bins and y axes bins to be passed and used independetly of one another
    binning = [logBins_a, nbins]
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(a, chirp_mass, binning)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, np.log10(Hmasked), cmap = 'summer')
    plt.xlabel(r'DCO separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\log_{10} \left( \mathrm{Number \, of \, DCOs} \right)$', 
                       fontsize = 16)
    plt.xscale('log', nonposy = 'clip') # Log scales the x axis
    
    plot_title = 'Chirp_mass_against_DCO_separation_distribution'
    plot_description = '2D density distribution of all merging binaries, of chirp mass against DCO separation'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def symmetrised_qDCO_chirp_densityPlot(initial_data, merging_data, BNS_data, 
                                       BHNS_data, BBH_data):
    """
    Plots 2D denisty distribution of the mass ratio and chirp mass for DCOs in 
    the merging data.
    """
    # Creating the symmetrised arrays for the DCO component masses
    m1DCO_sym, m2DCO_sym = m1DCOm2DCO_symmetrisation(merging_data)
    
    q = m2DCO_sym / m1DCO_sym
    
    chirp_mass = (((m1DCO_sym * m2DCO_sym) ** (3.0 / 5.0)) / 
                  ((m1DCO_sym + m2DCO_sym) ** (1.0 / 5.0)))
    
    # Obtains the masked 2D histogram data, and the x and y bin edges
    Hmasked, xedges, yedges = histogram2d_mask(chirp_mass, q, nbins)
    
    # Plots 2D density distribution using the summer colour map
    plt.figure()
    plt.pcolormesh(xedges, yedges, Hmasked, cmap = 'summer')
    plt.xlabel(r'Chrip mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'DCO mass ratio, $q_{\rm DCO}$', fontsize = 16)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'$\mathrm{Number \, of \, DCOs}$', 
                       fontsize = 16)
    
    plot_title = 'DCO_mass_ratio_against_chirp_mass_symmetrised_distribution'
    plot_description = 'Symmetrised 2D density distribution of all merging binaries, of DCO mass ratio against chirp mass'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BNSPops_m1ZAMSm2ZAMS_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                              BBH_data):
    """
    Plots comparitive distributions of the merging BNS population's 
    component ZAMS masses.
    """
    plt.figure(figsize = (8, 6))
    plt.hist(BNS_data['M1ZAMS'].as_matrix(), histtype = 'step', color = 'b',  bins = nbins, 
             range = [0, 150], label = r'$M_{1, {\rm ZAMS}}$')
    plt.hist(BNS_data['M2ZAMS'].as_matrix(), histtype = 'step', color = 'r',  bins = nbins, 
             range=[0, 150], label = r'$M_{2, {\rm ZAMS}}$')

    plt.xlabel(r'ZAMS mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of NSs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'BNS_ZAMS_mass_distributions'
    plot_description = 'Comparison of m1ZAMS and m2ZAMS distributions for the BNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BHNSPops_m1ZAMSm2ZAMS_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                               BBH_data):
    """
    Plots comparitive distributions of the merging BHNS population's 
    component ZAMS masses.
    """
    plt.figure(figsize = (8, 6))
    plt.hist(BHNS_data['M1ZAMS'].as_matrix(), histtype = 'step', color = 'b',  bins = nbins, 
             range = [0, 150], label = r'$M_{1, {\rm ZAMS}}$')
    plt.hist(BHNS_data['M2ZAMS'].as_matrix(), histtype = 'step', color = 'r',  bins = nbins, 
             range=[0, 150], label = r'$M_{2, {\rm ZAMS}}$')

    plt.xlabel(r'ZAMS mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of NSs or BHs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'BHNS_ZAMS_mass_distributions'
    plot_description = 'Comparison of m1ZAMS and m2ZAMS distributions for the BHNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BBHPops_m1ZAMSm2ZAMS_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                              BBH_data):
    """
    Plots comparitive distributions of the merging BBH population's 
    component ZAMS masses.
    """
    plt.figure(figsize = (8, 6))
    plt.hist(BBH_data['M1ZAMS'].as_matrix(), histtype = 'step', color = 'b',  bins = nbins, 
             range = [0, 150], label = r'$M_{1, {\rm ZAMS}}$')
    plt.hist(BBH_data['M2ZAMS'].as_matrix(), histtype = 'step', color = 'r',  bins = nbins, 
             range=[0, 150], label = r'$M_{2, {\rm ZAMS}}$')

    plt.xlabel(r'ZAMS mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of BHs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'BBH_ZAMS_mass_distributions'
    plot_description = 'Comparison of m1ZAMS and m2ZAMS distributions for the BBH population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BNSPops_m1DCOm2DCO_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                            BBH_data):
    """
    Plots comparitive distributions of the merging BNS population's component 
    DCO masses.
    """
    plt.figure()
    plt.hist(BNS_data['M1'].as_matrix(), histtype = 'step', color = 'b',  bins = nbins, 
            range = [0, 100], label = r'$M_1$')
    plt.hist(BNS_data['M2'].as_matrix(), histtype = 'step', color = 'r',  bins = nbins, 
            range = [0, 100], label = r'$M_2$')
    
    plt.xlabel(r'Mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of NSs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'BNS_mass_distributions'
    plot_description = 'Comparison of m1DCO and m2DCO distributions for the BNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BHNSPops_m1DCOm2DCO_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                             BBH_data):
    """
    Plots comparitive distributions of the merging BHNS population's component 
    DCO masses.
    """
    plt.figure()
    plt.hist(BHNS_data['M1'].as_matrix(), histtype = 'step', color = 'b',  bins = nbins, 
            range = [0, 100], label = r'$M_1$')
    plt.hist(BHNS_data['M2'].as_matrix(), histtype = 'step', color = 'r',  bins = nbins, 
            range = [0, 100], label = r'$M_2$')
    
    plt.xlabel(r'Mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of NSs or BHs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'BHNS_mass_distributions'
    plot_description = 'Comparison of m1DCO and m2DCO distributions for the BHNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BBHPops_m1DCOm2DCO_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                            BBH_data):
    """
    Plots comparitive distributions of the merging BBH population's component 
    DCO masses.
    """
    plt.figure()
    plt.hist(BBH_data['M1'].as_matrix(), histtype = 'step', color = 'b',  bins = nbins, 
            range = [0, 100], label = r'$M_1$')
    plt.hist(BBH_data['M2'].as_matrix(), histtype = 'step', color = 'r',  bins = nbins, 
            range = [0, 100], label = r'$M_2$')
    
    plt.xlabel(r'Mass of component, $M$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of BHs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'BBH_mass_distributions'
    plot_description = 'Comparison of m1DCO and m2DCO distributions for the BBH population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BNSPops_aInitial_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                          BBH_data):
    """
    Plots the distribution of the initial separation for the merging BNS 
    population.
    """
    a = BNS_data['separationInitial']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(a)))) to 1e3, with a number of bins equal
    # to nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(a))), 3, nbins)
    
    plt.figure(figsize = (7, 5))
    plt.hist(a.as_matrix(), histtype = 'step', bins = logBins, color = 'b')
    plt.xlabel(r'Initial separation, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of BNSs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin = 0.5)
    
    plot_title = 'BNS_initial_separation_distribution'
    plot_description = 'Initial separation distribution for the BNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BHNSPops_aInitial_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                           BBH_data):
    """
    Plots the distribution of the initial separation for the merging BHNS 
    population.
    """
    a = BHNS_data['separationInitial']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(a)))) to 1e3, with a number of bins equal
    # to nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(a))), 3, nbins)
    
    plt.figure(figsize = (7, 5))
    plt.hist(a.as_matrix(), histtype = 'step', bins = logBins, color = 'b')
    plt.xlabel(r'Initial separation, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of BHNSs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin = 0.5)
    
    plot_title = 'BHNS_initial_separation_distribution'
    plot_description = 'Initial separation distribution for the BHNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BBHPops_aInitial_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                          BBH_data):
    """
    Plots the distribution of the initial separation for the merging BBH 
    population.
    """
    a = BBH_data['separationInitial']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(a)))) to 1e3, with a number of bins equal
    # to nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(a))), 3, nbins)
    
    plt.figure(figsize = (7, 5))
    plt.hist(a.as_matrix(), histtype = 'step', bins = logBins, color = 'b')
    plt.xlabel(r'Initial separation, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of BBHs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin = 0.5)
    
    plot_title = 'BBH_initial_separation_distribution'
    plot_description = 'Initial separation distribution for the BBH population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BNSPops_aDCO_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                      BBH_data):
    """
    Plots the distribution of the separation at DCO formation for the merging 
    BNS population.
    """
    aDCO = BNS_data['separationDCOFormation']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(a)))) to 1e3, with a number of bins equal
    # to nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(aDCO))), 
                          np.ceil(np.log10(np.amax(aDCO))), nbins)  
    
    plt.figure(figsize = (8, 5))
    plt.hist(aDCO.as_matrix(), histtype = 'step', bins = logBins, color = 'b')
    plt.xlabel(r'Separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of BNSs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin=0.5)
    
    plot_title = 'BNS_separation_at_DCO_formation_distribution'
    plot_description = 'DCO separation distribution for the BNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BHNSPops_aDCO_plot(initial_data, merging_data, BNS_data, BHNS_data, 
                       BBH_data):
    """
    Plots the distribution of the separation at DCO formation for the merging 
    BHNS population.
    """
    aDCO = BHNS_data['separationDCOFormation']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(a)))) to 1e3, with a number of bins equal
    # to nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(aDCO))), 
                          np.ceil(np.log10(np.amax(aDCO))), nbins)  
    
    plt.figure(figsize = (7, 5))
    plt.hist(aDCO.as_matrix(), histtype = 'step', bins = logBins, color = 'b')
    plt.xlabel(r'Separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of BHNSs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin=0.5)
    
    plot_title = 'BHNS_separation_at_DCO_formation_distribution'
    plot_description = 'DCO separation distribution for the BHNS population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def BBHPops_aDCO_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots the distribution of the separation at DCO formation for the merging 
    BBH population.
    """
    aDCO = BBH_data['separationDCOFormation']
    
    # Creates a set of bin edges that look equal sized on a log scale, from
    # 1e(np.floor(np.log10(np.amin(a)))) to 1e3, with a number of bins equal
    # to nbins
    logBins = np.logspace(np.floor(np.log10(np.amin(aDCO))), 
                          np.ceil(np.log10(np.amax(aDCO))), nbins)  
    
    plt.figure(figsize = (7, 5))
    plt.hist(aDCO.as_matrix(), histtype = 'step', bins = logBins, color = 'b')
    plt.xlabel(r'Separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Number of BBHs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(ymin=0.5)
    
    plot_title = 'BBH_separation_at_DCO_formation_distribution'
    plot_description = 'DCO separation distribution for the BBH population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def DCO_tc_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots distribution of coalescence times for each DCO population.
    """
    Myrs = 60 * 60 * 24 * 365 * 1000000
    
    t_BNS = (BNS_data['tc'] + 1) / Myrs # + 1 to avoid 0 in logaritmic tc scale. This is negligible as tc is recorded in seconds in the merging data
    t_BHNS = (BHNS_data['tc'] + 1) / Myrs
    t_BBH = (BBH_data['tc'] + 1) / Myrs
    # Creates a set of bin edges, such that the histogram will appear to have
    # "equal" bin lengths on a logarithmic scale. Bins start at minimum value of
    # t, and finish at the maximum value of t, split into nbins bins
    logBins_BNS = np.logspace(np.floor(np.log10(np.amin(t_BNS))), 
                          np.ceil(np.log10(np.amax(t_BNS))), nbins)
    logBins_BHNS = np.logspace(np.floor(np.log10(np.amin(t_BHNS))), 
                          np.ceil(np.log10(np.amax(t_BHNS))), nbins)
    logBins_BBH = np.logspace(np.floor(np.log10(np.amin(t_BBH))), 
                          np.ceil(np.log10(np.amax(t_BBH))), nbins)
    
    tc_range = [-20, 5]
    
    plt.figure(figsize = (8, 6))
    plt.hist(t_BNS.as_matrix(), histtype = 'step', bins = logBins_BNS, range = tc_range, 
             color = 'r', label = r'BNS')
    plt.hist(t_BHNS.as_matrix(), histtype = 'step', bins = logBins_BHNS, range = tc_range, 
             color = 'y', label = r'BHNS')
    plt.hist(t_BBH.as_matrix(), histtype = 'step', bins = logBins_BBH, range = tc_range, 
             color = 'b', label = r'BBH')
    # plt.xlim(1,1e5) # Can be used to set x axis limits, normally used to ignore DCOs with tc = 0
    plt.legend(loc = 'best', shadow = False)
    
    plt.xlabel(r'Coalescence time, $t_{\rm c}$ (Myr)', fontsize = 16)
    plt.ylabel(r'Number of DCOs', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    plt.ylim(0.1, 1e6) # Lower limit set to 0.1 so bins with a value of 1 can be identified with ease
    
    plot_title = 'Coalesence_time_distribution_for_each_DCO_population'
    plot_description = 'Comparison of the tc distributions for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def DCO_chirp_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots distribution of chirp mass for each DCO population.
    """
    m1_BNS = BNS_data['M1']
    m2_BNS = BNS_data['M2']
    chirp_mass_BNS = (((m1_BNS * m2_BNS) ** (3.0 / 5.0)) / 
                      ((m1_BNS + m2_BNS) ** (1.0 / 5.0)))
    
    x_max_BNS = round(max(chirp_mass_BNS)) + 5.0
    
    binning_BNS = np.linspace(0, x_max_BNS, x_max_BNS + 1)
    
    m1_BHNS = BHNS_data['M1']
    m2_BHNS = BHNS_data['M2']
    chirp_mass_BHNS = (((m1_BHNS * m2_BHNS) ** (3.0 / 5.0)) / 
                       ((m1_BHNS + m2_BHNS) ** (1.0 / 5.0)))
    
    x_max_BHNS = round(max(chirp_mass_BHNS)) + 5.0
    
    binning_BHNS = np.linspace(0, x_max_BHNS, x_max_BHNS + 1)
    
    m1_BBH = BBH_data['M1']
    m2_BBH = BBH_data['M2']
    chirp_mass_BBH = (((m1_BBH * m2_BBH) ** (3.0 / 5.0)) / 
                      ((m1_BBH + m2_BBH) ** (1.0 / 5.0)))
    
    x_max_BBH = round(max(chirp_mass_BBH)) + 5.0
    
    binning_BBH = np.linspace(0, x_max_BBH, x_max_BBH + 1)
    
    
    plt.figure()
    plt.hist(chirp_mass_BNS.as_matrix(), histtype = 'step', bins = binning_BNS, 
             range = [0, x_max_BNS], color = 'r', label = r'BNS')
    plt.hist(chirp_mass_BHNS.as_matrix(), histtype = 'step', bins = binning_BHNS, 
             range = [0, x_max_BHNS], color = 'y', label = r'BHNS')
    plt.hist(chirp_mass_BBH.as_matrix(), histtype = 'step', bins = binning_BBH, 
             range = [0, x_max_BBH], color = 'b', label = r'BBH')
    plt.xlim(0, x_max_BBH)
    plt.xlabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Number of DCOs',fontsize = 16)
    plt.yscale('log', nonposy = 'clip')
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'Chirp_mass_distribution_for_each_DCO_population'
    plot_description = 'Comparison of the chirp mass distributions for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def DCO_qDCO_plot(initial_data, merging_data, BNS_data, BHNS_data, BBH_data):
    """
    Plots distribution of DCO mass ratio for each DCO population.
    """
    q_BNS = BNS_data['M2'] / BNS_data['M1']
    
    q_BHNS = BHNS_data['M2'] / BHNS_data['M1']
    
    q_BBH = BBH_data['M2'] / BBH_data['M1']
    
    plt.figure(figsize = (8, 6))
    plt.hist(q_BNS.as_matrix(), histtype = 'step', bins = nbins, range = [0, 10], 
             color = 'r', label = r'BNS')
    plt.hist(q_BHNS.as_matrix(), histtype = 'step', bins = nbins, range = [0, 10], 
             color = 'y', label = r'BHNS')
    plt.hist(q_BBH.as_matrix(), histtype = 'step', bins = nbins, range = [0, 10], 
             color = 'b', label = r'BBH')
    plt.xlabel(r'DCO mass ratio, $q_{\rm DCO}$', fontsize = 16)
    plt.ylabel(r'Number of DCOs', fontsize = 16)
    plt.legend(loc = 'best', shadow = False)
    
    plot_title = 'DCO_mass_ratio_distribution_for_each_DCO_population'
    plot_description = 'Comparison of the unsymmetrised DCO mass ratios for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def etaChirp_systemPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                        BBH_data):
    """
    Plots the 2D distribution of the symmetric mass ratio against chirp mass 
    for all merging DCOs.
    """
    m1_BNS = BNS_data['M1']
    m2_BNS = BNS_data['M2']
    
    eta_BNS = m1_BNS * m2_BNS * (m1_BNS + m2_BNS) ** (-2)
    chirp_mass_BNS = (((m1_BNS * m2_BNS) ** (3.0 / 5.0)) / 
                      ((m1_BNS + m2_BNS) ** (1.0 / 5.0)))
    
    m1_BHNS = BHNS_data['M1']
    m2_BHNS = BHNS_data['M2']
    
    eta_BHNS = m1_BHNS * m2_BHNS * (m1_BHNS + m2_BHNS) ** (-2)
    chirp_mass_BHNS = (((m1_BHNS * m2_BHNS) ** (3.0 / 5.0)) / 
                       ((m1_BHNS + m2_BHNS) ** (1.0 / 5.0)))
    
    m1_BBH = BBH_data['M1']
    m2_BBH = BBH_data['M2']
    
    eta_BBH = m1_BBH * m2_BBH * (m1_BBH + m2_BBH) ** (-2)
    chirp_mass_BBH = (((m1_BBH * m2_BBH) ** (3.0 / 5.0)) / 
                      ((m1_BBH + m2_BBH) ** (1.0 / 5.0)))
    
    plt.figure()
    plt.scatter(chirp_mass_BNS, eta_BNS, color = 'r', s = 5, alpha = 0.5, 
                label = r'BNS')
    plt.scatter(chirp_mass_BHNS, eta_BHNS, color = 'y', s = 5, alpha = 0.5, 
                label = r'BHNS')
    plt.scatter(chirp_mass_BBH, eta_BBH, color = 'b', s = 5, alpha = 0.5, 
                label = r'BBH')
    plt.legend(loc = 'best', shadow = False)
    
    plt.xlabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    plt.ylabel(r'Symmetric mass ratio, $\eta$', fontsize = 16)
    
    plot_title = 'Symmetric_mass_ratio_against_chirp_mass_distribution_for_each_DCO_population'
    plot_description = '2D distribution of the symmetric mass ratio against chirp mass for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def aDCOEccentricityDCO_systemPlot(initial_data, merging_data, BNS_data, 
                                   BHNS_data, BBH_data):
    """
    Plots the distribution of semi-major axes and eccentricity for each DCO
    population.
    """
    aDCO_BNS = BNS_data['separationDCOFormation']
    eccentricityDCO_BNS = BNS_data['eccentricityDCOFormation'] + 1e-15
    
    aDCO_BHNS = BHNS_data['separationDCOFormation']
    eccentricityDCO_BHNS = BHNS_data['eccentricityDCOFormation'] + 1e-15
    
    aDCO_BBH = BBH_data['separationDCOFormation']
    eccentricityDCO_BBH = BBH_data['eccentricityDCOFormation'] + 1e-15
    
    plt.figure()
    plt.scatter(aDCO_BNS, eccentricityDCO_BNS, color = 'r', s = 5, alpha = 0.5, 
               label = r'BNS')
    plt.scatter(aDCO_BHNS, eccentricityDCO_BHNS, color = 'y', s = 5, 
                alpha = 0.5, label = r'BHNS')
    plt.scatter(aDCO_BBH, eccentricityDCO_BBH, color = 'b', s = 5, alpha = 0.5, 
               label = r'BBH')
    plt.legend(loc = 'best', shadow = False)
    
    plt.xlabel(r'Separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Eccentricity, $e$', fontsize = 16)
    plt.xscale('log', nonposy = 'clip')
    plt.yscale('log', nonposy = 'clip')
    
    plot_title = 'Eccentricity_at_DCO_formation_against_separation_at_DCO_formation_distribution_for_each_DCO_population'
    plot_description = '2D distribution of DCO separation against DCO eccentricity for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def aPreSN_ePreSN_systemPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                             BBH_data):
    """
    Plots the 2D distribution of the orbital separation and eccentricity, just
    before the second supernova, for each DCO population.
    """
    a_BNS = BNS_data['separationPrior2ndSN']
    e_BNS = BNS_data['eccentricityPrior2ndSN']
    
    a_BHNS = BHNS_data['separationPrior2ndSN']
    e_BHNS = BHNS_data['eccentricityPrior2ndSN']
    
    a_BBH = BBH_data['separationPrior2ndSN']
    e_BBH = BBH_data['eccentricityPrior2ndSN']
    
    plt.figure()
    plt.scatter(a_BNS, e_BNS, color = 'r', s = 5, alpha = 0.5, label = r'BNS')
    plt.scatter(a_BHNS, e_BHNS, color = 'b', s = 5, alpha = 0.5, 
                label = r'BHNS')
    plt.scatter(a_BBH, e_BBH, color = 'b', s = 5, alpha = 0.5, label = r'BBH')
    plt.legend(loc = 'best', shadow = False)
      
    plt.xlabel(r'Separation, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'Eccentricity, $e$', fontsize = 16)
    plt.xscale('log', nonposy = 'clip') # Log scale on the x-axis
    
    plot_title = 'Eccentricity_just_before_second_SN_against_separation_just_before_second_SN_distribution_for_each_DCO_population'
    plot_description = '2D distribution of binary separation prior to the second SN against binary eccentricity prior to the second SN, for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def VrelVkick_systemPlot(initial_data, merging_data, BNS_data, BHNS_data, 
                         BBH_data):
    """
    Plots the distribution of the relative component's velocity and the 
    supernova kick velocity for each DCO population.
    """
    Vrel_BNS = BNS_data['relativeVelocity2ndSN']
    Vkick_BNS = BNS_data['drawnKick2']
    
    Vrel_BHNS = BHNS_data['relativeVelocity2ndSN']
    Vkick_BHNS = BHNS_data['drawnKick2']
    
    Vrel_BBH = BBH_data['relativeVelocity2ndSN']
    Vkick_BBH = BBH_data['drawnKick2']
    
    plt.figure()
    plt.scatter(Vkick_BNS, Vrel_BNS, color = 'r', s = 5, alpha = 0.5, 
                label = r'BNS')
    plt.scatter(Vkick_BHNS, Vrel_BHNS, color = 'y', s = 5, alpha = 0.5, 
                label = r'BHNS')
    plt.scatter(Vkick_BBH, Vrel_BBH, color = 'b', s = 5, alpha = 0.5, 
                label = r'BBH')
    plt.legend(loc = 'best', shadow = False)
      
    plt.xlabel(r'Kick velocity, $V_{{\rm kick}, 2}$', fontsize = 16)
    plt.ylabel(r'Relative velocity, $V_{\rm rel}$', fontsize = 16)
    
    plot_title = 'Relative_velocity_against_kick_velocity_distribution_for_each_DCO_population'
    plot_description = '2D distribution of the binary components relative velocity against the second SN kick velocity, for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def chirpSeparationZAMS_systemPlot(initial_data, merging_data, BNS_data, 
                                   BHNS_data, BBH_data):
    """
    Plots the distribution of chirp mass against initial separation, for each 
    DCO population.
    """
    m1_BNS = BNS_data['M1']
    m2_BNS = BNS_data['M2']
    a_BNS = BNS_data['separationInitial']
    
    chirp_mass_BNS = ((m1_BNS * m2_BNS) ** 0.6) / ((m1_BNS + m2_BNS) ** 0.2)
    
    m1_BHNS = BHNS_data['M1']
    m2_BHNS = BHNS_data['M2']
    a_BHNS = BHNS_data['separationInitial']
    
    chirp_mass_BHNS = (((m1_BHNS * m2_BHNS) ** 0.6) / 
                       ((m1_BHNS + m2_BHNS) ** 0.2))
    
    m1_BBH = BBH_data['M1']
    m2_BBH = BBH_data['M2']
    a_BBH = BBH_data['separationInitial']
    
    chirp_mass_BBH = ((m1_BBH * m2_BBH) ** 0.6) / ((m1_BBH + m2_BBH) ** 0.2)
    
    plt.figure()
    plt.scatter(a_BNS, chirp_mass_BNS, color = 'r', s = 5, alpha = 0.5, 
                label = r'BNS')
    plt.scatter(a_BHNS, chirp_mass_BHNS, color = 'y', s = 5, alpha = 0.5, 
                label = r'BHNS')
    plt.scatter(a_BBH, chirp_mass_BBH, color = 'b', s = 5, alpha = 0.5, 
                label = r'BBH')
    plt.legend(loc = 'best', shadow = False)
      
    plt.xlabel(r'Initial separation, $a$ (AU)', fontsize = 16)
    plt.ylabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    plt.xscale('log', nonposy = 'clip') # Log scales the x-axis
    
    plot_title = 'Chirp_mass_against_initial_separation_distribution_for_each_DCO_population'
    plot_description = '2D distribution of chirp mass against initial separation for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()

def chirpSeparationDCO_systemPlot(initial_data, merging_data, BNS_data, 
                                  BHNS_data, BBH_data):
    """
    Plots the distribution of chirp mass against DCO separation, for each DCO
    population.
    """
    m1_BNS = BNS_data['M1']
    m2_BNS = BNS_data['M2']
    a_BNS = BNS_data['separationDCOFormation']
    
    chirp_mass_BNS = ((m1_BNS * m2_BNS) ** 0.6) / ((m1_BNS + m2_BNS) ** 0.2)
    
    m1_BHNS = BHNS_data['M1']
    m2_BHNS = BHNS_data['M2']
    a_BHNS = BHNS_data['separationDCOFormation']
    
    chirp_mass_BHNS = (((m1_BHNS * m2_BHNS) ** 0.6) / 
                       ((m1_BHNS + m2_BHNS) ** 0.2))
    
    m1_BBH = BBH_data['M1']
    m2_BBH = BBH_data['M2']
    a_BBH = BBH_data['separationDCOFormation']
    
    chirp_mass_BBH = ((m1_BBH * m2_BBH) ** 0.6) / ((m1_BBH + m2_BBH) ** 0.2)
    
    plt.figure()
    plt.scatter(a_BNS, chirp_mass_BNS, color = 'r', s = 5, alpha = 0.5, 
                label = r'BNS')
    plt.scatter(a_BHNS, chirp_mass_BHNS, color = 'y', s = 5, alpha = 0.5, 
                label = r'BHNS')
    plt.scatter(a_BBH, chirp_mass_BBH, color = 'b', s = 5, alpha = 0.5, 
                label = r'BBH')
    plt.legend(loc = 'best', shadow = False)
      
    plt.xlabel(r'Separation, $a_{\rm DCO}$ (AU)', fontsize = 16)
    plt.ylabel(r'Chirp mass, $M_{\rm c}$ ($M_{\odot}$)', fontsize = 16)
    plt.xscale('log', nonposy = 'clip') # Log scales the x-axis
    
    plot_title = 'Chirp_mass_against_separation_at_DCO_formation_distribution_for_each_DCO_population'
    plot_description = '2D distribution of chirp mass against DCO separation for each DCO population'
    save_and_add_to_html(plot_title, plot_description)
    
    plt.close()


"""
Choosing plots
"""

plots = {'a_plot' : a_plot, 'm1_plot' :  m1_plot, 'q_plot' : q_plot, 
         'eccentricity_plot' : eccentricity_plot, 'chirp_plot' : chirp_plot, 
         'qDCO_plot' : qDCO_plot, 'tc_plot' : tc_plot, 
         'compare_m1DCOm2DCO_plot' : compare_m1DCOm2DCO_plot,
         'compare_m1ZAMSm2ZAMS_DCOplot' : compare_m1ZAMSm2ZAMS_DCOplot,
         'aInitial_DCOplot' : aInitial_DCOplot, 'aDCO_plot' : aDCO_plot,
         'periPreSN_plot' : periPreSN_plot, 
         'periPostSN_plot' : periPostSN_plot, 'periZAMS_plot' : periZAMS_plot,
         'm1ZAMSm2ZAMS_densityPlot_initial' : m1ZAMSm2ZAMS_densityPlot_initial,
         'totalZAMS_q_densityPlot' : totalZAMS_q_densityPlot,
         'm1ZAMS_q_densityPlot' : m1ZAMS_q_densityPlot,
         'm1ZAMSm2ZAMS_densityPlot_merging' : m1ZAMSm2ZAMS_densityPlot_merging,
         'm1DCOm2DCO_densityPlot' : m1DCOm2DCO_densityPlot,
         'chirp_qDCO_densityPlot' : chirp_qDCO_densityPlot,
         'totalDCO_qDCO_densityPlot' : totalDCO_qDCO_densityPlot,
         'chirp_q_densityPlot' : chirp_q_densityPlot,
         'totalDCO_q_densityPlot' : totalDCO_q_densityPlot,
         'totalZAMS_q_DCOdensityPlot' : totalZAMS_q_DCOdensityPlot,
         'm1ZAMS_qDCO_densityPlot' : m1ZAMS_qDCO_densityPlot,
         'm1ZAMS_q_DCOdensityPlot' : m1ZAMS_q_DCOdensityPlot,
         'm1DCO_qDCO_densityPlot' : m1DCO_qDCO_densityPlot,
         'etaChirp_densityPlot' : etaChirp_densityPlot,
         'aDCOEccentricityDCO_densityPlot' : aDCOEccentricityDCO_densityPlot,
         'symmetric_m1DCOm2DCO_densityPlot' : symmetric_m1DCOm2DCO_densityPlot,
         'm1ZAMSSemi_densityPlot' : m1ZAMSSemi_densityPlot,
         'totalDCOmass_chirp_densityPlot' : totalDCOmass_chirp_densityPlot,
         'aPreSN_ePreSN_densityPlot' : aPreSN_ePreSN_densityPlot,
         'VrelVkick_densityPlot' : VrelVkick_densityPlot,
         'chirp_separationZAMS_densityPlot' : chirp_separationZAMS_densityPlot,
         'chirp_separationDCO_densityPlot' : chirp_separationDCO_densityPlot,
         'symmetrised_qDCO_chirp_densityPlot' : symmetrised_qDCO_chirp_densityPlot,
         'BNSPops_m1ZAMSm2ZAMS_plot' : BNSPops_m1ZAMSm2ZAMS_plot,
         'BHNSPops_m1ZAMSm2ZAMS_plot' : BHNSPops_m1ZAMSm2ZAMS_plot,
         'BBHPops_m1ZAMSm2ZAMS_plot' : BBHPops_m1ZAMSm2ZAMS_plot,
         'BNSPops_m1DCOm2DCO_plot' : BNSPops_m1DCOm2DCO_plot,
         'BHNSPops_m1DCOm2DCO_plot' : BHNSPops_m1DCOm2DCO_plot,
         'BBHPops_m1DCOm2DCO_plot' : BBHPops_m1DCOm2DCO_plot,
         'BNSPops_aInitial_plot' : BNSPops_aInitial_plot,
         'BHNSPops_aInitial_plot' : BHNSPops_aInitial_plot,
         'BBHPops_aInitial_plot' : BBHPops_aInitial_plot,
         'BNSPops_aDCO_plot' : BNSPops_aDCO_plot,
         'BHNSPops_aDCO_plot' : BHNSPops_aDCO_plot,
         'BBHPops_aDCO_plot' : BBHPops_aDCO_plot,
         'DCO_tc_plot' : DCO_tc_plot,
         'DCO_chirp_plot' : DCO_chirp_plot,
         'DCO_qDCO_plot' : DCO_qDCO_plot,
         'etaChirp_systemPlot' : etaChirp_systemPlot,
         'aDCOEccentricityDCO_systemPlot' : aDCOEccentricityDCO_systemPlot,
         'aPreSN_ePreSN_systemPlot' : aPreSN_ePreSN_systemPlot,
         'VrelVkick_systemPlot' : VrelVkick_systemPlot,
         'chirpSeparationZAMS_systemPlot' : chirpSeparationZAMS_systemPlot,
         'chirpSeparationDCO_systemPlot' : chirpSeparationDCO_systemPlot}

plots_names = ['a_plot', 
               'm1_plot', 
               'q_plot',
               'eccentricity_plot',
               'chirp_plot',
               'qDCO_plot',
               'tc_plot',
               'compare_m1DCOm2DCO_plot',
               'compare_m1ZAMSm2ZAMS_DCOplot',
               'aInitial_DCOplot',
               'aDCO_plot',
               'periPreSN_plot',
               'periPostSN_plot',
               'periZAMS_plot',
               'm1ZAMSm2ZAMS_densityPlot_initial',
               'totalZAMS_q_densityPlot',
               'm1ZAMS_q_densityPlot',
               'm1ZAMSm2ZAMS_densityPlot_merging',
               'm1DCOm2DCO_densityPlot',
               'chirp_qDCO_densityPlot',
               'totalDCO_qDCO_densityPlot',
               'chirp_q_densityPlot',
               'totalDCO_q_densityPlot',
               'totalZAMS_q_DCOdensityPlot',
               'm1ZAMS_qDCO_densityPlot',
               'm1ZAMS_q_DCOdensityPlot',
               'm1DCO_qDCO_densityPlot',
               'etaChirp_densityPlot',
               'aDCOEccentricityDCO_densityPlot',
               'symmetric_m1DCOm2DCO_densityPlot',
               'm1ZAMSSemi_densityPlot',
               'totalDCOmass_chirp_densityPlot',
               'aPreSN_ePreSN_densityPlot',
               'VrelVkick_densityPlot',
               'chirp_separationZAMS_densityPlot',
               'chirp_separationDCO_densityPlot',
               'symmetrised_qDCO_chirp_densityPlot',
               'BNSPops_m1ZAMSm2ZAMS_plot',
               'BHNSPops_m1ZAMSm2ZAMS_plot',
               'BBHPops_m1ZAMSm2ZAMS_plot',
               'BNSPops_m1DCOm2DCO_plot',
               'BHNSPops_m1DCOm2DCO_plot',
               'BBHPops_m1DCOm2DCO_plot',
               'BNSPops_aInitial_plot',
               'BHNSPops_aInitial_plot',
               'BBHPops_aInitial_plot',
               'BNSPops_aDCO_plot',
               'BHNSPops_aDCO_plot',
               'BBHPops_aDCO_plot',
               'DCO_tc_plot',
               'DCO_chirp_plot',
               'DCO_qDCO_plot',
               'etaChirp_systemPlot',
               'aDCOEccentricityDCO_systemPlot',
               'aPreSN_ePreSN_systemPlot',
               'VrelVkick_systemPlot',
               'chirpSeparationZAMS_systemPlot',
               'chirpSeparationDCO_systemPlot']

def plottingRoutines():
    """""""""
    Data loading
    """""""""
    
    # Loading in initial data
    initial_data = pandas.read_table(initialdata_name, header = 1, 
                                     delim_whitespace = True)
    
    # Loads in merging data
    mergingData_raw = pandas.read_table(mergingdata_name, header = 1, 
                                        delim_whitespace = True)
    
    """
    Exclusion checks for Hubble time, RLOF and pessimistic prescription
    """
    
    # If statement to exclude/include binaries merging outside of a Hubble time, 
    # as per the user's request
    if Hubble == True:
        # Checks to see if binaries merge in a Hubble time
        HubbleFlag_data = mergingData_raw['mergesInHubbleTimeFlag']
        # Creates a truth table of events merging in less than a Hubble time
        HubbleFlag_truth = np.logical_not(HubbleFlag_data == 0)
        # Pulls out data that has mergers in less than a Hubble time
        merging_data1 = mergingData_raw[HubbleFlag_truth]
        print 'Binaries merging in over a Hubble time have been excluded.'
    else:
        merging_data1 = mergingData_raw # All merging data is kept
        print 'Binaries merging in over a Hubble time have been included.'
    
    # If statement to exclude/include binaries whose secondary star undergoes RLOF 
    # immediately after a CEE
    if RLOF == True:
        RLOFFlag_data = merging_data1['doubleCommonEnvelopeFlag']
        # Creates truth table of events which don't undergo RLOF immediately after 
        # a CEE
        RLOFFlag_truth = np.logical_not(RLOFFlag_data == 1)
        # Pulls out data which don't undergo RLOF immediately after a CEE
        merging_data2 = merging_data1[RLOFFlag_truth]
        print 'Binaries whose secondary star undergoes RLOF immediately after a CEE have been exlcuded.'
    else:
        merging_data2 = merging_data1 # Keeps all merging data after Hubble time 
                                        # exclusion check
        print 'Binaries whose secondary star undergoes RLOF immediately after a CEE have been included.'
    
    # If statement to exclude/include binaries which use the pessimistic CE 
    # prescription
    if Pessimistic == True:
        OptimisticFlag_data = merging_data2['optimisticCEFlag']
        # Creates truth table of events which don't use the optimistic CE 
        # prescription
        PessimisticFlag_truth = np.logical_not(OptimisticFlag_data == 1)
        # Pulls out data which don't use the optimistic CE prescription
        merging_data = merging_data2[PessimisticFlag_truth]
        print 'Binaries which use the optimistic CE prescription have been exlcuded.'
    else:
        merging_data = merging_data2
        print 'Binaries which use the optimistic CE prescription have been included.'
    
    """
    Output path checks
    """
    
    # True if output_path directory exists, False if not
    output_path_existenceCheck = os.path.isdir(output_path)
    # If statement to create output_path directory if the directory doesn't exist
    if output_path_existenceCheck == False:
        print 'The directory \'' + output_path + '\' does not exist, therefore creating this directory.'
        os.makedirs(output_path)
        print output_path, 'has been created!'
    
    # If statement to perform check relevant to the webpage alone
    if html_choice == True:
        # True if html_path directory exists, False if not
        htmlPath_existenceCheck = os.path.isdir(html_path)
        # If statement to create html_path directory if the directory doesn't exist
        if htmlPath_existenceCheck == False:
            print 'The directory \'' + html_path + '\' does not exist, therefore creating this directory.'
            os.makedirs(html_path)
            print html_path, 'has been created!'
    
    """""""""
    Data separation
    """""""""
    
    # Creating truth tables for DCOs which belong to a specific population
    BNS_truth = np.logical_and(merging_data['stellarType1'] == 13, 
                               merging_data['stellarType2'] == 13)
    BHNS_truth = (np.logical_and(merging_data['stellarType1'] == 14, 
                                (merging_data['stellarType2'] == 13)) + 
                  np.logical_and(merging_data['stellarType1'] == 13, 
                                merging_data['stellarType2'] == 14))
    BBH_truth = np.logical_and(merging_data['stellarType1'] == 14, 
                               merging_data['stellarType2'] == 14)
    
    # Pulling out the specific DCO population data
    BNS_data = merging_data[BNS_truth]
    BHNS_data = merging_data[BHNS_truth]
    BBH_data = merging_data[BBH_truth]
    
    """""""""
    Summary information
    """""""""
    
    table_header = ["All DCO", "BNS", "BHNS", "BBH"]
    exclusions = ["All", "Hubble", 
                  "RLOF", "Pessimistic"]
    summary = summaryTable(mergingData_raw)
    
    row_format ="{:>15}" * (len(table_header) + 1)
    print row_format.format("", *table_header)
    for exclusion, row in zip(exclusions, summary):
        print row_format.format(exclusion, *row)
    
    # Creating columns for HTML webpage
    col1 = ["", "All", "Mergers in a Hubble time", "and no post-CE RLOF", 
            "and pessimistic CE prescription"]
    col2 = ["All DCO"] + [row[0] for row in summary]
    col3 = ["BNS"] + [row[1] for row in summary]
    col4 = ["BHNS"] + [row[2] for row in summary]
    col5 = ["BBH"] + [row[3] for row in summary]
    table = zip(col1, col2, col3, col4, col5)
    
    
    f = open(programOptions_name)
    programOptions = f.readline()
    f.close()
    
    print programOptions

    """""""""
    Cycling through plots
    """""""""
    
    for i, name in enumerate(plots_names):
        plots[name](initial_data, merging_data, BNS_data, BHNS_data, BBH_data)
    
    """""""""
    HTML webpage
    """""""""
    
    webpage_numbers = np.arange(len(plots_names)) + 1
    
    webpage_figures = zip(webpage_titles, webpage_descriptions, webpage_imagepaths,
                          webpage_numbers)
    
    env = jj2.Environment(loader = jj2.FileSystemLoader('.'))
    temp = env.get_template('template.html')
    templateVars = {"title" : webpage_name,
                    "pagetitle" : webpage_title,
                    "webpage_figures" : webpage_figures,
                    "programOptions" : programOptions,
                    "summaryTable" : table,
                    "Hubble" : str(Hubble),
                    "RLOF" : str(RLOF),
                    "Pessimistic" : str(Pessimistic)}
    r = temp.render(templateVars)
    
    webpage_namePath = html_path + webpage_name
    f = open(webpage_namePath, 'w')
    f.write(r)
    f.close()
    
    webpage_commandLine = 'chmod -R 755 ' + html_path   # Bash command
    os.system(webpage_commandLine)                      # Executes bash command