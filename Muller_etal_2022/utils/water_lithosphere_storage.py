# -*- coding: utf-8 -*-
import pygplates
import numpy as np
import math
from scipy.special import erfinv
import time as tme
#functions to measure water content in ocean crust
#water calcs. taken from Jarrard (2003) and Merdith et al. (2019)
#lithospheric thickness taken from East et al. (2019)

T1 = 1150.
To = 0.
Tm = 1350.
kappa = 0.804e-6
Myr2sec=1e6*365*24*60*60

def plate_isotherm_depth(age, temp, *vartuple) :
    '''
    Computes the depth to the temp - isotherm in a cooling plate mode.\
    Solution by iteration. By default the plate thickness is 125 km as\
    in Parsons/Sclater.  Change given a 3rd parameter. Originally written
    en by Simon Williams(?), taken from East et al. (2019).
    '''

    if len(vartuple) != 0 :
        PLATE_THICKNESS_KM = vartuple[0]
    else :
        PLATE_THICKNESS_KM = 125

    PLATE_THICKNESS = PLATE_THICKNESS_KM * 1000


    # default depth is 0
    z = 0

    if age <= 0.0 :
        z_try = 0
        done = 1
    else :
        z_too_small = 0.0
        z_too_big = PLATE_THICKNESS
        done = 0
        n_try = 0

    while done != 1 and n_try < 20 :
        n_try += 1
        z_try = 0.5 * (z_too_small + z_too_big)
        t_try = plate_temp (age, z_try, PLATE_THICKNESS)
        t_wrong = temp - t_try

        if t_wrong < -0.001 :
            z_too_big = z_try
        elif t_wrong > 0.001 :
            z_too_small = z_try
        else :
            done = 1

        z = z_try
    return z

def plate_temp(age, z, PLATE_THICKNESS) :
    "Computes the temperature in a cooling plate for age = t\
    and at a depth = z."

    KAPPA = 0.804E-6
    T_MANTLE = 1350.0
    T_SURFACE = 0.0
    SEC_PR_MY = 3.15576e13

    t = T_SURFACE

    sum = 0
    sine_arg = math.pi * z / PLATE_THICKNESS
    exp_arg = -KAPPA * math.pi * math.pi * age * SEC_PR_MY / (PLATE_THICKNESS * PLATE_THICKNESS)
    for k in range(1, 20) :
        sum = sum + np.sin(k * sine_arg) * np.exp(k*k*exp_arg)/k

    if age <= 0.0 :
        t = T_MANTLE * np.ones(z.shape)
    else :
        t = t + 2.0 * sum * (T_MANTLE - T_SURFACE)/math.pi + (T_MANTLE - T_SURFACE) * z/PLATE_THICKNESS

    return t

def porosity(layer,age):

    if layer == 'Upper':
        microporosity = 7.8
        macroporosity = 13.01-5.625*np.log10(age)

    if layer == 'Lower':
        microporosity = 5.1
        macroporosity = (13.01-5.625*np.log10(age))/2.

    if layer == 'Dykes':
        microporosity = 2.2
        macroporosity = 0.84

    if layer == 'Gabbros':
        microporosity = 0.7
        macroporosity = 0

    total_porosity = microporosity + macroporosity

    return total_porosity, microporosity, macroporosity

def density(layer,age, total_porosity):

    fluid_density = 1.02

    if layer == 'Upper':
        matrix_density = 3.01-0.0631*np.log10(age)
    if layer == 'Lower':
        matrix_density = 3.01-0.5*0.0631*np.log10(age)
    if layer == 'Dykes':
        matrix_density = 2.98
    if layer == 'Gabbros':
        matrix_density = 2.99

    bulk_density = (0.01*total_porosity)*fluid_density + (1 - 0.01 * total_porosity) * matrix_density
    return bulk_density, matrix_density

def structural_H2O(layer, macroporosity, matrix_density):
    fluid_density = 1.02

    if layer == 'Upper' or layer == 'Lower':
        H2OS = (103.1 -34.27*matrix_density) + (0.17 * (13.01-macroporosity))
    if layer == 'Dykes':
        H2OS = 1.76
    if layer == 'Gabbros':
        H2OS = 0.79
    return H2OS

def pore_H2O(total_porosity, bulk_density):
    fluid_density = 1.02
    #print total_porosity, bulk_density
    H2OP = total_porosity * fluid_density/bulk_density

    return H2OP

def H2O_thickness(H2OS, H2OP, thickness, bulk_density):
    #print H2OS, H2OP, thickness, bulk_density
    H2OS_thick = H2OS * (1 - 0.01*H2OP)*thickness*bulk_density * 1/100
    H2OP_thick = H2OP * thickness*bulk_density * 1/100

    total_H2O_thick = H2OS_thick + H2OP_thick

    return total_H2O_thick, H2OS_thick, H2OP_thick

def get_total_H2O_thickness(layer, thick, ages):

    #layer_type = #name of layer
    total_porosity, microporosity, macroporosity = porosity(layer,ages)
    bulk_density, matrix_density = density(layer, ages, total_porosity)
    H2OS = structural_H2O(layer,macroporosity, matrix_density)
    H2OP = pore_H2O(total_porosity, bulk_density)
    total_H2O_thick, H2OS_thick, H2OP_thick = H2O_thickness(H2OS, H2OP, thick, bulk_density)

    return total_H2O_thick, H2OS_thick, H2OP_thick

def get_lithosphere_thickness(ages):
    #print np.shape(ages)
    lithos_thickos = np.zeros_like(ages)
    for ind, age in enumerate(ages) :
        #print i
    #for i in range(0, 100) :
        lithos_thickos[ind] = plate_isotherm_depth(age, T1)

    return lithos_thickos
