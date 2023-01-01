import datetime
from dateutil.relativedelta import relativedelta
import numpy as np
import os
import sys
import cv2
import math
import matplotlib.pyplot as plt
from numba import jit, njit
import pandas as pd

"""
This python script aims to calculate Boltzmann entropy of numerical data organized in raster data
"""

"""
Reference:
Gao, P., Zhang, H., & Li, Z. (2017). A hierarchy-based solution to calculate the configurational entropy of landscape gradients. Landscape ecology, 32(6), 1133-1146.
Gao, P., & Li, Z. (2019). Aggregation-based method for computing absolute Boltzmann entropy of landscape gradient with full thermodynamic consistency. Landscape Ecology, 34(8), 1837-1847.
Gao, P., & Li, Z. (2019). Computation of the Boltzmann entropy of a landscape: A review and a generalization. Landscape Ecology, 34(9), 2183-2196.
#regarding the problems of using this python script, please feel free to send an email to cxh9791156936@gmail.com
"""

# calculate Boltzmann entropy of nominal raster data
@jit(nopython=True, parallel=True, fastmath=True)
def One_to_four_nominal(pixel_data_window):
    """
    this function is the key one part of calculating Boltzmann entropy of nominal raster data (see Gao and Li (2019)).
    :param pixel_data_window: pixels included by a 2*2 window. Date type is nominal, e.g., cropland, water,impervious surface.
    :return: Boltzmann entropy for categorical raster data
    """
    res = 0
    # statistics of the number of kinds of particles

    arrWindow = pixel_data_window

    list_widonw = list(arrWindow)

    list_set = set(list_widonw)

    num_nominal = len(list_set)

    if num_nominal == 4:

        res = 24

    elif num_nominal == 3:

        res = 8

    elif num_nominal == 2:

        res = 6

    else:
        res = 1

    return res

@jit(nopython=True, parallel=True, fastmath=True)
def One_to_four(pixel_data_window):
    """
    this function is the key one part of calculating Boltzmann entropy and utilizes a 2*2 sliding window
    :param pixel_data_window: pixels included by a 2*2 window
    :return: the number of microstates, W.
    """
    res = 0
    Mean = np.around(np.mean(pixel_data_window))
    Sum = np.sum(pixel_data_window)
    Max = np.max(pixel_data_window)
    Min = np.min(pixel_data_window)

    Tempvalue = (Sum - Max - Min) / 2
    xa = np.floor(Tempvalue)
    xb = np.ceil(Tempvalue)
    da = xa - Min
    db = Max - xb
    d = min(da, db)

    if d == 0:
        if da == db:
            if xa == xb:
                res = 1
            else:
                res = 6
        else:
            if xa == xb:
                res = 4
            else:
                res = 12
    else:
        if da == db:
            if xa == xb:
                res = 12 + 24 * (d - 1) + 6
            else:
                res = 24 + 24 * (d - 1) + 6
        else:
            if xa == xb:
                res = 12 + 24 * (d - 1) + 12
            else:
                res = 24 + 24 * (d - 1) + 12

    return res, Mean

@jit(nopython=True,parallel=True,fastmath=True)
def Iteration_Resampling(banddata,mybase):
    """
    using a 2*2 window that slides with overlap to calculate Boltzmann entropy
    :param banddata: a two-dimensional array
    :param mybase: the base for calculating Boltzmann entropy
    :return:
    """
    rows_n, cols_n = banddata.shape
    rows_n, cols_n = rows_n - 1, cols_n - 1
    # data_n = np.array(rows_n,cols_n)
    new_grid = np.zeros((rows_n, cols_n), dtype=np.int16)
    result = np.zeros((rows_n, cols_n), dtype=np.float32)
    for rn in range(rows_n):
        for cn in range(cols_n):
            pixel_data_window = np.array(
                [banddata[rn, cn], banddata[rn + 1, cn], banddata[rn, cn + 1], banddata[rn + 1, cn + 1]])

            # pixel_data_window = pixel_data_window.flatten()
            temp, Mean = One_to_four(pixel_data_window)
            new_grid[rn, cn] = Mean
            # log_a = np.array([temp])
            # log_value = np.log2(np.array([temp]))
            result[rn, cn] = np.log(np.array([temp]))[0]/np.log(np.array([mybase]))[0]  # np.log([temp], mybase)  #math.log(temp, mybase)
    res_value = np.sum(result.flatten())
    return res_value,new_grid

@jit(nopython=True,parallel=True,fastmath=True)
def Iteration_Aggragation(banddata,mybase):
    """
    using a 2*2 window that slides without overlap to calculate Boltzmann entropy
    :param banddata: a two-dimensional array
    :param mybase: the base for calculating Boltzmann entropy
    :return:
    """
    RBE = 0
    rows_n, cols_n = banddata.shape
    rows_n, cols_n = int(rows_n / 2), int(cols_n / 2)
    new_grid2 = np.zeros((rows_n, cols_n), dtype=np.int16)

    for rn in range(rows_n):
        for cn in range(cols_n):
            pixel_data_window = np.array(
                [banddata[2 * rn, 2 * cn], banddata[2 * rn + 1, 2 * cn], banddata[2 * rn, 2 * cn + 1],
                 banddata[2 * rn + 1, 2 * cn + 1]])

            # pixel_data_window = pixel_data_window.flatten()
            temp, Mean = One_to_four(pixel_data_window)

            # log_a = np.array([temp])
            s = np.log(np.array([temp]))[0]/np.log(np.array([mybase]))[0]

            RBE += s
            new_grid2[rn, cn] = Mean

    return RBE,new_grid2

def GetBoltzmann_Resampling(banddata,mybase=2,relative = False):
    """
    calculating Boltzmann entropy based on a resampling-based method (see Gao and Li (2019))
    :param banddata: a two-dimensional array
    :param mybase: the base for calculating Boltzmann entropy
    :param relative: when the relative is "False", you will get the absolute Boltzmann entropy value that captures the multiscale information.
                     when the relative is "True", the relative Boltzmann entropy will be calculated.
    :return:
    """
    res_value = 0
    while(banddata.shape[0]>1 and banddata.shape[1]>1):

        rows_n, cols_n = banddata.shape
        rows_n, cols_n = rows_n - 1, cols_n - 1
        # data_n = np.array(rows_n,cols_n)
        # new_grid = np.zeros((rows_n, cols_n), dtype=np.int)
        res_value,new_grid2 = Iteration_Resampling(banddata, mybase)
        # print(res_value)
        # for rn in range(rows_n):
        #     for cn in range(cols_n):
        #         res_value += result[rn, cn]

        if (relative == True):
            break
        else:
            banddata = new_grid2

    return res_value

#calculating Boltzmann entropy based on an aggregation-based methond (see Gao and Li (2019))
def GetBoltzmann_Aggragation(banddata,mybase=2,relative = False):
    ## bandata is a two dimensional array
    ABE =0
    RBEs = []
    while(banddata.shape[0]>1 and banddata.shape[1]>1):
        RBE,new_grid = Iteration_Aggragation(banddata,mybase)

        RBEs.append(RBE) # appending Boltzmann entropy of raster data for each spatial transformation from fine to coarse scale.
        banddata = new_grid #recording the raster data at the next spatial scale.

    ABE = sum(RBEs) #calculating absolute Boltzmann entropy by summing relative Boltzmann entropy of all spatial transformation.

    if (relative == True):
        return RBEs[0] #relative Boltzmann entropy
    else:
        return ABE  #absolute Boltzmann entropy

    return ABE, RBEs[0]


