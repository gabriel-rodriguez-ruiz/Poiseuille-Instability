#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 10:42:33 2019

@author: gabriel
"""

import numpy as np


def basic_solution(y):
    """
    Dimentionless basic solution for an undisturbed plane Poiseuille flow.
    
    Parameters
    ----------
    y : float
        Coordinate pependicular to the flow.
   
    Returns
    -------
    float
        The value of the flow profile at height y.
    """
    return 1 - y**2

def H_cubic(i, y, h):
    """
    Hermite interpolation function for a cubic element.
    
    Parameters
    ----------
    i : int
        The order of the Hermite function between 1 and 4.
    y : float
        Coordinate pependicular to the flow.
    h : float
        Element size.
    
    Returns
    -------
    float
        The value of the ith Hermite function at height y.
    """
    xi_contour = [-1, 1]
    xi = 2*y/h
    if i==1 or i==2:
        return (xi + xi_contour[i-1])**2 * (xi*xi_contour[i-1] - 2)/4
    elif i==3 or i==4:
        return -xi_contour[i-1]*(xi + xi_contour[i-1])**2 * (xi*xi_contour[i-1] - 1)*h/8            
