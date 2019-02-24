#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 23 10:42:33 2019

@author: gabriel
"""

import numpy as np
from scipy.integrate import quad
import scipy

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

def H_cubic(y, i):
    """
    Hermite interpolation function for a cubic element.
    
    Parameters
    ----------
    i : int
        Number of the hermite funcion between 0 and 3.
    y : float
        Coordinate pependicular to the flow.
    
    Returns
    -------
    float
        The value of the ith Hermite function at height y.
    """
    if i==0:
        return (y + 2)*(1 - y)**2/4
    elif i==1:
        return (y + 1)*(1 - y)/4
    elif i==2:
        return (y + 1)**2 * (2 - y)/4
    elif i==3:
        return (y - 1)*(1 + y)**2/4
        
def H_first_derivative(y, i):
    """
    First derivative of the Hermite interpolation function for a cubic element.
    
    Parameters
    ----------
    i : int
        Number of the hermite funcion between 0 and 3.
    y : float
        Coordinate pependicular to the flow.
    
    Returns
    -------
    float
        The value of the first derivative of the ith Hermite function at height y.
    """        
    if i==0:
        return (-y + 1)**2/4 + (y + 2)*(2*y - 2)/4
    elif i==1:
        return -y/2
    elif i==2:
        return (-y + 2)*(2*y + 2)/4 - (y + 1)**2/4
    elif i==3:
        return (y - 1)*(2*y + 2)/4 + (y + 1)**2/4
    
def H_second_derivative(y, i):
    """
    Second derivative of the Hermite interpolation function for a cubic element.
    
    Parameters
    ----------
    i : int
        Number of the hermite funcion between 0 and 3.
    y : float
        Coordinate pependicular to the flow.
    
    Returns
    -------
    float
        The value of the second derivative of the ith Hermite function at height y.
    """        
    if i==0:
        return 3*y/2
    elif i==1:
        return -1/2
    elif i==2:
        return -3*y/2
    elif i==3:
        return 3*y/2 + 1/2


def A(N, h, k, alpha, Re):
    """
    The matrix A.
    
    Parameters
    ----------
    N : int
        Number of finite elements.
    h : float
        Element step.
    k : int
        Element number between 0 and N-1.
    Re : float
        Reynolds number.
    alpha : float
        Wavenumber.
    Returns
    -------
    A : array
        The matrix A for the element k.
    """
    M = 2*N
    A = np.zeros((4, 4), dtype=complex)
    for i in range(4):
        for j in range(4):
                
                integrand = lambda y,s,p: H_second_derivative(y,s)*H_second_derivative(y,p)\
                            + 2*alpha**2*H_first_derivative(y,s)*H_first_derivative(y,p)\
                            + alpha**4*H_cubic(y,s)*H_cubic(y,p)\
                            + 1j*alpha*Re*basic_solution(y)*(H_second_derivative(y,s) - alpha**2*H_cubic(y,s))*H_cubic(y,p)\
                            - 2*H_cubic(y,s)*H_cubic(y,p)
                integrand_real = lambda y,s,p: scipy.real(integrand(y,s,p))
                integrand_imaginary = lambda y,s,p: scipy.imag(integrand(y,s,p))               
                A[i, j] = quad(integrand_real, k*h, (k+1)*h, args=(j,i))[0] + 1j*quad(integrand_imaginary, k*h, (k+1)*h, args=(j,i))[0]
    return A
