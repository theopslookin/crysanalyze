"""
Data processing module for XRD analysis
"""

import numpy as np
from scipy.signal import savgol_filter, find_peaks
from scipy.ndimage import gaussian_filter1d

def subtract_background(data, window_size=51, poly_order=3):
    """
    Subtract background from XRD data using Savitzky-Golay filter.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        window_size (int): Size of the smoothing window
        poly_order (int): Order of the polynomial for fitting
        
    Returns:
        numpy.ndarray: Background-subtracted data
    """
    x, y = data.T
    background = savgol_filter(y, window_size, poly_order)
    return np.column_stack((x, y - background))

def find_peaks(data, min_height=100, min_distance=5):
    """
    Find peaks in XRD data.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        min_height (float): Minimum peak height
        min_distance (float): Minimum distance between peaks
        
    Returns:
        numpy.ndarray: Array of (2θ, intensity) pairs for each peak
    """
    x, y = data.T
    peaks, properties = find_peaks(y, height=min_height, distance=min_distance)
    return np.column_stack((x[peaks], y[peaks]))

def smooth_data(data, sigma=1):
    """
    Smooth XRD data using Gaussian filter.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        sigma (float): Standard deviation for Gaussian kernel
        
    Returns:
        numpy.ndarray: Smoothed data
    """
    x, y = data.T
    smoothed_y = gaussian_filter1d(y, sigma)
    return np.column_stack((x, smoothed_y))

def normalize_intensity(data):
    """
    Normalize intensity values to [0, 1] range.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        
    Returns:
        numpy.ndarray: Normalized data
    """
    x, y = data.T
    y_norm = (y - y.min()) / (y.max() - y.min())
    return np.column_stack((x, y_norm)) 