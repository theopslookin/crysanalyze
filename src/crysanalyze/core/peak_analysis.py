"""
Advanced peak analysis module for XRD data
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from scipy.special import voigt_profile

def fit_peak(x, y, peak_type='pseudo_voigt', background='linear'):
    """
    Fit a single peak with specified function.
    
    Args:
        x (numpy.ndarray): x values
        y (numpy.ndarray): y values
        peak_type (str): Type of peak function ('gaussian', 'lorentzian', 'pseudo_voigt')
        background (str): Type of background ('linear', 'quadratic', 'none')
        
    Returns:
        dict: Peak parameters and fit information
    """
    def gaussian(x, a, x0, sigma):
        return a * np.exp(-((x - x0) / sigma)**2)
    
    def lorentzian(x, a, x0, gamma):
        return a / (1 + ((x - x0) / gamma)**2)
    
    def pseudo_voigt(x, a, x0, sigma, eta):
        g = gaussian(x, a, x0, sigma)
        l = lorentzian(x, a, x0, sigma)
        return eta * l + (1 - eta) * g
    
    # Select peak function
    if peak_type == 'gaussian':
        func = gaussian
        p0 = [max(y), x[np.argmax(y)], 0.1]
    elif peak_type == 'lorentzian':
        func = lorentzian
        p0 = [max(y), x[np.argmax(y)], 0.1]
    elif peak_type == 'pseudo_voigt':
        func = pseudo_voigt
        p0 = [max(y), x[np.argmax(y)], 0.1, 0.5]
    else:
        raise ValueError(f"Unknown peak type: {peak_type}")
    
    # Add background if specified
    if background == 'linear':
        def model(x, *params):
            peak_params = params[:-2]
            m, c = params[-2:]
            return func(x, *peak_params) + m * x + c
        p0.extend([0, min(y)])
    elif background == 'quadratic':
        def model(x, *params):
            peak_params = params[:-3]
            a, b, c = params[-3:]
            return func(x, *peak_params) + a * x**2 + b * x + c
        p0.extend([0, 0, min(y)])
    else:
        model = func
    
    # Fit peak
    popt, pcov = curve_fit(model, x, y, p0=p0)
    
    # Calculate FWHM
    if peak_type == 'gaussian':
        fwhm = 2 * np.sqrt(2 * np.log(2)) * popt[2]
    elif peak_type == 'lorentzian':
        fwhm = 2 * popt[2]
    else:  # pseudo_voigt
        fwhm_g = 2 * np.sqrt(2 * np.log(2)) * popt[2]
        fwhm_l = 2 * popt[2]
        fwhm = popt[3] * fwhm_l + (1 - popt[3]) * fwhm_g
    
    return {
        'position': popt[1],
        'amplitude': popt[0],
        'fwhm': fwhm,
        'parameters': popt,
        'covariance': pcov,
        'type': peak_type
    }

def deconvolute_peaks(x, y, num_peaks=2, peak_type='pseudo_voigt'):
    """
    Deconvolute overlapping peaks.
    
    Args:
        x (numpy.ndarray): x values
        y (numpy.ndarray): y values
        num_peaks (int): Number of peaks to fit
        peak_type (str): Type of peak function
        
    Returns:
        list: List of peak parameters
    """
    # Find initial peak positions
    peaks, _ = find_peaks(y, height=max(y)*0.1, distance=len(x)//10)
    if len(peaks) < num_peaks:
        raise ValueError(f"Found only {len(peaks)} peaks, but requested {num_peaks}")
    
    # Sort peaks by intensity
    peak_intensities = y[peaks]
    peak_order = np.argsort(peak_intensities)[::-1]
    peaks = peaks[peak_order[:num_peaks]]
    
    # Create combined model
    def model(x, *params):
        result = np.zeros_like(x)
        for i in range(num_peaks):
            if peak_type == 'gaussian':
                result += gaussian(x, *params[i*3:(i+1)*3])
            elif peak_type == 'lorentzian':
                result += lorentzian(x, *params[i*3:(i+1)*3])
            else:  # pseudo_voigt
                result += pseudo_voigt(x, *params[i*4:(i+1)*4])
        return result
    
    # Initial parameters
    p0 = []
    for peak in peaks:
        if peak_type == 'pseudo_voigt':
            p0.extend([y[peak], x[peak], 0.1, 0.5])
        else:
            p0.extend([y[peak], x[peak], 0.1])
    
    # Fit peaks
    popt, pcov = curve_fit(model, x, y, p0=p0)
    
    # Extract individual peak parameters
    results = []
    for i in range(num_peaks):
        if peak_type == 'pseudo_voigt':
            params = popt[i*4:(i+1)*4]
            fwhm_g = 2 * np.sqrt(2 * np.log(2)) * params[2]
            fwhm_l = 2 * params[2]
            fwhm = params[3] * fwhm_l + (1 - params[3]) * fwhm_g
        else:
            params = popt[i*3:(i+1)*3]
            fwhm = 2 * np.sqrt(2 * np.log(2)) * params[2]
        
        results.append({
            'position': params[1],
            'amplitude': params[0],
            'fwhm': fwhm,
            'parameters': params,
            'type': peak_type
        })
    
    return results

def calculate_asymmetry(x, y, peak_position):
    """
    Calculate peak asymmetry.
    
    Args:
        x (numpy.ndarray): x values
        y (numpy.ndarray): y values
        peak_position (float): Peak position
        
    Returns:
        float: Asymmetry factor
    """
    # Find peak maximum
    peak_idx = np.argmin(np.abs(x - peak_position))
    
    # Calculate FWHM
    half_max = y[peak_idx] / 2
    left_idx = np.argmin(np.abs(y[:peak_idx] - half_max))
    right_idx = peak_idx + np.argmin(np.abs(y[peak_idx:] - half_max))
    
    # Calculate asymmetry
    left_width = peak_position - x[left_idx]
    right_width = x[right_idx] - peak_position
    
    return right_width / left_width

def calculate_integral_breadth(x, y):
    """
    Calculate integral breadth of a peak.
    
    Args:
        x (numpy.ndarray): x values
        y (numpy.ndarray): y values
        
    Returns:
        float: Integral breadth
    """
    # Calculate peak area
    area = np.trapz(y, x)
    
    # Calculate maximum intensity
    max_intensity = np.max(y)
    
    # Calculate integral breadth
    return area / max_intensity

def correct_instrumental_broadening(peak_fwhm, standard_fwhm):
    """
    Correct for instrumental broadening using a standard.
    
    Args:
        peak_fwhm (float): FWHM of the peak
        standard_fwhm (float): FWHM of the standard
        
    Returns:
        float: Corrected FWHM
    """
    # Use quadratic subtraction
    return np.sqrt(peak_fwhm**2 - standard_fwhm**2)

def calculate_crystallite_size(fwhm, wavelength, theta, shape_factor=0.9):
    """
    Calculate crystallite size using Scherrer equation.
    
    Args:
        fwhm (float): FWHM in radians
        wavelength (float): X-ray wavelength in Angstroms
        theta (float): Bragg angle in radians
        shape_factor (float): Shape factor (typically 0.9)
        
    Returns:
        float: Crystallite size in Angstroms
    """
    return shape_factor * wavelength / (fwhm * np.cos(theta)) 