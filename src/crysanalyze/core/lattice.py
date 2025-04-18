"""
Lattice parameter calculations and WH plot analysis
"""

import numpy as np
from scipy.optimize import curve_fit
from .indexing import calculate_theoretical_d

def calculate_lattice_parameters(indexed_peaks, crystal_system):
    """
    Calculate refined lattice parameters from indexed peaks.
    
    Args:
        indexed_peaks (list): List of indexed peaks
        crystal_system (str): Crystal system
        
    Returns:
        dict: Refined lattice parameters
    """
    def error_function(params, peaks, crystal_system):
        error = 0
        a, b, c = params[:3]
        angles = [90, 90, 90]  # Default angles
        
        if crystal_system == 'triclinic':
            angles = params[3:]
        elif crystal_system == 'monoclinic':
            angles[1] = params[3]
        
        for peak in peaks:
            hkl = peak['hkl']
            d_exp = peak['d_exp']
            d_calc = calculate_theoretical_d(hkl, a, b, c, *angles)
            error += (d_exp - d_calc)**2
        
        return np.sqrt(error)
    
    # Initial guess
    if crystal_system == 'cubic':
        params = [4.0]  # Single parameter
    elif crystal_system == 'tetragonal':
        params = [4.0, 4.0]  # a and c
    elif crystal_system == 'hexagonal':
        params = [4.0, 4.0]  # a and c
    elif crystal_system == 'orthorhombic':
        params = [4.0, 4.0, 4.0]  # a, b, c
    elif crystal_system == 'monoclinic':
        params = [4.0, 4.0, 4.0, 90]  # a, b, c, beta
    else:  # triclinic
        params = [4.0, 4.0, 4.0, 90, 90, 90]  # a, b, c, alpha, beta, gamma
    
    # Refine parameters
    result = minimize(
        error_function,
        params,
        args=(indexed_peaks, crystal_system),
        method='L-BFGS-B'
    )
    
    # Format results
    if crystal_system == 'cubic':
        return {'a': result.x[0], 'b': result.x[0], 'c': result.x[0],
                'alpha': 90, 'beta': 90, 'gamma': 90}
    elif crystal_system == 'tetragonal':
        return {'a': result.x[0], 'b': result.x[0], 'c': result.x[1],
                'alpha': 90, 'beta': 90, 'gamma': 90}
    elif crystal_system == 'hexagonal':
        return {'a': result.x[0], 'b': result.x[0], 'c': result.x[1],
                'alpha': 90, 'beta': 90, 'gamma': 120}
    elif crystal_system == 'orthorhombic':
        return {'a': result.x[0], 'b': result.x[1], 'c': result.x[2],
                'alpha': 90, 'beta': 90, 'gamma': 90}
    elif crystal_system == 'monoclinic':
        return {'a': result.x[0], 'b': result.x[1], 'c': result.x[2],
                'alpha': 90, 'beta': result.x[3], 'gamma': 90}
    else:  # triclinic
        return {'a': result.x[0], 'b': result.x[1], 'c': result.x[2],
                'alpha': result.x[3], 'beta': result.x[4], 'gamma': result.x[5]}

def calculate_wh_plot(indexed_peaks, wavelength):
    """
    Calculate Williamson-Hall plot for strain and crystallite size analysis.
    
    Args:
        indexed_peaks (list): List of indexed peaks
        wavelength (float): X-ray wavelength in Angstroms
        
    Returns:
        dict: WH plot parameters and data
    """
    # Calculate 1/d and βcosθ
    d_inv = []
    beta_cos_theta = []
    
    for peak in indexed_peaks:
        hkl = peak['hkl']
        two_theta = peak['position']
        fwhm = peak.get('fwhm', 0.1)  # Default FWHM if not provided
        
        theta = np.radians(two_theta / 2)
        d = 1 / peak['d_exp']
        
        d_inv.append(d)
        beta_cos_theta.append(fwhm * np.cos(theta))
    
    d_inv = np.array(d_inv)
    beta_cos_theta = np.array(beta_cos_theta)
    
    # Fit linear function
    def linear_func(x, a, b):
        return a * x + b
    
    popt, pcov = curve_fit(linear_func, d_inv, beta_cos_theta)
    
    # Calculate parameters
    strain = popt[0]  # Slope = strain
    crystallite_size = wavelength / (popt[1] * np.pi)  # Intercept = λ/(π*size)
    
    return {
        'strain': strain,
        'crystallite_size': crystallite_size,
        'd_inv': d_inv,
        'beta_cos_theta': beta_cos_theta,
        'fit_params': popt,
        'fit_cov': pcov
    }

def calculate_crystallinity(data, amorphous_region=(10, 20)):
    """
    Calculate crystallinity from XRD data.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        amorphous_region (tuple): 2θ range for amorphous contribution
        
    Returns:
        float: Crystallinity percentage
    """
    x, y = data.T
    
    # Find amorphous region
    mask = (x >= amorphous_region[0]) & (x <= amorphous_region[1])
    amorphous_intensity = np.mean(y[mask])
    
    # Subtract amorphous background
    y_corrected = y - amorphous_intensity
    y_corrected[y_corrected < 0] = 0
    
    # Calculate areas
    total_area = np.trapz(y_corrected, x)
    crystalline_area = np.trapz(y_corrected[y_corrected > 0], x[y_corrected > 0])
    
    # Calculate crystallinity
    crystallinity = (crystalline_area / total_area) * 100
    
    return crystallinity 