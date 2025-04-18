"""
Peak indexing module for XRD analysis
"""

import numpy as np
from itertools import product
from scipy.optimize import minimize
from .processing import calculate_d_spacing

# Crystal system parameters
CRYSTAL_SYSTEMS = {
    'cubic': {'constraints': lambda a, b, c, alpha, beta, gamma: (a == b == c) and (alpha == beta == gamma == 90)},
    'tetragonal': {'constraints': lambda a, b, c, alpha, beta, gamma: (a == b) and (alpha == beta == gamma == 90)},
    'hexagonal': {'constraints': lambda a, b, c, alpha, beta, gamma: (a == b) and (alpha == beta == 90) and (gamma == 120)},
    'orthorhombic': {'constraints': lambda a, b, c, alpha, beta, gamma: (alpha == beta == gamma == 90)},
    'monoclinic': {'constraints': lambda a, b, c, alpha, beta, gamma: (alpha == gamma == 90)},
    'triclinic': {'constraints': lambda a, b, c, alpha, beta, gamma: True}
}

def index_peaks(peaks, wavelength, crystal_system='cubic', max_hkl=5, tolerance=0.05):
    """
    Index XRD peaks to determine Miller indices.
    
    Args:
        peaks (numpy.ndarray): Array of (2θ, intensity) pairs
        wavelength (float): X-ray wavelength in Angstroms
        crystal_system (str): Crystal system to test
        max_hkl (int): Maximum Miller index to consider
        tolerance (float): Positional tolerance in degrees
        
    Returns:
        dict: Dictionary containing indexed peaks and lattice parameters
    """
    # Calculate d-spacings
    d_spacings = calculate_d_spacing(peaks[:, 0], wavelength)
    
    # Generate possible hkl combinations
    hkl_combinations = list(product(range(-max_hkl, max_hkl+1), repeat=3))
    hkl_combinations = [hkl for hkl in hkl_combinations if sum(hkl) > 0]  # Remove 000 and equivalent combinations
    
    # Initialize results
    indexed_peaks = []
    best_fit = None
    best_error = float('inf')
    
    # Try different lattice parameters
    for a in np.linspace(2, 20, 50):  # Adjust range based on expected lattice parameters
        for b in np.linspace(2, 20, 50):
            for c in np.linspace(2, 20, 50):
                if not CRYSTAL_SYSTEMS[crystal_system]['constraints'](a, b, c, 90, 90, 90):
                    continue
                
                # Calculate theoretical d-spacings
                theoretical_d = []
                for hkl in hkl_combinations:
                    d = calculate_theoretical_d(hkl, a, b, c, 90, 90, 90)
                    theoretical_d.append((hkl, d))
                
                # Match experimental and theoretical peaks
                matches = match_peaks(d_spacings, theoretical_d, tolerance)
                
                if len(matches) > len(indexed_peaks):
                    indexed_peaks = matches
                    best_fit = {'a': a, 'b': b, 'c': c, 'alpha': 90, 'beta': 90, 'gamma': 90}
                    best_error = calculate_fit_error(matches)
    
    return {
        'indexed_peaks': indexed_peaks,
        'lattice_parameters': best_fit,
        'error': best_error,
        'crystal_system': crystal_system
    }

def calculate_d_spacing(two_theta, wavelength):
    """Calculate d-spacing from 2θ using Bragg's law."""
    return wavelength / (2 * np.sin(np.radians(two_theta / 2)))

def calculate_theoretical_d(hkl, a, b, c, alpha, beta, gamma):
    """Calculate theoretical d-spacing for given hkl indices."""
    h, k, l = hkl
    alpha, beta, gamma = np.radians([alpha, beta, gamma])
    
    # Calculate reciprocal lattice parameters
    V = a * b * c * np.sqrt(1 - np.cos(alpha)**2 - np.cos(beta)**2 - np.cos(gamma)**2 + 
                           2 * np.cos(alpha) * np.cos(beta) * np.cos(gamma))
    
    a_star = b * c * np.sin(alpha) / V
    b_star = a * c * np.sin(beta) / V
    c_star = a * b * np.sin(gamma) / V
    
    alpha_star = np.arccos((np.cos(beta) * np.cos(gamma) - np.cos(alpha)) / 
                          (np.sin(beta) * np.sin(gamma)))
    beta_star = np.arccos((np.cos(alpha) * np.cos(gamma) - np.cos(beta)) / 
                         (np.sin(alpha) * np.sin(gamma)))
    gamma_star = np.arccos((np.cos(alpha) * np.cos(beta) - np.cos(gamma)) / 
                          (np.sin(alpha) * np.sin(beta)))
    
    # Calculate d-spacing
    d2 = (h**2 * a_star**2 + k**2 * b_star**2 + l**2 * c_star**2 +
          2 * h * k * a_star * b_star * np.cos(gamma_star) +
          2 * h * l * a_star * c_star * np.cos(beta_star) +
          2 * k * l * b_star * c_star * np.cos(alpha_star))
    
    return 1 / np.sqrt(d2)

def match_peaks(experimental_d, theoretical_d, tolerance):
    """Match experimental and theoretical peaks within tolerance."""
    matches = []
    used_theoretical = set()
    
    for exp_d in experimental_d:
        best_match = None
        best_diff = float('inf')
        
        for hkl, theo_d in theoretical_d:
            if hkl in used_theoretical:
                continue
                
            diff = abs(exp_d - theo_d)
            if diff < tolerance and diff < best_diff:
                best_match = (hkl, theo_d)
                best_diff = diff
        
        if best_match:
            matches.append({
                'hkl': best_match[0],
                'd_exp': exp_d,
                'd_theo': best_match[1],
                'error': best_diff
            })
            used_theoretical.add(best_match[0])
    
    return matches

def calculate_fit_error(matches):
    """Calculate the root mean square error of the fit."""
    if not matches:
        return float('inf')
    errors = [m['error'] for m in matches]
    return np.sqrt(np.mean(np.array(errors)**2))

def refine_lattice_parameters(indexed_peaks, crystal_system):
    """
    Refine lattice parameters using indexed peaks.
    
    Args:
        indexed_peaks (list): List of indexed peaks
        crystal_system (str): Crystal system
        
    Returns:
        dict: Refined lattice parameters
    """
    def objective(params):
        """Objective function for minimization."""
        if crystal_system == 'cubic':
            a = params[0]
            error = 0
            for peak in indexed_peaks:
                h, k, l = peak['h'], peak['k'], peak['l']
                d_calc = a / np.sqrt(h**2 + k**2 + l**2)
                error += (peak['d_spacing'] - d_calc)**2
        elif crystal_system == 'tetragonal':
            a, c = params
            error = 0
            for peak in indexed_peaks:
                h, k, l = peak['h'], peak['k'], peak['l']
                d_calc = 1 / np.sqrt((h**2 + k**2)/a**2 + l**2/c**2)
                error += (peak['d_spacing'] - d_calc)**2
        else:
            raise ValueError(f"Refinement not implemented for {crystal_system}")
        return error
    
    # Initial guess
    if crystal_system == 'cubic':
        x0 = [3.0]  # Typical lattice parameter
        bounds = [(0.1, 10.0)]
    elif crystal_system == 'tetragonal':
        x0 = [3.0, 3.0]
        bounds = [(0.1, 10.0), (0.1, 10.0)]
    
    # Minimize
    result = minimize(objective, x0, bounds=bounds)
    
    if crystal_system == 'cubic':
        return {'a': result.x[0]}
    elif crystal_system == 'tetragonal':
        return {'a': result.x[0], 'c': result.x[1]} 