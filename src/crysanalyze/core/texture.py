"""
Texture analysis module for XRD data
"""

import numpy as np
from pathlib import Path
from scipy.stats import norm

def analyze_texture(two_theta, intensity, reference_cif=None):
    """
    Perform basic texture analysis.
    
    Args:
        two_theta (np.ndarray): 2-theta values
        intensity (np.ndarray): Intensity values
        reference_cif (str, optional): Path to reference CIF file
        
    Returns:
        list: List of texture analysis results
    """
    results = []
    
    # If reference CIF is provided, try to parse it
    if reference_cif:
        try:
            reference_intensities = _parse_cif_intensities(reference_cif)
            results.extend(_compare_intensities(intensity, reference_intensities))
        except Exception as e:
            results.append(f"Warning: Could not parse reference CIF: {e}")
    
    # Basic texture analysis based on relative intensities
    results.extend(_analyze_relative_intensities(two_theta, intensity))
    
    return results

def _parse_cif_intensities(cif_path):
    """
    Parse intensities from a CIF file.
    
    Args:
        cif_path (str): Path to CIF file
        
    Returns:
        dict: Dictionary of (h,k,l): intensity
    """
    intensities = {}
    
    with open(cif_path, 'r') as f:
        for line in f:
            if line.startswith('_refln_F_squared_meas'):
                # Parse reflection data
                h = int(line.split()[1])
                k = int(line.split()[2])
                l = int(line.split()[3])
                intensity = float(line.split()[4])
                intensities[(h, k, l)] = intensity
    
    return intensities

def _compare_intensities(observed_intensities, reference_intensities):
    """
    Compare observed and reference intensities.
    
    Args:
        observed_intensities (np.ndarray): Observed intensities
        reference_intensities (dict): Reference intensities
        
    Returns:
        list: List of comparison results
    """
    results = []
    
    # Normalize observed intensities
    observed_max = np.max(observed_intensities)
    observed_norm = observed_intensities / observed_max
    
    # Normalize reference intensities
    reference_max = max(reference_intensities.values())
    reference_norm = {hkl: i/reference_max for hkl, i in reference_intensities.items()}
    
    # Compare relative intensities
    for hkl, ref_int in reference_norm.items():
        # Find closest observed peak
        obs_int = np.max(observed_norm)  # Simplified - should match peaks properly
        ratio = obs_int / ref_int
        
        if ratio > 1.5:
            results.append(f"Peak {hkl} shows intensity enhancement (ratio: {ratio:.2f})")
        elif ratio < 0.5:
            results.append(f"Peak {hkl} shows intensity reduction (ratio: {ratio:.2f})")
    
    return results

def _analyze_relative_intensities(two_theta, intensity):
    """
    Analyze relative intensities for texture.
    
    Args:
        two_theta (np.ndarray): 2-theta values
        intensity (np.ndarray): Intensity values
        
    Returns:
        list: List of analysis results
    """
    results = []
    
    # Normalize intensities
    intensity_norm = intensity / np.max(intensity)
    
    # Find major peaks
    peak_indices = np.where(intensity_norm > 0.5)[0]
    
    if len(peak_indices) > 0:
        # Calculate average intensity
        avg_intensity = np.mean(intensity_norm[peak_indices])
        
        # Check for significant deviations
        for idx in peak_indices:
            rel_intensity = intensity_norm[idx] / avg_intensity
            if rel_intensity > 1.5:
                results.append(
                    f"Peak at {two_theta[idx]:.2f}° shows significant intensity enhancement "
                    f"(relative intensity: {rel_intensity:.2f})"
                )
            elif rel_intensity < 0.5:
                results.append(
                    f"Peak at {two_theta[idx]:.2f}° shows significant intensity reduction "
                    f"(relative intensity: {rel_intensity:.2f})"
                )
    
    return results

def calculate_pole_figure(indexed_peaks, hkl, wavelength, num_points=100):
    """
    Calculate pole figure for a specific hkl reflection.
    
    Args:
        indexed_peaks (list): List of indexed peaks
        hkl (tuple): Miller indices (h,k,l)
        wavelength (float): X-ray wavelength in Angstroms
        num_points (int): Number of points in each direction
        
    Returns:
        dict: Pole figure data
    """
    # Generate grid of sample directions
    phi = np.linspace(0, 360, num_points)
    psi = np.linspace(0, 90, num_points)
    phi_grid, psi_grid = np.meshgrid(phi, psi)
    
    # Calculate pole figure intensity
    intensity = np.zeros_like(phi_grid)
    
    for peak in indexed_peaks:
        if peak['hkl'] == hkl:
            # Calculate pole figure coordinates
            two_theta = peak['position']
            theta = np.radians(two_theta / 2)
            
            # Calculate intensity distribution
            for i in range(num_points):
                for j in range(num_points):
                    # Calculate angle between sample direction and scattering vector
                    angle = np.arccos(np.sin(np.radians(psi_grid[i,j])) * 
                                    np.cos(np.radians(phi_grid[i,j] - peak['phi'])))
                    
                    # Add intensity contribution
                    intensity[i,j] += peak['intensity'] * norm.pdf(angle, 0, 0.1)
    
    return {
        'phi': phi_grid,
        'psi': psi_grid,
        'intensity': intensity,
        'hkl': hkl
    }

def calculate_texture_coefficients(pole_figure, max_l=4):
    """
    Calculate texture coefficients from pole figure.
    
    Args:
        pole_figure (dict): Pole figure data
        max_l (int): Maximum order of texture coefficients
        
    Returns:
        dict: Texture coefficients
    """
    phi = pole_figure['phi']
    psi = pole_figure['psi']
    intensity = pole_figure['intensity']
    
    # Normalize intensity
    intensity = intensity / np.max(intensity)
    
    # Calculate texture coefficients
    coefficients = {}
    for l in range(0, max_l + 1, 2):
        for m in range(-l, l + 1):
            # Calculate spherical harmonic
            Y = _spherical_harmonic(l, m, phi, psi)
            
            # Calculate coefficient
            coeff = np.sum(intensity * Y) / np.sum(Y**2)
            coefficients[f'C_{l}_{m}'] = coeff
    
    return coefficients

def _spherical_harmonic(l, m, phi, psi):
    """Calculate spherical harmonic function."""
    from scipy.special import sph_harm
    
    # Convert to spherical coordinates
    theta = np.radians(90 - psi)
    phi_rad = np.radians(phi)
    
    # Calculate spherical harmonic
    Y = sph_harm(m, l, phi_rad, theta)
    
    return Y.real

def calculate_orientation_distribution(pole_figures, resolution=5):
    """
    Calculate orientation distribution function (ODF) from multiple pole figures.
    
    Args:
        pole_figures (list): List of pole figure data
        resolution (float): Angular resolution in degrees
        
    Returns:
        dict: ODF data
    """
    # Generate Euler angle grid
    phi1 = np.arange(0, 360, resolution)
    Phi = np.arange(0, 90, resolution)
    phi2 = np.arange(0, 360, resolution)
    
    # Initialize ODF
    odf = np.zeros((len(phi1), len(Phi), len(phi2)))
    
    # Calculate ODF from pole figures
    for pf in pole_figures:
        # Calculate contribution to ODF
        for i, p1 in enumerate(phi1):
            for j, P in enumerate(Phi):
                for k, p2 in enumerate(phi2):
                    # Calculate corresponding pole figure coordinates
                    psi, phi = _euler_to_pole(p1, P, p2, pf['hkl'])
                    
                    # Add intensity contribution
                    odf[i,j,k] += _interpolate_pole_figure(pf, phi, psi)
    
    # Normalize ODF
    odf = odf / np.max(odf)
    
    return {
        'phi1': phi1,
        'Phi': Phi,
        'phi2': phi2,
        'odf': odf
    }

def _euler_to_pole(phi1, Phi, phi2, hkl):
    """Convert Euler angles to pole figure coordinates."""
    # Convert to radians
    phi1, Phi, phi2 = np.radians([phi1, Phi, phi2])
    
    # Calculate direction cosines
    h, k, l = hkl
    g = np.array([
        [np.cos(phi1)*np.cos(phi2) - np.sin(phi1)*np.sin(phi2)*np.cos(Phi),
         -np.cos(phi1)*np.sin(phi2) - np.sin(phi1)*np.cos(phi2)*np.cos(Phi),
         np.sin(phi1)*np.sin(Phi)],
        [np.sin(phi1)*np.cos(phi2) + np.cos(phi1)*np.sin(phi2)*np.cos(Phi),
         -np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2)*np.cos(Phi),
         -np.cos(phi1)*np.sin(Phi)],
        [np.sin(phi2)*np.sin(Phi),
         np.cos(phi2)*np.sin(Phi),
         np.cos(Phi)]
    ])
    
    # Calculate pole figure coordinates
    direction = g @ np.array([h, k, l])
    psi = np.degrees(np.arccos(direction[2]))
    phi = np.degrees(np.arctan2(direction[1], direction[0]))
    
    return psi, phi

def _interpolate_pole_figure(pf, phi, psi):
    """Interpolate pole figure intensity at given coordinates."""
    from scipy.interpolate import griddata
    
    # Ensure angles are in correct range
    phi = phi % 360
    psi = np.clip(psi, 0, 90)
    
    # Interpolate intensity
    points = np.column_stack((pf['phi'].ravel(), pf['psi'].ravel()))
    values = pf['intensity'].ravel()
    return griddata(points, values, (phi, psi), method='linear') 