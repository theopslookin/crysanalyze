"""
Physical constants and formulas for XRD analysis.
"""

import numpy as np

# Physical constants
PLANCK_CONSTANT = 6.62607015e-34  # J·s
SPEED_OF_LIGHT = 299792458  # m/s
ELECTRON_CHARGE = 1.602176634e-19  # C
ELECTRON_MASS = 9.1093837015e-31  # kg
AVOGADRO_NUMBER = 6.02214076e23  # mol^-1

# Common X-ray wavelengths (in Angstroms)
WAVELENGTHS = {
    'Cu-Kα1': 1.5406,
    'Cu-Kα2': 1.5444,
    'Cu-Kβ': 1.3922,
    'Mo-Kα1': 0.7093,
    'Mo-Kα2': 0.7136,
    'Co-Kα1': 1.7890,
    'Co-Kα2': 1.7929,
    'Fe-Kα1': 1.9360,
    'Fe-Kα2': 1.9400,
    'Cr-Kα1': 2.2897,
    'Cr-Kα2': 2.2936
}

def calculate_d_spacing(two_theta, wavelength):
    """
    Calculate d-spacing from 2-theta using Bragg's law.
    
    Args:
        two_theta (float): 2-theta angle in degrees
        wavelength (float): X-ray wavelength in Angstroms
        
    Returns:
        float: d-spacing in Angstroms
    """
    return wavelength / (2 * np.sin(np.radians(two_theta / 2)))

def calculate_two_theta(d_spacing, wavelength):
    """
    Calculate 2-theta from d-spacing using Bragg's law.
    
    Args:
        d_spacing (float): d-spacing in Angstroms
        wavelength (float): X-ray wavelength in Angstroms
        
    Returns:
        float: 2-theta angle in degrees
    """
    return 2 * np.degrees(np.arcsin(wavelength / (2 * d_spacing)))

def calculate_1_over_d_squared(h, k, l, crystal_system, a=None, b=None, c=None):
    """
    Calculate 1/d^2 for a given reflection and crystal system.
    
    Args:
        h, k, l (int): Miller indices
        crystal_system (str): Crystal system
        a, b, c (float, optional): Lattice parameters in Angstroms
        
    Returns:
        float: 1/d^2 in Angstroms^-2
    """
    if crystal_system == 'cubic':
        if a is None:
            raise ValueError("Lattice parameter 'a' required for cubic system")
        return (h**2 + k**2 + l**2) / a**2
    
    elif crystal_system == 'tetragonal':
        if a is None or c is None:
            raise ValueError("Lattice parameters 'a' and 'c' required for tetragonal system")
        return (h**2 + k**2)/a**2 + l**2/c**2
    
    elif crystal_system == 'hexagonal':
        if a is None or c is None:
            raise ValueError("Lattice parameters 'a' and 'c' required for hexagonal system")
        return 4/3 * (h**2 + h*k + k**2)/a**2 + l**2/c**2
    
    elif crystal_system == 'orthorhombic':
        if a is None or b is None or c is None:
            raise ValueError("Lattice parameters 'a', 'b', and 'c' required for orthorhombic system")
        return h**2/a**2 + k**2/b**2 + l**2/c**2
    
    else:
        raise ValueError(f"Unsupported crystal system: {crystal_system}")

def calculate_unit_cell_volume(crystal_system, a=None, b=None, c=None):
    """
    Calculate unit cell volume.
    
    Args:
        crystal_system (str): Crystal system
        a, b, c (float, optional): Lattice parameters in Angstroms
        
    Returns:
        float: Unit cell volume in cubic Angstroms
    """
    if crystal_system == 'cubic':
        if a is None:
            raise ValueError("Lattice parameter 'a' required for cubic system")
        return a**3
    
    elif crystal_system == 'tetragonal':
        if a is None or c is None:
            raise ValueError("Lattice parameters 'a' and 'c' required for tetragonal system")
        return a**2 * c
    
    elif crystal_system == 'hexagonal':
        if a is None or c is None:
            raise ValueError("Lattice parameters 'a' and 'c' required for hexagonal system")
        return np.sqrt(3)/2 * a**2 * c
    
    elif crystal_system == 'orthorhombic':
        if a is None or b is None or c is None:
            raise ValueError("Lattice parameters 'a', 'b', and 'c' required for orthorhombic system")
        return a * b * c
    
    else:
        raise ValueError(f"Unsupported crystal system: {crystal_system}")

def calculate_scherrer_crystallite_size(beta, wavelength, theta):
    """
    Calculate crystallite size using Scherrer equation.
    
    Args:
        beta (float): FWHM in radians
        wavelength (float): X-ray wavelength in Angstroms
        theta (float): Bragg angle in degrees
        
    Returns:
        float: Crystallite size in Angstroms
    """
    K = 0.9  # Scherrer constant
    return K * wavelength / (beta * np.cos(np.radians(theta))) 