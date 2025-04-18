"""
Plotting module for XRD data visualization
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D

def plot_data(two_theta, intensity, plot_type='raw', save_path=None):
    """
    Generate XRD plots.
    
    Args:
        two_theta (np.ndarray): 2-theta values
        intensity (np.ndarray): Intensity values
        plot_type (str): Type of plot to generate
        save_path (str, optional): Path to save the plot
    """
    plt.figure(figsize=(10, 6))
    
    if plot_type == 'raw':
        _plot_raw_data(two_theta, intensity)
    elif plot_type == 'bgsub':
        _plot_background_subtracted(two_theta, intensity)
    elif plot_type == 'peaks':
        _plot_peaks(two_theta, intensity)
    elif plot_type == 'indexed':
        _plot_indexed(two_theta, intensity)
    elif plot_type == 'scherrer':
        _plot_scherrer(two_theta, intensity)
    elif plot_type == 'whplot':
        _plot_williamson_hall(two_theta, intensity)
    else:
        raise ValueError(f"Unknown plot type: {plot_type}")
    
    # Customize plot
    plt.xlabel('2θ (degrees)')
    plt.ylabel('Intensity (a.u.)')
    plt.grid(True, alpha=0.3)
    plt.gca().xaxis.set_minor_locator(AutoMinorLocator())
    plt.gca().yaxis.set_minor_locator(AutoMinorLocator())
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()

def _plot_raw_data(two_theta, intensity):
    """Plot raw XRD data."""
    plt.plot(two_theta, intensity, 'k-', linewidth=1)
    plt.title('Raw XRD Data')

def _plot_background_subtracted(two_theta, intensity):
    """Plot background-subtracted XRD data."""
    from .processing import subtract_background
    bg_sub = subtract_background(two_theta, intensity)
    plt.plot(two_theta, bg_sub, 'b-', linewidth=1)
    plt.title('Background-Subtracted XRD Data')

def _plot_peaks(two_theta, intensity):
    """Plot XRD data with peaks marked."""
    from .processing import find_peaks
    
    # Plot data
    plt.plot(two_theta, intensity, 'k-', linewidth=1)
    
    # Find and mark peaks
    peaks = find_peaks(two_theta, intensity)
    for peak in peaks:
        plt.axvline(x=peak['position'], color='r', linestyle='--', alpha=0.5)
        plt.text(
            peak['position'],
            peak['intensity'],
            f"{peak['position']:.2f}°",
            rotation=90,
            va='bottom'
        )
    
    plt.title('XRD Data with Peaks Marked')

def _plot_indexed(two_theta, intensity):
    """Plot XRD data with indexed peaks labeled."""
    from .processing import find_peaks
    from .indexing import index_peaks
    
    # Plot data
    plt.plot(two_theta, intensity, 'k-', linewidth=1)
    
    # Find and index peaks
    peaks = find_peaks(two_theta, intensity)
    indexed_peaks = index_peaks(peaks, wavelength=1.5406)  # Default Cu Kα
    
    # Mark and label indexed peaks
    for peak in indexed_peaks:
        plt.axvline(x=peak['position'], color='r', linestyle='--', alpha=0.5)
        plt.text(
            peak['position'],
            peak['intensity'],
            f"({peak['h']}{peak['k']}{peak['l']})",
            rotation=90,
            va='bottom'
        )
    
    plt.title('XRD Data with Indexed Peaks')

def _plot_scherrer(two_theta, intensity):
    """Plot Scherrer analysis data."""
    from .processing import find_peaks
    
    # Find peaks and get FWHM
    peaks = find_peaks(two_theta, intensity, fit_peaks=True)
    
    if not peaks or not any(p['fwhm'] for p in peaks):
        raise ValueError("No valid peak fits found for Scherrer analysis")
    
    # Calculate Scherrer parameters
    theta = np.array([p['position']/2 for p in peaks if p['fwhm']])
    beta = np.array([np.radians(p['fwhm']) for p in peaks if p['fwhm']])
    beta_cos_theta = beta * np.cos(np.radians(theta))
    
    # Plot
    plt.plot(theta, beta_cos_theta, 'bo')
    plt.xlabel('θ (degrees)')
    plt.ylabel('βcosθ')
    plt.title('Scherrer Analysis')

def _plot_williamson_hall(two_theta, intensity):
    """Plot Williamson-Hall analysis data."""
    from .processing import find_peaks
    
    # Find peaks and get FWHM
    peaks = find_peaks(two_theta, intensity, fit_peaks=True)
    
    if not peaks or not any(p['fwhm'] for p in peaks):
        raise ValueError("No valid peak fits found for Williamson-Hall analysis")
    
    # Calculate Williamson-Hall parameters
    theta = np.array([p['position']/2 for p in peaks if p['fwhm']])
    beta = np.array([np.radians(p['fwhm']) for p in peaks if p['fwhm']])
    beta_cos_theta = beta * np.cos(np.radians(theta))
    sin_theta = np.sin(np.radians(theta))
    
    # Plot
    plt.plot(4*sin_theta, beta_cos_theta, 'bo')
    plt.xlabel('4sinθ')
    plt.ylabel('βcosθ')
    plt.title('Williamson-Hall Plot')

def plot_xrd(data, title='XRD Pattern', save_path=None):
    """
    Plot XRD pattern.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        title (str): Plot title
        save_path (str, optional): Path to save the plot
    """
    plt.figure(figsize=(10, 6))
    plt.plot(data[:, 0], data[:, 1], 'k-', linewidth=1)
    plt.xlabel('2θ (degrees)')
    plt.ylabel('Intensity (a.u.)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()

def plot_peaks(data, peaks, title='XRD Pattern with Peaks', save_path=None):
    """
    Plot XRD pattern with marked peaks.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        peaks (numpy.ndarray): Array of (2θ, intensity) pairs for peaks
        title (str): Plot title
        save_path (str, optional): Path to save the plot
    """
    plt.figure(figsize=(10, 6))
    plt.plot(data[:, 0], data[:, 1], 'k-', linewidth=1)
    
    # Mark peaks
    for peak in peaks:
        plt.axvline(x=peak[0], color='r', linestyle='--', alpha=0.5)
        plt.text(peak[0], peak[1], f"{peak[0]:.2f}°", rotation=90, va='bottom')
    
    plt.xlabel('2θ (degrees)')
    plt.ylabel('Intensity (a.u.)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()

def plot_wh_plot(wh_data, title='Williamson-Hall Plot', save_path=None):
    """
    Plot Williamson-Hall analysis.
    
    Args:
        wh_data (dict): WH plot data from calculate_wh_plot
        title (str): Plot title
        save_path (str, optional): Path to save the plot
    """
    plt.figure(figsize=(8, 6))
    
    # Plot data points
    plt.scatter(wh_data['d_inv'], wh_data['beta_cos_theta'], 
                c='k', marker='o', label='Data')
    
    # Plot fit line
    x_fit = np.linspace(min(wh_data['d_inv']), max(wh_data['d_inv']), 100)
    y_fit = wh_data['fit_params'][0] * x_fit + wh_data['fit_params'][1]
    plt.plot(x_fit, y_fit, 'r-', label='Fit')
    
    plt.xlabel('1/d (Å⁻¹)')
    plt.ylabel('βcosθ')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    # Add strain and size information
    text = f"Strain: {wh_data['strain']:.4f}\n"
    text += f"Crystallite Size: {wh_data['crystallite_size']:.1f} nm"
    plt.text(0.05, 0.95, text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()

def plot_pole_figure(pole_figure, title='Pole Figure', save_path=None):
    """
    Plot pole figure.
    
    Args:
        pole_figure (dict): Pole figure data
        title (str): Plot title
        save_path (str, optional): Path to save the plot
    """
    plt.figure(figsize=(8, 8))
    
    # Create polar plot
    ax = plt.subplot(111, projection='polar')
    
    # Plot intensity
    phi = pole_figure['phi']
    psi = pole_figure['psi']
    intensity = pole_figure['intensity']
    
    # Create custom colormap
    cmap = LinearSegmentedColormap.from_list('custom', ['white', 'blue', 'red'])
    
    # Plot contour
    c = ax.contourf(np.radians(phi), psi, intensity, cmap=cmap)
    plt.colorbar(c, label='Intensity')
    
    # Add hkl information
    h, k, l = pole_figure['hkl']
    plt.title(f"{title}\n({h}{k}{l})")
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()

def plot_odf(odf_data, title='Orientation Distribution Function', save_path=None):
    """
    Plot orientation distribution function.
    
    Args:
        odf_data (dict): ODF data
        title (str): Plot title
        save_path (str, optional): Path to save the plot
    """
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create meshgrid
    phi1, Phi, phi2 = np.meshgrid(odf_data['phi1'], odf_data['Phi'], odf_data['phi2'])
    
    # Plot isosurface
    ax.scatter(phi1.ravel(), Phi.ravel(), phi2.ravel(),
               c=odf_data['odf'].ravel(), cmap='viridis')
    
    ax.set_xlabel('φ₁')
    ax.set_ylabel('Φ')
    ax.set_zlabel('φ₂')
    plt.title(title)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close()

def plot_crystallinity(data, crystallinity, title='Crystallinity Analysis', save_path=None):
    """
    Plot crystallinity analysis.
    
    Args:
        data (numpy.ndarray): Array of (2θ, intensity) pairs
        crystallinity (float): Crystallinity percentage
        title (str): Plot title
        save_path (str, optional): Path to save the plot
    """
    plt.figure(figsize=(10, 6))
    
    # Plot XRD pattern
    plt.plot(data[:, 0], data[:, 1], 'k-', linewidth=1, label='Raw Data')
    
    # Add crystallinity information
    text = f"Crystallinity: {crystallinity:.1f}%"
    plt.text(0.05, 0.95, text, transform=plt.gca().transAxes,
             verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.xlabel('2θ (degrees)')
    plt.ylabel('Intensity (a.u.)')
    plt.title(title)
    plt.grid(True, alpha=0.3)
    plt.legend()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        plt.show()
    plt.close() 