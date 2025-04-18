"""
Input/Output module for handling XRD data files
"""

import numpy as np
from pathlib import Path

def read_xy_data(filename):
    """
    Read XRD data from a file.
    
    Args:
        filename (str): Path to the data file
        
    Returns:
        numpy.ndarray: Array of (2θ, intensity) pairs
    """
    filepath = Path(filename)
    if not filepath.exists():
        raise FileNotFoundError(f"File not found: {filename}")
        
    # Try to determine file format from extension
    if filepath.suffix.lower() in ['.xy', '.dat', '.txt']:
        # Simple XY format
        data = np.loadtxt(filename)
    elif filepath.suffix.lower() == '.csv':
        # CSV format
        data = np.loadtxt(filename, delimiter=',')
    else:
        raise ValueError(f"Unsupported file format: {filepath.suffix}")
        
    if data.shape[1] != 2:
        raise ValueError("Data must contain exactly two columns (2θ and intensity)")
        
    return data

def write_xy_data(filename, data):
    """
    Write XRD data to a file.
    
    Args:
        filename (str): Path to save the data
        data (numpy.ndarray): Array of (2θ, intensity) pairs
    """
    filepath = Path(filename)
    
    # Determine format from extension
    if filepath.suffix.lower() in ['.xy', '.dat', '.txt']:
        # Simple XY format
        np.savetxt(filename, data, fmt='%.6f')
    elif filepath.suffix.lower() == '.csv':
        # CSV format
        np.savetxt(filename, data, delimiter=',', fmt='%.6f')
    else:
        raise ValueError(f"Unsupported file format: {filepath.suffix}")

def read_xrd_data(file_path):
    """
    Read XRD data from various file formats.
    
    Args:
        file_path (str): Path to the XRD data file
        
    Returns:
        tuple: (two_theta, intensity) as numpy arrays
        
    Raises:
        ValueError: If file format is not supported or data cannot be parsed
        FileNotFoundError: If file does not exist
    """
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")
    
    # Try to determine file format from extension
    if file_path.suffix.lower() in ['.xy', '.dat']:
        return _read_xy_format(file_path)
    else:
        # Try to read as space/tab/comma delimited
        return _read_delimited_format(file_path)

def _read_xy_format(file_path):
    """Read .xy or .dat format files."""
    try:
        data = np.loadtxt(file_path)
        if data.shape[1] != 2:
            raise ValueError("File must contain exactly 2 columns")
        return data[:, 0], data[:, 1]
    except Exception as e:
        raise ValueError(f"Error reading .xy/.dat file: {e}")

def _read_delimited_format(file_path):
    """Read space/tab/comma delimited files."""
    try:
        # Try different delimiters
        for delimiter in [' ', '\t', ',']:
            try:
                data = np.loadtxt(file_path, delimiter=delimiter)
                if data.shape[1] == 2:
                    return data[:, 0], data[:, 1]
            except:
                continue
        raise ValueError("Could not parse file with any standard delimiter")
    except Exception as e:
        raise ValueError(f"Error reading delimited file: {e}")

def write_results(file_path, results, header=None):
    """
    Write analysis results to a file.
    
    Args:
        file_path (str): Path to save results
        results (list): List of results to write
        header (str, optional): Header text to write at the top of the file
    """
    with open(file_path, 'w') as f:
        if header:
            f.write(header + '\n')
        for result in results:
            f.write(str(result) + '\n')

def read_indexed_peaks(file_path):
    """
    Read a file containing indexed peaks.
    
    Args:
        file_path (str): Path to the indexed peaks file
        
    Returns:
        list: List of dictionaries containing peak information
    """
    peaks = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith('#'):
                try:
                    h, k, l, two_theta, d_spacing = map(float, line.strip().split())
                    peaks.append({
                        'h': int(h),
                        'k': int(k),
                        'l': int(l),
                        'position': two_theta,
                        'd_spacing': d_spacing
                    })
                except ValueError:
                    continue
    return peaks 