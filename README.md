# crysanalyze

A command-line interface (CLI) tool for rapid preliminary analysis of powder X-ray Diffraction (XRD) data. This tool is designed for quick initial screening of synthesized materials and provides essential XRD analysis capabilities.

## Features

- Read common XRD data formats (.xy, .dat)
- Background subtraction using polynomial fitting
- Peak finding with customizable parameters
- Multiple visualization options
- 2θ range selection
- Command-line interface for automation and scripting

## Installation

### From Source

```bash
git clone https://github.com/raymsm/crysanalyze.git
cd crysanalyze
pip install -r requirements.txt
```

### Dependencies

- Python 3.6+
- NumPy
- SciPy
- Matplotlib

## Usage

The basic syntax for using crysanalyze is:

```bash
python crysanalyze.py --input <data_file> --wavelength <wavelength> <command> [options]
```

### Required Arguments

- `--input`: Path to the XRD data file (.xy, .dat)
- `--wavelength`: X-ray wavelength in Angstroms (e.g., 1.5406 for Cu Kα1)

### Global Optional Arguments

- `--output`: Path to save text-based results
- `--range MIN MAX`: Process only data within this 2θ range

### Commands

#### Peak Finding

```bash
python crysanalyze.py --input sample.xy --wavelength 1.5406 peaks [options]
```

Options:
- `--min-height`: Minimum peak height for detection
- `--min-dist`: Minimum horizontal distance between peaks (in degrees)
- `--bg-poly`: Order of polynomial for background subtraction
- `--fit-peaks`: Enable peak fitting (if implemented)

Example:
```bash
python crysanalyze.py --input sample.xy --wavelength 1.5406 peaks --bg-poly 3 --min-height 100 --min-dist 0.5
```

#### Plotting

```bash
python crysanalyze.py --input sample.xy --wavelength 1.5406 plot [options]
```

Options:
- `--plot-type`: Type of plot to generate
  - `raw`: Raw XRD data
  - `bgsub`: Background-subtracted data
  - `peaks`: Data with marked peaks
- `--bg-poly`: Order of polynomial for background subtraction
- `--save-plot`: Path to save the plot

Example:
```bash
python crysanalyze.py --input sample.xy --wavelength 1.5406 plot --plot-type peaks --save-plot peaks.png
```

## Examples

1. Find peaks with background subtraction:
```bash
python crysanalyze.py --input sample.xy --wavelength 1.5406 peaks --bg-poly 3 --min-height 100
```

2. Plot data with marked peaks:
```bash
python crysanalyze.py --input sample.xy --wavelength 1.5406 plot --plot-type peaks --save-plot peaks.png
```

3. Analyze specific 2θ range:
```bash
python crysanalyze.py --input sample.xy --wavelength 1.5406 --range 20 60 peaks
```

## Input Format

The input file should be a two-column text file containing:
- Column 1: 2θ values in degrees
- Column 2: Intensity values

## Output Format

### Peak Finding
```
Peak Analysis Results
2-theta  Intensity
30.700   1200.0
40.700   1500.0
```

### Plots
- Raw data plot
- Background-subtracted plot
- Peak-marked plot with labels

## Limitations

- Currently supports only two-column data files
- Basic polynomial background subtraction
- Simple peak finding algorithm
- Limited peak fitting capabilities

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This tool uses the following libraries:
- NumPy for numerical computations
- SciPy for peak finding and optimization
- Matplotlib for plotting

## Citation

If you use this tool in your research, please cite:

```
crysanalyze: A CLI Tool for Rapid XRD Analysis
https://github.com/yourusername/crysanalyze
```