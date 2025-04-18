from setuptools import setup, find_packages

setup(
    name="crysanalyze",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "click",
    ],
    entry_points={
        "console_scripts": [
            "crysanalyze=crysanalyze.cli:cli",
        ],
    },
) 