"""Setup script for AAV Nanobody Display Pipeline."""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_path = Path(__file__).parent / "README.md"
long_description = readme_path.read_text() if readme_path.exists() else ""

setup(
    name="aav_nanobody_display",
    version="0.1.0",
    description="Visualization pipeline for AAV capsids with nanobody insertions",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Voyager Therapeutics",
    python_requires=">=3.8",
    packages=find_packages(),
    include_package_data=True,
    package_data={
        "aav_nanobody_display": ["config/*.yaml"],
    },
    install_requires=[
        "numpy>=1.21.0",
        "biopython>=1.79",
        "pyyaml>=6.0",
        "py3dmol>=2.0.0",
    ],
    extras_require={
        "esmfold": [
            "torch>=1.12.0",
            "fair-esm>=2.0.0",
        ],
        "dev": [
            "pytest>=7.0.0",
            "black>=22.0.0",
            "mypy>=0.990",
        ],
    },
    entry_points={
        "console_scripts": [
            "aav-viz=aav_nanobody_display.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
)
