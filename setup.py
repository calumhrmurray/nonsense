#!/usr/bin/env python3
"""
Setup script for AbacusSummit Intrinsic Alignment Analysis
"""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="abacus-ia-analysis",
    version="1.0.0",
    author="LRG HOD Pipeline",
    description="LRG HOD Pipeline for AbacusSummit with Intrinsic Alignment Analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "run-lrg-hod=src.test_lrg_hod:main",
            "compute-tidal-shear=src.efficient_tidal_shear:test_tidal_shear",
            "validate-catalog=src.validate_lrg_catalog:main",
        ],
    },
)