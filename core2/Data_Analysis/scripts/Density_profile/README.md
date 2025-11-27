# 🧪 3D Lipid Membrane Density Analysis

A Python toolkit for **3D grid-based analysis of lipid membrane systems** from molecular dynamics simulations. This project identifies constricted regions in membrane structures and quantifies lipid density distributions in three-dimensional space using parallel processing.

![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)
![MDAnalysis](https://img.shields.io/badge/MDAnalysis-2.0+-orange.svg)
![GROMACS](https://img.shields.io/badge/GROMACS-Compatible-green.svg)
![License](https://img.shields.io/badge/License-MIT-yellow.svg)

---

## ✨ Features

- ✅ **Automatic constriction detection** - Finds the most narrow regions in cylindrical membrane systems
- ✅ **3D volumetric grid analysis** - Divides membrane volume into customizable 3D grid cells
- ✅ **Lipid-specific analysis** - Focuses on NC3 headgroups to track lipid distributions
- ✅ **Parallel processing** - Analyzes multiple trajectory frames simultaneously for efficiency
- ✅ **VTK format output** - Generates 3D density data compatible with visualization software
- ✅ **Comprehensive CSV logging** - Records grid occupancy data for statistical analysis
- ✅ **GROMACS compatibility** - Works with standard .gro and .xtc trajectory files

---

## 🎯 Scientific Applications

- **Nanopore characterization** - Analyze lipid packing around membrane constrictions
- **Membrane deformation studies** - Quantify 3D density variations in deformed bilayers
- **Lipid redistribution analysis** - Track how lipids reorganize around structural features
- **Volumetric density mapping** - Create 3D density profiles of membrane systems

---

## 📦 Installation

### Prerequisites
```bash
python >= 3.7
MDAnalysis >= 2.0.0
numpy >= 1.21.0
pandas >= 1.3.0
matplotlib >= 3.5.0

# Configure analysis parameters
Max_Mesh_Num = 300    # 3D grid resolution
NUMBER_FRAMES = 64    # Number of frames to analyze

# Run the analysis
python DP_3D.py
