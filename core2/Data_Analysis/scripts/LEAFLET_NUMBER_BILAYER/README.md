# 🔬 Cylindrical Membrane Lipid Analysis

A Python tool for **classifying and tracking inner/outer leaflet lipids** in cylindrical membrane systems from molecular dynamics simulations. This project analyzes lipid distribution across membrane leaflets and tracks their dynamics over simulation time.

![Python](https://img.shields.io/badge/Python-3.7+-blue.svg)
![MDAnalysis](https://img.shields.io/badge/MDAnalysis-2.0+-orange.svg)
![GROMACS](https://img.shields.io/badge/GROMACS-Compatible-green.svg)
![Analysis](https://img.shields.io/badge/Analysis-Leaflet%20Classification-green.svg)

---

## ✨ Features

- ✅ **Automatic leaflet classification** - Identifies inner vs. outer leaflet lipids based on radial position
- ✅ **Segmented cylindrical analysis** - Divides membrane tube into customizable segments along the axis
- ✅ **PO4 headgroup tracking** - Uses phosphate groups for accurate lipid positioning
- ✅ **3D & 2D visualization** - Interactive 3D scatter plots and 2D top-down views of lipid distribution
- ✅ **Frame-by-frame dynamics** - Tracks lipid counts and redistribution over simulation time
- ✅ **Radius calculation** - Computes membrane radius from tailgroup (C3B) median distances
- ✅ **Statistical time series** - Generates plots and data files for temporal analysis

---

## 🧪 Scientific Applications

- **Leaflet composition analysis** - Quantify asymmetric distribution in cylindrical membranes
- **Lipid flipping studies** - Track transitions between inner and outer leaflets
- **Membrane curvature effects** - Analyze how curvature influences lipid packing
- **Nanotube characterization** - Study lipid organization in membrane nanotubes
- **Dynamic redistribution** - Monitor lipid movement during simulation

---

## 📦 Input Requirements

To use this code, you need:
- **Trajectory files**: `.xtc` format from GROMACS simulations
- **Topology files**: `.tpr` system topology
- **Lipid selection**: DLPC lipids with PO4 headgroups and C3B tailgroups

---

## 🔧 Core Components

### 1. Lipid Classification (`count_lipids_cylindrical_x`)
- Calculates center of mass for each cylindrical segment
- Determines membrane radius from tailgroup median distances
- Classifies lipids as inner (≤ radius) or outer (> radius) leaflet
- Generates 3D and 2D visualization plots

### 2. Multi-frame Analysis (`analyze_first_10_frames`)
- Processes trajectory frames with progress tracking
- Divides cylinder into configurable segments
- Computes statistics and averages across segments
- Exports comprehensive data tables and trend plots

---

## ⚙️ Configuration Parameters

```python
lipid_resnames = ["DLPC"]     # Lipid types to analyze
segment_width = 100           # Segment size in Angstroms (default: 100Å)
headgroup_bead = "PO4"        # Headgroup atom for positioning
tailgroup_bead = "C3B"        # Tailgroup atom for radius calculation
