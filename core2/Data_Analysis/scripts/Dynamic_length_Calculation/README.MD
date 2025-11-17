
# 🧪 Tube Membrane Radius Analysis 

A Python script to **dynamically calculate the radius of a cylindrical membrane** from molecular dynamics (MD) simulations using **MDAnalysis**.  
This project analyzes per-segment membrane radii, computes statistical properties, and visualizes Gaussian distributions and mean radii across segments.

---

## ✨ Features
- ✅ To use this code, you need the .tpr and .xtc files from your molecular dynamics trajectory.
- ✅ Splits cylinder into **dynamic segments** along the z-axis. You can choose number of segments according to your need  
- ✅ Computes **center-of-mass** and **mean radial distance** for each segment  
- ✅ Computes **radius** by measuring the average distance between enter-of-mass and lipid C3B particles
- ✅ Handles frames with missing atoms gracefully by assigning `NaN`  
- ✅ Saves **per-segment mean radii**, **averages**, and **standard deviations** to CSV/TXT  
- ✅ Generates **histograms with Gaussian fits** for each segment  
- ✅ Plots **mean radius vs. segment index** with error bars  
- ✅ Fully configurable: number of segments, atom selection, and number of frames analyzed  

---

## 🛠️ Dependencies

- [MDAnalysis](https://www.mdanalysis.org/)  
- `numpy`  
- `pandas`  
- `matplotlib`  
- `scipy`  

Install via pip:

```bash
pip install MDAnalysis numpy pandas matplotlib scipy


python3 Dynamic_Length_Calculation1.py --traj traj_comp.xtc --tpr nvt_auto.tpr --frames 100 --num_segments_global 37
