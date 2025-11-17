# 📊 FFT & Autocorrelation Analysis for tube membranes in molecular dynamics

A Python script to **analyze fluctuation modes and relaxation times of a tube membrane in Fourier space** from molecular dynamics (MD) simulation data.  

---

## ✨ Features

- ✅ To use this code, you must first generate the file "Radius_mean_segment.txt" by running Dynamic_Length_Calculation.py on your MD trajectory (*.tpr and *.xtc required).
- ✅ "Radius_mean_segment.txt" shows how the radius fluctuates over time for each tube membrane configuration. Each row in this file represents the average radius of a segment of a certain size (R(z)), and each column corresponds to the time evolution of the simulation.
- ✅ FFT3_Magnitude.py extracts each row from "Radius_mean_segment.txt" and sees that as **R(z)** and applies a Fast Fourier transform to analyze its frequency components. The resulting frequency magnitudes are then averaged to obtain the ensemble average of fluctuation modes.
- 
- ## ✨ Output
  
- ✅ Computes **mode fluctuation amplitudes** as a function of wave vector \(q\)  
- ✅ Calculates **autocorrelation function** for each mode  
- ✅ Fits autocorrelation to **exponential decay** to extract **relaxation times**  
- ✅ Saves results to CSV/TXT files (`fft_magnitudes_all.csv`, `relaxation_times.txt`)

- ## ✨ plots
  - **Fluctuation spectrum as a function of wave vector**  
  - **Autocorrelation** of each mode and the corresponding relation time. 
  - **Fluctuation amplitude as a function of time** for each mode
  - **Gaussian distribution of each mode** 

---

## 🛠️ Dependencies

- `numpy`  
- `pandas`  
- `matplotlib`  
- `scipy`  
- [MDAnalysis](https://www.mdanalysis.org/)  

Install via pip:

```bash
pip install numpy pandas matplotlib scipy MDAnalysis


Flow:






