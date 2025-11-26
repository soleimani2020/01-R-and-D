import numpy as np
import matplotlib.pyplot as plt
import os


skip=7

def plot_radius_vs_position(num_bins=50, prefix="average_radius_vs_position_bin_", suffix=".txt"):
    plt.figure(figsize=(10, 6))
    
    # Loop from 3 up to num_bins-3
    for i in range(skip, num_bins-skip):
        filename = f"{prefix}{i}{suffix}"
        if not os.path.exists(filename):
            continue
        try:
            data = np.loadtxt(filename, skiprows=1)  # skip header row
            x, y = data[:,0], data[:,1]
            plt.plot(x, y, linestyle=":",  alpha=0.9)
        except Exception as e:
            print(f"Could not read {filename}: {e}")
    
    plt.xlabel("z / Å")
    plt.ylabel("R(z) / Å")
    plt.ylim(30, 60)  # set y-axis limits
    plt.legend(fontsize=8, ncol=2, bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("Constrict.png")
    plt.show()

# Example usage
plot_radius_vs_position()

