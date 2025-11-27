import pandas as pd
import numpy as np

mesh_resolution = 200
num_z_cells = 201


# Read the CSV file
df = pd.read_csv("grid_summary_all.csv")

# Convert grid_key from string to tuple if it's stored as string
if isinstance(df['grid_key'].iloc[0], str):
    df['grid_key'] = df['grid_key'].apply(lambda x: eval(x) if isinstance(x, str) else x)

# Create a dictionary to store num_atoms for each grid cell
cell_atoms = {}

# Group by grid_key and collect all num_atoms values
for grid_key, group in df.groupby('grid_key'):
    num_atoms_list = group['num_atoms'].tolist()
    cell_atoms[grid_key] = num_atoms_list

# Calculate averages for each grid cell
cell_averages = {}
for grid_key, atoms_list in cell_atoms.items():
    average_atoms = np.mean(atoms_list)
    cell_averages[grid_key] = average_atoms
    print(f"Grid cell {grid_key}: Average atoms = {average_atoms:.2f}")

# If you want to create a 3D array of averages
average_matrix = np.zeros((mesh_resolution, mesh_resolution, num_z_cells))

for (i, j, k), avg in cell_averages.items():
    average_matrix[i, j, k] = avg

print(f"Overall average across all cells: {np.mean(list(cell_averages.values())):.2f}")













print("End of the code . good luck")




import pandas as pd
import numpy as np


# Read the CSV file
df = pd.read_csv("grid_summary_all.csv")

# Convert grid_key from string to tuple if it's stored as string
if isinstance(df['grid_key'].iloc[0], str):
    df['grid_key'] = df['grid_key'].apply(lambda x: eval(x) if isinstance(x, str) else x)

# Create a dictionary to store num_atoms for each grid cell
cell_atoms = {}

# Group by grid_key and collect all num_atoms values
for grid_key, group in df.groupby('grid_key'):
    num_atoms_list = group['num_atoms'].tolist()
    cell_atoms[grid_key] = num_atoms_list

# Calculate averages for each grid cell
cell_averages = {}
for grid_key, atoms_list in cell_atoms.items():
    average_atoms = np.mean(atoms_list)
    cell_averages[grid_key] = average_atoms
    print(f"Grid cell {grid_key}: Average atoms = {average_atoms:.2f}")

# If you want to create a 3D array of averages
average_matrix = np.zeros((mesh_resolution, mesh_resolution, num_z_cells))

for (i, j, k), avg in cell_averages.items():
    average_matrix[i, j, k] = avg

print(f"Overall average across all cells: {np.mean(list(cell_averages.values())):.2f}")





# Read the CSV file
df = pd.read_csv("grid_summary_all.csv")

# Convert grid_key to tuple if needed
if isinstance(df['grid_key'].iloc[0], str):
    df['grid_key'] = df['grid_key'].apply(lambda x: eval(x))

# Calculate average for each grid cell
cell_averages = df.groupby('grid_key')['num_atoms'].mean().to_dict()

# Write averages to file in the specified order
with open("normalized_values.txt", "w") as f:
    for k in range(num_z_cells):          # z loop
        for j in range(mesh_resolution):  # y loop
            for i in range(mesh_resolution):  # x loop
                grid_key = (i, j, k)
                
                # Get the average value for this grid cell
                avg_value = cell_averages.get(grid_key, 0.0)
                
                # Write to file
                f.write(f"{avg_value}\n")

print(f"Total values written: {mesh_resolution * mesh_resolution * num_z_cells}")
