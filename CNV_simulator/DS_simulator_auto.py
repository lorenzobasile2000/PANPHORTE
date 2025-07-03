import pandas as pd
import random
import numpy as np

# Function to generate uniformly distributed start positions
def generate_positions_optimized(total_length, num_positions, min_distance):
    positions = np.linspace(1, total_length, num_positions, dtype=int)
    jitter = np.random.randint(-min_distance // 2, min_distance // 2, size=num_positions)
    positions += jitter
    positions = np.clip(positions, 1, total_length - 1) 
    return sorted(set(positions)) 

# Function to generate CNV data and save to CSV
def generate_cnv_csv(file_name, total_length, num_cnvs=10000):
    # Define proportions for the categories and corresponding factors (F)
    num_cnvs_small = int(num_cnvs * 0.40)  # Length < 100, F = 100
    num_cnvs_medium = int(num_cnvs * 0.30)  # Length 100-1000, F = 100
    num_cnvs_large = int(num_cnvs * 0.20)  # Length 1000-10000, F = 10
    num_cnvs_very_large = int(num_cnvs * 0.10)  # Length 10000-100000, F = 10

    # Minimum distance for uniform distribution
    min_distance = total_length // num_cnvs

    # Generate uniformly distributed start positions
    positions = generate_positions_optimized(total_length, num_cnvs, min_distance)

    # Prepare the dataset
    data_with_factors = []

    # Assign CNVs with corresponding lengths and factors
    for i, start_position in enumerate(positions):
        if i < num_cnvs_small:
            length = random.randint(1, 99)
            F = 100
        elif i < num_cnvs_small + num_cnvs_medium:
            length = random.randint(100, 999)
            F = 100
        elif i < num_cnvs_small + num_cnvs_medium + num_cnvs_large:
            length = random.randint(1000, 9999)
            F = 10
        else:
            length = random.randint(10000, 99999)
            F = 10

        scaled_length = length * F
        data_with_factors.append([start_position, length, scaled_length])

    # Create a DataFrame
    cnv_with_factors_df = pd.DataFrame(data_with_factors, columns=["R", "x", "y"])

    # Save to CSV
    csv_with_factors_path = f"t2t/csv/{file_name}.csv"
    cnv_with_factors_df.to_csv(csv_with_factors_path, index=False)
    print(f"Dataset saved to {csv_with_factors_path}")

# Read the input file and process each line
input_file = "t2t/lengths.txt"
with open(input_file, "r") as f:
    for line in f:
        name, length = line.strip().split(": ")
        total_length = int(length)
        generate_cnv_csv(name, total_length)