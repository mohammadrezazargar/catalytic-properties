import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

def michaelis_menten(S, Vmax, Km):
    """Michaelis-Menten equation"""
    return (Vmax * S) / (Km + S)

def analyze_enzyme_kinetics(file_path, enzyme_concentration=20e-12):
    """
    Reads an Excel file, extracts kinetic data, and computes Km, Vmax, kcat, and catalytic efficiency.

    Parameters:
    - file_path (str): Path to the Excel file.
    - enzyme_concentration (float): Enzyme concentration in M (default: 20 pM = 20e-12 M).

    Returns:
    - dict: Contains Km, Vmax, kcat, and catalytic efficiency.
    """

    # Load the Excel file (assuming the first sheet contains the data)
    df = pd.read_excel(file_path, sheet_name=0, header=None)

    # Extract substrate concentrations from the first row (B to Q)
    substrate_concentrations = df.iloc[0, 1:-1].astype(float).values

    # Extract time points from the last column (starting from row 2)
    time_points = df.iloc[2:, -1].astype(float).values

    # Extract CO2 disappearance data (columns B to Q, starting from row 2)
    reaction_data = df.iloc[2:, 1:-1].astype(float).values

    # Calculate initial reaction velocities (v0) using linear regression on first few points
    num_points_for_slope = 5
    initial_rates = []

    for i in range(reaction_data.shape[1]):
        slope, _ = np.polyfit(time_points[:num_points_for_slope], reaction_data[:num_points_for_slope, i], 1)
        initial_rates.append(-slope)  # Negative because CO2 is disappearing

    initial_rates = np.array(initial_rates)

    # Fit the Michaelis-Menten equation
    popt, _ = curve_fit(michaelis_menten, substrate_concentrations, initial_rates, bounds=(0, np.inf))
    Vmax, Km = popt

    # Calculate kcat and catalytic efficiency
    kcat = Vmax / enzyme_concentration
    catalytic_efficiency = kcat / Km

    # Store results in a dictionary
    results = {
        "Km (µM)": Km,
        "Vmax (µM/s)": Vmax,
        "kcat (s^-1)": kcat,
        "Catalytic Efficiency (M^-1s^-1)": catalytic_efficiency
    }

    return results

# Example Usage:
file_path = "data.xlsx"  # Change this to your actual file path
kinetics_results = analyze_enzyme_kinetics(file_path)

# Print results
for key, value in kinetics_results.items():
    print(f"{key}: {value:.3e}")
