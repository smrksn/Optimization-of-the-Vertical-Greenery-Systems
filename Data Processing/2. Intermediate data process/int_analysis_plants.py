import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from collections import defaultdict

################################################################################

# Directories

filename = '1. Preliminary data process\Output\plants.xlsx'
sheet_name = 'Impedance Ratio'

figures_folderpath = r'2. Intermediate data process\Figures'
if not os.path.exists(figures_folderpath):
    os.makedirs(figures_folderpath)

output_folderpath = r'2. Intermediate data process\Outputs'
if not os.path.exists(output_folderpath):
    os.makedirs(output_folderpath)

################################################################################

# Iteration conditions

# 1: 'Ajuga',
# 2: 'Bergenia',
# 3: 'Grass',
# 4: 'Heuchera',
# 5: 'Waldestania'
plant_numbers = [1,2,3,4,5]

# Cavity combinations (modified two-cavity method)
combinations = [[0, 0.01], [0, 0.02], [0, 0.03], [0, 0.04], [0, 0.05]]

################################################################################

# Other variables

z_0 = 412 # rayl
c = 343 # m/s, speed of sound in the air

################################################################################

# MAIN functions

def select(plant_num, L1, L2):
    plants = {
        1: {
            'name': 'Ajuga',
            'columns': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
            'cavities': [0, 0.01, 0.02, 0.03, 0.04, 0.05], #meters
            'thickness': [0.08] #meters
        },
        2: {
            'name': 'Bergenia',
            'columns': [13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
            'cavities': [0, 0.01, 0.02, 0.03, 0.04, 0.05], #meters
            'thickness': [0.05] #meters
        },
        3: {
            'name': 'Grass',
            'columns': [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36],
            'cavities': [0, 0.01, 0.02, 0.03, 0.04, 0.05], #meters
            'thickness': [0.06] #meters
        },
        4: {
            'name': 'Heuchera',
            'columns': [37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48],
            'cavities': [0, 0.01, 0.02, 0.03, 0.04, 0.05], #meters
            'thickness': [0.15] #meters
        },
        5: {
            'name': 'Waldestania',
            'columns': [49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60],
            'cavities': [0, 0.01, 0.02, 0.03, 0.04, 0.05], #meters
            'thickness': [0.08] #meters
        }
    }
    # Accessing the plant
    plant = plants.get(plant_num)
    # Retrieving the corresponding columns
    columns = []
    if L1 in plant['cavities'] and L2 in plant['cavities']:
        index1 = plant['cavities'].index(L1)
        index2 = plant['cavities'].index(L2) 
        # Getting the corresponding columns for each cavity depth
        columns1 = plant['columns'][index1 * 2: index1 * 2 + 2]
        columns2 = plant['columns'][index2 * 2: index2 * 2 + 2]
        # Combining the columns
        columns = columns1 + columns2
    # Retrieving the corresponding sample thickness
    thickness = plant['thickness'][0]
    name = plant['name']
    return name, columns, thickness



def load_excel(filename, sheet_name, columns, z_0):
    df = pd.read_excel(filename, sheet_name=sheet_name, usecols=[0]+columns, skiprows=1)
    frequencies = df.iloc[:, 0].values
    z_s1_values = (df.iloc[:, 1].values + 1j * df.iloc[:, 2].values) * z_0
    z_s2_values = (df.iloc[:, 3].values + 1j * df.iloc[:, 4].values) * z_0
    return frequencies, z_s1_values, z_s2_values         


def check(frequencies, L1, L2, c):
    # According to the Utsuno 1989, the appropriate set of air space depths must be selected so as to not satisfy:
    # Equation: f(L-L') = nc/2
    # Storing all depth differences (L-L')
    cavities = [L1, L2] # m
    L_diffs = []
    for i in range(len(cavities)):
        for j in range(i + 1, len(cavities)):
            diff = round(np.abs(cavities[i] - cavities[j]), 2) # Round to 2 decimal points
            # Check if diff already exists in L_diffs
            if diff not in L_diffs:
                L_diffs.append(diff)
    n_list = np.arange(0, 10) # List of integers from 0 to 10
    for f in frequencies:
        for L in L_diffs:
            for n in n_list:
                right = f * L
                left = n * c/2
                if abs(right - left) <= 0.5: # Tolerance (chosen so to meet the outcome of Utsuno)
                    print(f"The cavity depths difference {L} m for frequency {f} Hz fails to meet the condition.")


def calc_z_l2(z_0,L2, frequencies, c):
    z_l2_values = []
    for f in frequencies:
        k = 2*np.pi*f/c
        z_l2 = -1j*z_0*(1/np.tan(k*L2))
        z_l2_values.append(z_l2)
    return z_l2_values


def calc_z_c(z_s1_values, z_s2_values, z_l2_values):
    z_c_values = []
    for z_s1, z_s2, z_l2 in zip(z_s1_values, z_s2_values, z_l2_values):
        z_c = np.sqrt(z_l2*(z_s2-z_s1)+(z_s2*z_s1))
        z_c_real = abs(z_c.real) 
        z_c_reconstructed = complex(z_c_real, z_c.imag)
        z_c_values.append(z_c_reconstructed)
    return z_c_values


def calc_k_c(z_s1_values, z_c_values, thickness):
    k_c_values = []
    numerator = -1j
    denominator = 2 * thickness
    for z_s1, z_c in zip(z_s1_values, z_c_values):
        argument = (z_s1+z_c)/(z_s1-z_c)
        k_c = (numerator / denominator) * np.log(argument)
        k_c_values.append(k_c)
    return k_c_values


################################################################################

# Figures functions

def fig_eff(z_c_values, k_c_values, name, L1, L2, figures_folderpath):
    rows = 5
    columns = 2

    fig2 = plt.figure(2, figsize=(14, 15)) # width, height in inches
    plt.rcParams.update({'font.size': 8})

    ax = fig2.add_subplot(rows, columns, idx*columns - 1)
    real_part_zc = [z.real for z in z_c_values]
    imag_part_zc = [z.imag for z in z_c_values]
    ax.plot(frequencies, real_part_zc, label='Real', color='red')
    ax.plot(frequencies, imag_part_zc, label='Imag', color='blue')
    # Adjusting the rest:
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel(r'$z_{c}$ [rayl]')
    ax.grid(True)
    ax.legend()
    ax.set_title(f'Characteristic impedance', y=1.1)
    ax.text(0, 1.05 ,f'{name} (cavities: {L1}m and {L2}m)', transform=ax.transAxes)

    ax = fig2.add_subplot(rows, columns, idx*columns)
    real_part_kc = [z.real for z in k_c_values]
    imag_part_kc = [z.imag for z in k_c_values]
    ax.plot(frequencies, real_part_kc, label='Real', color='red')
    ax.plot(frequencies, imag_part_kc, label='Imag', color='blue')
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel('k [rad/m]')
    ax.grid(True)
    ax.legend()
    ax.set_title(f'Wavenumber', y=1.1)
    ax.text(0, 1.05 ,f'{name} (cavities: {L1}m and {L2}m)', transform=ax.transAxes)

    return fig2


def save_figures(figures_folderpath, fig2):
    fig2.tight_layout()

    output_file2 = f'untreated_eff_plant_{name}.png'
    fig2.savefig(os.path.join(figures_folderpath, output_file2), dpi=200, bbox_inches='tight')
    plt.close(fig2)

################################################################################

# DATA MANIPULATION functions

def z_c_mean_values(z_c_values_per_plant):
    z_c_mean_values_per_plant = {}
    z_c_std_values_per_plant = {}
    for name, z_c_values_list in z_c_values_per_plant.items():
        # Stack the list of arrays to a 2D array
        z_c_values_array = np.vstack(z_c_values_list)
        # Compute the mean
        mean_z_c_values = np.mean(z_c_values_array, axis=0)
        # Compute the standard deviation - split the complex array into real and imaginary parts
        real_parts = z_c_values_array.real
        imag_parts = z_c_values_array.imag
        std_z_c_values_real = np.std(real_parts, axis=0)
        std_z_c_values_imag = np.std(imag_parts, axis=0)
        # Store values
        z_c_mean_values_per_plant[name] = mean_z_c_values
        z_c_std_values_per_plant[name] = std_z_c_values_real  + 1j * std_z_c_values_imag
    return z_c_mean_values_per_plant, z_c_std_values_per_plant


def fig_z_c_mean(z_c_mean_values_per_plant, frequencies, figures_folderpath, z_c_std_values_per_plant):
    plt.rcParams.update({'font.size': 8})
    for name, z_c_mean_values in z_c_mean_values_per_plant.items():
        std_z_c_values = z_c_std_values_per_plant[name]
        plt.figure(figsize=(5, 5))
        plt.plot(frequencies, z_c_mean_values.real, label='Real', color='red')
        plt.plot(frequencies, z_c_mean_values.imag, label='Imag', color='blue')
        # Plotting the shaded area for standard deviation
        plt.fill_between(frequencies, 
                         (z_c_mean_values.real - std_z_c_values.real), 
                         (z_c_mean_values.real + std_z_c_values.real), 
                         color='red', alpha=0.2, label='std')
        plt.fill_between(frequencies, 
                         (z_c_mean_values.imag - std_z_c_values.imag), 
                         (z_c_mean_values.imag + std_z_c_values.imag), 
                         color='blue', alpha=0.2, label='std')
        # Rest
        plt.xlabel('f [Hz]')
        plt.ylim(-600, 2000)
        plt.ylabel(r'$z_{c}$ [rayl]')
        plt.legend()
        plt.grid(True)
        plt.text(0, 1.05 ,f'{name}', transform=plt.gca().transAxes)
        plt.title(f'Characteristic impedance', y=1.1)
        plt.tight_layout()
        # Save
        output_file = f'zc_mean_plant_{name}.png'
        plt.savefig(os.path.join(figures_folderpath, output_file), dpi=200, bbox_inches='tight')
        plt.close()

################################################################################

# Storing functions

def store_zs1_values(z_s1_values, name, output_folderpath):
    filtered_zs1_values = [z_s1 for z_s1, f in zip(z_s1_values, frequencies) if f >= 200] # Store above 200 Hz (data cleaning)
    zs1_filepath = os.path.join(output_folderpath, f'zs1_{name}.npy')
    np.save(zs1_filepath, filtered_zs1_values)


def store_z_c_mean(z_c_mean_values_per_plant, frequencies, output_folderpath):
    for name, z_c_values in z_c_mean_values_per_plant.items():
        filtered_z_c_values = [z_c for z_c, f in zip(z_c_values, frequencies) if f >= 200] # Store above 200 Hz (data cleaning)
        zc_filepath = os.path.join(output_folderpath, f'zc_mean_{name}.npy')
        np.save(zc_filepath, filtered_z_c_values)


################################################################################

# PROCESS

# Initialize a dictionary to collect z_c_values (for all the combinations) for the specific plant
z_c_values_per_plant = defaultdict(list)

for plant_num in plant_numbers:

    print(f"Processing with plant_num = {plant_num}")

    # Initialize a list to collect z_c_values per combination
    z_c_values_per_combination = []

    for idx, (L1, L2) in enumerate(combinations, start=1):

        # MAIN functions
        name, columns, thickness = select(plant_num, L1, L2)
        frequencies, z_s1_values, z_s2_values = load_excel(filename, sheet_name, columns, z_0)
        check(frequencies, L1, L2, c)
        z_l2_values = calc_z_l2(z_0, L2, frequencies, c)
        z_c_values = calc_z_c(z_s1_values, z_s2_values, z_l2_values)
        k_c_values = calc_k_c(z_s1_values, z_c_values, thickness)

        ### Figure
        fig2 = fig_eff(z_c_values, k_c_values, name, L1, L2, figures_folderpath)

        ### Storing
        if idx == 1: 
            store_zs1_values(z_s1_values, name, output_folderpath) # Store z_s1_values (0 cm cavity) for the subsequent data analysis
        
        ###
        z_c_values_per_combination.append(z_c_values) # Store the z_c_values for the current combination

        # END
        print(f"Loop {idx} is done")

    save_figures(figures_folderpath, fig2)

    # Store the z_c_values for the current plant
    z_c_values_per_plant[name] = z_c_values_per_combination

# DATA MANIPULATION functions
# Calculate mean z_c_values (from all the combinations) for the specific plant
z_c_mean_values_per_plant, z_c_std_values_per_plant = z_c_mean_values(z_c_values_per_plant)
# Store mean z_c_values for the subsequent data analysis
store_z_c_mean(z_c_mean_values_per_plant, frequencies, output_folderpath)

# Display mean z_c_values with standard deviation for a specific plant
fig_z_c_mean(z_c_mean_values_per_plant, frequencies, figures_folderpath, z_c_std_values_per_plant)