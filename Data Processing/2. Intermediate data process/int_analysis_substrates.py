import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

################################################################################

# Directories

filenames = ['1. Preliminary data process\Output\substrates_dry.xlsx', 
             '1. Preliminary data process\Output\substrates_sat.xlsx']
sheet_name = 'Impedance Ratio'

figures_folderpath = r'2. Intermediate data process\Figures'
if not os.path.exists(figures_folderpath):
    os.makedirs(figures_folderpath)

output_folderpath = r'2. Intermediate data process\Outputs'
if not os.path.exists(output_folderpath):
    os.makedirs(output_folderpath)

################################################################################

# Iteration conditions

# Clay balls: 1
# Coco coir: 2
# Coco husk: 3
# Moss: 4
# Perlite: 5
# Pumice: 6
# Vermiculite: 7
substrate_numbers = [1,2,3,4,5,6,7]

################################################################################

# Other variables

# Air impedance
z_0 = 412 # rayl

################################################################################

# MAIN functions

def select(substrate_num):
    substrates = {
        1: {
            'name': 'Clay balls',
            'columns': [1, 2, 3, 4],
            'thickness': [0.04, 0.08]
        },
        2: {
            'name': 'Coco coir',
            'columns': [5, 6, 7, 8],
            'thickness': [0.03, 0.06]
        },
        3: {
            'name': 'Coco husk',
            'columns': [9, 10, 11, 12],
            'thickness': [0.04, 0.08]
        },
        4: {
            'name': 'Moss',
            'columns': [13, 14, 15, 16],
            'thickness': [0.03, 0.06]
        },
        5: {
            'name': 'Perlite',
            'columns': [17, 18, 19, 20],
            'thickness': [0.03, 0.06]
        },
        6: {
            'name': 'Pumice',
            'columns': [21, 22, 23, 24],
            'thickness': [0.04, 0.08]
        },
        7: {
            'name': 'Vermiculite',
            'columns': [25, 26, 27, 28],
            'thickness': [0.03, 0.06]
        }
    }
    # Accessing the info associated with a key
    substrate = substrates.get(substrate_num)
    name = substrate['name']
    columns = substrate['columns']
    thicknesses = substrate['thickness']
    return name, columns, thicknesses


def load_excel(filename, sheet_name, columns, z_0):
    df = pd.read_excel(filename, sheet_name=sheet_name, usecols=[0]+columns, skiprows=1)
    frequencies = df.iloc[:, 0].values
    z_s1_values = (df.iloc[:, 1].values + 1j * df.iloc[:, 2].values) * z_0
    z_s2_values = (df.iloc[:, 3].values + 1j * df.iloc[:, 4].values) * z_0
    return frequencies, z_s1_values, z_s2_values


def calc_z_c(z_s1_values, z_s2_values):
    z_c_values = []
    for z_s1, z_s2 in zip(z_s1_values, z_s2_values):
        z_c = np.sqrt(z_s1 * (2 * z_s2 - z_s1))
        z_c_values.append(z_c)
    return z_c_values


def calc_k_c(z_s1_values, z_s2_values, thicknesses):
    k_c_values = []
    numerator = -1j
    denominator = thicknesses[1]
    for z_s1, z_s2 in zip(z_s1_values, z_s2_values):
        argument = np.sqrt((2*z_s2-z_s1)/z_s1)
        k_c=(numerator/denominator)*np.log((1+argument)/(1-argument))
        k_c_values.append(k_c)
    return k_c_values


def check(k_c_values, thicknesses, frequencies):
    n_list = np.arange(-10, 11)
    for i, k_c in enumerate(k_c_values):
        for n in n_list:
            left = k_c.real * thicknesses[0] # left side (kd1)
            right = k_c.real * thicknesses[1] + n * np.pi # right side (kd2 ± nπ)
            nominator = np.abs(right-left)
            denominator = np.abs((left+right)/2)
            percentage_difference = (nominator / denominator) * 100
            tolerance_percentage = 5 # right and left must differ at least 5%
            if percentage_difference < tolerance_percentage:
                print(f"For frequency {frequencies[i]} and n = {n}, the condition of z1 ≠ z2 is not satisfied.")

################################################################################

# Figures functions

def fig_eff(z_c_values, k_c_values, state, name, figures_folderpath):
    rows = 2
    columns = 2

    fig2 = plt.figure(2, figsize=(10, 10)) # width, height in inches
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
    ax.text(0, 1.05 ,f'{name} {state}', transform=ax.transAxes)

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
    ax.text(0, 1.05 ,f'{name} {state}', transform=ax.transAxes)

    return fig2


def save_figures(figures_folderpath, fig2):
    fig2.tight_layout()

    output_file2 = f'untreated_eff_substrate_{name}.png'
    fig2.savefig(os.path.join(figures_folderpath, output_file2), dpi=200, bbox_inches='tight')
    plt.close(fig2)

    

################################################################################

# Storing functions

def store_zc_values(z_c_values, name, output_folderpath, state):
    filtered_zc_values = [z_c for z_c, f in zip(z_c_values, frequencies) if f >= 200] # Store above 200 Hz (data cleaning)
    zc_filepath = os.path.join(output_folderpath, f'zc_{name}_{state}.npy')
    np.save(zc_filepath, filtered_zc_values)


def store_kc_values(k_c_values, name, output_folderpath, state):
    filtered_kc_values = [k_c for k_c, f in zip(k_c_values, frequencies) if f >= 200] # Store above 200 Hz (data cleaning)
    zc_filepath = os.path.join(output_folderpath, f'kc_{name}_{state}.npy')
    np.save(zc_filepath, filtered_kc_values)

################################################################################

# PROCESS

for substrate_num in substrate_numbers:
    
    print(f"Processing with substrate_num={substrate_num}")
    
    for idx, filename in enumerate(filenames, 1): # Start from iteration from 1
        
        # MAIN Functions
        name, columns, thicknesses = select(substrate_num)
        frequencies, z_s1_values, z_s2_values = load_excel(filename, sheet_name, columns, z_0)
        z_c_values = calc_z_c(z_s1_values, z_s2_values)
        k_c_values = calc_k_c(z_s1_values, z_s2_values, thicknesses)
        check(k_c_values, thicknesses, frequencies)

        # Storing Functions
        if 'dry' in filename:
            state = 'dry'
        else:
            state = 'saturated'
        store_zc_values(z_c_values, name, output_folderpath, state)
        store_kc_values(k_c_values, name, output_folderpath, state)

        # Figures Function
        fig2 = fig_eff(z_c_values, k_c_values, state, name, figures_folderpath)

        # END
        print(f"Loop {idx} is done")
    
    # Figures Function
    save_figures(figures_folderpath, fig2)