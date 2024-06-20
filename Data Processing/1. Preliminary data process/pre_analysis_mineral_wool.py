import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import os

################################################################################

# Directories

filenames = ['1. Preliminary data process\Input\mineral_wool_dry_v1.xlsx', 
             '1. Preliminary data process\Input\mineral_wool_dry_v2.xlsx', 
             '1. Preliminary data process\Input\mineral_wool_sat_v1.xlsx', 
             '1. Preliminary data process\Input\mineral_wool_sat_v2.xlsx']
sheet_name = 'Impedance Ratio'

figures_folderpath = r'1. Preliminary data process\Figures'
if not os.path.exists(figures_folderpath):
    os.makedirs(figures_folderpath)

################################################################################

# Iteration conditions

# Mineral wool + fleece: 1
# Mineral wool: 2
substrate_numbers = [2]

################################################################################

# Other variables

z_0 = 412 # rayl

################################################################################

# MAIN functions

def sub_select(substrate_num):
    substrates = {
        1: {
            'name': 'Mineral wool + fleece',
            'columns': [1, 2, 3, 4],
            'thickness': [0.025, 0.05]
        },
        2: {
            'name': 'Mineral wool',
            'columns': [5, 6, 7, 8],
            'thickness': [0.025, 0.05]
        }
    }
    # Accessing the info associated with a key
    substrate = substrates.get(substrate_num)
    name = substrate['name']
    columns = substrate['columns']
    thickness = substrate['thickness']
    return name, columns, thickness


def load_excel(filename, sheet_name, columns, z_0):
    # Loading the Excel file and creating the dataframe with nessessary data
    df = pd.read_excel(filename, sheet_name=sheet_name, usecols=[0]+columns, skiprows=1)
    # Extracting values from the dataframe
    frequencies = df.iloc[:, 0].values
    z_s1_values = (df.iloc[:, 1].values + 1j * df.iloc[:, 2].values) * z_0
    z_s2_values = (df.iloc[:, 3].values + 1j * df.iloc[:, 4].values) * z_0
    return frequencies, z_s1_values, z_s2_values


################################################################################

# Figures functions

def fig_zs(frequencies, z_s1_values, z_s2_values, name, state, side):
    fig1 = plt.figure(2, figsize=(10, 10)) # width, height in inches
    plt.rcParams.update({'font.size': 8})

    # Plotting the subplots with the same y-axis limits
    y_limits = [
    min([min([z.real for z in z_s1_values]), min([z.imag for z in z_s1_values])]),
    max([max([z.real for z in z_s1_values]), max([z.imag for z in z_s1_values])]),
    min([min([z.real for z in z_s2_values]), min([z.imag for z in z_s2_values])]),
    max([max([z.real for z in z_s2_values]), max([z.imag for z in z_s2_values])])
    ]
    y_min = min(y_limits)
    y_max = max(y_limits)
    
    ax = fig1.add_subplot(4, 2, idx*2 - 1)
    real_part_zs1 = [z.real for z in z_s1_values]
    imag_part_zs1 = [z.imag for z in z_s1_values]
    ax.plot(frequencies, real_part_zs1, label='Real', color='red')
    ax.plot(frequencies, imag_part_zs1, label='Imag', color='blue')
    # Adjusting the rest:
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel(r'$z_{s1}$ [rayl]')
    ax.grid(True)
    ax.legend()
    ax.set_ylim(y_min, y_max)
    ax.text(0, 1.05 ,f'{name} ({state}) - side:{side} - (sample thickness: {thickness[0]} m)', transform=ax.transAxes)

    ax = fig1.add_subplot(4, 2, idx*2)
    real_part_zs2 = [z.real for z in z_s2_values]
    imag_part_zs2 = [z.imag for z in z_s2_values]
    ax.plot(frequencies, real_part_zs2, label='Real', color='red')
    ax.plot(frequencies, imag_part_zs2, label='Imag', color='blue')
    # Adjusting the rest:
    ax.set_xlabel('f [Hz]')
    ax.set_ylabel(r'$z_{s2}$ [rayl]')
    ax.grid(True)
    ax.legend()
    ax.set_ylim(y_min, y_max)
    ax.text(0, 1.05 ,f'{name} ({state}) - side:{side} - (sample thickness: {thickness[1]} m)', transform=ax.transAxes)
    
    return fig1


def save_figures(figures_folderpath, fig1):
    fig1.tight_layout()

    output_file1 = f'zs_substrate_{name}.png'
    fig1.savefig(os.path.join(figures_folderpath, output_file1), dpi=400, bbox_inches='tight')
    plt.close(fig1)



################################################################################

# DATA MANIPULATION

def z_s_mean_values(z_s1_values_per_dry_substrate, z_s2_values_per_dry_substrate, z_s1_values_per_saturated_substrate, z_s2_values_per_saturated_substrate):
    z_s1_mean_values_per_dry_substrate = {}
    z_s2_mean_values_per_dry_substrate = {}
    z_s1_mean_values_per_saturated_substrate = {}
    z_s2_mean_values_per_saturated_substrate = {}
    z_s1_std_values_per_dry_substrate = {}
    z_s2_std_values_per_dry_substrate = {}
    z_s1_std_values_per_saturated_substrate = {}
    z_s2_std_values_per_saturated_substrate = {}
    
    for (name, state), z_s1_values_list in z_s1_values_per_dry_substrate.items():
        # Stack the list of arrays to a 2D array
        z_s1_values_array = np.vstack(z_s1_values_list)
        # Compute the mean
        mean_z_s1_values = np.mean(z_s1_values_array, axis=0)
        # Compute the standard deviation - split the complex array into real and imaginary parts
        real_parts = z_s1_values_array.real
        imag_parts = z_s1_values_array.imag
        std_z_s1_values_real = np.std(real_parts, axis=0)
        std_z_s1_values_imag = np.std(imag_parts, axis=0)
        # Store values
        z_s1_mean_values_per_dry_substrate[name, state] = mean_z_s1_values
        z_s1_std_values_per_dry_substrate[name, state] = std_z_s1_values_real  + 1j * std_z_s1_values_imag

    for (name, state), z_s1_values_list in z_s1_values_per_saturated_substrate.items():
        # Stack the list of arrays to a 2D array
        z_s1_values_array = np.vstack(z_s1_values_list)
        # Compute the mean
        mean_z_s1_values = np.mean(z_s1_values_array, axis=0)
        # Compute the standard deviation - split the complex array into real and imaginary parts
        real_parts = z_s1_values_array.real
        imag_parts = z_s1_values_array.imag
        std_z_s1_values_real = np.std(real_parts, axis=0)
        std_z_s1_values_imag = np.std(imag_parts, axis=0)
        # Store values
        z_s1_mean_values_per_saturated_substrate[name, state] = mean_z_s1_values
        z_s1_std_values_per_saturated_substrate[name, state] = std_z_s1_values_real  + 1j * std_z_s1_values_imag
    
    for (name, state), z_s2_values_list in z_s2_values_per_dry_substrate.items():
        # Stack the list of arrays to a 2D array
        z_s2_values_array = np.vstack(z_s2_values_list)
        # Compute the mean
        mean_z_s2_values = np.mean(z_s2_values_array, axis=0)
        # Compute the standard deviation - split the complex array into real and imaginary parts
        real_parts = z_s2_values_array.real
        imag_parts = z_s2_values_array.imag
        std_z_s2_values_real = np.std(real_parts, axis=0)
        std_z_s2_values_imag = np.std(imag_parts, axis=0)
        # Store values
        z_s2_mean_values_per_dry_substrate[name, state] = mean_z_s2_values
        z_s2_std_values_per_dry_substrate[name, state] = std_z_s2_values_real  + 1j * std_z_s2_values_imag

    for (name, state), z_s2_values_list in z_s2_values_per_saturated_substrate.items():
        # Stack the list of arrays to a 2D array
        z_s2_values_array = np.vstack(z_s2_values_list)
        # Compute the mean
        mean_z_s2_values = np.mean(z_s2_values_array, axis=0)
        # Compute the standard deviation - split the complex array into real and imaginary parts
        real_parts = z_s2_values_array.real
        imag_parts = z_s2_values_array.imag
        std_z_s2_values_real = np.std(real_parts, axis=0)
        std_z_s2_values_imag = np.std(imag_parts, axis=0)
        # Store values
        z_s2_mean_values_per_saturated_substrate[name, state] = mean_z_s2_values
        z_s2_std_values_per_saturated_substrate[name, state] = std_z_s2_values_real  + 1j * std_z_s2_values_imag

    return z_s1_mean_values_per_dry_substrate, z_s2_mean_values_per_dry_substrate, z_s1_mean_values_per_saturated_substrate, z_s2_mean_values_per_saturated_substrate, z_s1_std_values_per_dry_substrate, z_s2_std_values_per_dry_substrate,z_s1_std_values_per_saturated_substrate, z_s2_std_values_per_saturated_substrate



def fig_zs_mean(z_s1_mean_values_per_dry_substrate, z_s2_mean_values_per_dry_substrate, z_s1_mean_values_per_saturated_substrate, z_s2_mean_values_per_saturated_substrate, z_s1_std_values_per_dry_substrate, z_s2_std_values_per_dry_substrate,z_s1_std_values_per_saturated_substrate, z_s2_std_values_per_saturated_substrate):
    rows = 2
    columns = 2

    fig2 = plt.figure(figsize=(10, 10))  # width, height in inches
    plt.rcParams.update({'font.size': 8})

    idx_s1 = 0
    for (name, state), z_s1_mean_values in z_s1_mean_values_per_dry_substrate.items():
        std_z_s1_values = z_s1_std_values_per_dry_substrate[name, state]
        idx_s1 += 1
        ax = fig2.add_subplot(rows, columns, idx_s1)  # Adjusted subplot index
        ax.plot(frequencies, z_s1_mean_values.real, label='Real', color='red')
        ax.plot(frequencies, z_s1_mean_values.imag, label='Imag', color='blue')
        # Plotting the shaded area for standard deviation
        ax.fill_between(frequencies, 
                        (z_s1_mean_values.real - std_z_s1_values.real), 
                        (z_s1_mean_values.real + std_z_s1_values.real), 
                        color='red', alpha=0.2, label='std')
        ax.fill_between(frequencies, 
                        (z_s1_mean_values.imag - std_z_s1_values.imag), 
                        (z_s1_mean_values.imag + std_z_s1_values.imag), 
                        color='blue', alpha=0.2, label='std')
        # Rest
        ax.set_xlabel('f [Hz]')
        ax.set_ylabel(r'$z_{s1}$ [rayl]')
        ax.legend()
        ax.grid(True)
        ax.text(0, 1.05 ,f'{name} ({state}) (sample thickness: {thickness[0]} m)', transform=ax.transAxes)
        #ax.set_title(f'Speciment thickness: {thickness[0]} m', y=1.1)
        fig1.tight_layout()

    idx_s11 = 0
    for (name, state), z_s1_mean_values in z_s1_mean_values_per_saturated_substrate.items():
        std_z_s1_values = z_s1_std_values_per_saturated_substrate[name, state]
        idx_s11 += 1
        ax = fig2.add_subplot(rows, columns, idx_s11+2)  # Adjusted subplot index
        ax.plot(frequencies, z_s1_mean_values.real, label='Real', color='red')
        ax.plot(frequencies, z_s1_mean_values.imag, label='Imag', color='blue')
        # Plotting the shaded area for standard deviation
        ax.fill_between(frequencies, 
                        (z_s1_mean_values.real - std_z_s1_values.real), 
                        (z_s1_mean_values.real + std_z_s1_values.real), 
                        color='red', alpha=0.2, label='std')
        ax.fill_between(frequencies, 
                        (z_s1_mean_values.imag - std_z_s1_values.imag), 
                        (z_s1_mean_values.imag + std_z_s1_values.imag), 
                        color='blue', alpha=0.2, label='std')
        # Rest
        ax.set_xlabel('f [Hz]')
        ax.set_ylabel(r'$z_{s1}$ [rayl]')
        ax.legend()
        ax.grid(True)
        ax.text(0, 1.05 ,f'{name} ({state}) (sample thickness: {thickness[0]} m)', transform=ax.transAxes)
        #ax.set_title(f'Speciment thickness: {thickness[0]} m', y=1.1)
        fig1.tight_layout()

    idx_s2 = 1
    for (name, state), z_s2_mean_values in z_s2_mean_values_per_dry_substrate.items():
        std_z_s2_values = z_s2_std_values_per_dry_substrate[name, state]
        idx_s2 += 1
        ax = fig2.add_subplot(rows, columns, idx_s2)  # Adjusted subplot index
        ax.plot(frequencies, z_s2_mean_values.real, label='Real', color='red')
        ax.plot(frequencies, z_s2_mean_values.imag, label='Imag', color='blue')
        # Plotting the shaded area for standard deviation
        ax.fill_between(frequencies, 
                        (z_s2_mean_values.real - std_z_s2_values.real), 
                        (z_s2_mean_values.real + std_z_s2_values.real), 
                        color='red', alpha=0.2, label='std')
        ax.fill_between(frequencies, 
                        (z_s2_mean_values.imag - std_z_s2_values.imag), 
                        (z_s2_mean_values.imag + std_z_s2_values.imag), 
                        color='blue', alpha=0.2, label='std')
        # Rest
        ax.set_xlabel('f [Hz]')
        ax.set_ylabel(r'$z_{s2}$ [rayl]')
        ax.legend()
        ax.grid(True)
        ax.text(0, 1.05 ,f'{name} ({state}) (sample thickness: {thickness[1]} m)', transform=ax.transAxes)
        #ax.set_title(f'Speciment thickness: {thickness[1]} m', y=1.1)
        fig1.tight_layout()

    idx_s22 = 1
    for (name, state), z_s2_mean_values in z_s2_mean_values_per_saturated_substrate.items():
        std_z_s2_values = z_s2_std_values_per_saturated_substrate[name, state]
        idx_s22 += 1
        ax = fig2.add_subplot(rows, columns, idx_s22 + 2)  # Adjusted subplot index
        ax.plot(frequencies, z_s2_mean_values.real, label='Real', color='red')
        ax.plot(frequencies, z_s2_mean_values.imag, label='Imag', color='blue')
        # Plotting the shaded area for standard deviation
        ax.fill_between(frequencies, 
                        (z_s2_mean_values.real - std_z_s2_values.real), 
                        (z_s2_mean_values.real + std_z_s2_values.real), 
                        color='red', alpha=0.2, label='std')
        ax.fill_between(frequencies, 
                        (z_s2_mean_values.imag - std_z_s2_values.imag), 
                        (z_s2_mean_values.imag + std_z_s2_values.imag), 
                        color='blue', alpha=0.2, label='std')
        # Rest
        ax.set_xlabel('f [Hz]')
        ax.set_ylabel(r'$z_{s2}$ [rayl]')
        ax.legend()
        ax.grid(True)
        ax.text(0, 1.05 ,f'{name} ({state}) (sample thickness: {thickness[1]} m)', transform=ax.transAxes)
        #ax.set_title(f'Speciment thickness: {thickness[1]} m', y=1.1)
        
        fig2.tight_layout()
        output_file2 = f'zs_mean_substrate_{name}.png'
        fig2.savefig(os.path.join(figures_folderpath, output_file2), dpi=400, bbox_inches='tight')
        plt.close(fig2)

    return fig2


################################################################################

# PROCESS

for substrate_num in substrate_numbers:
    
    print(f"Processing with substrate_num={substrate_num}")

    z_s1_values_per_dry_substrate = defaultdict(list)
    z_s2_values_per_dry_substrate = defaultdict(list)
    z_s1_values_per_saturated_substrate = defaultdict(list)
    z_s2_values_per_saturated_substrate = defaultdict(list)

    for idx,filename in enumerate(filenames,1):
        
        # MAIN functions
        name, columns, thickness = sub_select(substrate_num)
        frequencies, z_s1_values, z_s2_values = load_excel(filename, sheet_name, columns, z_0)
        
        # Conditions for sample sides and states
        state = 'dry' if 'dry' in filename else 'sat'
        side = 'v1' if 'v1' in filename else 'v2'

        # Storing surface impedance values for the current condition and substrate
        if state == 'dry':
            key = (name, 'dry')
            z_s1_values_per_dry_substrate[key].append(z_s1_values)
            z_s2_values_per_dry_substrate[key].append(z_s2_values)
        else:
            key = (name, 'saturated')
            z_s1_values_per_saturated_substrate[key].append(z_s1_values)
            z_s2_values_per_saturated_substrate[key].append(z_s2_values)

        # Figures Function
        fig1 = fig_zs(frequencies, z_s1_values, z_s2_values, name, state, side)
    
    # Figures Function
    save_figures(figures_folderpath, fig1)

# DATA MANIPULATION functions
z_s1_mean_values_per_dry_substrate, z_s2_mean_values_per_dry_substrate, z_s1_mean_values_per_saturated_substrate, z_s2_mean_values_per_saturated_substrate, z_s1_std_values_per_dry_substrate, z_s2_std_values_per_dry_substrate,z_s1_std_values_per_saturated_substrate, z_s2_std_values_per_saturated_substrate = z_s_mean_values(z_s1_values_per_dry_substrate, z_s2_values_per_dry_substrate, z_s1_values_per_saturated_substrate, z_s2_values_per_saturated_substrate)
fig2 = fig_zs_mean(z_s1_mean_values_per_dry_substrate, z_s2_mean_values_per_dry_substrate, z_s1_mean_values_per_saturated_substrate, z_s2_mean_values_per_saturated_substrate, z_s1_std_values_per_dry_substrate, z_s2_std_values_per_dry_substrate,z_s1_std_values_per_saturated_substrate, z_s2_std_values_per_saturated_substrate)