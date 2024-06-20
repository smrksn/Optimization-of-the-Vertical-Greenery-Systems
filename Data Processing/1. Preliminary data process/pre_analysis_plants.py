import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
import os

################################################################################

# Directories

filename = '1. Preliminary data process\Input\plants.xlsx'
sheet_name = 'Impedance Ratio'

figures_folderpath = r'1. Preliminary data process\Figures'
if not os.path.exists(figures_folderpath):
        os.makedirs(figures_folderpath)

################################################################################

# Iteration conditions

# 1: 'Ajuga',
# 2: 'Bergenia',
# 3: 'Grass',
# 4: 'Heuchera',
# 5: 'Waldestania'
plant_numbers = [1,2,3,4,5]

# Cavity combinations
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



# According to the Utsuno 1989, the appropriate set of air space depths must be selected so as to not satisfy:
# Equation: f(L-L') = nc/2
def check(frequencies, L1, L2, c):
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
    return True    

################################################################################

# Figures functions

def fig_zs(z_s1_values, z_s2_values, name, L1, L2, figures_folderpath):
        rows = 6
        columns = 1

        fig1 = plt.figure(1, figsize=(8, 15)) # width, height in inches
        plt.rcParams.update({'font.size': 8})
        
        y_limits = [
        min([min([z.real for z in z_s1_values]), min([z.imag for z in z_s1_values])]),
        max([max([z.real for z in z_s1_values]), max([z.imag for z in z_s1_values])]),
        min([min([z.real for z in z_s2_values]), min([z.imag for z in z_s2_values])]),
        max([max([z.real for z in z_s2_values]), max([z.imag for z in z_s2_values])])
        ]
        y_min = min(y_limits)
        y_max = max(y_limits)
        
        if (idx == 1 and L1 == 0):
            ax = fig1.add_subplot(rows, columns, 1)
            real_part_zs1 = [z.real for z in z_s1_values]
            imag_part_zs1 = [z.imag for z in z_s1_values]
            ax.plot(frequencies, real_part_zs1, label='Real', color='red')
            ax.plot(frequencies, imag_part_zs1, label='Imag', color='blue')
            ax.set_xlabel('f [Hz]')
            ax.set_ylabel(r'$z_{s1}$ [rayl]')
            ax.grid(True)
            ax.legend()
            ax.set_ylim(y_min, y_max)
            #ax.set_title(f'Surface impedance', y=1.1)
            ax.text(0, 1.05 ,f'{name} (cavity: {L1}m)', transform=ax.transAxes)

        if (idx == 1 and L2 == 0.01):
            index = 2
        if (idx == 2 and L2 == 0.02):
            index = 3
        if (idx == 3 and L2 == 0.03):
            index = 4
        if (idx == 4 and L2 == 0.04):
            index = 5
        if (idx == 5 and L2 == 0.05):
            index = 6

        ax = fig1.add_subplot(rows, columns, index)
        real_part_zs2 = [z.real for z in z_s2_values]
        imag_part_zs2 = [z.imag for z in z_s2_values]
        ax.plot(frequencies, real_part_zs2, label='Real', color='red')
        ax.plot(frequencies, imag_part_zs2, label='Imag', color='blue')
        ax.set_xlabel('f [Hz]')
        ax.set_ylabel(r'$z_{s2}$ [rayl]')
        ax.grid(True)
        ax.legend()
        ax.set_ylim(y_min, y_max)
        #ax.set_title(f'Surface impedance', y=1.1)
        ax.text(0, 1.05 ,f'{name} (cavity: {L2}m)', transform=ax.transAxes)

        return fig1


def save_figures(figures_folderpath, fig1):
    fig1.tight_layout()

    output_file1 = f'zs_plant_{name}.png'
    fig1.savefig(os.path.join(figures_folderpath, output_file1), dpi=200, bbox_inches='tight')
    plt.close(fig1)

################################################################################

# PROCESS

for plant_num in plant_numbers:

    print(f"Processing with plant_num={plant_num}")

    for idx, (L1, L2) in enumerate(combinations, start=1):

        # MAIN functions
        name, columns, thickness = select(plant_num, L1, L2)
        frequencies, z_s1_values, z_s2_values = load_excel(filename, sheet_name, columns, z_0)
        result = check(frequencies, L1, L2, c)

        # Figures functions
        fig1 = fig_zs(z_s1_values, z_s2_values, name, L1, L2, figures_folderpath)
    
    save_figures(figures_folderpath, fig1)