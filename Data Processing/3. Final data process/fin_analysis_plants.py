import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

################################################################################

# Directories

input_folderpath = r'2. Intermediate data process\Outputs'

figures_folderpath = r'3. Final data process\Figures'
if not os.path.exists(figures_folderpath):
    os.makedirs(figures_folderpath)

output_folderpath = r'3. Final data process\Output'
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

################################################################################

# Other variables

# Min and Max substrate thickness (RL graph)
start_thickness = 0.1
end_thickness = 0.5
layer_thicknesses_range = np.linspace(start_thickness, end_thickness, 100)

z_0 = 413 # rayl
c_0 = 343 # m/s, speed of sound in the air
frequencies = np.arange(200, 1601, 2)

################################################################################

# MAIN functions

def select(plant_num):
    plants = {
        1: {
            'name': 'Ajuga',
            'thickness': [0.08] #meters
        },
        2: {
            'name': 'Bergenia',
            'thickness': [0.05] #meters
        },
        3: {
            'name': 'Grass',
            'thickness': [0.06] #meters
        },
        4: {
            'name': 'Heuchera',
            'thickness': [0.15] #meters
        },
        5: {
            'name': 'Waldestania',
            'thickness': [0.08] #meters
        }
    }
    # Accessing the plant
    plant = plants.get(plant_num)
    name = plant['name']
    thickness = plant['thickness'][0]
    return name, thickness
   

def load_zs1_values(name, input_folderpath):
    files = os.listdir(input_folderpath)
    for file in files:
        if name in file and "zs1" in file:
            filepath= os.path.join(input_folderpath, file)
            z_s1_values = np.load(filepath)
    return z_s1_values


def load_z_c_values(name, input_folderpath):
    files = os.listdir(input_folderpath)
    for file in files:
        if name in file and "zc" in file:
            filepath = (os.path.join(input_folderpath, file))
            z_c_values = np.load(filepath)
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


def unwrapphase(k_c_values, frequencies):
    k_c_real_values = [k.real for k in k_c_values]
    k_c_imag_values = [k.imag for k in k_c_values]
    k_c_real_values_above = []
    k_c_real_values_below = []
    for i, (k, f) in enumerate(zip(k_c_real_values, frequencies)):
        if f>300:
            k_c_real_values_above.append(k)
        else:
            k_c_real_values_below.append(k)
    k_c_real_values_above_unwrapped = np.unwrap(k_c_real_values_above)
    k_c_real_values = k_c_real_values_below + list(k_c_real_values_above_unwrapped)
    k_c_values = np.array(k_c_real_values) + 1j * np.array(k_c_imag_values)
    return k_c_values


def cleaning(k_c_values):
    # Create a boolean mask for elements that meet the condition
    mask = (k_c_values.real > 0) & (k_c_values.imag < 0)
    # Replace elements that do not meet the condition with NaN + 1j * NaN
    k_c_values[~mask] = np.nan + 1j * np.nan
    return k_c_values


def calc_D(frequencies, k_c_values):
    c_values=[]
    D1_values=[]
    D2_values=[]
    D3_values=[]
    D4_values=[]
    for i, (f, k) in enumerate(zip(frequencies, k_c_values)):
        c = (2 * np.pi * f) / k.real
        c_values.append(c)
        wavelength = c/f
        D1 = 1*(wavelength/4) # m=1; particle pressure minima
        D2 = 2*(wavelength/4) # n=2; particle pressure maxima
        D3 = 3*(wavelength/4) # m=3; particle pressure minima
        D4 = 4*(wavelength/4) # n=4; particle pressure maxima
        D1_values.append(D1)
        D2_values.append(D2)
        D3_values.append(D3)
        D4_values.append(D4)
    return c_values, D1_values, D2_values, D3_values, D4_values


def calc_T(frequencies, layer_thicknesses_range, k_c_values, z_c_values):
    T_values = np.zeros((len(layer_thicknesses_range), len(frequencies), 2, 2), dtype=np.complex128)
    for i, d in enumerate(layer_thicknesses_range):
        for j, frequency in enumerate(frequencies):
            T = np.array([[np.cos(k_c_values[j]*d), 1j*z_c_values[j]*np.sin(k_c_values[j]*d)],
                          [1j/z_c_values[j]*np.sin(k_c_values[j]*d), np.cos(k_c_values[j]*d)]], dtype=np.complex128)
            T_values[i, j] = T
    return T_values


def calc_RL(frequencies, layer_thicknesses_range, T_values, z_0):
    RL_values = np.zeros((len(layer_thicknesses_range), len(frequencies)))
    for i, d in enumerate(layer_thicknesses_range):
        for j, frequency in enumerate(frequencies):
            T = T_values[i, j]
            numerator_r = T[0,0] + T[0,1]/z_0  - T[1,0]*z_0 - T[1,1]
            denominator_r = T[0,0] + T[0,1]/z_0 + T[1,0]*z_0 + T[1,1]
            if denominator_r != 0 and not np.isnan(denominator_r):
                reflection_coefficient = numerator_r / denominator_r
                reflection_power = np.abs(reflection_coefficient)**2
                RL = -10*np.log10(reflection_power)
                RL_values[i, j] = RL
            else:
                RL_values[i, j] = np.nan
    return RL_values

################################################################################

# Figures functions

def fig_eff(z_c_values, k_c_values, name):
    rows = 1
    columns = 2

    fig2 = plt.figure(2, figsize=(8, 4)) # width, height in inches
    plt.rcParams.update({'font.size': 8})

    idx = 0
    idx += 1

    ax = fig2.add_subplot(rows, columns, idx*columns - 1)
    real_part_zc = [z.real for z in z_c_values]
    imag_part_zc = [z.imag for z in z_c_values]
    ax.plot(frequencies, real_part_zc, label='Real', color='red')
    ax.plot(frequencies, imag_part_zc, label='Imag', color='blue')
    # Adjusting the rest:
    ax.set_xlabel('f [Hz]')
    ax.set_xlim(200, 1600)
    ax.set_ylabel(r'$z_{c}$ [rayl]')
    ax.set_ylim(-100, 1750)
    ax.grid(True)
    ax.legend()
    ax.set_title(f'Characteristic impedance', y=1.1)
    ax.text(0, 1.05 ,f'{name}', transform=ax.transAxes)

    ax = fig2.add_subplot(rows, columns, idx*columns)
    real_part_kc = [z.real for z in k_c_values]
    imag_part_kc = [z.imag for z in k_c_values]
    ax.plot(frequencies, real_part_kc, label='Real', color='red')
    ax.plot(frequencies, imag_part_kc, label='Imag', color='blue')
    ax.set_xlabel('f [Hz]')
    ax.set_xlim(200, 1600)
    ax.set_ylabel('k [rad/m]')
    ax.set_ylim(-30, 50)
    ax.grid(True)
    ax.legend()
    ax.set_title(f'Wavenumber', y=1.1)
    ax.text(0, 1.05 ,f'{name}', transform=ax.transAxes)

    return fig2


def fig_analysis(z_c_values, name):
    rows = 1
    columns = 1

    fig3 = plt.figure(3, figsize=(6, 5)) # width, height in inches
    plt.rcParams.update({'font.size': 8})

    idx = 0
    idx += 1

    ax = fig3.add_subplot(rows, columns, idx*columns)
    z_real = [z.real for z in z_c_values]
    z_imag = [z.imag for z in z_c_values]
    # Separate positive and negative values of z_imag
    ax.plot(frequencies, z_real, label='Resistance, R', color='red')
    ax.axhline(y=z_0, color='black', linestyle='--', label=r'Impedance of air, $z_0$')
    for i in range(len(frequencies)):
        if z_imag[i] >= 0:
            ax.plot(frequencies[i], z_imag[i], marker='o', markersize = 2, color='green')
        else:
            ax.plot(frequencies[i], abs(z_imag[i]), marker='o', markersize = 2, color='purple')
    # Manually adding custom legend entries
    custom_line1 = plt.Line2D([0], [0], color='green', linestyle='dotted')
    custom_line2 = plt.Line2D([0], [0], color='purple', linestyle='dotted')
    handles, labels = plt.gca().get_legend_handles_labels()
    handles.extend([custom_line1, custom_line2])
    labels.extend([r'$Mass-reactance, X_m$', r'$Stiffness-reactance, X_s$'])
    plt.legend(handles=handles, labels=labels)
    # Adjusting the rest:
    ax.set_xlabel('f [Hz]')
    ax.set_xlim(200, 1600)
    ax.set_ylabel('z [rayl]')
    ax.set_ylim(0, 1750)
    ax.grid(True)
    #ax.set_title(f'Stiffness- and mass-reactance', y=1.1)
    ax.text(0, 1.05 ,f'{name}', transform=ax.transAxes)

    return fig3


def fig_result(frequencies, layer_thicknesses_range, RL_values, name):
    rows = 1
    columns = 1
    
    fig4 = plt.figure(4, figsize=(6, 5))
    plt.rcParams.update({'font.size': 8}) 

    idx = 0
    idx += 1

    ax = fig4.add_subplot(rows,columns,idx*columns)

    im = ax.imshow(RL_values, extent=(min(frequencies), max(frequencies), min(layer_thicknesses_range), max(layer_thicknesses_range)),
                aspect='auto', origin='lower', cmap='cubehelix', vmax=30)
    plt.colorbar(im, ax=ax, label='RL [dB]')
    # Adding quarter wavelengths:
    ax.plot(frequencies, D1_values, label='m=1', color='deepskyblue', linestyle='--')
    ax.plot(frequencies, D2_values, label='n=2', color='orangered', linestyle='--')
    ax.plot(frequencies, D3_values, label='m=3', color='lightskyblue', linestyle='--')
    ax.plot(frequencies, D4_values, label='n=4', color='tomato', linestyle='--')
    # Rest
    ax.set_xlabel('f [Hz]')
    ax.set_xlim(200, 1600)
    ax.set_ylabel('d [m]')
    ax.set_ylim(start_thickness, end_thickness)
    #ax.set_title('Return Loss', y=1.1)
    ax.text(0, 1.05 ,f'{name}', transform=ax.transAxes)
    ax.legend(loc='upper right')

    return fig4


def save_figures(figures_folderpath, fig2, fig3, fig4):
    fig2.tight_layout()
    fig3.tight_layout()
    fig4.tight_layout()

    output_file2 = f'eff_plant_{name}.png'
    fig2.savefig(os.path.join(figures_folderpath, output_file2), dpi=200, bbox_inches='tight')
    plt.close(fig2)

    output_file3 = f'ana_plant_{name}.png'
    fig3.savefig(os.path.join(figures_folderpath, output_file3), dpi=200, bbox_inches='tight')
    plt.close(fig3)

    output_file4 = f'res_plant_{name}.png'
    fig4.savefig(os.path.join(figures_folderpath, output_file4), dpi=200, bbox_inches='tight')
    plt.close(fig4)

################################################################################

# Storing functions

def store_eff_values(name, z_c_values, k_c_values, output_folderpath):
    zc_filepath = os.path.join(output_folderpath, f'zc_{name}.npy')
    kc_filepath = os.path.join(output_folderpath, f'kc_{name}.npy')
    np.save(zc_filepath, z_c_values)
    np.save(kc_filepath, k_c_values)

################################################################################

# PROCESS

for plant_num in plant_numbers:
    
    print(f"Processing with substrate_num = {plant_num}")

    # MAIN Functions
    name, thickness = select(plant_num)
    z_s1_values = load_zs1_values(name, input_folderpath)
    z_c_values = load_z_c_values(name, input_folderpath)
    k_c_values = calc_k_c(z_s1_values, z_c_values, thickness)
    k_c_values = unwrapphase(k_c_values, frequencies)
    k_c_values = cleaning(k_c_values)
    c_values, D1_values, D2_values, D3_values, D4_values = calc_D(frequencies, k_c_values)
    T_values = calc_T(frequencies, layer_thicknesses_range, k_c_values, z_c_values)
    RL_values = calc_RL(frequencies, layer_thicknesses_range, T_values, z_0)

    # Storing Functions
    store_eff_values(name, z_c_values, k_c_values, output_folderpath)

    # Figures Functions
    fig2 = fig_eff(z_c_values, k_c_values, name)
    fig3 = fig_analysis(z_c_values, name)
    fig4 = fig_result(frequencies, layer_thicknesses_range, RL_values, name)
    save_figures(figures_folderpath, fig2, fig3, fig4)

    # END
    print(f"Loop for {name} is done")