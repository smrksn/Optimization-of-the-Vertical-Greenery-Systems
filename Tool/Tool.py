import numpy as np
import warnings
import matplotlib.pyplot as plt
import os
import tkinter as tk
from tkinter import ttk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

##################################################################################################################

folderpath = r'Effective Parameters'

##################################################################################################################

# Variables
z_0 = 413 # air impedance
c_0 = 343 # speed of sound in air
frequencies = np.arange(200, 1601, 2) # working frequency range

##################################################################################################################

# Global variables

# Default values (backup)
substrate_active = True
plant_active = True
cavity_active = True
substrate = "Clay balls"
state = "dry"  # dry/saturated
plant = "Ajuga"
d_substrate = 0.1  # meters
d_plant = 0.1  # meters
d_cavity = 0.1 # meters

def save_variables():
    np.savez("cache.npz", 
             substrate_active=substrate_active,
             plant_active=plant_active,
             cavity_active=cavity_active,
             substrate=substrate,
             state=state,
             plant=plant,
             d_substrate=d_substrate,
             d_plant=d_plant,
             d_cavity=d_cavity)

def load_variables():
    global substrate_active, plant_active, cavity_active, substrate, state, plant, d_substrate, d_plant, d_cavity
    try:
        data = np.load("cache.npz")
        substrate_active = bool(data["substrate_active"])
        plant_active = bool(data["plant_active"])
        cavity_active = bool(data["cavity_active"])
        substrate = str(data["substrate"])
        state = str(data["state"])
        plant = str(data["plant"])
        d_substrate = float(data["d_substrate"])
        d_plant = float(data["d_plant"])
        d_cavity = float(data["d_cavity"])
    except FileNotFoundError:
        print("No saved variables found. Using default values.")
    except KeyError:
        print("The saved variables file exists but is empty. Using default values.")

##################################################################################################################

# INPUT GUI

def update_layer_activness():
    global plant_active, substrate_active, cavity_active 
    substrate_active = substrate_active_boolean.get()
    plant_active = plant_active_boolean.get()
    cavity_active = cavity_active_boolean.get()
    save_variables()

def update_layer_selection(event):
    global substrate, plant
    substrate = substrate_options.get()
    plant = plant_options.get()
    save_variables()

def update_substrate_state():
    global state
    state = state_var.get()
    save_variables()

def update_substrate_thickness(value):
    global d_substrate
    d_substrate = round(float(value) / 100, 2)
    save_variables()

def update_plant_thickness(value):
    global d_plant
    d_plant = round(float(value) / 100, 2)
    save_variables()

def update_cavity_thickness(value):
    global d_cavity
    d_cavity = round(float(value) / 100, 2)
    save_variables()



def input_interface():
    global substrate_active_boolean, plant_active_boolean, cavity_active_boolean, substrate_options, plant_options, state_var, substrate_slider, plant_slider, cavity_slider 

    load_variables()  # Loading variables when GUI starts

    # Creating main window
    root = tk.Tk()
    root.title("Performance Assessment Tool")
    # Configuring the style
    style = ttk.Style()
    style.configure("TLabel", padding=5)
    style.configure("TCheckbutton", padding=5)
    style.configure("TRadiobutton", padding=5)
    style.configure("TButton", padding=5)


    # Adding title
    intro_label = ttk.Label(root, text="Configure Vertical Greenery System", font=("Helvetica", 12, "bold"))
    intro_label.grid(row=0, column=0, columnspan=4, padx=5, pady=5, sticky="nsew")


    # Including scheme
    fig = fig_scheme_opening()
    fig_frame = ttk.Frame(root, relief=tk.SUNKEN, borderwidth=2)
    fig_frame.grid(row=1, column=0, columnspan=4, padx=5, pady=5, sticky="nsew")
    # Adding figure to the frame
    canvas = FigureCanvasTkAgg(fig, master=fig_frame)
    canvas.draw()
    canvas.get_tk_widget().grid(row=0, column=0, columnspan=4, sticky="nsew")
    # Adding explanatory text below the figure
    explanatory_text = ttk.Label(fig_frame, text="Rigidly-Backed System - Transfer Matrix Scheme", font=("Helvetica", 10, "italic"), background="white")
    explanatory_text.grid(row=1, column=0, columnspan=4, sticky="nsew")


    # Layers activation/deactivation
    active_frame = ttk.LabelFrame(root, text="Activate/Deactivate Layers")
    active_frame.grid(row=2, column=0, columnspan=4, padx=10, pady=10, sticky="ew")
    
    plant_active_boolean = tk.BooleanVar(value=plant_active)
    plant_checkbutton = ttk.Checkbutton(active_frame, text="Foliage", variable=plant_active_boolean, command=update_layer_activness)
    plant_checkbutton.grid(row=0, column=1, padx=5, pady=5)
    
    substrate_active_boolean = tk.BooleanVar(value=substrate_active)
    substrate_checkbutton = ttk.Checkbutton(active_frame, text="Substrate", variable=substrate_active_boolean, command=update_layer_activness)
    substrate_checkbutton.grid(row=0, column=2, padx=5, pady=5)

    cavity_active_boolean = tk.BooleanVar(value=cavity_active)
    cavity_checkbutton = ttk.Checkbutton(active_frame, text="Cavity", variable=cavity_active_boolean, command=update_layer_activness)
    cavity_checkbutton.grid(row=0, column=3, padx=5, pady=5)


    # Plant selection
    plant_frame = ttk.LabelFrame(root, text="Plant Selection")
    plant_frame.grid(row=3, column=0, columnspan=4, padx=10, pady=10, sticky="ew")

    plant_options = tk.StringVar(value=plant)
    plant_menu = ttk.OptionMenu(plant_frame, plant_options, plant, "Ajuga", "Bergenia", "Grass", "Heuchera", "Waldestania", command=update_layer_selection)
    plant_menu.grid(row=0, column=0, padx=5, pady=5)


    # Substrate and corresponding state selection
    substrate_frame = ttk.LabelFrame(root, text="Substrate Selection and State")
    substrate_frame.grid(row=4, column=0, columnspan=4, padx=10, pady=10, sticky="ew")

    substrate_options = tk.StringVar(value=substrate)
    substrate_menu = ttk.OptionMenu(substrate_frame, substrate_options, substrate, "Clay balls", "Coco coir", "Coco husk", "Moss", "Perlite", "Pumice", "Vermiculite", "Mineral wool", command=update_layer_selection)
    substrate_menu.grid(row=0, column=0, padx=5, pady=5)

    state_var = tk.StringVar(value=state)
    dry_button = ttk.Radiobutton(substrate_frame, text="Dry", variable=state_var, value="dry", command=update_substrate_state)
    dry_button.grid(row=0, column=1, padx=5, pady=5)
    
    saturated_button = ttk.Radiobutton(substrate_frame, text="Saturated", variable=state_var, value="saturated", command=update_substrate_state)
    saturated_button.grid(row=0, column=2, padx=5, pady=5)


    # Sliders for foliage layer thickness
    plant_slider_frame = ttk.LabelFrame(root, text="Foliage Layer Thickness")
    plant_slider_frame.grid(row=5, column=0, columnspan=4, padx=10, pady=10, sticky="ew")
    plant_slider = ttk.Scale(plant_slider_frame, from_=0, to=20, orient="horizontal", length=200)
    plant_slider.set(d_plant * 100)
    plant_slider.grid(row=0, column=0, padx=5, pady=5)
    plant_slider_value = ttk.Label(plant_slider_frame, text=f"d (plant) = {int(d_plant * 100)} cm")
    plant_slider_value.grid(row=0, column=1, padx=5, pady=5)
    # Updating the plant slider value label dynamically
    def update_plant_slider_label(value):
        plant_slider_value.config(text=f"d (plant) = {int(float(value))} cm")
        update_plant_thickness(value)
    plant_slider.config(command=update_plant_slider_label)


    # Sliders for substrate layer thickness
    substrate_slider_frame = ttk.LabelFrame(root, text="Substrate Layer Thickness")
    substrate_slider_frame.grid(row=6, column=0, columnspan=4, padx=10, pady=10, sticky="ew")
    substrate_slider = ttk.Scale(substrate_slider_frame, from_=0, to=20, orient="horizontal", length=200)
    substrate_slider.set(d_substrate * 100)
    substrate_slider.grid(row=0, column=0, padx=5, pady=5)
    substrate_slider_value = ttk.Label(substrate_slider_frame, text=f"d (substrate) = {int(d_substrate * 100)} cm")
    substrate_slider_value.grid(row=0, column=1, padx=5, pady=5)
    # Updating the substrate slider value label dynamically
    def update_substrate_slider_label(value):
        substrate_slider_value.config(text=f"d (substrate) = {int(float(value))} cm")
        update_substrate_thickness(value)
    substrate_slider.config(command=update_substrate_slider_label)


    # Sliders for cavity thickness
    cavity_slider_frame = ttk.LabelFrame(root, text="Cavity Thickness")
    cavity_slider_frame.grid(row=7, column=0, columnspan=4, padx=10, pady=10, sticky="ew")
    cavity_slider = ttk.Scale(cavity_slider_frame, from_=0, to=100, orient="horizontal", length=200)
    cavity_slider.set(d_cavity * 100)
    cavity_slider.grid(row=0, column=0, padx=5, pady=5)
    cavity_slider_value = ttk.Label(cavity_slider_frame, text=f"d (cavity) = {int(d_cavity * 100)} cm")
    cavity_slider_value.grid(row=0, column=1, padx=5, pady=5)
    # Updating the cavity slider value label dynamically
    def update_cavity_slider_label(value):
        cavity_slider_value.config(text=f"d (cavity) = {int(float(value))} cm")
        update_cavity_thickness(value)
    cavity_slider.config(command=update_cavity_slider_label)


    # Button to update selections and print (check) assigned variables
    def update(root):
        root.destroy()
    # Update and Print button
    update_button = ttk.Button(root, text="Compile", command=lambda: update(root))
    update_button.grid(row=8, column=0, columnspan=4, padx=10, pady=10, sticky="ew")

    root.mainloop()

##################################################################################################################

# MAIN FUNCTIONS

def check():
    print("Plant layer active:", plant_active)
    print("Substrate layer active:", substrate_active)
    print("Cavity active:", cavity_active)
    print("Chosen plant:", plant)
    print(f"Chosen substrate: {substrate}, state: {state}")
    print("Plant layer thickness:", d_plant)
    print("Substrate layer thickness:", d_substrate)
    print("Cavity thickness:", d_cavity)


def extract_filepaths(folderpath, substrate, plant):
    kc_filepaths = []
    zc_filepaths = []
    files = os.listdir(folderpath)
    for file in files:
        if any(word in file for word in [substrate, plant]):
            if "kc" in file:
                kc_filepaths.append(os.path.join(folderpath, file))
            elif "zc" in file:
                zc_filepaths.append(os.path.join(folderpath, file))
    return kc_filepaths, zc_filepaths


def load_eff_values(kc_filepaths, zc_filepaths, state):
    if state == "dry":
        kc_values_substrate = np.load(kc_filepaths[0])
        zc_values_substrate = np.load(zc_filepaths[0])
    if state == "saturated":
        kc_values_substrate = np.load(kc_filepaths[1])
        zc_values_substrate = np.load(zc_filepaths[1])
    kc_values_plant = np.load(kc_filepaths[2])
    zc_values_plant = np.load(zc_filepaths[2])
    return kc_values_substrate, kc_values_plant, zc_values_substrate, zc_values_plant


def calc_T_for_d(z_0, c_0, substrate_active, plant_active, frequencies, d_substrate, d_plant, kc_values_substrate, kc_values_plant, zc_values_substrate, zc_values_plant):
    
    T_values_cavity = []
    T_values_substrate = []
    T_values_plant = []

    for j, f in enumerate(frequencies):

        # Transfer Matrix of the Cavity
        if cavity_active == True:
            omega_cavity = 2*np.pi*f
            T_cavity = np.array([[np.cos(omega_cavity * d_cavity / c_0), 1j * z_0 * np.sin(omega_cavity * d_cavity / c_0)],
                        [1j * np.sin(omega_cavity * d_cavity / c_0) / z_0, np.cos(omega_cavity * d_cavity / c_0)]], dtype=np.complex128)
            T_values_cavity.append(T_cavity)
        elif cavity_active == False:
            T_values_cavity.append(np.eye(2)) # If cavity is not active, appending an identity matrix [[1, 0],[0, 1]]
        
        # Transfer Matrix of the Substrate
        if substrate_active == True:
            T_substrate = np.array([[np.cos(kc_values_substrate[j]*d_substrate),  1j*zc_values_substrate[j]*np.sin(kc_values_substrate[j]*d_substrate)],
                                    [1j/zc_values_substrate[j]*np.sin(kc_values_substrate[j]*d_substrate),  np.cos(kc_values_substrate[j]*d_substrate)]], dtype=np.complex128)
            T_values_substrate.append(T_substrate)
        elif substrate_active == False:
            T_values_substrate.append(np.eye(2)) # If substrate is not active, appending an identity matrix [[1, 0],[0, 1]]
        
        # Transfer Matrix of the Plant
        if plant_active == True:
            T_plant = np.array([[np.cos(kc_values_plant[j]*d_plant), 1j*zc_values_plant[j]*np.sin(kc_values_plant[j]*d_plant)],
                          [1j/zc_values_plant[j]*np.sin(kc_values_plant[j]*d_plant), np.cos(kc_values_plant[j]*d_plant)]], dtype=np.complex128)
            T_values_plant.append(T_plant)
        elif plant_active == False:
            T_values_plant.append(np.eye(2)) # If plant is not active, appending an identity matrix [[1, 0],[0, 1]]
    
    T_values_cavity = np.array(T_values_cavity)
    T_values_substrate = np.array(T_values_substrate)
    T_values_plant = np.array(T_values_plant)

    return T_values_cavity, T_values_substrate, T_values_plant


def T_total(T_values_cavity, T_values_substrate, T_values_plant):  
    # Performing transfer matrix multiplication (sequance as sound wave encounters the layers)
    T_values = T_values_plant @ T_values_substrate @ T_values_cavity
    return T_values


def RL_total(T_values, z_0):
    RL_values = []
    R_values = []
    for T in T_values:
        # T for the rigidly backed system
        numerator_r =   T[0,0] - T[1,0]*z_0
        denominator_r = T[0,0] + T[1,0]*z_0
        if denominator_r != 0 and not np.isnan(denominator_r):
            # Pressure reflection coefficient
            r = numerator_r / denominator_r
            # Reflection power
            R = np.abs(r)**2 # R = |r|^2
            R_values.append(R)
            # Return Loss
            RL = -10*np.log10(R)
            RL_values.append(RL)
        else:
            R_values.append(np.nan)
            RL_values.append(np.nan)
    return R_values, RL_values


def abs_total(R_values):
    abs_values = []
    for R in R_values:
        # Absorption coefficient
        abs = 1 - R
        abs_values.append(abs)
    return abs_values


def abs_to_third_octave_bands(frequencies, abs_values):
    # List of tuples containing frequency-abs pairs.
    abs_values = list(zip(frequencies, abs_values))
    # Dictionary to store third-octave band absorption coefficients
    abs_third_octave = {}
    # List of center frequencies for third-octave bands
    center_frequencies = [200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600]
    for center_frequency in center_frequencies:
        # The third-octave band limits
        k = np.sqrt(2) ** (1/3)
        lower_limit = center_frequency / k
        upper_limit = center_frequency * k
        # Filtering frequencies within the current third-octave band
        filtered_abs = [abs_val for freq, abs_val in abs_values if lower_limit <= freq <= upper_limit]
        # The average absorption coefficient for the band, ignoring nan values
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            average_abs = np.nanmean(filtered_abs)
        abs_third_octave[center_frequency] = average_abs
    return abs_third_octave    


def practical_abs(abs_third_octave):

    # Rounding function as specificed in ISO 11654
    def custom_round(value):
        if np.isnan(value):
            return np.nan
        value = round(value * 100) / 100  # Round to the second decimal place
        fractional_part = value * 100 % 10
        if fractional_part < 3:
            rounded_value = np.floor(value * 10) / 10
        elif fractional_part < 8:
            rounded_value = np.floor(value * 10) / 10 + 0.05
        else:
            rounded_value = np.ceil(value * 10) / 10
        return min(rounded_value, 1.0)  # Capping the value at 1.0
    
    octave_bands = [250, 500, 1000, 1600]
    band_mappings = {250: [200, 250, 315],
                    500: [400, 500, 630],
                    1000: [800, 1000, 1250],
                    1600: [1600]}
    abs_full_octave = {}
    for octave in octave_bands:
        one_third_values = [abs_third_octave.get(band) for band in band_mappings[octave] if abs_third_octave.get(band) is not None]
        # The arithmetic mean value of the three 1/3 octave abs as sepcified in ISO 11654
        mean_value = np.mean(one_third_values)
        abs_full_octave[octave] = custom_round(mean_value)
    
    return abs_full_octave


def abs_weighted(abs_full_octave):
    # Reference curve based on ISO 11654
    reference_curve = {250: 0.8, 500: 1, 1000: 1, 1600: 1}
    # Initialized shifted reference curve
    shifted_reference_curve = reference_curve.copy()
    # Initialized step for shifting the reference curve towards the measured value 
    step = 0.05

    # Iteratively shifting the reference curve and calculate the sum of unfavorable deviations
    while True:
        unfavorable_deviations_sum = 0
        for octave, abs_value in abs_full_octave.items():
            if not np.isnan(abs_value):
                measured_value = shifted_reference_curve.get(octave, np.nan)
                if not np.isnan(measured_value):
                    deviation = max(0, measured_value - abs_value)
                    unfavorable_deviations_sum += deviation

        # Checking if sum of unfavorable deviations is less than or equal to 0.1
        if unfavorable_deviations_sum <= 0.1:
            break

        # Shifting the entire reference curve down uniformly by the step size
        for octave in shifted_reference_curve.keys():
            if not np.isnan(shifted_reference_curve[octave]):
                shifted_reference_curve[octave] -= step

    # Weighted sound absorption
    abs_w = shifted_reference_curve[500]
    abs_w = round(abs_w,2)

    # Shape indicators
    shape_indicator = ""
    if abs_full_octave.get(250, np.nan) - shifted_reference_curve[250] > 0.25:
        shape_indicator += "(L)"
    if abs_full_octave.get(500, np.nan) - shifted_reference_curve[500] > 0.25 or \
       abs_full_octave.get(1000, np.nan) - shifted_reference_curve[1000] > 0.25:
        shape_indicator += "(M)"

    if shape_indicator:
        additional_note = "Additional note: It is strongly recommended to use this single-number rating in combination with the complete sound absorption coefficient curve."
    else:
        additional_note = ""
    
    return abs_w, shape_indicator, shifted_reference_curve, additional_note


def determine_absorption_class(abs_w):
    if 0.90 <= abs_w <= 1.00:
        return "A"
    elif 0.80 <= abs_w < 0.90:
        return "B"
    elif 0.60 <= abs_w < 0.80:
        return "C"
    elif 0.30 <= abs_w < 0.60:
        return "D"
    elif 0.00 <= abs_w < 0.30:
        return "E"

##################################################################################################################

# FIGURES

def fig_scheme_opening():
    # Scheme Attributes
    image_paths = {
        1: "Scheme Attributes\scheme_frontwaves.png",
        2: "Scheme Attributes\scheme_foliage.png",
        3: "Scheme Attributes\scheme_substrate.png",
        4: "Scheme Attributes\scheme_cavity.png",
        5: "Scheme Attributes\scheme_interimwaves.png",
        6: "Scheme Attributes\scheme_rigidbacking.png"
    }
    # Loading attributes
    images = {num: plt.imread(path) for num, path in image_paths.items()}
    # Defining the sequence
    sequence = [1, 2, 5, 3, 5, 4, 6]
    fig, axes = plt.subplots(1, len(sequence), figsize=(4, 2))
    fig.subplots_adjust(wspace=0.01)
    for ax, num in zip(axes, sequence):
        ax.imshow(images[num])
        ax.axis('off')
    return fig


def fig_scheme():
    # Scheme Attributes
    image_paths = {
        1: "Scheme Attributes\scheme_frontwaves.png",
        2: "Scheme Attributes\scheme_foliage.png",
        3: "Scheme Attributes\scheme_substrate.png",
        4: "Scheme Attributes\scheme_cavity.png",
        5: "Scheme Attributes\scheme_interimwaves.png",
        6: "Scheme Attributes\scheme_rigidbacking.png"
    }
    # Loading attributes
    images = {num: plt.imread(path) for num, path in image_paths.items()}
    # Image sequences for different conditions:
    # plant_active, substrate_active, cavity_active
    sequences = {
        (True, True, True): [1, 2, 5, 3, 5, 4, 6],
        (False, True, True): [1, 3, 5, 4, 6],
        (True, False, True): [1, 2, 5, 4, 6],
        (True, True, False): [1, 2, 5, 3, 6],
        (True, False, False): [1, 2, 6],
        (False, False, True): [1, 4, 6],
        (False, True, False): [1, 3, 6],
    }
    # Checking if condition exists in sequences
    condition = (plant_active, substrate_active, cavity_active)
    if condition in sequences:
        # Sequence for the current condition
        sequence = sequences[condition]
        fig, axes = plt.subplots(1, len(sequence), figsize=(4, 2))
        fig.subplots_adjust(wspace=0.01)
        for ax, num in zip(axes, sequence):
            ax.imshow(images[num])
            ax.axis('off')
    return fig


def fig_RL(RL_values, frequencies):
    fig1 = plt.figure(figsize=(6, 4))  # Standard ISO size adjustment
    ax = fig1.add_subplot(1, 1, 1)
    ax.plot(frequencies, RL_values, linestyle='-', color='black')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('RL [dB]')
    ax.set_title('Frequency-dependent return loss')
    ax.grid(True)
    return fig1


def fig_abs(abs_third_octave):
    frequencies = list(abs_third_octave.keys())
    abs_values = list(abs_third_octave.values())
    fig2 = plt.figure(figsize=(6, 4))  # Standard ISO size adjustment
    ax = fig2.add_subplot(1, 1, 1)
    ax.plot(frequencies, abs_values, linestyle='-', color='black', marker='o')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(r'$\alpha_s$ [-]')
    ax.set_ylim(0,1)
    ax.set_title('One-third-octave band sound absorption coefficients')
    ax.grid(True)
    return fig2


def fig_ISO11654(abs_full_octave, shifted_reference_curve):
    octaves = sorted(abs_full_octave.keys())
    abs_values = [abs_full_octave[octave] for octave in octaves]
    ref_values = [shifted_reference_curve[octave] if not np.isnan(shifted_reference_curve[octave]) else None for octave in octaves]
    fig3, ax = plt.subplots(figsize=(5, 5))
    ax.plot(octaves, abs_values, 'bo-', label='Derived Absorption')
    ax.plot(octaves, ref_values, 'ro-', label='Shifted Reference Curve')
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel(r'$\alpha_p$ [-]')
    ax.set_title('Practical sound absorption coefficients and the reference curve according to ISO 11654')
    ax.grid(True)
    ax.legend()
    plt.close(fig3)
    return fig3


##################################################################################################################

# OUTPUT GUI

def output_interface(fig1, fig2, abs_w, shape_indicator, additional_note, absorption_class):
    # New top-level GUI window for displaying simulation output
    root = tk.Tk()
    root.title("Simulation Output")
    root.geometry("1200x800")
    root.configure(bg='#f0f0f0')
    # Grid layout with specific thicknesses for columns
    root.grid_columnconfigure(0, weight=1)  # Set thickness for column 0
    root.grid_columnconfigure(1, weight=1)  # Set thickness for column 1
    root.grid_rowconfigure(0, weight=1)       # Configure row 0 to expand
    root.grid_rowconfigure(1, weight=1)       # Configure row 1 to expand

    # Working with the first frame
    fig = fig_scheme()
    frame1 = ttk.Frame(root, relief=tk.SUNKEN, borderwidth=2)
    frame1.grid(row=0, column=0, padx=5, pady=5, sticky="nsew")
    frame1.grid_columnconfigure(0, weight=1)
    frame1.grid_columnconfigure(1, weight=1)
    frame1.grid_rowconfigure(0, weight=1)
    frame1.grid_rowconfigure(1, weight=0)     
    # Adding figure to the frame
    canvas_fig = FigureCanvasTkAgg(fig, master=frame1)
    canvas_fig.draw()
    canvas_fig.get_tk_widget().grid(row=0, column=0, sticky='nsew')
    # Adding heading below the figure
    heading = ttk.Label(frame1, text="Configured Vertical Greenery System - Transfer Matrix Scheme", font=("Helvetica", 10, "italic"), background="white")
    heading.grid(row=1, column=0, columnspan=2, sticky="nsew")
    # Adding parameters
    text_1 = ""
    text1 = ""
    text_2 = ""
    text_22 = ""
    text2 = ""
    text3 = ""
    if plant_active == True:
        text_1 = f"Foliage Layer: {plant}"
        text1 = f"d (plant) = {int(d_plant*100)} cm"
    if substrate_active == True:
        text_2 = f"Substrate Layer: {substrate}"
        text_22 = f"Substrate State: {state}"
        text2 = f"d (substrate) = {int(d_substrate*100)} cm"
    if cavity_active == True:
        text3 = f"d (cavity) = {int(d_cavity*100)} cm"
    components_text =  text_1 + "\n" + text_2 + "\n" + text_22 + "\n" + "\n"
    parameters_text = text1 + "\n" + text2 + "\n" + text3
    text = components_text + parameters_text 
    paramters = ttk.Label(frame1, text=text, font=("Helvetica", 10, "italic"), background="white")
    paramters.grid(row=0, column=1, rowspan=2, sticky="nsew")

   # Working with second frame
    frame2 = ttk.Frame(root, relief=tk.SUNKEN, borderwidth=2)
    frame2.grid(row=0, column=1, padx=5, pady=5, rowspan=2, sticky='nsew')
    frame2.grid_columnconfigure(0, weight=1)
    frame2.grid_rowconfigure(0, weight=1) 
    frame2.grid_rowconfigure(2, weight=2) 
    # Adding a heading
    heading = ttk.Label(frame2, text="Parameters According to ISO 11654", font=("Helvetica", 14, "bold"), background="white", style="Heading.TLabel")
    heading.grid(row=0, column=0, sticky="nsew")
    # Adding explanatory text with some formatting
    text4 = f"Weighted Sound Absorption Coefficient: {abs_w}{shape_indicator}"
    text5 = f"Absorption Class: {absorption_class}"
    text6 = f"{additional_note}"
    # Combining text into one block
    results_text = f"{text4}\n\n{text5}\n\n{text6}"
    # Creating text widget for formatted text
    text_widget = tk.Text(frame2, wrap="word", width=50, height=10, background="white", relief=tk.FLAT)
    text_widget.insert(tk.END, results_text)
    text_widget.config(state="disabled", font=("Helvetica", 12), padx=10, pady=10)
    text_widget.grid(row=1, column=0, sticky='nsew')
    # Adding figure 2
    canvas_fig2 = FigureCanvasTkAgg(fig2, master=frame2)
    canvas_fig2.draw()
    canvas_fig2.get_tk_widget().grid(row=2, column=0, sticky='nsew')

    # Working with third frame
    frame3 = ttk.Frame(root, relief=tk.SUNKEN, borderwidth=2)
    frame3.grid(row=1, column=0, padx=5, pady=5, sticky='nsew')
    canvas_fig1 = FigureCanvasTkAgg(fig1, master=frame3)
    canvas_fig1.draw()
    canvas_fig1.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Running the application
    root.mainloop()

##################################################################################################################

# CALLING THE FUNCTIONS

# Input GUI
input_interface()
check()

# Loading efficient values
kc_filepaths, zc_filepaths = extract_filepaths(folderpath, substrate, plant)
kc_values_substrate, kc_values_plant, zc_values_substrate, zc_values_plant = load_eff_values(kc_filepaths, zc_filepaths, state)

# Deriving total transfer matrix of the system
T_values_cavity, T_values_substrate, T_values_plant = calc_T_for_d(z_0, c_0, substrate_active, plant_active, frequencies, d_substrate, d_plant, kc_values_substrate, kc_values_plant, zc_values_substrate, zc_values_plant)
T_values = T_total(T_values_cavity, T_values_substrate, T_values_plant)

# Calcualting return loss of the system
R_values, RL_values = RL_total(T_values, z_0)

# Calcualting absoprtion of the system
abs_values = abs_total(R_values)
abs_third_octave = abs_to_third_octave_bands(frequencies, abs_values)
abs_full_octave = practical_abs(abs_third_octave)
abs_w, shape_indicator, shifted_reference_curve, additional_note = abs_weighted(abs_full_octave)
absorption_class = determine_absorption_class(abs_w)

# Output GUI
fig1 = fig_RL(RL_values, frequencies)
fig2 = fig_abs(abs_third_octave)
output_interface(fig1, fig2, abs_w, shape_indicator, additional_note, absorption_class)