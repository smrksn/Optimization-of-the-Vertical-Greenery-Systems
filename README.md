# Optimization of the Vertical Greenery Systems for Noise Mitigation in Urban Environments
This repository contains products developed as part of the Master Thesis at TU Delft.

<img src="https://github.com/smrksn/Optimization-of-the-Vertical-Greenery-Systems/assets/144154829/273f9bb7-db9d-4a66-aee3-9aae7c7b3af0" width=15% height=15%>

Copyright © Sofia Markson, 2024

## Overview:

Vertical Greenery Systems (VGS) offer effective solutions for reducing urban noise pollution through strategic implementation. This project focuses on optimizing multi-layered VGS to maximize their noise mitigation capabilities. A pivotal advancement involves developing a tool based on the Transfer Matrix Method (TMM), a mathematical technique for analyzing wave transmission and reflection in layered structures. In this approach, VGS layers are modeled as equivalent one-dimensional fluid-like layers with complex and frequency-dependent effective parameters.

The tool includes a materials library that allows the configuration of various VGS setups, incorporating components such as foliage and substrate. The main aim of the tool is to provide information on the noise reduction potential of configured VGS designs. Key performance metrics that aid stakeholders in decision-making include frequency-dependent Return Loss, weighted absorption coefficient with sound absorption class. Moreover, the tool provides frequency-dependent absorption coefficients that can be used as inputs for further urban simulations. These metrics provide information for normal incidence performance only!

![Tool execution](https://github.com/smrksn/Optimization-of-the-Vertical-Greenery-Systems/assets/144154829/9a10b7e8-3606-4b10-a099-7e9107cee726)

## Contents
### Tool
- `Tool.py`: Main script for running the tool.
- `VGS_layers.txt`: Notes on materials and corresponding suppliers.
- `Tool_workflow.html`: Detailed description of tool functions, functionalities, and execution sequence.
- `cache.npz`: Binary file storing variables from the input GUI.
- **Effective_Parameters/**: Folder containing binary files of complex effective parameters for each material.
- **Scheme_Attributes/**: Folder with attributes for compiling the Transfer Matrix Scheme within the tool GUI.

### Data Processing
The scripts provided facilitate the derivation of complex effective parameters from impedance tube measurements (ISO 10534-2). These scripts guide through a structured data processing procedure divided into three phases: preliminary, intermediate, and final. Each phase of the workflow is clearly defined to ensure a systematic approach to deriving accurate results from the measurement data.

## How to Use
### Prerequisites
- Anaconda3-2022.10
- VS Code version 1.90.2 (May 2024 update) or compatible IDE
- Python 3.9.13

### Dependencies
- **numpy**: Numerical computing library for handling arrays and mathematical operations.
- **matplotlib**: Plotting library for creating visualizations like graphs, charts, and plots.
- **pandas**: Data manipulation and analysis library.
- **tkinter**: GUI toolkit for creating user interfaces in Python.

1. **Setup**:
   - Clone the repository and install necessary dependencies.
2. **Navigate**:
   - Open the "Tool\" folder in VS Code or your preferred IDE.
3. **Run Tool.py**:
   - Execute the script `Tool.py`.
   - Adjust the relevant variables as required within the interactive GUI.
   - Wait for the final pop-up GUI displaying performance metrics.
4. **Documentation**:
   - For detailed methodology, refer to the TU Delft repository [link](http://resolver.tudelft.nl/uuid:165f7d97-f6ac-4bc4-be74-eb534aaea012).
