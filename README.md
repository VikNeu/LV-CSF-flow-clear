# Cerebrospinal fluid flow and clearance driven by lateral ventricle volume oscillations

This repository stands for this project:

Cerebrospinal fluid flow and clearance driven by lateral ventricle volume oscillations
Viktor Neumaier, Moritz Bonhoeffer, Melissa Thalhammer, Julia Schulz, Gabriel Hoffmann, Juliana Zimmermann, Felix Brandl, Matthias Brendel, Igor Yakushev, Josef Priller, Jan Kirschke, Christine Preibisch, Benedikt Zott, Christian Sorg
bioRxiv 2025.10.20.683529; doi: https://doi.org/10.1101/2025.10.20.683529

## Abstract

In the human brain, substance clearance is intimately connected to the cerebrospinal fluid (CSF) and its flow. CSF extends from the lateral ventricles (LVs) to the parenchyma’s perivascular spaces. Macroscopic undulating CSF flow is present during both wakefulness and sleep and can be experimentally induced. However, the mechanisms generating this flow and its contribution to brain clearance remain unclear. Using fMRI and PET across various conditions, we demonstrate that LV-volume oscillations drive undulating CSF flow in the ventricles and subarachnoid basal cisternae. LV oscillations are driven by cortical blood-volume changes induced by neuronal activity, heartbeat and respiration. LV oscillations’ amplitudes determine PET-tracer clearance from the LVs. Conclusively, induced by extra- and intracranial physiological drivers and mediated by cortical blood-volume changes, LV-volume oscillations drive macroscopic CSF flow and clearance.

## Installation

This repository contains the code to reproduce the analyses and figures of the project.  
The stats code is mainly written in Python (and one script in MATLAB).
The analysis code starting from the preprocessed data is written in MATLAB

The instructions below assume a Unix-like system (macOS / Linux).  
On Windows, commands may need minor adjustments.

### Software requirements

- Python >= 3.13.2  
- MATLAB >= R2021b (for mediation analysis)

### Installation steps

1. **Create a directory and clone the repository**

   ```bash
   mkdir /path/to/your/setup && cd /path/to/your/setup
   git clone --recursive https://github.com/VikNeu/LV-CSF-flow-clear


## Environment setup

### Python environment

2. **Create a virtual environment with **Python 3.12.3** and activate it**

   ```bash
   python3.12 -m venv .venv 
   source .venv/bin/activate

3. **Install Python dependencies**

   ```bash
   pip install -r requirements.txt

### Install Toolbox for MATLAB
4. **Mediation Toolbox**
  ```bash
  git clone https://github.com/canlab/MediationToolbox code/MediationToolbox
  ```


## Usage (Note: Still WIP)

### General

The `/analysis` scripts do contain the MATLAB scripts to obtain the metrics analysed in `/stats`, starting with the preprocessed images.
`/data` contains all the relevant data tables used in `/stats`.
`/stats` contains all the jupyter notebooks used for creating plots and statistics - The respective scripts are organized figurewise and should be ready to run.

