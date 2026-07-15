# Motor-Imagery-BCI-Analysis

Analysis of EEG recordings acquired during a multi-day Motor Imagery (MI) Brain–Computer Interface (BCI) experiment.

## Project Overview

This repository contains the MATLAB code developed for the Neurorobotics course assignment. The project focuses on the analysis and decoding of EEG signals recorded from healthy participants performing motor imagery tasks.

The implemented pipeline includes:

- EEG preprocessing
- Grand-average analysis across subjects
- Feature extraction and selection
- Subject-specific classifier training
- Offline and online decoding evaluation
- Evidence accumulation for decision making
- Performance assessment and visualization

The goal is to investigate the decoding of three mental states (both hands, both feet, and rest) and evaluate the performance of a subject-specific BCI system on both calibration and online recordings.

## Dataset

The dataset consists of EEG recordings collected from healthy participants during a three-day Motor Imagery BCI experiment.

Each participant completed:

- **Day 1**
  - 3 offline (calibration) runs
  - 2 online runs
- **Days 2–3**
  - 2 online runs per day

EEG was recorded from 16 electrodes at a sampling frequency of **512 Hz**.

> **Note:** The dataset is **not included** in this repository. If interested, please contact me for more info.

## Repository Structure

```text
project/
├── dataset/               # EEG recordings (not included)
├── functions/             # Utility functions
├── *.m                    # Main MATLAB scripts
└── README.md
```

## Instructions for running the scripts
- place the 'dataset' and the (provided) 'functions' folders in the same parent folder

- when in Matlab, enter in the parent folder where you have placed 'dataset' and 'functions'

- add the paths for biosig and eeglab packages

- now you can run the 'main_data_processing.m' script, which exploits the customized function 'process_all_gdf.m' function, which is inside the 'functions' folder; it creates a new folder, 'processed_dataset' that will be used in the following analysis

- run 'main_grand_average.m' to see the results of the Grand Average Analysis

- run 'main_calibration.m' to perform the training on the offline files for each subject; the calibrated decoder and the selected feature of each subject are stored inside its sub-folder in the 'processed_dataset' folder

- run 'main_evaluation.m' to test the model of each subject on the online files; a display classification report is printed in the Command Window

## Requirements

- MATLAB
- Required toolboxes (e.g., Signal Processing Toolbox)

## Reference:

Tonin L et al. The role of the control framework for continuous tele-operation of a BMI driven mobile robot. IEEE Transactions on Robotics, 36(1):78-91, 2020. doi: 10.1109/TRO.2019.2943072

Pfurtscheller G et al. Motor imagery and direct brain-computer communication. Proceedings of the IEEE, 89(7):1123-34, 2001. doi: 10.1109/5.939829

Wolpaw JR et al. Control of a two-dimensional movement signal by a noninvasive brain-computer interface in humans. Proc Natl Acad Sci USA, 101(51):17849-54, 2004.
doi: 10.1073/pnas.0403504101

Leeb R et al. Transferring brain–computer interfaces beyond the laboratory: Successful application control for motor-disabled users. Artificial Intelligence in Medicine, 59(2):121-32, 2013.
doi: 10.1016/j.artmed.2013.08.004

Perdikis S et al. The Cybathlon BCI race: Successful longitudinal mutual learning with two tetraplegic users. PLOS Biology 16(5):e2003787, 2018. doi: 10.1371/journal.pbio.2003787


## Authors

Michele Simoncelli

University of Padova
