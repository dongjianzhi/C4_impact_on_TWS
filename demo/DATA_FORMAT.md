# Data and code for reproducing the main figures

This folder contains the processed data and MATLAB plotting scripts required to reproduce Figures 1–3 of the manuscript on C4 vegetation and terrestrial water storage.

## Original data

The analyses integrate satellite-based terrestrial water-storage observations, vegetation and land-cover information, and monthly hydroclimatic variables. The original datasets were processed onto a common 0.5° global land grid and aligned to a common monthly time axis before analysis.

The processed data needed to reproduce the main figures are provided in a single MATLAB data file. Users do not need to download or reprocess the original large datasets to run the three plotting scripts in this folder.

Details of the original datasets, preprocessing procedures, analysis periods, and quality-control criteria are provided in the Methods and Data Availability sections of the manuscript.

## Missing values

Missing observations, excluded grid cells, ocean pixels, and locations outside the analysis domain are stored as `NaN`.

## Input data format

The plotting programs require the MATLAB v7.3 file:

- `Main_Figures_Data.mat`

This file contains three structures:

- `Figure1` — processed data used for Figure 1
- `Figure2` — processed data used for Figure 2
- `Figure3` — processed data used for Figure 3

The public variable and field names are source-neutral. Each structure contains the mapped variables, land indices, statistical results, and other derived quantities required by its corresponding plotting program.

## Required files

Keep the following four files in the same folder:

- `Main_Figures_Data.mat`
- `Plot_Figure1_From_Main_Data.m`
- `Plot_Figure2_From_Main_Data.m`
- `Plot_Figure3_From_Main_Data.m`

## File descriptions

### `Main_Figures_Data.mat`

Combined processed input data for all three main figures. The file contains only the structures `Figure1`, `Figure2`, and `Figure3`.

### `Plot_Figure1_From_Main_Data.m`

Reproduces Figure 1, including the spatial pattern of terrestrial water-storage depletion during drought, mean responses for the vegetation groups, and their contributions to global depletion.

### `Plot_Figure2_From_Main_Data.m`

Reproduces Figure 2, including the conditional terrestrial water-storage response, C4 vegetation fraction, regression analysis, predictor importance, and bootstrapped slope distribution.

### `Plot_Figure3_From_Main_Data.m`

Reproduces Figure 3, including changes in C4 vegetation fraction, changes in terrestrial water-storage response, comparisons between decreasing and increasing C4 areas, and the relative C4 response to atmospheric CO2.

## Software requirements

- MATLAB R2022b or later
- Mapping Toolbox
- Statistics and Machine Learning Toolbox
- Image Processing Toolbox
- Optimization Toolbox

The plotting scripts also use the MATLAB colormap functions called by the original figure programs. These functions must be available on the MATLAB path.

## Usage

Place all four required files in the same directory. Open MATLAB, change the current folder to that directory, and run:

```matlab
Plot_Figure1_From_Main_Data
Plot_Figure2_From_Main_Data
Plot_Figure3_From_Main_Data
```

The resulting PNG files are written to the automatically created `Generated_Figures` folder at 300 dpi.

