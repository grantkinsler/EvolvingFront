This directory contains the code and data for Kinsler et al 2024:
"A high-resolution two-step evolution experiment in yeast reveals a shift from pleiotropic to modular adaptation"

## INSTALLATION

To run the scripts and notebooks contained here, it's best to use a virtual environment. 

You can install all the same packages/versions used here via:

pip3 install -r venv_requirements.txt

## PROCESSING DATA

All code used to process data is contained in the directory code/processing.

Code for processing data Fitness Measurement experiments is contained in the directory code/processing/FitnessMeasurements. It can be run using the following workflow:

(1) Map raw fastq files to barcode counts using BarcodeCounter2 the files and scripts used are contained in the SequenceToCount directory.
(2) Process barcode counts to fitness estimates

Code for processing mutation data is in the directory code/processing/WGS

Post-processing was done by running the wrapper frequency_trajectories.ipynb which also contains some analysis as to quality of fitness measurment experiments.

fitness_analysis.ipynb contains analysis of replicate-replicate correlations and frequency-dependent fitness effects. 

organizing_mutations_and_fitness.ipynb combines the fitness and mutation analysis and calculates performances to generate the key data table (data/fitness_withMutations.csv) from which most of the analysis is derived.

## DATA ANALYSIS

Code for the majority of the data analysis is contained in the directory code/analysis.

1. 2D fitness vs stationary performance.ipynb
Contains code to generate Figure 4 and analyses involving stationary phase performance.

2. Basic Properties of Isolated Mutants.ipynb
Contains code to generate panels in Figure 2 that depict fitness effects of identified mutations.
This also contains basic analyses and calculations of mutational properties.

3. calculating performances.ipynb
Contains code to generate the performance calculation example in Figure 3A.

4. Evo1D_mutant_analysis.ipynb
Contains code for generating Evo1D analysis figure S2. 

5. Exploring Tradeoffs Between Growth Phases.ipynb
Contains code for generating the majority of the tradeoff figures, including Figures 3, 5, S3, S4, S5, S7, S8, S9, and S10. 


### figures/
	analysis/ 
		Contains figures from various analysis, including panels used in main figures
	main_text/
		Contains figures in the paper and keynote file used to assemble panels.


### data/ 
	Contains data tables used in analysis as well as supplemental data files uploaded to the journal.


### EvolutionTracking/
	
	This directory contains code for processing of evolution tracking data, calculating DFEs, and comparing between fitness measurement experiments and evolution fitness.

	processing/
		Contains code used to process count data into barcode counts

	FitMut2_processing/
		Contains FitMut output

	DFE inference.ipynb
		This is the script used to calculate the DFE shown in Figure 2. 

	remeasurement vs evolution fitness.ipynb
		This is the analysis script for comparing between evolution and fitness measurement fitness.
		This code generates Figure S6.