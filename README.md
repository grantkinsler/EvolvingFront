EvolvingFront Project

INSTALLATION

To run the scripts and notebooks contained here, it's best to use a virtual environment. 

You can install all the same packages/versions used here via:

pip3 install -r venv_requirements.txt

PROCESSING DATA

All code used to process data is contained in the directory code/processing.

Code for processing data Fitness Measurement experiments is contained in the directory code/processingFitnessMeasurements. It can be run using the following workflow:

(1) Map raw fastq files to barcode counts using BarcodeCounter2 the files and scripts used are contained in the SequenceToCount directory.
(2) Process barcode counts to fitness estimates

Code for processing mutation data is in the directory code/processing/WGS

Post-processing was done by running the wrapper frequency_trajectories.ipynb which also contains some analysis as to quality of fitness measurment experiments.











DiversityCheck
	To check diversity of transformed barcode pools.

EvolutionTracking
	Barcode sequencing and analysis of the evolution trajectories.

Metagrid
	Barcode and sanger sequencing of locations of barcodes in sorted 96 well plates.

FitnessMeasurements
	Barcode sequencing and analysis of fitness measurement experiments.
