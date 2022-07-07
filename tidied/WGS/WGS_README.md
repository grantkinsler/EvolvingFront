The general procedure is laid out in the iPython notebook, but I’ll also describe it here:
(1) Process files using gatk
    1.1 Run yeast_alignment.sbatch. Align the samples to the reference genome with bwa and so some standard filtering (mark duplicates, read groups, etc.). Note: you need a yeast_alignment_samples.inp with all the sequencing files you want to process and to change the —array=1-X where X is the number of files you are processing
    1.2 Run yeast_mergeGVCFs+callgenotypes.sbatch. This merges all the files together and uses GATK to call genotypes on the full set.
    1.3 Calculate average coverage for each sample, keeping only samples above your favorite sample depth (I kept everything above 20X coverage on average)

(2) Filter the data using “hard filters” in the iPython notebook.
You might want to play with these a bit, but I’ve mostly deviated from GATK default hard filters with mapping quality (which was causing a lot of SNPs to show up that really didn’t need to be there. You can see how many variants get filtered out, etc. in the notebook.

(3) Call “variable sites”. This gets rid of all sites that are the same for all your samples (likely represent differences between your ancestral strain and the reference), sites where there are lots of low quality calls, and mitochondrial SNPs. I then output these to a .csv file.

(4) Manual curation with IGV. Now that we’ve done some filtering, hopefully we’ve reduced our SNPs to a much more reasonable number (I got ~50). Now, I look at IGV to see if these SNPs are real (or not). Most of mine were clearly in repeat regions (might think about trying to figure out a way to do this programmatically rather than manually in the future). I then mark which SNPs are true calls remaining (reduces to ~20 SNPs). 

(5) Annotate SNPs to effects with snpEff. This gives us the SNPs in more interpretable effects (e.g. IRA1 nonsense_variant). Note that this gives only the most likely strongest effect associated with each of these SNPs instead of a more complete list.

(6) Output to a .csv file with SNPs called per sample. 

The rest that follows is a start on my copy number variant copying. Essentially, I’m looking at deviations in coverage for regions across the genome. I’ll try to remember to let you know once I’m done completing this stuff.



I’ve also attached a couple scripts that are used to extract barcodes from the WGS data. 

extract_barcodes.sbatch - This runs extract_barcodes.py using the same yeast_alignment_samples.inp file as the WGS scripts.
extract_barcodes.py - This extracts the barcode region using grep and uses regex to extract the barcodes, removing reads that might be amplicon contamination based on the location of the barcode in the read
extract_barcodes_assembly.py - This compiles all the extracted barcodes into a single file for ease of use. 
