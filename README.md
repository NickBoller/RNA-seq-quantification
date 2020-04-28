# CMSC423_rna_quant
Implemented two versions of the Expectation Maximization algorithm for transcript-level quantification from RNA-seq data.

## How to run
**Run Full Model EM**
- python squant.py --out PATH-TO-SAVE-OUTPUT --i PATH-OF-ALIGNMENTS.CS423
	- ie: python squant.py --out quant.tsv --i data/alignments.cs423

**Run Equivalence Class EM**
- python squant.py --eqc True --out PATH-TO-SAVE-OUTPUT --i PATH-OF-ALIGNMENTS.CS423
	- ie: python squant.py --eqc True --out eqc_quant.tsv --i data/alignments.cs423

**Run Calculate Statistics**
- python squant.py --stats True --em_estimation_file PATH-TO-EM-OUTPUT-TSV-FILE --true_counts_file PATH-OF-TRUE_COUNTS.TSV
	- ie: python squant.py --stats True --em_estimation_file quant.tsv --true_counts_file data/true_counts.tsv

**Libraries Used**
- argparse
	- pip install argparse
- numpy
	- pip install numpy
- scipy
	- pip install scipy
- numba
	- pip install numba
