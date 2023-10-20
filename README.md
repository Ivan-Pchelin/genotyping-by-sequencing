Here one can find python scripts, used for sorting the world sample of *Trichophyton mentagrophytes* and *T. interdigitale* ITS region sequences.

The script **puttorights.py** selects nucleotide sequences from infile.txt (GenBank format), containing exact matches with test element. The element should be modified directly in the code.

The script **compactor.py** merges identical sequences into groups (genotypes). Its input is in FASTA format. The sequences in the input must be cropped.

The script **attributetogenotypes.py** compares two fasta files: the dataset "sequence.fasta" and the file "referenceset.fasta" with reference sequences. The result will appear in the file "attributed.txt", where the first column contains names from "sample.fasta" and the second one contains a name from the reference dataset, if found. 
The code and the *T. mentagrophytes / T. interdigitale* reference dataset of ribosomal ITS region sequences were described by [Taghipour et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31444823/), the corrected genotyping approach and the updated dataset were published by [Nikkholgh et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37429606/).
