Here one can find python scripts, used for sorting the world sample of Trichophyton mentagrophytes and T. interdigitale ITS region sequences.

The script puttorights.py selects nucleotide sequences from infile.txt (Genbank format), containing exact matches with test element. The element should be modified directly in the code.

The script compactor.py merges identical sequences into groups (genotypes). Its input is in FASTA format. The sequences in the input must be cropped.

The script attributetogenotypes.py compares two fasta files: the dataset "sample.fasta" and file "referenceset.fasta" with reference sequences. The result is the file "attributed.txt", where the first column contains names from "sample.fasta" and the second -- a name from reference dataset, if found.
