Thank you for visiting this page. Here one can find python scripts, used for sorting the world sample of *Trichophyton mentagrophytes* and *T. interdigitale* ITS region sequences.

The single most important file here may be **referenceset.fasta** with a collection of published genotypes. If you found even another brand new *Trichophyton* ITS region sequence and have an intent to publish it with a fancy Latin number, **please do** accomplish the following:
(1) Ensure that the differentiating nucleotide substitution(s) are within the boundaries of ITS region. See for example [Nikkholgh et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37429606/) or another typing paper.
(2) Manually check the substitution(s) on both forward and reverse electrophoregrams.
(3) Perform a literature search to see if someone have already reserved that Latin number for another sequence.

The script **puttorights.py** selects nucleotide sequences from infile.txt (GenBank format), containing exact matches with test element. The element should be modified directly in the code.

The script **compactor.py** merges identical sequences into groups (genotypes). Its input is in FASTA format. The sequences in the input must be cropped.

The script **attributetogenotypes.py** compares two fasta files: the dataset "sequence.fasta" and the file "referenceset.fasta" with reference sequences. The result will appear in the file "attributed.txt", where the first column contains names from "sample.fasta" and the second one contains a name from the reference dataset, if found. 
The code and the *T. mentagrophytes / T. interdigitale* reference dataset of ribosomal ITS region sequences were described by [Taghipour et al. 2019](https://pubmed.ncbi.nlm.nih.gov/31444823/), the corrected genotyping approach and the updated dataset were published by Nikkholgh et al.

I wrote the scripts during my early days in science, and there are far better ways of accomplishing the mentioned tasks. For example, extracting non-typable sequences from a dataset "sequence.fasta" can be done with [SeqKit](https://github.com/shenwei356/seqkit.git) using the following commands:

```
conda activate seqkit
grep -v ">" referenceset.fasta > clean_sequences.txt
seqkit grep -s -v -f clean_sequences.txt sequence.fasta -o unmatched.fasta
```
If you have some suggestions or questions, please feel free to contact me. Good luck in research!
