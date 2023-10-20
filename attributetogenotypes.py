'''
This script compares two FASTA files: a dataset "sequence.fasta"
and the file "referenceset.fasta" with reference sequences.
The result is the file "attributed.txt", where the names from "sequence.fasta"
are supplemented with relevant accessions from "referenceset.fasta".
If there is no perfect match between analyzed sequence and one
of the reference sequences, no result is indicated.
To launch the script, one should have installed Python interpreter
from www.python.org/downloads/. Windows users can place the script
with the files sequence.fasta and referenceset.fasta in the same folder,
then type cmd in the address bar. In the black window appeared one
should type "python attributetogenotypes.py", without the quotes.
Version 2023-10-20 by Ivan Pchelin, arcella.oraia@gmail.com.
If I can be of assistance, please do not hesitate to contact me.
'''

import re
import os
import sys

temporaryfiles=set()

# this re-writes FASTA sequences in single lines, capitalizes them
# and removes spaces
def irontonew (file):
    with open(file) as inf:
        with open('neat_'+file, 'w') as ouf:
            applicant = ''
            k = 0
            line = inf.readline()
            while line:
                if '>' in line:
                    applicant = applicant.replace('-', '').upper()
                    ouf.write(applicant)
                    if k != 0:
                        ouf.write('\n')
                    k += 1
                    applicant = ''
                    ouf.write(line)
                else:
                    line = line.replace('\n', '')
                    applicant = applicant + line
                line = inf.readline()
            if not line:
                applicant = applicant.replace('-', '').upper()
                ouf.write(applicant)
    temporaryfiles.add('neat_'+file)

irontonew('referenceset.fasta')
irontonew('sequence.fasta')

analyzednames = [] # a list with names from analyzed file to check whether for
# all sequences types were available 
analyzednamesSET = set() # used in the check for duplicated entries
sequences = []
with open ('neat_sequence.fasta') as inf:
    lines = inf.readlines()
    runner = 0
    while runner < len(lines):
        c = []
        analyzednames.append([re.findall(r'\S+', lines[runner])[0]][0])
        analyzednamesSET.add([re.findall(r'\S+', lines[runner])[0]][0])
        c.append([re.findall(r'\S+', lines[runner])[0]])
        c.append([re.findall(r'\S+', lines[runner+1])[0]])
        sequences.append(c)
        runner += 2

nnames = len(analyzednamesSET)
qwerty = list(analyzednames)

# check for duplicated entries
for i in analyzednamesSET:
    qwerty.remove(i)
if len(qwerty) > 0:
    qwerty = set(qwerty)
    print('\n', 'In the analyzed file, there are multiple entries of:')
    for i in qwerty:
        print('', i[1:])
    print ('\n', '------  Please fix it  ------')
    for i in temporaryfiles:
        os.remove(i)
    sys.exit()


# create dictionaries with reference and analyzed sequences
references = []
with open ('neat_referenceset.fasta') as inf:
    lines = inf.readlines()
    runner = 0
    while runner < len(lines):
        c = []
        try:
            c.append([re.findall(r'\S+', lines[runner])[0]])
            c.append([re.findall(r'\S+', lines[runner+1])[0]])
            references.append(c)
            runner += 2
        except:
            pass
            runner += 1

listwithvariants = []
generalcount = 0

# write the file with genotyping results
with open ('attributed.txt', 'w') as ouf:
    for i in references:
        variantsofcurrentgenotype = set()
        for j in sequences:
            if i[1][0] in j[1][0]:
                generalcount += 1
                variantsofcurrentgenotype.add(j[1][0])
                ouf.write (j[0][0][1:])
                ouf.write ('\t')
                ouf.write (i[0][0][1:])
                ouf.write ('\n')
        if len(variantsofcurrentgenotype) > 1:
            number = 1
            for k in variantsofcurrentgenotype:
                listwithvariants.append(i[0][0]+'-'+'Variant'+str(number))
                listwithvariants.append(k)
                number+=1

# write the file with sequence variants, if any
if listwithvariants != []:
    with open('variants.fasta', 'w') as ouf:
        for i in listwithvariants:
            ouf.write(i)
            ouf.write('\n')

print('\n', nnames, 'sequences analyzed')
print('', generalcount, 'sequences with assigned genotypes')

# delete temporary files
for i in temporaryfiles:
    os.remove(i)
