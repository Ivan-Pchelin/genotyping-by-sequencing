# This script compares two fasta files: the dataset "sequence.fasta" and file
# "referenceset.fasta" with reference sequences. The result is the file
# "attributed.txt", where the first column contains names
# from "sample.fasta" and the second -- a name from reference dataset, if found

# Version 2020-01-09

import re
import os
current_dir = os.getcwd()
files = [f.lower() for f in os.listdir(current_dir)]
if 'sequence.fasta' in files:
    sourcefile = 'sequence.fasta'
elif 'sequence.fst' in files:
    sourcefile = 'sequence.fst'
elif 'sequence.fa' in files:
    sourcefile = 'sequence.fst'
elif 'sequence.txt' in files:
    sourcefile = 'sequence.txt'

with open(sourcefile) as inf:
    with open('neat.fasta', 'w') as ouf:
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
                line = line.replace('\n', '').replace(' ', '')
                applicant = applicant + line
            line = inf.readline()
        if not line:
            applicant = applicant.replace('-', '').upper()
            ouf.write(applicant)

# Этот блок переписывает референсные последовательности в одну линию
with open('referenceset.fasta') as inf:
    with open('referenceset_plain.fasta', 'w') as ouf:
        applicant = ''
        k = 0
        line = inf.readline()
        while line:
            if '>' in line:
                applicant = applicant.replace('-','').upper()
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
            applicant = applicant.replace('-','').upper()
            ouf.write(applicant)

# Этот блок ищет референсные последовательности
# из referenceset_plain.fasta в файле neat.fasta.
with open('referenceset_plain.fasta') as inf:
    references = inf.readlines()

generalcount = 0
poolofnames = set()
accnos = set()
with open('neat.fasta') as inf:
    with open('attributed.txt', 'w') as ouf:
        line = inf.readline().strip()
        while line:
            if line[0] == '>':
                generalcount += 1
                elements = re.findall(r'\S+', line)
                if generalcount > 1:
                    ouf.write('\n')
                if (elements[0])[1:] in poolofnames:
                    print('Multiple entries of', (elements[0])[1:])
                ouf.write(elements[0][1:]) # номер изолята
                poolofnames.add(elements[0][1:])
            else:
                runner = 0
                while runner < len(references):
                    if references[runner].strip() in line:
                        ouf.write('\t')
                        match = re.search(r'^\S+', references[runner-1])#
                        accno = references[runner-1][1:match.end()]
                        ouf.write(accno)
                        accnos.add(accno)
                    runner += 1
            line = inf.readline().strip()

# Найдем случаи, в которых разные последовательности, приписанные
# одному и тому же референсу, различаются.
# Варианты выпишем в variants.fasta
with open('neat.fasta') as inf:
    rawseqs = inf.readlines()

switch = 0
with open('attributed.txt') as inf:
    with open('variants.fasta', 'w') as ouf:
        line = inf.readline().strip()
        while line:
            accnu = re.findall(r'\S+', line)
            accno = accnu[-1]
            if accno in accnos:
                caption = '>' + accno + '-'
                accnos.remove(accno)
                runner = 0
                while runner < len(references):
                    if accno in references[runner]:
                        currentseq = references[runner+1].strip()
                    runner += 1
                runner = 0
                currentvariants = set()
                while runner < len(rawseqs):
                    if currentseq in rawseqs[runner]:
                        currentvariants.add(rawseqs[runner])
                    runner += 1
                runner = 1
                if len(currentvariants) > 1:
                    switch = 1
                    for i in currentvariants:
                        ouf.write(caption)
                        ouf.write(str(runner))
                        ouf.write('\n')
                        ouf.write(i)
                        runner += 1
            line = inf.readline()

print('\n', generalcount, 'sequences analyzed')

# Удалим временные файлы
os.remove('referenceset_plain.fasta')
#os.remove('neat.fasta')
if switch == 0:
    os.remove('variants.fasta')