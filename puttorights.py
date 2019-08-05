# The script selects sequences from infile.txt (Genbank format), containing
# exact matches with test element. The element should be modified directly
# in the code.

# Version 2019-01-19

# This defines test element to select sequences from the infile 
testseq = 'TTCCGGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCCGTGAATCATCGAATCTTTGAACGCACATTGCGCC'

import re
import string
import sys

def complement(sequence):
	sequence = sequence.upper()
	c_sequence = sequence[::-1]
	c_sequence = c_sequence.replace('A', 'F').replace('T', 'A').replace('F', 'T')
	c_sequence = c_sequence.replace('G', 'F').replace('C', 'G').replace('F', 'C')
	c_sequence = c_sequence.replace('M', 'F').replace('K', 'M').replace('F', 'K')
	c_sequence = c_sequence.replace('R', 'F').replace('Y', 'R').replace('F', 'Y')
	c_sequence = c_sequence.replace('V', 'F').replace('B', 'V').replace('F', 'B')
	c_sequence = c_sequence.replace('H', 'F').replace('D', 'H').replace('F', 'D')
	return (c_sequence)
c_testseq = complement(testseq)
testseq = testseq.upper()

# Этот элемент принимает соглашение с условием отбора последовательностей
# или останавливает выполнение скрипта
print('\n')
print('Sequences with the following element will be included in the outfile:')
print('\n', testseq)
print('\n', 'y/n')
reply = str(input())
if reply.strip() == 'n':
	sys.exit()
elif 'y' not in reply:
	while 'n' or 'y' not in reply:
		print(' y or n can be read')
		reply = str(input())
		if reply.strip() == 'n':
			sys.exit()
		if reply.strip() == 'y':
			break

# Этот блок считает количество записей в исходном файле и проверяет формат Genbank
n_Genbank_w = {} # инициация пустого словаря
text = open('infile.txt', 'r')
for word in re.findall(r'LOCUS', text.read()):
	count = n_Genbank_w.get('LOCUS', 0)
	n_Genbank_w[word] = count + 1
if n_Genbank_w != {}:
	ne = n_Genbank_w['LOCUS']
else:
	print('Error: cannot recognize file format')
	sys.exit()

# Эта часть выписывает в множество accnos все номера из infile
with open('infile.txt') as inf:
	accnos = set()
	lines = inf.readlines()
	i = 0
	while i < len(lines):
		if 'LOCUS' in lines[i]:
			accno = re.findall(r'\S+', lines[i])[1]
			accnos.add(accno)
		i += 1

# Выпишем в файл clean.fasta все уникальные последовательности,
# попутно исправляя комплементарные.
transferredseqs = 0
excess = 0
waisted = 0
trash = set()
emptyaccnos = []
with open('infile.txt') as inf:
	with open('clean.fasta', 'w') as ouf:
		line = inf.readline()
		while line:
			if 'LOCUS' in line[:5]:
				accno = re.findall(r'\S+', line)[1]
				skip = 0
				if accno in accnos:
					accnos.remove(accno)
					applicant = ''
					while 'ORIGIN' not in line[:6] and '//' not in line[:2]:
						line = inf.readline()
					if '//' in line:
						print('Empty entry:', accno)
						emptyaccnos.append(accno)
						skip = 1
					if skip == 0:
						while line.strip() != '//':
							line = inf.readline()
							letts = line[10:100]
							letts = letts.replace(' ', '')
							letts = letts.strip()
							applicant = applicant + letts
						applicant = applicant.upper()
						if testseq in applicant:
							ouf.write('>')
							ouf.write(accno)
							ouf.write('\n')
							ouf.write(applicant)
							ouf.write('\n')
							transferredseqs += 1
						elif c_testseq in applicant:
							ouf.write('>')
							ouf.write(accno)
							ouf.write('\n')
							applicant = complement(applicant)
							ouf.write(applicant)
							ouf.write('\n')
							transferredseqs += 1
						else:
							waisted += 1
							trash.add(accno)
				else:
					excess += 1
			line = inf.readline()

# Создадим файл-помойку и выпишем туда отброшенные последовательности
if len(trash) != 0:
	with open('infile.txt') as inf:
		with open('abandoned.fasta', 'w') as ouf:
			for line in inf:
				if 'LOCUS' in line[:5]:
					accno = line[12:21].strip()
					if accno in trash:
						trash.remove(accno)
						ouf.write('>')
						ouf.write(accno)
						ouf.write('\n')
						while 'ORIGIN' not in line[:6]:
							line = inf.readline()
						while line.strip() != '//':
							line = inf.readline()
							letts = line[10:100]
							letts = letts.replace(' ', '')
							ouf.write(letts)
				line = line.strip()

# Выведем результаты подсчетов
print ('\n', 'General number of entries =', ne)
if excess == 0:
	print(' Single copies throughout the sample', '\n')
else:
	print(' Initial number of unique entries =', ne - excess)
	print(' A total of', excess, 'excess sequences', '\n')
	# В исходнике ожидается не больше двух копий каждой записи
	if transferredseqs + waisted + excess != ne:
		print(' Unexpected number of copies in original file')
print(' Transferred to clean.fasta', transferredseqs, 'sequences')
print ('', waisted, 'sequences waisted')
if len(accnos) > 0:
	print(accnos, 'wrong count of accessions')
print('', len(emptyaccnos), 'empty entries')