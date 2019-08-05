# Version 2019-01-23

import re
import os.path
import sys

# Проверка дублирования номеров доступа
accnos = set()
excess = []
ne = 0
with open ('single.fasta') as inf:
	line = inf.readline()
	while line:
		if '>' in line:
			accno = re.findall(r'^\S+', line)[0][1:]
			if accno in accnos:
				print('Multiple entries of', accno)
				excess.append(accno)
			else:
				ne += 1
			accnos.add(accno)
		line = inf.readline()
print('\n A total of', ne, 'sequences', end='')
if excess != []:
	print(', excluding', len(excess), 'repeats')

# Этот блок спросит, отсевать ли последовательности
# с вырожденными нуклеотидами
print('\n')
print(' Remove sequences with non-standard bases? y/n')
reply = str(input())
if reply.strip() == 'n':
	decision_ns = 0
elif reply.strip() == 'y':
	decision_ns = 1
elif 'y' or 'n' not in reply:
	while 'n' or 'y' not in reply:
		print(' y or n can be read')
		reply = str(input())
		if reply.strip() == 'n':
			decision_ns = 0
			break
		if reply.strip() == 'y':
			decision_ns = 1
			break

# Этот блок спросит, отсевать ли последовательности, депонированные
# только из одной страны
print('\n Remove sequences, deposited less than three times ', end='')
print('and/or from only one country? y/n')
reply = str(input())
if reply.strip() == 'n':
	decision_oc = 0
elif reply.strip() == 'y':
	decision_oc = 1
elif 'y' or 'n' not in reply:
	while 'n' or 'y' not in reply:
		print(' y or n can be read')
		reply = str(input())
		if reply.strip() == 'n':
			decision_oc = 0
			break
		if reply.strip() == 'y':
			decision_oc = 1
			break

# Если выбрано удаление последовательностей, найденных только в одной стране,
# этот блок проверит, существует ли файл entry_origin.txt со странами
# происхождения записей Генбанка, и если его нет, то создает заново,
# открывая infile.txt с форматом GenBank
countries = []
if decision_oc == 1:
	useq = 0
	listed = set()
	if not os.path.exists('entry_origin.txt'):
		print(' The file infile.txt with Genbank entries is needed ', end='')
		print('to extract labs locations')
		with open ('infile.txt') as inf:
			with open ('entry_origin.txt', 'w') as ouf:
				line = inf.readline()
				while line:
					if 'LOCUS' in line[:5]:
						accno = re.findall(r'\S+', line)[1]
						if accno not in listed and accno in accnos:
							ouf.write('>')
							ouf.write(accno)
							ouf.write('\n')
							useq += 1
							listed.add(accno)
							while line.strip() != '//':
								if 'Direct Submission' in line:
									line = inf.readline()
									line = inf.readline()
									while not re.findall('\S', line[:5]):
										preline = line
										line = inf.readline()
									country = re.findall(r'\S+$', preline)[0]
								line = inf.readline()
							ouf.write(country)
							ouf.write('\n')
					line = inf.readline()

# Создадим список со строками из выровненного и обрезанного файла
# с последовательностями. Он пригодится для проверки файла
# с географией лабораторий и потом, в основном блоке. Перевод в 1 линию!
allseqs = []
lines = []
applicant = ''
letts = ''
d = 0
with open('single.fasta') as inf:
	rawlines = inf.readlines()
	for i in rawlines:
		if '>' not in i:
			qwerty = i.replace('-', '')
		else:
			qwerty = i
		lines.append(qwerty)
	for i in lines:
		if re.findall('\w', i):
			if '>' in i:
				if d != 0:
					allseqs.append(applicant.strip().upper())
					allseqs.append(letts)
					letts = ''
				applicant = i
				d += 1
			else:
				letts += i.strip().upper()
allseqs.append(applicant.strip().upper())
allseqs.append(letts)

# Если надо, заберем страны из файла entry_origin.txt в список, проверив,
# что в нем строки с угольными скобками чередуются со строками
# с непробельными символами
a = 0
check1 = 0
origins = []
if decision_oc == 1:
	with open('entry_origin.txt') as inf:
		line = inf.readline()
		while line:
			if '>' in line:
				a += -1
				origins.append(line.strip().upper())
			elif re.findall('\S', line):
				a += 1
				origins.append(line.strip().upper())
			line = inf.readline()
	for m in allseqs: # все строки из single
		if '>' in m:
			t = 0
			for n in origins: # все строки из entry_origin
				if m[:m.find('.')] in n:
					t = 1
			if t == 0:
				print ('', 'no country for', m[1:])
				a += 1
	if a != 0:
		print (' Please correct the file entry_origin.txt')
		sys.exit()
	else:
		print (' List of country assignments successfully checked')

# Это основной блок, который выпишет последовательности по группам,
# укажет объем каждой группы в ее названии, укажет расшифровку групп
# после их названия. Здесь реализован выбор, включать или не включать
# в итоговый файл последовательности с вырожденными нуклеотидами
accnos = set()
group_origin = set()
group = 1
allrewseqs = 0
nonstandard = 0
runner = 0
with open('grouped.fasta', 'w') as ouf:
	while runner <= len(allseqs)-1:
		line = allseqs[runner]
		if '>' in line:
			accno = line.strip().upper()
			if accno not in accnos:
				sample = allseqs[runner+1].strip().upper()
				allrewseqs += 1 # Это счетчик генотипов
				if re.findall('[^ATGC]', sample):
					nonstandard += 1 # Счетчик нест. генотипов
				group_origin.add(accno)
				#print(accno)
				i = 0
				while i < len(allseqs)-1: # Этот цикл соберет список с номерами доступа
					counterpart_no = allseqs[i].strip()
					counterpart = allseqs[i + 1]
					if sample in counterpart and counterpart_no not in accnos:
						accnos.add(counterpart_no)
						group_origin.add(counterpart_no)
					i += 2
				if len(group_origin) > 1: # Этот цикл выпишет группы
					if decision_ns == 0:
						ouf.write('>')
						ouf.write('Group')
						ouf.write(str(group))
						group += 1
						ouf.write('-')
						ouf.write(str(len(group_origin)))
						ouf.write(' ')
						for element in group_origin:
							element = element.replace('>', '')
							ouf.write(element.upper())
							ouf.write(' ')
					elif decision_ns == 1:
						if not re.findall('[^ATGC]', sample):
							ouf.write('>')
							ouf.write('Group')
							ouf.write(str(group))
							group += 1
							ouf.write('-')
							ouf.write(str(len(group_origin)))
							ouf.write(' ')
							for element in group_origin:
								element = element.replace('>', '')
								ouf.write(element.upper())
								ouf.write(' ')
				else: # Этот цикл выпишет уникальные последовательности
					if decision_ns == 0:
						smth = str(re.findall(r'\S+', str(group_origin))[0][2:])
						ouf.write(smth)
					elif decision_ns == 1:
						if not re.findall('[^ATGC]', sample):
							smth = str(re.findall(r'\S+', str(group_origin))[0][2:])
							ouf.write(smth)
				if decision_ns == 0:
					ouf.write('\n')
					ouf.write(sample)
					ouf.write('\n')
				elif decision_ns == 1:
					if not re.findall('[^ATGC]', sample):
						ouf.write('\n')
						ouf.write(sample)
						ouf.write('\n')
			group_origin.clear()
		runner += 1
print('', len(allseqs)//2, 'sequences in total')
print('', allrewseqs, 'genotypes in the sample')
print('', nonstandard, 'genotypes with non-standard bases in the sample')

# Если указано принимать во внимание последовательности из не менее,
# чем двух стран, то этот блок обратится к файлу entry_origin.txt
# и перепишет в новый файл соответствующую часть содержания grouped.fasta
origins = set()
accorigins = []
if decision_oc == 1:
	with open ('entry_origin.txt') as inf:
		accorigins_raw = inf.readlines()
		for k in accorigins_raw:
			k = k.upper()
			accorigins.append(k)
	print(' The list of countries is following: ', end='')
	k = 1
	displaycountries = set()
	while k < len(accorigins):
		displaycountries.add(accorigins[k])
		k+=2
	for q in displaycountries:
		print(q, end='')
	with open ('grouped.fasta') as inf:
		with open ('grouped_nloc.fasta', 'w') as ouf:
			line = inf.readline()
			while line:
				if '>' in line:
					accessions = re.split(' ', line.strip())
					if len (accessions) > 2:
						caption = accessions[0]
						del accessions[0]
						for unit in accessions:
							i = 0
							while i < len (accorigins):
								if unit.upper() in accorigins[i]:
									origins.add(accorigins[i+1].strip())
								i += 1
						if len (origins) > 1:
							ouf.write(caption)
							ouf.write(' ')
							smth = str(re.findall(r'\w+', str(accessions))).replace('\'', '')
							smth = smth.replace('[', '').replace(']', '')
							ouf.write(smth)
							ouf.write('\n')
							line = inf.readline()
							ouf.write(line)
						origins.clear()
				line = inf.readline()

# Это элемент выписывает номера доступа по группам
if decision_oc == 1:
	with open ('grouped_nloc.fasta') as inf:
		with open ('groups.txt', 'w') as ouf:
			line = inf.readline()
			while line:
				if 'Group' in line:
					elements = re.split(' ', line.strip())
					ouf.write(elements[0].replace('>', ''))
					ouf.write('\n')
					k = 1
					while k < len(elements):
						ouf.write(elements[k])
						ouf.write(' ')
						k += 1
					ouf.write('\n')
				line = inf.readline()

	# Этот скрипт посчитает объем выборки, сложив окончания заголовков
	# в формате "GroupN-quant" например, "Group7-237" (how_many.py)
	import re
	groups = 0
	quantity = 0
	with open ('groups.txt') as inf:
		line = inf.readline()
		while line:
			if 'Group' in line:
				groups += 1
				fragments = re.findall(r'\w+', line)
				quantity += int(fragments[1])
			line = inf.readline()
	print (quantity, 'accessions in', groups, 'groups.txt')