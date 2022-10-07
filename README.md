## Домашнее задание №1

Создаем символическую ссылку на каждый из файлов:
```
ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq
```

С помощью команды seqtk выбираем случайно 5 миллионов чтений типа paired-end и 1.5 миллиона чтений типа mate-pairs (мой seed = 920)
```
seqtk sample -s920 oil_R1.fastq 5000000 > sub1.fq
seqtk sample -s920 oil_R2.fastq 5000000 > sub2.fq
seqtk sample -s920 oilMP_S4_L001_R1_001.fastq 1500000 > msub1.fq
seqtk sample -s920 oilMP_S4_L001_R2_001.fastq 1500000 > msub2.fq
```

С помощью программ fastQC и multiQC оцениваем качество исходных чтений и получаем по ним общую **статистику**.
```
fastqc sub1.fq
fastqc sub2.fq
fastqc msub1.fq
fastqc msub2.fq
multiqc .
```

С помощью программ platanus_trim и platanus_internal_trim подрезаем чтения по качеству и удалить адаптеры
```
platanus_trim sub1.fq sub2.fq
platanus_internal_trim msub1.fq msub2.fq
```

Перенесем файлы в отдельную папку trimmed.
```
mkdir trimmed
mv *.trimmed trimmed
mv *.int_trimmed trimmed
```

После подрезания чтений удалим исходные .fastq файлы, полученные с помощью программы seqtk.
```
rm sub1.fq
rm sub2.fq
rm msub2.fq
rm msub1.fq
```

С помощью программы fastQC и multiQC оценить качество подрезанных чтений и получаем по ним общую **статистику**.
```
fastqc *trimmed
multiqc .
```

С помощью программы “platanus assemble” собраем контиги из подрезанных чтений.
```
platanus assemble -f sub1.fq.trimmed sub2.fq.trimmed 2> assemble.log
```
**Анализ полученных контигов.**

Ссылка на collab: https://colab.research.google.com/drive/1kQxTdcArcjqGPdhR3LVAd5iLewEFFa2i#scrollTo=DcbVqy5Ta8cP
```
contigs_number = 0
contigs_len = 0
max_contig_len = 0
contigs_lens = []

# Filling contigs info list.
with open("out_contig.fa") as file:
  for line in file:
    if line.startswith(">") :
      len = line.split('len')[1].split('_')[0]
      len = int(len)

      contigs_len += len
      contigs_number += 1
      contigs_lens.append(len)

max_contig_len = max(contigs_lens)
contigs_lens.sort(reverse = True)
print(contigs_lens)

print('общее кол-во контигов', contigs_number)
print('их общая длина', contigs_len)
print('длина самого длинного контига', max_contig_len)

sum = 0
n50 = 0
for len in contigs_lens:
  sum += len
  if sum * 2 >= contigs_len : 
    n50 = len
    break

print('N50', n50)
```
- общее кол-во контигов 636
- их общая длина 3927148
- длина самого длинного контига 179307
- N50 47993
 
С помощью программы “ platanus scaffold” собрараем скаффолды из контигов, а также из подрезанных чтений.
```
platanus gap_close -o GAP_CLOSE -c SCAF_scaffold.fa ../sub1.fq.trimmed ../sub2.fq.trimmed -OP2 ../msub1.fq.int_trimmed ../msub2.fq.int_trimmed 2> gap_close.log
```
**Анализ полученных скафолдов.**
```
scaffold_number = 0
scaffold_len = 0
max_scaffold_len = 0
scaffolds_lens = []

with open("SCAF_scaffold.fa") as file:
  for line in file:
    if line.startswith(">") :
      len = line.split('len')[1].split('_')[0]
      len = int(len)

      scaffold_len += len
      scaffold_number += 1
      scaffolds_lens.append(len)

max_scaffold_len = max(scaffolds_lens)
scaffolds_lens.sort(reverse = True)
print(scaffolds_lens)

print('общее кол-во скаффолдов', scaffold_number)
print('их общая длина', scaffold_len)
print('длина самого длинного скаффолда', max_scaffold_len)

sum = 0
n50 = 0
for len in scaffolds_lens:
  sum += len
  if sum * 2 >= scaffold_len : 
    n50 = len
    break

print('N50', n50)
```
- общее кол-во скаффолдов 71
- их общая длина 3874005
- длина самого длинного скаффолда 3832100

Для самого длинного скаффолда посчитаем количество гэпов (участков, состоящих из букв NNNN) и их общую длину.
```

```

С помощью программы “ platanus gap_close” уменьшим кол-во гэпов с помощью подрезанных чтений.
```
platanus gap_close -o GAP_CLOSE -c SCAF_scaffold.fa -IP1 ../sub1.fq.trimmed ../sub2.fq.trimmed -OP2 ../msub1.fq.int_trimmed ../msub2.fq.int_trimmed 2> gap_close.log
```
Для самого длинного скаффолда посчитаем количество гэпов (участков, состоящих из букв NNNN) и их общую длину.
```

```
