fastqs = []

for i in range(1, 13):
    if len(str(i)) == 1:
        fastqs.append("barcode0" + str(i) + ".fastq")
    else:
        fastqs.append("barcode" + str(i) + ".fastq")
fastqs.append("unclassified.fastq")

for fastq_f in fastqs:
    with open(fastq_f, 'r') as fastq:
        with open(fastq_f[:-5] + 'txt', 'w') as txt:
            line = fastq.readline()
            while line:
                if line.startswith('@'):
                    line = fastq.readline() # go to seq
                    txt.write(line)
                    line = fastq.readline() # go to '+'
                    line = fastq.readline() # go to PHRED
                line = fastq.readline() # go to next sequence