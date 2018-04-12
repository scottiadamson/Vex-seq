import csv

a=open('68_97_final_sequences_fixed.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
#chromosome	position	reference	alternative	Upstream Int Adj	Down Int Adj	Exon up	Exon down	strand	Sequence	Final Sequence	Barcode
#chr1	16069286	Control	Original	16069268	16069433	16069330	16069402	+	ATGAAGCCTGCCCTCCTAACACCCCTGCACGCTGCCTCGCTGCCTCCCTCCCGCCCCCTCAGGCCTGATCGGGGCCCTTGGAGTCTTGGTCCTGAACAGCCTCCTGAAAGTTTACTTCTTCGTGGGCTGTGCCAAGTGAGTGCCCACCTCATCCCCCCAGGGAATG	TAGCGTCTGTCCGTCTGCAGATGAAGCCTGCCCTCCTAACACCCCTGCACGCTGCCTCGCTGCCTCCCTCCCGCCCCCTCAGGCCTGATCGGGGCCCTTGGAGTCTTGGTCCTGAACAGCCTCCTGAAAGTTTACTTCTTCGTGGGCTGTGCCAAGTGAGTGCCCACCTCATCCCCCCAGGGAATGCAATTGACTACTAGTTCTTAGGCTCTAGACAACTACTACTACAG	TCTTAGGC
next(reader)

base ='/home/CAM/adamson/Vex-seq/single_way/68_97_pool/161212_AV7M1_Vex-seq_Delta/Refs_GTFs/'

def make_1o(line, tot):
    seq = str(line[10])[20:-29]
    FWD_primer = 'CCACTGACTCTCTCTGCCTCTGCAG'
    REV_primer = 'TCTAGAGGGCCCGTTTAAACCCGCTA'
    Illumina_REV = 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
    Illumina_FWD = 'ACACTCTTTCCCTACACGACGCTCTTCCGATCT'
    if tot == True:
	return Illumina_FWD + FWD_primer + seq + 'NNNNNNNN' + REV_primer + Illumina_REV
    else:
        return FWD_primer + seq + str(line[11]) + REV_primer
seen = {}
d=open('1o_ref_indels.fa', 'wb')
for line in reader:
    ID = '>' + '_'.join(line[0:4])
    if ID not in seen:
	if 'Control' in str(ID):
	    d.write(ID +'\n' + make_1o(line, True)+'\n')
	elif len(str(line[2])) != len(str(line[3])):
	    d.write(ID +'\n' + make_1o(line, True)+'\n')
    seen [ID] = 1
d.close()



