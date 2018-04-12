import csv, itertools
#chromosome[0]      position[1]        reference[2]       alternative[3]     Upstream Int Adj[4]        Down Int Adj[5]    Exon up[6] Exon down[7]       strand[8]  Sequence[9]        Final Sequence[10]  Barcode[11]
base_dir = '/home/CAM/adamson/Vex-seq/single_way/68_97_pool/'
a= open(base_dir + '68_97_final_sequences_fixed.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
BCs = {};seen_randomer = {};read_BCsK= {}; read_BCsH={}

i=0
for line in reader:
    if i ==0:
	i=1
    else:
	BC = str(line[11])
	ID = '_'.join(line[0:4])
	BCs[BC] = ID
#	read_BCsH[BC] = 0; read_BCsK[BC] = 0
	seen_randomer[BC] = {}
	samples = ['K562-1_S1_L001_R1_001.fastq','K562-2_S2_L001_R1_001.fastq', 'K562-3_S3_L001_R1_001.fastq', 'HepG2-1_S4_L001_R1_001.fastq', 'HepG2-2_S5_L001_R1_001.fastq','HepG2-3_S6_L001_R1_001.fastq']
	for sample in samples:
	    open(base_dir +'170411_B3PGD_solid_amp/Refs_GTFs/' + ID + '/' + sample, 'w').close()
	    open(base_dir +'170411_B3PGD_solid_amp/Refs_GTFs/' + ID + '/' + sample.replace('R1', 'R2'), 'w').close()
a.close()

print "PHASE 1 complete"
def RC(string):
    comp = {'N': 'N', 'A' : 'T', 'C':'G', 'G': 'C', 'T' : 'A'}
    new_string = string[::-1]
    rev_comp = []
    for nuc in new_string:
	rev_comp.append(comp[nuc])
    return ''.join(rev_comp)

BC_reads = {}
samples = ['K562-1_S1', 'K562-2_S2', 'K562-3_S3', 'HepG2-1_S4', 'HepG2-2_S5', 'HepG2-3_S6']
read_fields = ['read', 'seq', 'strand', 'qual']
for sample in samples:
    PCR_dup = 0; no_E3 = 0; BC_foul = 0; bad_primer=0; bp1 =0; pass_=0
    b= open(base_dir + '170411_B3PGD_solid_amp/'+ sample + '_L001_R1_001.fastq', 'r')
    c= open(base_dir + '170411_B3PGD_solid_amp/'+ sample + '_L001_R2_001.fastq', 'r')
    i=0;j=0;F_read_info = {}; R_read_info = {}
    for F_line, R_line in itertools.izip(b, c):
	F_line = F_line.replace('\n', ''); R_line = R_line.replace('\n', '')
	if i ==4:
	    j=j+1; i=0
	    if RC(R_read_info['seq'])[-64:-58] == 'TCTAGA':#TCTAGAGGACCCGTTTAAACCCACTGATCAGCCTCGACTGTGCCTTCTAGTTGCNNNNNNNNNN
		BC = RC(R_read_info['seq'])[-72:-64]
		pass_ +=1#pass_ with TCTAGA site
		if BC in BCs:
		    Randomer = RC(R_read_info['seq'])[-11:-1]
		    if Randomer not in seen_randomer[BC]:
			if RC(R_read_info['seq'])[-82:-72] == 'TCGCTCTAGT':#need exon 3
			    if BC not in BC_reads:
			    	BC_reads[BC] = str(F_read_info['read']).replace('>', '')
			    else:
				BC_reads[BC] = BC_reads[BC]  + ',' + str(F_read_info['read'])
			    d = open(base_dir +'170411_B3PGD_solid_amp/Refs_GTFs/' + BCs[BC] + '/' + sample + '_L001_R1_001.fastq', 'a')
 			    e = open(base_dir +'170411_B3PGD_solid_amp/Refs_GTFs/' + BCs[BC] + '/' + sample + '_L001_R2_001.fastq', 'a')
		    	    for field in read_fields:
			    	d.write(F_read_info[field] + '\n')
			    	e.write(R_read_info[field] + '\n')
		    	    d.close()
		    	    e.close()
			    seen_randomer[BC][Randomer] = True
			else:
			    no_E3 +=1
			    #print RC(R_read_info['seq'])[-82:-72]
			    #print RC(R_read_info['seq'])[-87:-65]
		    else:
			PCR_dup +=1
			#print RC(R_read_info['seq'])[-11:-1]
			#print RC(R_read_info['seq'])[-15:-1]
		else:
		    #print RC(R_read_info['seq'])[-72:-64]
		    #print RC(R_read_info['seq'])[-76:-60]
		    BC_foul +=1
	    else:
		#print RC(R_read_info['seq'])[-64:-58]
		#print RC(R_read_info['seq'])[-70:-50]
		if RC(R_read_info['seq'])[-63:-57] == 'TCTAGA'or RC(R_read_info['seq'])[-65:-59] == 'TCTAGA':
		    bp1+=1
		bad_primer +=1
	    F_read_info = {}; R_read_info ={}
	if i ==0:
	    F_read_info['read'] = F_line; R_read_info['read'] = R_line
	if i ==1:
	    F_read_info['seq'] = F_line; R_read_info['seq'] = R_line
	if i== 2:
	    F_read_info['strand'] = F_line; R_read_info['strand'] = R_line
	if i==3:
	    F_read_info['qual'] = F_line; R_read_info['qual'] = R_line
	i=i+1
    b.close()
    c.close()
    #print seen_randomer
    print sample
    print 'total reads: ' + str(j)
    print 'BC off by a nuc: ' + str(bp1)
    print 'bad primer near BC: ' + str(bad_primer)
    print 'barcode not in pool: ' + str(BC_foul)
    print 'PCR duplicate: ' +str(PCR_dup)
    print 'no exon 3: ' + str(no_E3)
    print 'good reads: ' + str(j-bp1-bad_primer-BC_foul-PCR_dup-no_E3)
    print 'percent happy reads: ' + str((float(j-bp1-bad_primer-BC_foul-PCR_dup-no_E3)/j) *100)

#d= open(base_dir + '161212_AV7M1_Vex-seq_Delta/barcode_info.tsv', 'wb')
#writer = csv.writer(d, delimiter = '\t')
#writer.writerow(['BC', 'Vexon', 'HepG2_reads', 'K562_reads', 'read_ids'])
#one=0;five=0;ten=0
#theres a key error in here somewhere, but it may not matter if you don't need the barcode info column
#for item in BCs:
#    writer.writerow([item, BCs[item], read_BCsH[item], read_BCsK[item], BC_reads[BC]])
#    tot = read_BCsH[item]+ read_BCsK[item]
#    if tot >0:
#	one +=1
#    if tot > 4:
#	five +=1
#    if tot > 9:
#	ten +=1
#d.close()
