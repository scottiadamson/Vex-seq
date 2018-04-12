import csv
variants = open('variants.txt', 'r').read().replace('\n', '').split(' ')
samples=['HepG2-1_S4', 'HepG2-2_S5', 'HepG2-3_S6','K562-1_S1', 'K562-2_S2', 'K562-3_S3']
assigned_BCs = {}
e=open('68_97_final_sequences_fixed.tsv', 'r')
reader = csv.reader(e ,delimiter = '\t')
next(reader)
for line in reader:
    assigned_BCs[str(line[-1])] = '_'.join(line[0:4])
e.close()
c=open('../processed_files/splice_results_by_BC.tsv', 'wb')
writer = csv.writer(c, delimiter = '\t')
writer.writerow(['variant', 'BC', 'H1_in', 'H1_out', 'H1_other', 'H2_in', 'H2_out', 'H2_other', 'H3_in', 'H3_out', 'H3_other', 'K1_in', 'K1_out', 'K1_other', 'K2_in', 'K2_out', 'K2_other', 'K3_in', 'K3_out', 'K3_other'])

unannotated = 0;annotated=0;seen = {}
variant_reads = dict(zip(samples, [{},{},{},{},{},{}]))
for variant in variants:
    variant_reads= dict(zip(samples, [{},{},{},{},{},{}]))
    for sample in samples:
	i=0
#chr22_17990875_T_C	.	transcript	1	1216	.	+	.	transcript_id "chr22_17990875_T_C"
#chr22_17990875_T_C	.	exon	1	44	.	+	.	transcript_id "chr22_17990875_T_C"; exon_id "exon_1"
#chr22_17990875_T_C	.	exon	216	299	.	+	.	transcript_id "chr22_17990875_T_C"; exon_id "cassette"
#chr22_17990875_T_C	.	exon	1003	1216	.	+	.	transcript_id "chr22_17990875_T_C"; exon_id "exon_3"
	b=open('Refs_GTFs/' +variant +'/'+ variant + '.gtf', 'r')
	gtf_reader = csv.reader(b, delimiter = '\t')
	for line in gtf_reader:
	    if i ==1:
		I1_start = int(line[4]) +1#	45
	    elif i==2:
		I1_end = int(line[3])-1#	215
		I2_start = int(line[4]) +1#	300
	    elif i==3:
		I2_end = int(line[3]) -1#	1002
	    i+=1
	b.close()
#M01212:141:000000000-B3PGD:1:1101:10644:5218	99	chr13_49775314_CAT_C	1	255	44M165N54M8D29M702N23M	=	1066	1205	GGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGAAACATCCAATCAAAAATAATGGGATGCATTCAACTTATAGATAACATTAGGAGATGATATTCAACTTCATTGAACAGTTCAACTCCTGGGCAACGTGCTGGTCTG	ABBBBBBFF5DFFGFFGGGGGGGFHHGGGGFGHHGCEHHHGHFGGGHHHHHHHFHGGFFGGHFGHHHHFFFFHHHHGHHDGHHHHHHHHHHHHHEHHHHFHHHFGHHHHHHGHHHHHHHHHGHHHHHHGHHHFCFGGFHGHHEFHGHHHB	jI:B:i,45,209,301,1002
#M01212:141:000000000-B3PGD:1:1101:10644:5218	147	chr13_49775314_CAT_C	1066	255	140M10S	=	1	-1205	GTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTCTAGTACGTCCACTCTAGAGGGCCCGTTTAAACCCGCTGATCAGCCTCGACTGTGCCTTCTAGTTGCTGATTTCGCG	GGHHGGGGHGHFHGFGHHGHFHHFHFFHFHDDFHGCHHGGHEGGHGGHGGGHHHHHHFGHHGFFHGGGGGHGHHHGGHGGHFGGFHHFHFHHEEFEGGHHHHGGGGGGHHHGHHFGEGGGGHHHHHGHGGGGFGGGFGFFFFBBDABBBB	jI:B:i,-1
#M01212:141:000000000-B3PGD:1:1101:13100:8208	99	chr13_49775314_CAT_C	1	255	44M958N106M	=	1066	1205	GGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATG	ABBCBCCFFFFFGGGGGGGGGGHHHHGGGGFHHHGGGHHHGHHGGGHHHHHGHHGHGHHHGHHHHHHHHHHHGHHGHHGHHGHHHHHHHHGHHHHHHHHHHHGGGGHHHHHHHHGHHGHHHHHHHHHHHHGHHGGHHGGHHHGGHHHHHG	jI:B:i,45,1002
#M01212:141:000000000-B3PGD:1:1101:13100:8208	147	chr13_49775314_CAT_C	1066	255	140M10S	=	1	-1205	GTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTCTAGTACGTCCACTCTAGAGGGCCCGTTTAAACCCGCTGATCAGCCTCGACTGTGCCTTCTAGTTGCTGTCTGACCT	FHHGHHGGFHHHHHDHHHHHHHHHHHEHFHGGHHHHHHGHHGGGHHGHGGGHHGHHHHHHHHHGHFGGGGHHHHFGHHGGGGGGBHHGBHHHGGGGGGHHHFEGGGGGHHHHHGHHFGGFGHHHHHHHGGGFGGGGGGBFFFFFFAABAA	jI:B:i,-1
	BCs = {}
	a=open('Refs_GTFs/' +variant + '/' +sample + '_' + variant + 'Aligned.out.sam', 'r')
	reader = csv.reader(a, delimiter = '\t')
	next(reader);next(reader);next(reader);next(reader);i=0
	for line in reader:
	    if i%2 ==0:
		inclusion =10
		jI_flag = str(line[-1]).replace('jI:B:i,', '').split(',')
		if len(jI_flag) ==2:
		    if int(jI_flag[0]) == I1_start and int(jI_flag[1]) == I2_end:
			inclusion = False
		elif len(jI_flag) ==4:
		    if int(jI_flag[0]) == I1_start and int(jI_flag[1]) == I1_end and int(jI_flag[2]) == I2_start and int(jI_flag[3]) == I2_end:
			inclusion = True
	    else:
		BC =str(str(line[9])[-72:-64])
		if BC not in BCs:
		    BCs[BC] = [0, 0, 0]#in, out, other
		if inclusion == True:
		    BCs[BC][0] +=1;annotated +=1
		elif inclusion == False:
		    BCs[BC][1] +=1; annotated +=1
		else:
		    BCs[BC][2] +=1; unannotated +=1
	    i+=1
	a.close()
	variant_reads[sample]= BCs
    for BC in BCs:
	seen[BC] =1
	for sample in samples:
	    if BC not in variant_reads[sample]:
		variant_reads[sample][BC] = [0,0,0]
	writer.writerow([variant, BC, variant_reads[samples[0]][BC][0], variant_reads[samples[0]][BC][1],variant_reads[samples[0]][BC][2],variant_reads[samples[1]][BC][0], variant_reads[samples[1]][BC][1], variant_reads[samples[1]][BC][2], variant_reads[samples[2]][BC][0], variant_reads[samples[2]][BC][1], variant_reads[samples[2]][BC][2], variant_reads[samples[3]][BC][0], variant_reads[samples[3]][BC][1], variant_reads[samples[3]][BC][2], variant_reads[samples[4]][BC][0], variant_reads[samples[4]][BC][1], variant_reads[samples[4]][BC][2], variant_reads[samples[5]][BC][0], variant_reads[samples[5]][BC][1], variant_reads[samples[5]][BC][2]])
    print variant
for BC in assigned_BCs:
    if BC not in seen:
	writer.writerow([assigned_BCs[BC], BC, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
#'HepG2-1_S33', 'HepG2-2_S34', 'HepG2-3_S35', 'K562-1_S30', 'K562-2_S31', 'K562-3_S32']
c.close()
print annotated
print unannotated

