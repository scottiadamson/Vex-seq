import csv, sys

d=open('1o_ref.fa', 'r')
i=0;BC_loc = {}
for line in d:
    if i%2 ==0:
        control_name = str(line).replace('>', '').replace('\n', '')
    else:
        seq = str(line).replace('\n', '')
        this_BC = [seq.find('TCTAGA')-7, seq.find('TCTAGA')+1]
        BC_loc[control_name] =this_BC
    i+=1

#chromosome	position	reference	alternative	Upstream Int Adj	Down Int Adj	Exon up	Exon down	strand	Sequence	Final Sequence	Barcode
#chr1	16069286	Control	Original	16069268	16069433	16069330	16069402	+	ATGAAGCCTGCCCTCCTAACACCCCTGCACGCTGCCTCGCTGCCTCCCTCCCGCCCCCTCAGGCCTGATCGGGGCCCTTGGAGTCTTGGTCCTGAACAGCCTCCTGAAAGTTTACTTCTTCGTGGGCTGTGCCAAGTGAGTGCCCACCTCATCCCCCCAGGGAATG	TAGCGTCTGTCCGTCTGCAGATGAAGCCTGCCCTCCTAACACCCCTGCACGCTGCCTCGCTGCCTCCCTCCCGCCCCCTCAGGCCTGATCGGGGCCCTTGGAGTCTTGGTCCTGAACAGCCTCCTGAAAGTTTACTTCTTCGTGGGCTGTGCCAAGTGAGTGCCCACCTCATCCCCCCAGGGAATGCAATTGACTACTAGTTCTTAGGCTCTAGACAACTACTACTACAG	TCTTAGGC
a=open('../../161212_AV7M1_Vex-seq_Delta_flash/68_97_final_sequences_fixed.tsv', 'r')
reader = csv.reader(a, delimiter = '\t')
next(reader)
controls = {};seen = {};var_order = [];coordinates = {};strand= {};assigned_BCs={}; BC_counts = {}
for line in reader:
    ID = '_'.join(line[0:4])
    BC_counts[ID] = {}
    assigned_BCs[str(line[-1])] = ID
    if 'Control' == str(line[2]):
	controls[ID] = {}
    if ID not in seen:
	seen[ID] = 1
	strand[ID] = str(line[8])
	var_order.append(ID)
	if strand[ID] == '+':
            if int(line[4]) not in coordinates:
	    	coordinates[int(line[4])] = [ID]
	    else:
	    	coordinates[int(line[4])].append(ID)
	else:
	    if int(line[5]) not in coordinates:
		coordinates[int(line[5])] = [ID]
	    else:
		coordinates[int(line[5])].append(ID)
a.close()

for coor in coordinates:
    relative_positions = {}
    for ID in coordinates[coor]:#list of IDs
	if 'Control' not in ID:
	    pos = int(ID.split('_')[1])
	    ref = str(ID.split('_')[2])
	    alt = str(ID.split('_')[3])	    
	    if strand[ID] == '+':
	        rel_pos = 58 + pos - coor+1
	    else:
		rel_pos = 58 - pos + coor+1
	    if len(ref) != len(alt) and strand[ID] == '-':
		rel_pos = 58 - (pos+len(ref)) + coor +1 
	    if rel_pos not in relative_positions:
	        relative_positions[rel_pos] = [ID]
	    else:
		relative_positions[rel_pos].append(ID)
	elif 'Original' in ID:
	    control = ID
    try:
	controls[control] = relative_positions#chr1_3241_Control_Original = {3: [chr1_3245_A_C]}
    except NameError:
	continue
print controls['chr3_50212479_Control_Original']
#print controls['chr2_191862921_Control_Original']

def RC(seq):
    comp = {'A':'T', 'C':'G', 'G' : 'C', 'T':'A'}
    new_seq = []
    for nuc in list(seq):
	new_seq.append(comp[nuc])
    return ''.join(new_seq[::-1])

#bash sam2tsv.sh | python parse_sam2tsv.py -
#read_name	read_flags	reference_name	readpos	readbase	readqual	refpos	refbase	cigarop
#M01212:132:000000000-AV7M1:1:1102:23271:14562	0	chr1_16069286_Control_Original	76	C	G	111	C	M
i=0;prev_read_ID=None;current_read_ID = None;  match_maker = {}; read_BCs = {};matches = []; BC =[]
m_num = 100000;del_nums = [100000, 100000];inserted_nucs = []; insertion_annotated = []; del_var = '';relevant_lines = []; prev_line = ['.'] *9
for line in sys.stdin:
    line = str(line).replace('\n', '').split('\t')
    if i ==0:
	i=1
    else:
	current_read_ID = str(line[0])
	if str(line[2]) != '.':
	    current_alignment = str(line[2])
	    these_vars = controls[current_alignment]
	    if prev_read_ID != current_read_ID and prev_read_ID != None and prev_alignment != '.':
		if len(matches) == 0:
		    matches.append(prev_alignment)
		match_maker[prev_read_ID] = matches #{M01212:132:000000000-AV7M1:1:1102:23271:14562: [chr1_123432_A_C, chr1_13242_G_A]}
		read_BCs[prev_read_ID] = ''.join(BC)
		BC=[]; matches = []
	    if str(line[6]) != '.':
		if int(line[6]) in these_vars:
		    for var in these_vars[int(line[6])]:
			ID = var
			ref = str(ID.split('_')[2])
			alt = str(ID.split('_')[3])
 			pos = int(line[6])
			if strand[ID] == '+':
			     new_alt = alt[0]
			else:
			     new_alt = RC(alt[0])
			if str(line[8]) == 'M' and str(line[4]) == str(new_alt) and len(ref) == len(alt):
			    if str(line[4]) == str(new_alt):
				matches.append(ID)
			elif str(line[8]) == 'D' and len(ref)> len(alt) and str(prev_line[8]) == 'M':
		    	    del_var = ID
			    if strand[del_var] == '+':
				del_nums = [pos, pos + pos + len(ref)]
				m_num = pos + len(ref) +1
		if int(line[6]) in range(int(del_nums[0]), int(del_nums[1])+1):
		    if str(line[8]) !='D':
			del_nums = [100000, 100000]; m_num = 10000
		elif str(line[8]) == 'M' and int(line[6]) == m_num:
		    matches.append(del_var)
		    del_nums =[100000, 100000]; m_num =10000
		if int(line[6]) > m_num:#reset info
		    del_nums = [100000, 100000]; m_num =10000
	    else:
		if str(prev_line[6]) != '.' and str(line[8]) == 'I':
		    if int(prev_line[6]) in these_vars:
#		    	for var in these_vars[int(prev_line[6])]:
#		    	    ref = str(ID.split('_')[2])
#		    	    alt = str(ID.split('_')[3])
#			    if len(alt) > len(ref):
# 		    	    	pos = int(prev_line[6])
#			    	insertion_annotated = list(alt)
#			    	insertion_ID = var
#			    	if str(prev_line[8]) == 'M':
#				    inserted_nucs.append(str(line[4]))
#	    if len(inserted_nucs) !=0 and str(line[6]) != '.':#line after insertion variaant is detected
#		if str(line[8]) == 'M' and inserted_nucs == insertion_annotated:#this will probably only work for insertion of length one right now, but that's okay 
#		    matches.append(insertion_ID)
#		if str(prev_line[8]) != 'I':#reset
#		    inserted_nucs = []; insertion_annotated = []; insertion_ID = ''
	    if str(line[6]) != '.':#get barcode
	    	if int(line[6]) in range(this_BC[0], this_BC[1]):
	    	    BC.append(str(line[4]))
	    prev_alignment = current_alignment
	    prev_read_ID = current_read_ID
	    prev_line = line
read_BCs[current_read_ID] = ''.join(BC)
if len(matches) ==0:
    matches.append(str(line[2]))
match_maker[current_read_ID] = matches

for read in match_maker:
   BC = read_BCs[read]
   if BC in assigned_BCs:
    	if len(match_maker[read]) ==1:
	    snp = str(match_maker[read][0])
	    if snp != '.':
	    	if BC not in BC_counts[snp]:
	            BC_counts[snp][BC] = 1
	    	else:
	    	    BC_counts[snp][BC] +=1
    	else:#assigns to the designed one if multiple snps discovered
	    for match in match_maker[read]:
	    	if str(assigned_BCs[BC]) == str(match):
		    if BC not in BC_counts[match]:
	    	    	BC_counts[match][BC] = 1
		    else:
		    	BC_counts[match][BC] += 1

e=open('correct_BC_custom_novo_final.tsv', 'wb')
writer = csv.writer(e, delimiter = '\t')
writer.writerow(['variant', 'correct_count', 'incorrect_count', 'correct_BCs', 'incorrect_BCs'])
for variant in var_order:
    correct_BCs = []; incorrect_BCs= [];corr_count=0; incorr_count=0
    for BC in BC_counts[variant]:
	if variant == assigned_BCs[BC]:
	    correct_BCs.append(BC + ':' +str(BC_counts[variant][BC]))
	    corr_count += int(BC_counts[variant][BC])
	else:
	    incorrect_BCs.append(BC + ':' +str(BC_counts[variant][BC]))
	    incorr_count += int(BC_counts[variant][BC])
    writer.writerow([variant, corr_count, incorr_count, ','.join(correct_BCs), ','.join(incorrect_BCs)]) 
e.close()
