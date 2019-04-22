#chromosome[0]      position[1]        reference[2]       alternative[3]     Upstream Int Adj[4]        Down Int Adj[5]    Exon up[6] Exon down[7]       strand[8]  Sequence[9]        Final Sequence[10]  Barcode[11]
import csv
base_dir = '/home/CAM/adamson/Vex-seq/single_way/68_97_pool/'
a = open(base_dir+'68_97_final_sequences_fixed.tsv', 'r')
reader = csv.reader(a ,delimiter = '\t')

rev_primer_XbaI = 'TCTAGAGGACCCGTTTAAACCCACTGATCAGCCTCGACTGTGCCTTCTAGTTGCNNNNNNNNNN'
E3_dig_bits = 'CTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTATCACTAAGCTCGCTCTAGT'
I2_on = 'CAATTGTTTTCTTTTGTTTAATTCTTGCTTTCTTTTTTTTTCTTCTCCGCAATTTTTACTATTATACTTAATGCCTTAACATTGTGTATAACAAAAGGAAATATCTCTGAGATACATTAAGTAACTTAAAAAAAAACTTTACACAGTCTGCCTAGTACATTACTATTTGGAATATATGTGTGCTTATTTGCATATTCATAATCTCCCTACTTTATTTTCTTTTATTTTTAATTGATACATAATCATTATACATATTTATGGGTTAAAGTGTAATGTTTTAATATCGATACACATATTGACCAAATCAGGGTAATTTTGCATTTGTAATTTTAAAAAATGCTTTCTTCTTTTAATATACTTTTTTGTTTATCTTATTTCTAATACTTTCCCTAATCTCTTTCTTTCAGGGCAATAATGATACAATGTATCATGCCTCTTTGCACCATTCTAAAGAATAACAGTGATAATTTCTGGGTTAAGGCAATAGCAATATTTCTGCATATAAATATTTCTGCATATAAATTGTAACTGATGTAAGAGGTTTCATATTGCTAATAGCAGCTACAATCCAGCTACCATTCTGCTTTTATTTTATGGTTGGGATAAGGCTGGATTATTCTGAGTCCAAGCTAGGCCCTTTTGCTAATCATGTTCATACCTCTTATCTTCCTCCCACAG'
#						     #INTRON													  #PstI
E1_on = 'GGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGTTGGTATCAAGGTTACAAGACAGGTTTAAGGAGACCAATAGAAACTGGGCATATGGAGACAGAGAAGACTCTTGGGTTTCTGATAGGGCCCACTGACTCTCTCTGCCTCTGCAG'
E1 = 'GGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAG'


def get_seqs(line_):
    new_seq = str(line[9])
    exons={}#cannibalized from other script; too lazy to rewrite better
    ID = '_'.join(line_[0:4])
    definition = line_[4:9]
    exons[ID] = definition
    for ID in exons:
    	if str(ID.split('_')[2]) != 'Control':
	    adjustment = len(ID.split('_')[3]) - len(ID.split('_')[2])#CCCT_CCT ---> -1
    	else:
	    adjustment = 0
        adjustment_acceptor = 0; adjustment_donor = 0
    	splice_acceptor = int(exons[ID][2]) - int(exons[ID][0])
    	splice_donor = int(exons[ID][3]) - int(exons[ID][0])
    	if str(exons[ID][4]) =='+':
	    if adjustment !=0:#[upstream_loc[0], downstream_loc[1], exon_up_loc[2], exon_down_loc[3], strand[4]]
	    	mut_loc = int(ID.split('_')[1])
	    	if mut_loc < int(exons[ID][3]):
		    if mut_loc < int(exons[ID][2]):
		    	adjustment_acceptor = adjustment; adjustment_donor = adjustment
		    elif mut_loc > int(exons[ID][2]):
                    	adjustment_acceptor = 0; adjustment_donor = adjustment
	    E2 = new_seq[splice_acceptor+adjustment_acceptor:splice_donor+adjustment_donor+1]
        else:
	    splice_acceptor = int(exons[ID][1]) - int(exons[ID][3])
	    splice_donor = int(exons[ID][1]) - int(exons[ID][2])
	    if adjustment != 0:
	    	mut_loc = int(ID.split('_')[1])
	    	if mut_loc >int(exons[ID][2]):
		    if mut_loc > int(exons[ID][3]):
		    	adjustment_acceptor = adjustment; adjustment_donor = adjustment
		    elif mut_loc > int(exons[ID][2]):
		    	adjustment_acceptor = 0; adjustment_donor = adjustment
	    E2 = new_seq[splice_acceptor+adjustment_acceptor:splice_donor+adjustment_donor+1]
    Intron = new_seq[splice_donor+adjustment_donor+1:]
    return [new_seq, E2, Intron]

i=0
spliced_out = E1 + E3_dig_bits +'NNNNNNNN'+rev_primer_XbaI
class_prev = 'chr1_16069286_Control_Original'; seqs = {} 
for line in reader:
    if i==0:
	i=1
    else:
	class_curr = '_'.join(line[0:4])
	if class_curr != class_prev:
	    b = open(base_dir + '161212_AV7M1_Vex-seq_Delta/Refs_GTFs/'+class_prev+'/' +class_prev+ '.gtf', 'wb')
	    c = open(base_dir + '161212_AV7M1_Vex-seq_Delta/Refs_GTFs/'+class_prev+'/' +class_prev+ '.fa', 'wb')
	    c.write('>' + class_prev + '\n')
	    full_seq = E1_on + seqs[BC][0] + I2_on + E3_dig_bits + 'NNNNNNNN' + rev_primer_XbaI
	    c.write(full_seq)
	    c.close()
	    attribute = 'transcript_id "%s"' % (class_prev)
            b.write('\t'.join([class_prev, '.', 'transcript', '1', str(len(full_seq)+1), '.', '+', '.', attribute+ '\n']))
            b.write('\t'.join([class_prev, '.', 'exon', '1', '44', '.', '+', '.', attribute + '; exon_id "exon_1"\n']))
	    b.write('\t'.join([class_prev, '.', 'exon',str(160+seqs[BC][0].find(seqs[BC][1])),  str(159+(len(seqs[BC][0]) - len(seqs[BC][2]))),'.', '+', '.', attribute + '; exon_id "cassette"\n']))
#								#new full sequence.find(E2)	#new full sequence - distal Intron length
	    seqs = {}
	    b.write('\t'.join([class_prev, '.', 'exon', str(1+ len(full_seq) - len(E3_dig_bits + 'BARCODE_' + rev_primer_XbaI)),str(1+ len(full_seq)), '.', '+', '.', attribute +'; exon_id "exon_3"']))
	    b.close()
	BC = str(line[11])
	seqs[BC] = get_seqs(line)
	class_prev = '_'.join(line[0:4])
a.close()
#for the last one
b = open(base_dir + '161212_AV7M1_Vex-seq_Delta/Refs_GTFs/' + class_curr + '/' +class_curr+ '.gtf', 'wb')
c = open(base_dir + '161212_AV7M1_Vex-seq_Delta/Refs_GTFs/' + class_curr + '/' +class_curr+ '.fa', 'wb')
c.write('>' + class_curr + '\n')
full_seq = E1_on + seqs[BC][0] + I2_on +E3_dig_bits + 'NNNNNNNN' + rev_primer_XbaI
c.write(full_seq)
c.close()
attribute = 'transcript_id "%s"' %(class_curr)
b.write('\t'.join([class_curr, '.', 'transcript', '1', str(len(full_seq)+1), '.', '+', '.', attribute+ '\n']))
b.write('\t'.join([class_curr, '.', 'Exon 1', '1', '44', '.', '+', '.', attribute +'; exon_id "exon_1"\n']))
b.write('\t'.join([class_curr, '.', 'exon',str(160+seqs[BC][0].find(seqs[BC][1])), str(159+(len(seqs[BC][0]) - len(seqs[BC][2]))),'.', '+', '.', attribute + '; exon_id "cassette"\n']))
b.write('\t'.join([class_curr, '.', 'Exon 3', str(len(full_seq) - len(E3_dig_bits + 'BARCODE_' + rev_primer_XbaI)+1),str(1+ len(full_seq)), '.', '+', '.', attribute +'; exon_id "exon_3"']))
b.close()

