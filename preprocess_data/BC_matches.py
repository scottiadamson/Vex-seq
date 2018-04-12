import csv

a=open('correct_BC_custom_novo_final.tsv', 'r')
reader =csv.reader(a, delimiter = '\t')
#chr1_16069287_C_T	613	11	ACCTACTT:109,CGTCCTCG:293,CTTCGGTC:211	ACCTTCTT:1,CATTCTCG:1,GTCCAACA:1,CTGTATTG:1,ACCTACTA:1,TCATCGAG:1,TGCGACTA:1,GGCACGCC:1,GGCACGCG:1,CGTCCTCC:1,GGGTAGCA:1
next(reader);correct_BC = {}; incorrect_BC = {}; all_BCs = {}
for line in reader:
    if str(line[3]) != '':
    	for BC in str(line[3]).split(','):
	    correct_BC[str(BC.split(':')[0])] = int(BC.split(':')[1])
	    all_BCs[str(BC.split(':')[0])] =True
    if str(line[4]) != '':
	for BC in str(line[4]).split(','):
	    all_BCs[str(BC.split(':')[0])] =True
	    if str(BC.split(':')[0]) not in incorrect_BC:
	    	incorrect_BC[str(BC.split(':')[0])] = int(BC.split(':')[1])
 	    else:
	    	incorrect_BC[str(BC.split(':')[0])] += int(BC.split(':')[1])
a.close()
b=open('BC_correctness_novo_test.tsv', 'wb')
writer=csv.writer(b, delimiter = '\t')
writer.writerow(['BC', 'correct', 'incorrect', 'total', 'correct_PCT'])
for BC in all_BCs:
    if BC not in correct_BC:
	correct_BC[BC] =0
    if BC not in incorrect_BC:
	incorrect_BC[BC] =0
    writer.writerow([BC, correct_BC[BC], incorrect_BC[BC],correct_BC[BC]+ incorrect_BC[BC], float(correct_BC[BC])/(correct_BC[BC] + incorrect_BC[BC]) *100])
b.close()

