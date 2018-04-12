novoindex -n /Users/sadamson/Desktop/scripts_and_such/Vex-seq/delta_data/170411_B3PGD_solid_amp/1o_ref_indels_novo 1o_ref_indels.fai 1o_ref_indels.fa
gunzip 1o_merged.extendedFrags.fastq.gz -c | novoalign -f /dev/stdin -F ILM1.8 -d ./1o_ref_indels_novo -o SAM |samtools view -b - >1o_indels_novo.bam
novoindex -n /Users/sadamson/Desktop/scripts_and_such/Vex-seq/delta_data/170411_B3PGD_solid_amp/1o_ref 1o_ref.fai 1o_ref.fa
gunzip 1o_merged.extendedFrags.fastq.gz -c | novoalign -f /dev/stdin -F ILM1.8 -d ./1o_ref_indels_novo -o SAM |samtools view -b - >1o_sorted.bam



