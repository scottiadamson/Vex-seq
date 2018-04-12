#PBS -l nodes=8:ppn=12,walltime=100:00:00
#PBS -m e
#PBS -M adamson@uchc.edu

base_dir="/home/CAM/adamson/Vex-seq/single_way/68_97_pool/170411_B3PGD_solid_amp/Refs_GTFs"

readarray a < /home/CAM/adamson/Vex-seq/single_way/68_97_pool/170411_B3PGD_solid_amp/variants.txt
samples="HepG2-1_S4 HepG2-2_S5 HepG2-3_S6 K562-1_S1 K562-2_S2 K562-3_S3"

for entry in $a; do
    echo $entry;
    for sample in $samples; do
	STAR --genomeDir $base_dir/$entry --runThreadN 8 --readFilesIn $base_dir/$entry/$sample"_L001_R1_001.fastq" $base_dir/$entry/$sample"_L001_R2_001.fastq" --sjdbGTFfile $base_dir/$entry/$entry.gtf --outFileNamePrefix $base_dir/$entry/$sample"_"$entry --outSAMattributes jI;
    done;
done

