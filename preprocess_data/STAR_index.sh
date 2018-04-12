base_dir="/home/CAM/adamson/Vex-seq/single_way/68_97_pool/170411_B3PGD_solid_amp/Refs_GTFs"

readarray a < /home/CAM/adamson/Vex-seq/single_way/68_97_pool/170411_B3PGD_solid_amp/variants.txt

for entry in $a; do
    echo $entry;
    STAR --runMode genomeGenerate --runThreadN 8 --sjdbGTFfile $base_dir/$entry/$entry".gtf" --sjdbOverhang 259 --genomeDir $base_dir/$entry --genomeSAindexNbases 4 --genomeFastaFiles $base_dir/$entry/$entry".fa";
done
