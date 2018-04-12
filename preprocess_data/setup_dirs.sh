dir="/home/CAM/adamson/Vex-seq/single_way/68_97_pool/161212_AV7M1_Vex-seq_Delta/Refs_GTFs"
readarray a < $dir/variants.txt
current_dir="/home/CAM/adamson/Vex-seq/single_way/68_97_pool/170411_B3PGD_solid_amp/Refs_GTFs/"

for entry in $a; do
    mkdir $current_dir/$entry;
    cp $dir/$entry/$entry".fa" $current_dir/$entry;
    cp $dir/$entry/$entry".gtf" $current_dir/$entry;
done


