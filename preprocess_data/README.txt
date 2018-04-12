These are scripts that were used to preprocess the data, these are not neccessarily functional in the context of this repository.

These are more an illustration of how preprocessing was done; if you have any questions, email adamson@uchc.edu



Example RNA analysis pipeline:

bash setup_dirs.sh
bash STAR_index.sh
python sort_barcodes.py
bash STAR_align.sh
python splicing_overview_by_BC.py

Example of plasmid quality pipeline:

python make_1o_plasmid_ref.py
bash novoalign.sh
java -jar ~/git_clones/jvarkit/dist/sam2tsv.jar -r ../../161212_AV7M1_Vex-seq_Delta_flash/1o_ref.fa < 1o_sorted.bam | python parse_sam2tsv2.py
java -jar ~/git_clones/jvarkit/dist/sam2tsv.jar -r 1o_ref_indels.fa < 1o_indels_novo_sorted.bam | python parse_indel_alignment_novo.py
python BC_matches.py


