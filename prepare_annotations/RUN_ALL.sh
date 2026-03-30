



python PT_bacterial_assembly_interogation/prepare_annotations/parse_eggnog_annotations.py \
    --eggnog_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation/eggnog_results \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/eggnog_parsed



python PT_bacterial_assembly_interogation/prepare_annotations/parse_dbcan_results.py \
    --dbcan_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation/dbcan_results \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/dbcan_parsed


python PT_bacterial_assembly_interogation/prepare_annotations/parse_antismash_results.py \
    --antismash_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/functional_annotation/antismash_results \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/antismash_parsed


python PT_bacterial_assembly_interogation/prepare_annotations/build_ko_matrices.py \
    --ko_long_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/eggnog_parsed/eggnog_ko_long.tsv \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/ko_matrices


 python PT_bacterial_assembly_interogation/prepare_annotations/screen_pgpr_genes.py    \
    --eggnog_master_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/eggnog_parsed/eggnog_master_annotations.tsv  \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/pgpr_screen


python PT_bacterial_assembly_interogation/prepare_annotations/plot_ko_pcoa_and_intersections.py \
    --ko_matrix_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/ko_matrices/ko_presence_absence_matrix.tsv \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/ko_pcoa_and_intersections \
    --top_n_intersections 20


python PT_bacterial_assembly_interogation/prepare_annotations/ko_group_enrichment.py \
    --ko_matrix_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/ko_matrices/ko_presence_absence_matrix.tsv \
    --metadata_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/metadata/isolate_groups.tsv \
    --group_column group \
    --group_a rhizosphere \
    --group_b bulk \
    --out_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/ko_enrichment/ko_group_enrichment.tsv



python PT_bacterial_assembly_interogation/prepare_annotations/plot_pgpr_heatmap.py \
    --pgpr_matrix_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/pgpr_screen/pgpr_presence_absence_matrix.tsv \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/pgpr_heatmap


###################
# GO

python PT_bacterial_assembly_interogation/prepare_annotations/build_go_matrices.py \
    --go_long_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/eggnog_parsed/eggnog_go_long.tsv \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/go_matrices


python plot_go_pcoa_and_intersections.py \
    --go_matrix_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/go_matrices/go_presence_absence_matrix.tsv \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/go_pcoa_and_intersections \
    --top_n_intersections 20

##################
# dbCAN
python summarise_dbcan_results.py \
    --cazy_matrix_tsv /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/dbcan_parsed/dbcan_cazy_family_presence_absence_matrix.tsv \
    --out_dir /home/pthorpe001/david_bul/pthorpe001/2026_20th_March_Microbes_NG/post_annotation/dbcan_summary \
    --top_n_families 40