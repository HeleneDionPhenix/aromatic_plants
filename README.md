# aromatic_plants

Effect of aromatic plants on nest and eggshell bacterial microbiota and nestling traits 
by Hélène Dion-Phénix & Gabrielle Gingras

This repository contains R scripts for testing the effect of aromatic plants on bacterial microbiota of nest samples in an experimental approach, and on bacterial microbiota of nest material and eggshell and on nestling traits and behavior in an observational approach. The code files are numbered and can be run sequentially, starting with the raw data, or any specific code can be run independently using intermediate *RData* objects. The general structure of the repository is given below. A README file is available for the data and the codes in the corresponding folders.

## Project Directory Structure

|─ aromatic_plants.Rproj <br> 
|─ LICENSE <br>
|─ README.md <br>
|─ code <br>
│&nbsp; &nbsp; &nbsp; |─ README.md <br>
│&nbsp; &nbsp; &nbsp; |─ 1.2-phyloseq_figures.R <br>
│&nbsp; &nbsp; &nbsp; |─ 5-repeatability_table.R <br>
│&nbsp; &nbsp; &nbsp; |─ exp<br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_formatting.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 1.1-phyloseq.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2.1-diversity_analysis.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2.2-diversity_figures.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2.2-diversity_tables.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3.1-community_analysis.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3.2-community_figures.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3.2-community_tables.R <br>
│&nbsp; &nbsp; &nbsp; |─ obs<br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_eggshell.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_nest.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_formatting.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 1.1-phyloseq_eggshell.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 1.1-phyloseq_nest.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2.1-diversity_analysis.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2.2-diversity_figures.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2.2-diversity_tables.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3.1-community_analysis.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3.2-community_figures.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3.3-community_tables.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 4.1-nestling_traits_analysis.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 4.2-nestling_traits_figures.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 4.3-nestling_traits_tables.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 5.1-complementary_analysis.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 5.2-complementary_figures.R <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 5.3-complementary_tables.R <br>
|─ data <br>
│&nbsp; &nbsp; &nbsp; |─ README.md <br>
│&nbsp; &nbsp; &nbsp; |─ exp<br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_nest.csv <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_sample.csv <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 1-phyloseq_objects.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2-output_models_diversity.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3-output_models_community.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- bioinfo <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_complete.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- fastq (not included) <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- seqtab.rds <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- taxa_dada2.txt <br>
│&nbsp; &nbsp; &nbsp; |─ obs<br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_eggshell_sample.csv <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_fledging.csv <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_nest_sample.csv <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data_nestlings.csv <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-data.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-nest_composition.csv <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 1-phyloseq_objects_eggshell.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 1-phyloseq_objects_nest.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 2-output_models_diversity.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 3-output_models_community.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 4-output_models_nestling_traits.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 5-output_models_complementary.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- bioinfo <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_eggshell_complete.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_eggshell.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_nest_complete.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_nest.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- fastq (not included) <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- seqtab_eggshell.rds <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- seqtab_nest.rds <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- taxa_eggshell_dada2.txt <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- taxa_nest_dada2.txt <br>
│&nbsp; &nbsp; &nbsp; |─ SILVA138.1 (not included) <br>
|─ figure <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_2.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_3.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_4.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_5.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_6.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S1.1.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S1.2.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S1.3.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S1.4.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S1.5.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S1.6.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S2.1.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S2.2.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S2.3.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S4.1.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S4.2.png <br>
│&nbsp; &nbsp; &nbsp; |─ Figure_S4.3.png <br>
|─ table <br>
│&nbsp; &nbsp; &nbsp; |─ Table_2_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_2.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_3_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_3.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_4_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_4.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_5_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_5.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S2.1_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S2.1.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.1_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.1.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.2_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.2.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.3_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.3.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.4_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S3.4.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S4.1_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S4.1.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S4.2_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S4.2.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S4.3_formatted.docx <br>
│&nbsp; &nbsp; &nbsp; |─ Table_S4.3.docx <br>

    
*The folder data/exp/bioinfo/fastq, data/obs/bioinfo/fastq, and data/SILVA138.1 need to be added to the directory where specified above to run codes code/exp/0-dada2.R, code/obs/0-dada2_eggshell.R, and code/obs/0-dada2_nest.R*

*To load the fastq files : https://figshare.com/projects/Aromatic_plants_decrease_nest_bacterial_diversity_and_improve_nestling_condition_in_Blue_Tits/230232*

*To load the SILVA Ribosomal RNA Gene Database Project, version 138.1: https://www.arb-silva.de/news/view/2024/07/11/silva-release-1382/*

