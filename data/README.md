# Aromatic plants - Data

## Metadata

### Experimental approach

The file *data/exp/0-data_nest.csv* give information about the natural nests that where predated or abandon by the parents and that where collected for the experiment.

| Variable name | Description | Values |
| :------------ | ----------- | ------ |
| nichoir | Nest box identifier | (ex: ari11) |
| population | Blue tit population | [DMuro, EMuro, or EPirio] D stand for Deciduous, E stand for Evergreen |
| station | Site name | [ari: Arinelle, ava: Avapessa, fil: Filagna, gra: Grassa, mur: Muro, pir: Pirio, tua: Tuarelli] |
| bloc | Temporal experimental block | All samples from a block processed at the same dates |
| stage_pheno | Phenological stage at collection | [prebrood: Pre-incubation, brood: Incubation, breed: Post-hatching] |
| stage_chick | Chick developmental stage at collection | [p: laying stage - cold eggs, w: incubated eggs, n: nestlings] |
| failure | Cause of the nest failure | [aba: abandonment by the parents, pre: predation of the clutch] |
| date_pret | Estimated date when the nest construction was achieve | (ex: 2023-04-05) |
| date_recup | Date of collection on the filed | (ex: 2023-05-12) |
| date_j1 | First day of the experiment | (ex: 2023-05-21) |

The file *data/exp/0-data_sample.csv* give information about each sample of the experiment.

| Variable name | Description | Values |
| :------------ | ----------- | ------ |
| n_ext | Sample identifier for the extraction | (ex: 1) |
| sample.names | Sample identifier | (ex. 1-pir54-lm) |
| type | Sample type | [nest: nest material sample, plant: fresh aromatic plant sample, ctrl_ext: negative extraction control, ctrl_PCR_neg: negative PCR control, ctrl_PCR_pos: positive PCR control - mock community] |
| nichoir | Nest box identifier | (ex: ari11) |
| pot | Fake nest box identifier | (ex: 15) unique in each temporal block |
| traitement | plant addition treatment | Single species treatment (ach: *Achillea ligustica*, pul: *Pulicaria odorata*, men: *Mentha suaveolens insularis*, lav: *Lavandula stoechas*, imm: *Helichrysum italicum*) Two species treatment: (lavmen: *Lavandula* & *Mentha*, lavimm: *Lavandula* & *Helichrysum*, immmen: *Helichrysum* & *Mentha*), control treatment: (ctrl: no plant addition) |
| mass_1sp | Mass of plant added of one species at day 1, 3, and 5 (g) | (ex: 0.408) same as mass_tot for single species treatment |
| mass_tot | Total mass of plant added at day 1, 3, and 5 (g) | (ex: 0.816) represented one third of the sample mass |
| category | Sample type (broader category) | [ech: nest or plant sample, ctrl: control samples] |


### Observational approach

The file *data/obs/0-nest_composition.csv* give the percent cover of the five aromatic plant species studied in this project.

| Variable name | Description | Values |
| :------------ | ----------- | ------ |
| id | Clutch and stage identifier | (ex: ARI2_NC) composed of clutch id and pehnological stage (NC: incubation, NE: post-hatching) |
| clutch | Clucth identifier | (ex: ARI2) composed of nest box id, and "-1" or "-2" when two clutch occupied the same nest box |
| nichoir | Nest box identifier | (ex: ARI2) |
| station | Site name | [ARI: Arinelle, AVA: Avapessa, FEL: Felicheto, FIL: Filagna, GRA: Grassa, MUR: Muro, PIR: Pirio, TUA: Tuarelli] |
| pop | Blue tit population | [DMuro, EMuro, or EPirio] D stand for Deciduous, E stand for Evergreen |
| date | Date of the observation | (ex: 2023-04-18) |
| stade | Phenological stage at observation | [couvaison: Incubation, elevage: Post-hatching] |
| obs | Initial of the observer | (ex: CD) |
| ACH | Percent cover of *Achillea ligustica* | [0, 0.5, 1, 2, 3, 4, 5] details below |
| PUL | Percent cover of *Pulicaria odorata* | [0, 0.5, 1, 2, 3, 4, 5] details below |
| MEN | Percent cover of *Mentha suaveolens insularis* | [0, 0.5, 1, 2, 3, 4, 5] details in table below |
| LAV | Percent cover of *Lavandula stoechas* | [0, 0.5, 1, 2, 3, 4, 5] details in table below |
| IMM | Percent cover of *Helichrysum italicum* | [0, 0.5, 1, 2, 3, 4, 5] details in table below |

Percent cover classes
| code | description | median |
| :--- | ----------- | ------ |
| 0.5 | 1-2 fragments, < 5% | 1% |
| 1 | 3-5 fragments, < 5%  | 3% |
| 2 | 5%-25% | 15% |
| 3 | 25%-50% | 37.5% |
| 4 | 50%-75% | 62.5% |
| 5 | 75%-100% | 87.5% |

The file *data/obs/0-data_eggshell_sample.csv* give information on eggshell swab samples.
The file *data/obs/0-data_nest_sample.csv* give information on nest material samples.

| Variable name | Description | Values |
| :------------ | ----------- | ------ |
| category | Sample type (broader category) | [nid: nest sample, swab: eggshell swab sample, ctrl: control samples] |
| num_ext | Clucth identifier | (ex: ARI2_NC) composed of nest box id and pehnological stage (NC: incubation, NE: post-hatching) |
| an | Year | (2023) |
| nichoir | Nest box identifier | (ex: ARI2) |
| date | Date of sampling | (ex: 2023-04-25) |
| type | Sample type | [nid_couv: nest material sampled during incubation, nid_elev: nest material sampled post-hatching, swab_w: eggshell swab, field_ctrl: Field negative control, ext_ctrl: negative extraction control, PCR_ctrl: PCR control] |
| nom_cermo | Sample identifier | (ex: 1S-4-2023-FIL12W) |
| clutch | Clucth identifier | (ex: ARI2) composed of nest box id, and "-1" or "-2" when two clutch occupied the same nest box |
| id | Clutch and stage identifier | (ex: ARI2_NC) composed of clutch id and pehnological stage (NC: incubation, NE: post-hatching) |

The file *data/obs/0-data_nestlings.csv* give information on nest material samples.

| Variable name | Description | Values |
| :------------ | ----------- | ------ |
| chick_id | Nestling identifier | (ex: chick1) |
| clutch | Clucth identifier | (ex: ARI2) composed of nest box id, and "-1" or "-2" when two clutch occupied the same nest box |
| poids | Body mass (g) | (ex: 10.4) |
| tarsed | Right tarsus length (mm) | (ex: 15.7) |
| date_mesure | Date of sampling | (ex: 2023-04-29) |
| age_plume | Feather development age (day) | (ex. 15) Aging based on typical feather development |
| docilite | Aggressivness in hand | Number of movements nestlings do in hand in 10 seconds |
| date_eclo | Hatching date | (ex: 2023-04-29) |

The file *data/obs/0-data_fledging.csv* give information on nest material samples.

| Variable name | Description | Values |
| :------------ | ----------- | ------ |
| nichoir | Nest box identifier | (ex: ARI6) |
| pulenv | Number of nestlings fledgling in the clutch | (ex. 1) |
| echec_fledg | Number of nestlings not fledgling in the clutch | (ex: 9) |
| pop | Blue tit population | [DMuro, EMuro, or EPirio] D stand for Deciduous, E stand for Evergreen |
| an | Year | Between 2004 and 2023 |

## Fastq files

The fastq files containing the 16S ARNr sequences to run these analysis are available here:
https://figshare.com/projects/Aromatic_plants_decrease_nest_bacterial_diversity_and_improve_nestling_condition_in_Blue_Tits/230232

To assign the taxonomy, use the SILVA Ribosomal RNA Gene Database Project,
version 138.1:
https://www.arb-silva.de/news/view/2024/07/11/silva-release-1382/

## Intermediate data files

Alternatively, you can use the .RData files to run some specific 
part of the code without needing to download the raw data.

## Data Directory Structure

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
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_eggshell.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- 0-dada2_nest.RData <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- fastq (not included) <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- seqtab_eggshell.rds <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- seqtab_nest.rds <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- taxa_eggshell_dada2.txt <br>
│&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; │&nbsp; &nbsp; &nbsp; |-- taxa_nest_dada2.txt <br>
│&nbsp; &nbsp; &nbsp; |─ SILVA138.1 (not included) <br>