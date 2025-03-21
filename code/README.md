# Aromatic plants - Code

The code files are numbered and can be run sequentially, starting with the raw data, or any specific code can be run independently using intermediate *RData* objects.

## Code file descriptions

#### 1.2-phyloseq_figures.R

Creation of rarefaction curves and visual comparison of bacterial composition between controls and samples, between PCR positive and negative controls, and between sample groups.

#### 5-repeatability_table.R

Creation of a table of nest repeatability for bacterial diversity and composition from the experimental and observational studies, and for the quantity of aromatic plants female bring to the nest. This table reported the adjusted and non-adjusted repeatability for the random effect "nest" from the models created in the following codes: exp/2.1-diversity_analysis.R, exp/3.1-community_analysis.R, obs/2.1-diversity_analysis.R, and obs/3.1-community_analysis.R.

### exp

#### 0-dada2.R

DADA2 analysis following Kembel lab script https://github.com/kembel-lab/scripts/tree/91225181acba420f2bb01acba1aa8fa4806255a1/dada2
based on dada2 tutorial https://benjjneb.github.io/dada2/tutorial.html

#### 0-data_formatting.R

R formatting of *.csv* file. This code generate the *data/exp/0-data.RData* file used in different subsequent scripts.

#### 1.1-phyloseq.R

Rarefaction and creation of physloseq object from ASVs tables and metadata.
For more details: https://github.com/joey711/phyloseq

#### 2.1-diversity_analyses.R

Comparison of bacterial Shannon diversity among treatments and phenological stages, controlling for the nest of origin for nest blue tit bacterial microbiota using a linear mixed model.

##### 2.2-diversity_figures.R

Visual representation of the model created in the exp/2.1-diversity_analysis.R script showing raw data and model estimates with 95% confidence intervals.

##### 2.3-diversity_tables.R

Model estimates and 95% confidence intervals of the model created in the exp/2.1-diversity_analysis.R script.

#### 3.1-community_analyses

Tests of multivariate homogeneity of dispersion to compare samples homogeneity between control and treated nest samples, among treatments, among phenological stages, and among nest of origin for bacterial community of nest samples.
Distance-based RDA (dbRDA) on the Hellinger distance to test if the bacterial composition of nest samples is explained by plant addition and plant species, both in interaction with phenological stage and among nest of origin.

##### 3.2-community_figures.R

Representation of the model created in the exp/3.1-community_analysis.R script visualised in a dbRDA ordination.

##### 3.3-community_tables.R

Proportion (adjusted) of the total variation in the composition among nest samples explained by plant addition and plant species, both in interaction with phenological stage and among nest of origin.

### obs

#### 0-dada2_eggshell.R
#### 0-dada2_nest.R

DADA2 analysis following Kembel lab script https://github.com/kembel-lab/scripts/tree/91225181acba420f2bb01acba1aa8fa4806255a1/dada2
based on dada2 tutorial https://benjjneb.github.io/dada2/tutorial.html

#### 0-data_formatting.R

R formatting of *.csv* file. This code generate the *data/obs/0-data.RData* file used in different subsequent scripts.

#### 1.1-phyloseq_eggshell.R
#### 1.1-phyloseq_nest.R

Rarefaction and creation of physloseq object from ASVs tables and metadata.
For more details: https://github.com/joey711/phyloseq

#### 2.1-diversity_analyses.R

*Nest:* Relation between the bacterial Shannon diversity in nest material and the quantity of aromatic plants added by the female blue tit in interaction with the phenological stage at collection and the blue tit population (triple interaction), controlling for the nest identity using a linear mixed model.
*Eggshell:* Relation between the bacterial Shannon diversity on eggshell in incubation and the quantity of aromatic plants added by the female blue tit in interaction with the blue tit population using a linear model.
*Aromatic plant quantity:* Comparison of the aromatic plant quantities female blue tits added to their nest between the three studied populations, controlling for nest identity using a linear mixed model.

##### 2.2-diversity_figures.R

Visual representation of the three models created in the obs/2.1-diversity_analysis.R script showing raw data and model estimates with 95% confidence intervals.

##### 2.3-diversity_tables.R

Model estimates and 95% confidence intervals of the three models created in the obs/2.1-diversity_analysis.R script.

#### 3.1-community_analyses

Tests of multivariate homogeneity of dispersion to compare samples homogeneity among population, among phenological stages, and among nest identity for bacterial community of nest and eggshell samples.
*Nest:* Distance-based RDA (dbRDA) on the Hellinger distance to test if the bacterial composition of nest samples is explained by the quantity of aromatic plants in interaction with the phenological stage and the blue tit population and by the nest identity.
*Eggshell:* Distance-based RDA (dbRDA) on the Hellinger distance to test if the bacterial composition of eggshell samples is explained by the quantity of aromatic plants in interaction with the the blue tit population.

##### 3.2-community_figures.R

Representation of the dbRDA model on nest bacterial community created in the exp/3.1-community_analysis.R script visualised in a ordination.

##### 3.3-community_tables.R

Proportion (adjusted) of the total variation in the composition among nest samples explained by the quantity of aromatic plants in interaction with the phenological stage and the blue tit population and by the nest identity.

#### 4.1-nestling_traits_analyses

Principal Coordinate Analysis (PCoA) using euclidean distance on the blue tit nestling traits (body mass, tarsus length, developmental index, and aggressiveness in hand).
Relation between the first and second axis of the PCoA and the quantity of aromatic plants in interaction with the blue tit population, controlling for the clutch identity using two linear mixed models (one for each axis).

##### 4.2-nestling_traits_figures.R

Visualization of the PCoA ordination and of the two models created in the obs/4.1-nestling_traits_analysis.R script showing raw data and model estimates with 95% confidence intervals.

##### 4.3-nestling_traits_tables.R

Model estimates and 95% confidence intervals of the two models created in the obs/4.1-nestling_traits_analysis.R script.

#### 5.1-complementary_analyses

*Developmental index:* Relation between developmental index and hatching date and population, controlling for the clutch identity using a linear mixed model.
*Aromatic plant quantity:* Relation between aromatic plant quantity and hatching date and population using a linear mixed model.
*Fledging success:* Mean fledging success of clutches between 2003 and 2023 as assessed by a binomial linear model of the number of fledging successes compared to fledging failures by population and year, and mean fledging success over 20 years as assessed by a binomial linear model of the number of fledging successes compared to fledging failures by population.

##### 5.2-complementary_figures.R

Visualization of the models created in the obs/5.1-complementary_analysis.R script showing raw data and model estimates with 95% confidence intervals.

##### 5.3-complementary_tables.R

Model estimates and 95% confidence intervals of the models created in the obs/5.1-complementary_analysis.R script.

## Code Directory Structure

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