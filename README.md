# variant-effect-predictors
Collection of Variant Effect prediction tools

| Tool         | Year         | Model        | Features     | Scope        | Predicts     | Impacts      | Training     |
| :---         | :---         | :---         | :---         | :---         | :---         | :---         | :---         |
| ENVISION | 2018 | Stochastic Gradient Boosting  | SEQ,STR | nsSNP | F | protein | DMS3 [21,026 variants in eight proteins] |
| MutPred2 | 2017 | Neural Networks (FF1) | SEQ | nsSNP | P* | protein | HGMD + SwissVar5 + dbSNP6 + inter-species pairwise alignments [53,180 pathogenic / 206,946 unlabeled (putatively neutral) variants] |
| REVEL | 2016 | RandomForest  | ENS | nsSNP | P | protein | HGMD + ESP + ARIC + dbNSFP [6,182 disease-related variants / 123,706 rare neutral exome sequencing variants] |
| PANTHER-PSEP | 2016 | Phylogenetic Analysis | SEQ | nsSNP | P | protein | HumVar (cite: doi:10.1093/bioinformatics/btl423) |
| SNAP2 | 2015 | Neural Networks (FF1) | SEQ | nsSNP | F | protein | PMD + Swiss-Prot + OMIM + HumVar [61,037 effect and 40,478 neutral in 9,744 proteins] |
| PON-P2 | 2015 | RandomForest  | SEQ,STR | nsSNP | P | protein | dbSNP, VariBench [14,086 pathogenic variants in 1,082 proteins / 14,848 neutral variants in 6,598 proteins] |
| wKinMut-2 | 2015 | Annotation Summary | ENS,LIT,KB | nsSNP | P | protein | UniProt/Swiss-Prot [disease dataset: 865 mutations in 65 proteins / neutral dataset: 2,627 mutations in 447 human proteins) |
| PredictSNP | 2014 | Consensus Scoring | ENS | nsSNP | F | protein | Swiss-Prot [SNPs&GO dataset; 58,057 mutations], Swiss-Prot + HGMD [MutPred dataset, 65,654 mutations], dbSNP + PhenCode + Idbases + 16 individual locus- specific databases [PON-P datase, 39,670 mutations], Humsavar [36,994 neutral and disease-related mutations] |
| FATHMM-DS | 2014 | Hidden Markov Models  | SEQ | nsSNP | F | protein | HGMD [damaging mutations], Swiss-Prot/TrEMBL [putative neutral polymorphisms] |
| PolyPhen-2 | 2013 | Naïve Bayes Classifier | SEQ,STR | nsSNP | S,F | protein | HumDiv or HumVar [7,070 neutral and 5,322 deleterious variants] |
| FATHMM | 2013 | Hidden Markov Models  | SEQ | nsSNP | F | protein | HGMD [filtered for "damaging mutations" annotation] |
| VEST | 2013 | RandomForest  | SEQ | nsSNP | F | protein | PolyPhen-2 (v2.2.2) training set [HumDiv or HumVar] |
| FATHMM-cancer | 2013 | Hidden Markov Models  | SEQ | nsSNP | F | protein | CanProVar [12,720 positives], UniProt [36,928 negatives] |
| Meta-SNP | 2013 | RandomForest  | ENS | nsSNP | P | protein | SwissVar [SV-2009,  35,766 nsSNVs (17,883 disease variants) from 8,667 proteins] |
| PON-P | 2012 | RandomForest  | ENS | nsSNP | P | protein | PhenCode (2009) + IDbases + 16 individual locus- specific databases (LSDB) [14,610 pathogenic, manual curation, SwissVar/LSDB disease annotations], dbSNP [Build 131, 17,393 neutral] |
| KinMut | 2012 | Support Vector Machine | SEQ | nsSNP | P | protein | Swiss-Prot |
| MutationAssessor | 2011 | Functional Impact Scoring | SEQ,STR | nsSNP | F | protein | UniProt |
| MutPred | 2009 | RandomForest  | SEQ | nsSNP | S,F | protein | HGMD |
| PoPMuSiC-2.0 | 2009 | Energy Function | SEQ,STR,KB | nsSNP | S | protein | ProTherm [2,648 different point mutations, 131 proteins] |
| iMutant3 | 2008 | Support Vector Machine | SEQ|STR | nsSNP | S,F | protein | ProTherm [6,398 mutations, 55 proteins] |
| SNAP | 2007 | Neural Networks (FF1) | SEQ | nsSNP | F | protein | Swiss-Prot, Protein Mutant Database (PMD) [40,830 neutral and 39,987 deleterious variants] |
| PhD-SNP | 2006 | Support Vector Machine | SEQ | nsSNP | P | protein | Swiss-Prot/HumVar [12,944 disease-related and 8241 neutral polymorphisms, 3587 proteins] |
| Align-GVGD | 2006 | Extended GD2 Scoring | SEQ,STR | nsSNP | F | protein | p53 protein MSA  |
| FoldX | 2005 | FoldX Force Field | STR | nsSNP | S | protein | not applicable |
| iMutant2 | 2005 | Support Vector Machine | SEQ/STR | nsSNP | S | protein | ProTherm [2087 single mutations in 65 proteins] |
| MAPP | 2005 | Functional Impact Scoring | SEQ | nsSNP | F | protein | not applicable |
| nsSNPAnalyzer | 2005 | RandomForest  | SEQ,STR | nsSNP | P | protein | Swiss-Prot |
| PolyPhen | 2002 | Rule-based Classifier | SEQ,STR | nsSNP | S,F | protein | No training dataset, NRDB for alignments |
| SIFT | 2001 | PSSM based Probabilities | SEQ | nsSNP | D | protein | No training dataset , Swiss-Prot for MSA |
| PoPMuSiC | 2000 | Energy Function | SEQ,STR,KB | nsSNP | S | protein | Set of 141 high-resolution (􏰹2.5 Å) protein X-ray structures (< 􏰺25% sequence identity) |
| IDSV | 2019 | RandomForest  | SEQ | sSNP | D | protein | dbDSM (ClinVar, PubMed, Web of Knowledge variants as disease-causing [300 variants]) and VariSNP [300 neutral varians] |
| ARVIN | 2018 | RandomForest  | SEQ,NET | SNV | P | regulatory | HGMD |
| LINSIGHT | 2017 | Linear + Probabilistic Model  | SEQ | SNV | P | regulatory | High-coverage genome sequences for 54 unrelated individuals from the "69 Genome" data set from Complete Genomics |
| DDIG-SN | 2017 | Support Vector Machine | SEQ | sSNP | P | protein | HGMD [592 disease-causing variants] and 1kGP [10,925 putatively benign variants] in 318 genes |
| regSNPs-splicing | 2017 | RandomForest  | SEQ | sSNP | F,SP | protein | HGMD [1373 disease-causing synonymous SNVs (sSNVs)] and 1000 Genomes Project [7231 neutral VIE variants,  329 VSS variants] |
| GWAVA | 2014 | RandomForest  | SEQ,KB | SNV | P | regulatory | HGMD [1,614 disease-implicated] and 1000 Genomes Project [5,027 variants] |
| SilVA | 2013 | RandomForest  | SEQ | sSNP | P | protein | Curated, literature based dataset of rare (allele frequency <5%) synonymous variants [33 variants] and 1000 Genomes Project [746 rare synonymous variants in one individual] |
| Cscape | 2017 | Multiple Kernel Learning | SEQ | SNV | P | protein + regulatory | COSMIC (pathogenic variants) and 1000 Genomes Project (control) [46,420 coding examples and 131,714 non-coding examples] |
| DANN | 2015 | Deep Neural Network | SEQ | SNV | D | protein + regulatory | 16,627,775 observed variants and 49,407,057 simulated variants  |
| FATHMM-MKL | 2015 | Multiple Kernel Learning | SEQ,KB | SNV | F | protein + regulatory | coding: HGMD [17,362 coding; 3063 non-coding] and 1000 Genomes Project [4853 coding; 5252 non-coding] |
| CADD | 2014 | Support Vector Machine | SEQ | SNV/Indel | D | protein + regulatory | observed (14,893,290 SNVs, 627,071 insertions and 1,107,414 deletions) and simulated variants |
| MutationTaster2 | 2014 | Naïve Bayes Classifier | SEQ,KB | SNV/Indel | P | protein + regulatory | HGMD + ClinVar [>100,000 disease-associated mutations] and 1000 Genomes Project [>6,000,000 single base exchanges and short indels] |
| MutationTaster | 2010 | Naïve Bayes Classifier | SEQ,KB | SNV/Indel | P | protein | Ensembl, dbSNP, HapMap, UniProt/Swiss-Prot  |
| PROVEAN | 2012 | Delta Alignments Scoring | SEQ | SNVs/Indels | F | protein + regulatory | NCBI NR for MSAs |
| SIFT Indel | 2012 | Decision Tree | SEQ,KB | Indel | D | protein | HGMD (2010.2) [1,292 disease indels] 1000 Genomes Project [2,602 neutral indels] |
