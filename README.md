# Collection of Variant Effect prediction tools
Last updated: 02/10/20 [#45]

| Tool         | Year         | Model        | Features     | Scope        | Predicts     | Impacts      | Training     | Availability |
| :---         | :---         | :---         | :---         | :---         | :---         | :---         | :---         | :---         |
| ENVISION | 2018 | Stochastic Gradient Boosting  | SEQ,STR | nsSNP | F | protein | DMS[3] [21,026 variants in eight proteins] | https://envision.gs.washington.edu/shiny/envision_new/ |
| MutPred2 | 2017 | Neural Networks (FF[1]) | SEQ | nsSNP | P* | protein | HGMD + SwissVar + dbSNP + inter-species pairwise alignments [53,180 pathogenic / 206,946 unlabeled (putatively neutral) variants] | http://mutpred.mutdb.org |
| REVEL | 2016 | RandomForest  | ENS | nsSNP | P | protein | HGMD + ESP + ARIC + dbNSFP [6,182 disease-related variants / 123,706 rare neutral exome sequencing variants] | https://sites.google.com/site/revelgenomics |
| PANTHER-PSEP | 2016 | Phylogenetic Analysis | SEQ | nsSNP | P | protein | HumVar (cite: doi:10.1093/bioinformatics/btl423) | http://pantherdb.org/tools/csnpScoreForm.jsp |
| SNAP2 | 2015 | Neural Networks (FF[1]) | SEQ | nsSNP | F | protein | PMD + Swiss-Prot + OMIM + HumVar [61,037 effect and 40,478 neutral in 9,744 proteins] | https://rostlab.org/services/snap2web |
| PON-P2 | 2015 | RandomForest  | SEQ,STR | nsSNP | P | protein | dbSNP, VariBench [14,086 pathogenic variants in 1,082 proteins / 14,848 neutral variants in 6,598 proteins] | http://structure.bmc.lu.se/PON-P2 |
| wKinMut-2 | 2015 | Annotation Summary | ENS,LIT,KB | nsSNP | P | protein | UniProt/Swiss-Prot [disease dataset: 865 mutations in 65 proteins / neutral dataset: 2,627 mutations in 447 human proteins) | https://github.com/Rbbt-Workflows/kin_mut2 |
| PredictSNP | 2014 | Consensus Scoring | ENS | nsSNP | F | protein | Swiss-Prot [SNPs&GO dataset; 58,057 mutations], Swiss-Prot + HGMD [MutPred dataset, 65,654 mutations], dbSNP + PhenCode + Idbases + 16 individual locus- specific databases [PON-P datase, 39,670 mutations], Humsavar [36,994 neutral and disease-related mutations] | http://loschmidt.chemi.muni.cz/predictsnp |
| FATHMM-DS | 2014 | Hidden Markov Models  | SEQ | nsSNP | F | protein | HGMD [damaging mutations], Swiss-Prot/TrEMBL [putative neutral polymorphisms] | http://fathmm.biocompute.org.uk/disease.html |
| PolyPhen-2 | 2013 | Naïve Bayes Classifier | SEQ,STR | nsSNP | S,F | protein | HumDiv or HumVar [7,070 neutral and 5,322 deleterious variants] | http://genetics.bwh.harvard.edu/pph2 |
| FATHMM | 2013 | Hidden Markov Models  | SEQ | nsSNP | F | protein | HGMD [filtered for "damaging mutations" annotation] | http://fathmm.biocompute.org.uk |
| VEST | 2013 | RandomForest  | SEQ | nsSNP | F | protein | PolyPhen-2 (v2.2.2) training set [HumDiv or HumVar] | http://www.cravat.us/CRAVAT |
| FATHMM-cancer | 2013 | Hidden Markov Models  | SEQ | nsSNP | F | protein | CanProVar [12,720 positives], UniProt [36,928 negatives] | http://fathmm.biocompute.org.uk/cancer.html |
| Meta-SNP | 2013 | RandomForest  | ENS | nsSNP | P | protein | SwissVar [SV-2009,  35,766 nsSNVs (17,883 disease variants) from 8,667 proteins] | http://snps.biofold.org/meta-snp |
| PON-P | 2012 | RandomForest  | ENS | nsSNP | P | protein | PhenCode (2009) + IDbases + 16 individual locus- specific databases (LSDB) [14,610 pathogenic, manual curation, SwissVar/LSDB disease annotations], dbSNP [Build 131, 17,393 neutral] | https://www.ncbi.nlm.nih.gov/pubmed/22505138 |
| KinMut | 2012 | Support Vector Machine | SEQ | nsSNP | P | protein | Swiss-Prot | http://kinmut.bioinfo.cnio.es |
| MutationAssessor | 2011 | Functional Impact Scoring | SEQ,STR | nsSNP | F | protein | UniProt | http://mutationassessor.org/r3 |
| MutPred | 2009 | RandomForest  | SEQ | nsSNP | S,F | protein | HGMD | http://mutpred1.mutdb.org/ |
| PoPMuSiC-2.0 | 2009 | Energy Function | SEQ,STR,KB | nsSNP | S | protein | ProTherm [2,648 different point mutations, 131 proteins] | http://babylone.ulb.ac.be/popmusic |
| iMutant3 | 2008 | Support Vector Machine | SEQ|STR | nsSNP | S,F | protein | ProTherm [6,398 mutations, 55 proteins] | http://gpcr2.biocomp.unibo.it/cgi/predictors/I-Mutant3.0/I-Mutant3.0.cgi |
| SNAP | 2007 | Neural Networks (FF[1]) | SEQ | nsSNP | F | protein | Swiss-Prot, Protein Mutant Database (PMD) [40,830 neutral and 39,987 deleterious variants] | https://hub.docker.com/r/bromberglab/snap |
| PhD-SNP | 2006 | Support Vector Machine | SEQ | nsSNP | P | protein | Swiss-Prot/HumVar [12,944 disease-related and 8241 neutral polymorphisms, 3587 proteins] | http://snps.biofold.org/phd-snp/phd-snp.html |
| Align-GVGD | 2006 | Extended GD[2] Scoring | SEQ,STR | nsSNP | F | protein | p53 protein MSA  | http://agvgd.iarc.fr |
| FoldX | 2005 | FoldX Force Field | STR | nsSNP | S | protein | not applicable | http://foldxsuite.crg.eu |
| iMutant2 | 2005 | Support Vector Machine | SEQ/STR | nsSNP | S | protein | ProTherm [2087 single mutations in 65 proteins] | http://gpcr.biocomp.unibo.it/cgi/predictors/I-Mutant2.0/I-Mutant2.0.cgi |
| MAPP | 2005 | Functional Impact Scoring | SEQ | nsSNP | F | protein | not applicable | http://mendel.stanford.edu/supplementarydata/stone_MAPP_2005 |
| nsSNPAnalyzer | 2005 | RandomForest  | SEQ,STR | nsSNP | P | protein | Swiss-Prot | http://snpanalyzer.uthsc.edu |
| PolyPhen | 2002 | Rule-based Classifier | SEQ,STR | nsSNP | S,F | protein | No training dataset, NRDB for alignments | http://genetics.bwh.harvard.edu/pph |
| SIFT | 2001 | PSSM based Probabilities | SEQ | nsSNP | D | protein | No training dataset , Swiss-Prot for MSA | https://sift.bii.a-star.edu.sg |
| PoPMuSiC | 2000 | Energy Function | SEQ,STR,KB | nsSNP | S | protein | Set of 141 high-resolution (2.5 Å) protein X-ray structures (< 25% sequence identity) | http://babylone.ulb.ac.be/popmusic |
| IDSV | 2019 | RandomForest  | SEQ | sSNP | D | protein | dbDSM (ClinVar, PubMed, Web of Knowledge variants as disease-causing [300 variants]) and VariSNP [300 neutral varians] | http://bioinfo.ahu.edu.cn:8080/IDSV |
| ARVIN | 2018 | RandomForest  | SEQ,NET | SNV | P | regulatory | HGMD | https://github.com/gaolong/arvin |
| LINSIGHT | 2017 | Linear + Probabilistic Model  | SEQ | SNV | P | regulatory | High-coverage genome sequences for 54 unrelated individuals from the "69 Genome" data set from Complete Genomics | http://compgen.cshl.edu/~yihuang/LINSIGHT |
| DDIG-SN | 2017 | Support Vector Machine | SEQ | sSNP | P | protein | HGMD [592 disease-causing variants] and 1kGP [10,925 putatively benign variants] in 318 genes | https://sparks-lab.org/server/ddig |
| regSNPs-splicing | 2017 | RandomForest  | SEQ | sSNP | F,SP | protein | HGMD [1373 disease-causing synonymous SNVs (sSNVs)] and 1000 Genomes Project [7231 neutral VIE variants,  329 VSS variants] | http://regsnps-splicing.ccbb.iupui.edu |
| GWAVA | 2014 | RandomForest  | SEQ,KB | SNV | P | regulatory | HGMD [1,614 disease-implicated] and 1000 Genomes Project [5,027 variants] | https://www.sanger.ac.uk/science/tools/gwava |
| SilVA | 2013 | RandomForest  | SEQ | sSNP | P | protein | Curated, literature based dataset of rare (allele frequency <5%) synonymous variants [33 variants] and 1000 Genomes Project [746 rare synonymous variants in one individual] | http://compbio.cs.toronto.edu/silva |
| Cscape | 2017 | Multiple Kernel Learning | SEQ | SNV | P | protein + regulatory | COSMIC (pathogenic variants) and 1000 Genomes Project (control) [46,420 coding examples and 131,714 non-coding examples] | http://CScape.biocompute.org.uk |
| DANN | 2015 | Deep Neural Network | SEQ | SNV | D | protein + regulatory | 16,627,775 observed variants and 49,407,057 simulated variants  | https://cbcl.ics.uci.edu/public_data/DANN |
| FATHMM-MKL | 2015 | Multiple Kernel Learning | SEQ,KB | SNV | F | protein + regulatory | coding: HGMD [17,362 coding; 3063 non-coding] and 1000 Genomes Project [4853 coding; 5252 non-coding] | http://fathmm.biocompute.org.uk/fathmmMKL.htm |
| CADD | 2014 | Support Vector Machine | SEQ | SNV/Indel | D | protein + regulatory | observed (14,893,290 SNVs, 627,071 insertions and 1,107,414 deletions) and simulated variants | https://cadd.gs.washington.edu |
| MutationTaster2 | 2014 | Naïve Bayes Classifier | SEQ,KB | SNV/Indel | P | protein + regulatory | HGMD + ClinVar [>100,000 disease-associated mutations] and 1000 Genomes Project [>6,000,000 single base exchanges and short indels] | http://www.mutationtaster.org |
| MutationTaster | 2010 | Naïve Bayes Classifier | SEQ,KB | SNV/Indel | P | protein | Ensembl, dbSNP, HapMap, UniProt/Swiss-Prot  | http://www.mutationtaster.org |
| PROVEAN | 2012 | Delta Alignments Scoring | SEQ | SNVs/Indels | F | protein + regulatory | NCBI NR for MSAs | http://provean.jcvi.org/index.php |
| SIFT Indel | 2012 | Decision Tree | SEQ,KB | Indel | D | protein | HGMD (2010.2) [1,292 disease indels] 1000 Genomes Project [2,602 neutral indels] | https://sift.bii.a-star.edu.sg/www/SIFT_indels2.html |

### Footnotes:
- [1] FF: Feed-Forward
- [2] GD: Grantham Difference
- [3] DMS: Deep Mutational Scanning

### Datasources:
- HGMD: Human Gene Mutation Database (http://www.hgmd.cf.ac.uk/)
- SwissVar: Swiss-Prot disease variants (https://swissvar.expasy.org/)
- dbSNP: dbSNP is world's largest database for nucleotide variations (https://www.ncbi.nlm.nih.gov/snp/)
- ESP: Exome Sequencing Project (https://science.sciencemag.org/content/337/6090/64)
- ARIC: Atherosclerosis Risk in Communities (https://academic.oup.com/aje/article/129/4/687/87924)
- dbNSFP: Annotation database for non-synonymous SNPs (https://sites.google.com/site/jpopgen/dbNSFP)
- HumVar: Curated set of human disease-causing and common (neutral) mutations (https://doi.org/10.1093/bioinformatics/btl423)
- HumDiv: Curated set of variants causing human Mendelian diseases from the UniProt Knowledgebase (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2855889/)
- PMD: Literature based database of effect variants (http://pmd.ddbj.nig.ac.jp/~pmd/)
- Swiss-Prot: Manually annotated and reviewed section of the UniProt Knowledgebase (https://www.uniprot.org/uniprot/?query=reviewed:yes)
- TrEMBL: Automatic annotation section of the UniProt Knowledgebase (https://www.uniprot.org/uniprot/?query=reviewed:no) 
- OMIM: Online Mendelian Inheritance in Man (https://omim.org/)
- VariBench: Benchmark database suite comprising variation datasets (http://structure.bmc.lu.se/VariBench/)
- PhenCode: Collaborative project to better understand the relationship between genotype and phenotype in humans (http://phencode.bx.psu.edu/)
- Idbases: Immunodeficiency mutation databases (https://www.ncbi.nlm.nih.gov/pubmed/17004234)
- CanProVar: Human Cancer Proteome Variation Database (http://canprovar.zhang-lab.org/)(https://github.com/bzhanglab/canprovar)
- CanProVar: Human Cancer Proteome Variation Database v2 (http://canprovar2.zhang-lab.org/)
- ProTherm: Thermodynamic database Database for Proteins and Mutants (https://www.iitm.ac.in/bioinfo/ProTherm/)
- dbDSM: Curated database of Deleterious Synonymous Mutations (http://bioinfo.ahu.edu.cn:8080/dbDSM/index.jsp)
- VariSNP: Benchmark database for variations from dbSNP (http://structure.bmc.lu.se/VariSNP)
- 1kGP: 1000 Genomes Project (https://www.internationalgenome.org/)
- COSMIC: Catalogue Of Somatic Mutations In Cancer (https://cancer.sanger.ac.uk/cosmic)
- ClinVar: Public archive of reports of the relationships among human variations and phenotypes (https://www.clinicalgenome.org/data-sharing/clinvar/)


### Abbreviations:
SEQ = sequence-derived, STR = structure-derived, LIT = extracted from literature, KB = extracted from knowledgebase, NET = extracted from regulatory network. F = effect on function, P = pathogenicity, S = effect on structure, D = deleteriousness, SP = effect on splicing
