# Collection of Variant Effect prediction tools
Last updated: 02/11/20 [#45]

| Tool         | Year         | Model        | Features     | Scope        | Predicts     | Impacts      | Training     | Availability |
| :---         | :---         | :---         | :---         | :---         | :---         | :---         | :---         | :---         |
| ENVISION | 2018 | Stochastic Gradient Boosting  | SEQ,STR | nsSNP | F | protein | DMS3 [21,026 variants in eight proteins] | https://envision.gs.washington.edu/shiny/envision_new/ |
| MutPred2 | 2017 | Neural Networks (FF1) | SEQ | nsSNP | P* | protein | HGMD + SwissVar + dbSNP + inter-species pairwise alignments [53,180 pathogenic and 206,946 unlabeled (putatively neutral) variants] | http://mutpred.mutdb.org |
| REVEL | 2016 | RandomForest  | ENS | nsSNP | P | protein | HGMD + ESP + ARIC + dbNSFP [6,182 disease-related variants and 123,706 rare neutral exome sequencing variants] | https://sites.google.com/site/revelgenomics |
| PANTHER-PSEP | 2016 | Phylogenetic Analysis | SEQ | nsSNP | P | protein | HumVar [12,944 disease-related variants and 8,241 neutral polymorphisms in 3,587 proteins] | http://pantherdb.org/tools/csnpScoreForm.jsp |
| SNAP2 | 2015 | Neural Networks (FF1) | SEQ | nsSNP | F | protein | PMD + Swiss-Prot + OMIM + HumVar [61,037 effect and 40,478 neutral variants in 9,744 proteins] | https://rostlab.org/services/snap2web |
| PON-P2 | 2015 | RandomForest  | SEQ,STR | nsSNP | P | protein | dbSNP + VariBench [14,086 pathogenic variants in 1,082 proteins and 14,848 neutral variants in 6,598 proteins] | http://structure.bmc.lu.se/PON-P2 |
| wKinMut-2 | 2015 | Annotation Summary | ENS,LIT,KB | nsSNP | P | protein | Swiss-Prot [865 disease-related variants in 65 proteins and 2,627 neutral variants in 447 proteins (human)] | https://github.com/Rbbt-Workflows/kin_mut2 |
| PredictSNP | 2014 | Consensus Scoring | ENS | nsSNP | F | protein | Swiss-Prot [SNPs&GO dataset; 58,057 mutations], Swiss-Prot + HGMD [MutPred dataset, 65,654 mutations], dbSNP + PhenCode + Idbases + 16 individual locus-specific databases (LSDB) [PON-P dataset, 39,670 mutations], Humsavar [36,994 neutral and disease-related mutations] | http://loschmidt.chemi.muni.cz/predictsnp |
| FATHMM-DS | 2014 | Hidden Markov Models  | SEQ | nsSNP | F | protein | HGMD [damaging variants] and Swiss-Prot/TrEMBL [putative neutral polymorphisms] | http://fathmm.biocompute.org.uk/disease.html |
| PolyPhen-2 | 2013 | Naïve Bayes Classifier | SEQ,STR | nsSNP | S,F | protein | HumDiv or HumVar [7,070 neutral and 5,322 deleterious variants] | http://genetics.bwh.harvard.edu/pph2 |
| FATHMM | 2013 | Hidden Markov Models  | SEQ | nsSNP | F | protein | HGMD [49,532 disease-causing variants] and UniProt [36,928 putatively neutral variants] | http://fathmm.biocompute.org.uk |
| VEST | 2013 | RandomForest  | SEQ | nsSNP | F | protein | HGMD Professional (v2012.2) [47,724 missense variants] and Exome Sequencing Project (ESP6500 accessed 07/2012) [45,818 likely neutral variants] | http://www.cravat.us/CRAVAT |
| FATHMM-cancer | 2013 | Hidden Markov Models  | SEQ | nsSNP | F | protein | CanProVar [12,720 positives] and UniProt [36,928 negatives] | http://fathmm.biocompute.org.uk/cancer.html |
| Meta-SNP | 2013 | RandomForest  | ENS | nsSNP | P | protein | SwissVar (SV-2009) [35,766 nsSNVs (17,883 disease variants) from 8,667 proteins] | http://snps.biofold.org/meta-snp |
| PON-P | 2012 | RandomForest  | ENS | nsSNP | P | protein | PhenCode (2009) + IDbases + 16 individual locus-specific databases (LSDB) [14,610 pathogenic variants (manual curation, SwissVar/LSDB disease annotations)] and dbSNP (Build 131) [17,393 neutral variants] | https://www.ncbi.nlm.nih.gov/pubmed/22505138 |
| KinMut | 2012 | Support Vector Machine | SEQ | nsSNP | P | protein | Swiss-Prot [865 disease and 2,627 neutral kinase mutations] | http://kinmut.bioinfo.cnio.es |
| MutationAssessor | 2011 | Functional Impact Scoring | SEQ,STR | nsSNP | F | protein | UniProt for MSA[4] [10,000 variants subject to MSA] | http://mutationassessor.org/r3 |
| MutPred | 2009 | RandomForest  | SEQ | nsSNP | S,F | protein | HGMD [34,336 disease-causing variants] and Swiss-Prot + dbSNP [23,426 putatively neutral variants] | http://mutpred1.mutdb.org/ |
| PoPMuSiC-2.0 | 2009 | Energy Function | SEQ,STR,KB | nsSNP | S | protein | ProTherm [2,648 different point mutations in 131 proteins] | http://babylone.ulb.ac.be/popmusic |
| iMutant3 | 2008 | Support Vector Machine | SEQ|STR | nsSNP | S,F | protein | ProTherm [6,398 mutations in 55 proteins] | http://gpcr2.biocomp.unibo.it/cgi/predictors/I-Mutant3.0/I-Mutant3.0.cgi |
| SNAP | 2007 | Neural Networks (FF1) | SEQ | nsSNP | F | protein | Swiss-Prot + PMD [40,830 neutral and 39,987 deleterious variants] | https://hub.docker.com/r/bromberglab/snap |
| PhD-SNP | 2006 | Support Vector Machine | SEQ | nsSNP | P | protein | HumVar [12,944 disease-related and 8241 neutral polymorphisms in 3587 proteins] | http://snps.biofold.org/phd-snp/phd-snp.html |
| Align-GVGD | 2006 | Extended GD2 Scoring | SEQ,STR | nsSNP | F | protein | p53 protein MSA[4] | http://agvgd.iarc.fr |
| FoldX | 2005 | FoldX Force Field | STR | nsSNP | S | protein | not applicable; evaluation of physicochemical constraints | http://foldxsuite.crg.eu |
| iMutant2 | 2005 | Support Vector Machine | SEQ/STR | nsSNP | S | protein | ProTherm [2087 single mutations in 65 proteins] | http://gpcr.biocomp.unibo.it/cgi/predictors/I-Mutant2.0/I-Mutant2.0.cgi |
| MAPP | 2005 | Functional Impact Scoring | SEQ | nsSNP | F | protein | not applicable; evaluation of physicochemical constraints | http://mendel.stanford.edu/supplementarydata/stone_MAPP_2005 |
| nsSNPAnalyzer | 2005 | RandomForest  | SEQ,STR | nsSNP | P | protein | ModSNP (curated dataset based on Swiss-Prot) | http://snpanalyzer.uthsc.edu |
| PolyPhen | 2002 | Rule-based Classifier | SEQ,STR | nsSNP | S,F | protein | NRDB for MSA[4] | http://genetics.bwh.harvard.edu/pph |
| SIFT | 2001 | PSSM based Probabilities | SEQ | nsSNP | D | protein | Swiss-Prot for MSA[4] | https://sift.bii.a-star.edu.sg |
| PoPMuSiC | 2000 | Energy Function | SEQ,STR,KB | nsSNP | S | protein | Set of 141 high-resolution (􏰹2.5 Å) protein X-ray structures (< 􏰺25% sequence identity) | http://babylone.ulb.ac.be/popmusic |
| IDSV | 2019 | RandomForest  | SEQ | sSNP | D | protein | dbDSM (ClinVar, PubMed, Web of Knowledge variants as disease-causing) [300 disease-related variants] and VariSNP [300 neutral varians] | http://bioinfo.ahu.edu.cn:8080/IDSV |
| ARVIN | 2018 | RandomForest  | SEQ,NET | SNV | P | regulatory | HGMD [233 disease-causing variants] and 1000 Genomes Project [2330 common variants] | https://github.com/gaolong/arvin |
| LINSIGHT | 2017 | Linear + Probabilistic Model  | SEQ | SNV | P | regulatory | High-coverage genome sequences for 54 unrelated individuals from the "69 Genome" data set from Complete Genomics | http://compgen.cshl.edu/~yihuang/LINSIGHT |
| DDIG-SN | 2017 | Support Vector Machine | SEQ | sSNP | P | protein | HGMD Professional (v2015.3) [592 disease-causing variants] and 1000 Genomes Project [10,925 putatively benign variants] in 318 genes | https://sparks-lab.org/server/ddig |
| regSNPs-splicing | 2017 | RandomForest  | SEQ | sSNP | F,SP | protein | HGMD [1,373 disease-causing synonymous SNVs (sSNVs)] and 1000 Genomes Project [7,231 neutral (in internal exon) variants, 329 (on consensus splice site) variants] | http://regsnps-splicing.ccbb.iupui.edu |
| GWAVA | 2014 | RandomForest  | SEQ,KB | SNV | P | regulatory | HGMD [1,614 disease-implicated] and 1000 Genomes Project [5,027 putatively benign variants] | https://www.sanger.ac.uk/science/tools/gwava |
| SilVA | 2013 | RandomForest  | SEQ | sSNP | P | protein | Curated, literature based dataset of rare (allele frequency <5%) synonymous variants [33 variants] and 1000 Genomes Project [746 rare synonymous variants in one individual] | http://compbio.cs.toronto.edu/silva |
| Cscape | 2017 | Multiple Kernel Learning | SEQ | SNV | P | protein + regulatory | COSMIC (pathogenic variants) and 1000 Genomes Project (control) [46,420 coding examples and 131,714 non-coding examples] | http://CScape.biocompute.org.uk |
| DANN | 2015 | Deep Neural Network | SEQ | SNV | D | protein + regulatory | 16,627,775 observed variants and 49,407,057 simulated variants  | https://cbcl.ics.uci.edu/public_data/DANN |
| FATHMM-MKL | 2015 | Multiple Kernel Learning | SEQ,KB | SNV | F | protein + regulatory | HGMD [17,362 coding; 3063 non-coding] and 1000 Genomes Project [4853 coding; 5252 non-coding] | http://fathmm.biocompute.org.uk/fathmmMKL.htm |
| CADD | 2014 | Support Vector Machine | SEQ | SNV/Indel | D | protein + regulatory | observed (14,893,290 SNVs, 627,071 insertions and 1,107,414 deletions) and simulated variants | https://cadd.gs.washington.edu |
| MutationTaster2 | 2014 | Naïve Bayes Classifier | SEQ,KB | SNV/Indel | P | protein + regulatory | HGMD Professional + ClinVar [>100,000 disease-associated mutations] and 1000 Genomes Project [>6,000,000 single base exchanges and short indels] | http://www.mutationtaster.org |
| MutationTaster | 2010 | Naïve Bayes Classifier | SEQ,KB | SNV/Indel | P | protein | HGMD + OMIM + literature [57,100 disease-causing variants] and dbSNP + HapMap [523,425 putatively neutral variants] | http://www.mutationtaster.org |
| PROVEAN | 2012 | Delta Alignments Scoring | SEQ | SNVs/Indels | F | protein + regulatory | NCBInr for MSA[4] | http://provean.jcvi.org/index.php |
| SIFT Indel | 2012 | Decision Tree | SEQ,KB | Indel | D | protein | HGMD (2010.2) [1,292 disease indels] 1000 Genomes Project [2,602 neutral indels] | https://sift.bii.a-star.edu.sg/www/SIFT_indels2.html |

### Footnotes:
- [1] FF: Feed-Forward
- [2] GD: Grantham Difference
- [3] DMS: Deep Mutational Scanning
- [4] MSA: Multipe Sequence Alignment

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
- NRDB: Non-redundant composite of the following sources: PDB sequences, SWISS-PROT, SWISS-PROTupdate, PIR, GenPept and GenPeptupdate (http://130.88.97.239/bioactivity/owlinfo.html)
- NCBInr: Non-redundant protein sequences from GenPept, Swissprot, PIR, PDF, PDB, and NCBI RefSeq (ftp://ftp.ncbi.nlm.nih.gov/blast/db/README)
- Humsavar: Human polymorphisms and disease mutations: index (https://www.uniprot.org/docs/humsavar)
- ModSNP: Database for Sequence and Structure Information on Human Protein Variants (https://www.ncbi.nlm.nih.gov/pubmed/15108278)


### Abbreviations:
SEQ = sequence-derived, STR = structure-derived, LIT = extracted from literature, KB = extracted from knowledgebase, NET = extracted from regulatory network. F = effect on function, P = pathogenicity, S = effect on structure, D = deleteriousness, SP = effect on splicing
