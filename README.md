This is a collection of bioinformatics tools I have sourced from recent literature, organized by topic. I have not used most of these tools.

Tools with dead links (404, DNS failure, etc.) have been moved to [DEAD.md](DEAD.md). Tools hosted in repositories with no activity in 2+ years have been moved to [UNMAINTAINED.md](UNMAINTAINED.md). If a tool has been moved to one of these files by mistake, please [open an issue](https://github.com/jdidion/biotools/issues/new).

**Table of Contents**

<!-- toc -->

- [Clinical data](#clinical-data)
  - [EHR](#ehr)
- [Data Management](#data-management)
- [Data Sets](#data-sets)
- [Discovery](#discovery)
- [AI](#ai)
- [Genomics](#genomics)
  - [General Information](#general-information)
  - [Algorithms](#algorithms)
    - [Pairwise Alignment](#pairwise-alignment)
  - [Assay Design](#assay-design)
  - [Candidate Prioritization](#candidate-prioritization)
  - [Databases](#databases)
  - [Data Formats](#data-formats)
  - [Functional Enrichment/Ontology](#functional-enrichmentontology)
  - [GWAS/QTL](#gwasqtl)
  - [Microarray](#microarray)
  - [Motif/TFBS](#motiftfbs)
  - [Network Analysis](#network-analysis)
  - [Population genetics](#population-genetics)
  - [Prediction](#prediction)
  - [Security](#security)
  - [Sequencing Protocols](#sequencing-protocols)
  - [Simulation](#simulation)
  - [Variant annotation](#variant-annotation)
  - [Sequence Analysis](#sequence-analysis)
    - [General-purpose](#general-purpose)
    - [QC](#qc)
    - [Peak-based](#peak-based)
      - [ChIP-seq](#chip-seq)
      - [Chromatin accessibility](#chromatin-accessibility)
    - [Chromatin Interactions](#chromatin-interactions)
    - [DNA](#dna)
    - [Footprinting](#footprinting)
    - [Gene editing](#gene-editing)
    - [Metagenomics](#metagenomics)
    - [Methylation](#methylation)
    - [MNase-seq](#mnase-seq)
    - [Nanopore/PacBio](#nanoporepacbio)
    - [RNA](#rna)
    - [Single-cell](#single-cell)
    - [Somatic](#somatic)
    - [Multi-omics/integration](#multi-omicsintegration)
- [Proteomics](#proteomics)
  - [MSA](#msa)
- [General Programming Resources](#general-programming-resources)
  - [C/C++](#cc)
  - [R](#r)
    - [Find packages](#find-packages)
    - [Database](#database)
    - [Data Cleaning](#data-cleaning)
    - [Reporting](#reporting)
    - [Misc](#misc)
  - [Python](#python)
  - [HPC](#hpc)
  - [Command Line (OSX/Linux)](#command-line-osxlinux)
  - [Databases](#databases-1)
  - [Reproducibility/Containerization](#reproducibilitycontainerization)
    - [Building Pipelines](#building-pipelines)
    - [Software Distribution](#software-distribution)
  - [Other](#other)
- [Statistics/Machine Learning](#statisticsmachine-learning)
  - [Methods/algorithms](#methodsalgorithms)
  - [Python](#python-1)
  - [Deep Learning](#deep-learning)
  - [Web APIs](#web-apis)
  - [Data Sets](#data-sets-1)
  - [Text classification](#text-classification)
- [Visualization](#visualization)
  - [Genome browsers](#genome-browsers)
  - [Networks](#networks)
  - [Phylogenetic trees](#phylogenetic-trees)
  - [R](#r-1)
    - [ggplot2](#ggplot2)
    - [Plot Types](#plot-types)
    - [Data Types](#data-types)
    - [Interactive](#interactive)
  - [Python](#python-2)
  - [Javascript](#javascript)
  - [Examples](#examples)
- [Publication/Archiving](#publicationarchiving)
  - [Code/Data sharing](#codedata-sharing)
  - [Journals](#journals)
  - [Writing](#writing)
  - [Posters](#posters)
  - [Slides](#slides)
  - [CV](#cv)
- [Promising methods without software implementation](#promising-methods-without-software-implementation)

<!-- tocstop -->

# Clinical data

## EHR

- BlueButton related tools https://github.com/amida-tech

# Data Management

- Data version control https://github.com/vsbuffalo/scidataflow

# Data Sets

- 1000 Genomes (RNA-Seq, ChIP-Seq):
  - http://archive.gersteinlab.org/docs/2015/06.04/1kg_fun_studies.htm
  - In S3: http://s3.amazonaws.com/1000genomes
  - 2500 Phase 3 samples at 30X: https://www.ebi.ac.uk/ena/data/view/PRJEB31736
- Reference panel of ~250 Dutch families http://biorxiv.org/content/early/2016/01/18/036897
- FANTOM consortium has CAGE (5' single molecule RNA counting) data from ~1000 human cell/tissue/cell-line samples from ~300 different cell/tissue types
- Full text of all PMC papers from 2008-present: ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/
- Phased genome sequences
  - &gt;100 fully-phased: http://gigascience.biomedcentral.com/articles/10.1186/s13742-016-0148-z
  - Statistically phased: http://www.haplotype-reference-consortium.org/
  - Phased variants: http://genome.cshlp.org/content/early/2016/11/25/gr.210500.116.full.pdf+html
- Clinical trials: http://vivli.org/
- Targets for drug discovery: https://www.targetvalidation.org/
- Species images: http://phylopic.org/image/browse/
- Public RNA-Seq data: 
  - https://jhubiostatistics.shinyapps.io/recount/ 
  - with phenotype predictions https://bioconductor.org/packages/release/bioc/html/recount.html
  - human and mouse: http://amp.pharm.mssm.edu/archs4
- Whole genomes of 150 Danish individuals: http://www.nature.com/nature/journal/v548/n7665/full/nature23264.html
- Simons Genome Diversity Project: https://www.simonsfoundation.org/simons-genome-diversity-project/
- BloodPAC: https://doi.bloodpac.org/BLOODPAC.0001/
- Reference haplotypes for STR phasing: https://www.biorxiv.org/content/early/2018/03/06/277673
- Database of Curated Mutations: http://www.docm.info/
- Normalized, batch-corrected RNA-seq for both normal (GTEx) and cancer (TCGA) data: https://www.nature.com/articles/sdata201861
- Orthogonally validated TP and FP CNVs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5473400/
- Microarray transcriptome profiles prepared for ML applications: https://www.biorxiv.org/content/biorxiv/early/2018/06/22/353698.full.pdf
- Paired mRNA and protein abundance in 26 human tissues: https://www.biorxiv.org/content/biorxiv/early/2018/06/27/357137.full.pdf
- Cancer cell lines: correlation with primary tumor data (analysis for use as models): http://comphealth.ucsf.edu/TCGA110/
- gnomAD: SNV, small InDels, MNV, SV, and LoF metrics from large numbers of exomes and whole-genomes: https://gnomad.broadinstitute.org/downloads
  - Blog post on lessons and tools for filtering of gnomAD annotations https://www.whiffinlab.org/blog/six-lessons-for-variant-interpretation
- 1000G Structural variation: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/working/20181025_EEE_SV-Pop_1/AltReference_EEE_SV-Pop_1/
- PGP UK: https://www.personalgenomes.org.uk/data/
- Platinum genomes phased variant calls: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs001224.v1.p1
- ICGC: https://dcc.icgc.org/
- Cancer Cell Lines
  - Encylopedia: https://gdc.cancer.gov/
  - Deep sequencing of 3 cancer cell lines (2 colorectal, one breast cancer) on 2 sequencing platforms: https://www.biorxiv.org/content/10.1101/623702v1
  - Deep sequencing of another breast cancer cell line (SEQC2 consoritum): https://www.biorxiv.org/content/biorxiv/early/2019/05/13/625624.full.pdf
- Variation benchmark datasets: http://structure.bmc.lu.se/VariBench/
- RNA Atlas https://hgserver1.amc.nl/cgi-bin/r2/main.cgi?&dscope=RNAATLAS&option=about_dscope
- COLO829 sequenced on multiple short- and long-read platforms https://www.biorxiv.org/content/10.1101/2020.10.15.340497v1?rss=1
- HGDP
  - SV calls ftp://ngs.sanger.ac.uk/scratch/project/team19/HGDP_SV/
- MANE
  - Predicted protein structures: isoform.io
- OMICSPRED: an atlas of genetic scores https://omicspred.org
- Exclusion regions: https://github.com/dozmorovlab/excluderanges/blob/master/man/figures/Table_S1.csv
- Benchmark data set based on Chinese quartet https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03109-2

# Discovery

When looking for a bioinformatics tool for a specific application:

- OmicTools: http://omictools.com/
- GitXiv: http://www.gitxiv.com/?cat%5B0%5D=bioinformatics
- Bio.Tools: https://bio.tools/
- Biosharing: https://biosharing.org

# AI

- Interface to Chat-GPT for bioinformatics applications https://github.com/JoshuaChou2018/AutoBA

# Genomics

## General Information

- Migrating to GRCh38: https://software.broadinstitute.org/gatk/blog?id=8180
- Command-line interface to download and manage reference genomes: https://github.com/simonvh/genomepy
- Useful cheat-sheets: https://github.com/NBISweden/GAAS/tree/master/annotation/CheatSheet

## Algorithms

- Suffix arrays: http://almob.biomedcentral.com/articles/10.1186/s13015-016-0068-6
- gapped k-mer SVM: classifiers for DNA and protein sequences http://www.beerlab.org/gkmsvm/
  - Version for large-scale data: https://github.com/Dongwon-Lee/lsgkm/
- Wrapper around AMQs to increase speed and decrease false-positive rate https://github.com/lrobidou/findere


### Pairwise Alignment

#### General

* Alignment-free tools for pairwise sequence comparison:
    * http://www.combio.pl/alfree/tools/
    * http://www.combio.pl/alfree
    * http://bioinformatics.org.au/tools/decaf+py/
    * Idea: Train a NN to approximate the pairwise alignment distance https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty887/5140215?redirectedFrom=fulltext
    * Optimized embedding function: https://www.biorxiv.org/content/10.1101/2020.05.24.113852v2
* EDLIB: edit distance-based alignment: https://github.com/Martinsos/edlib
* The SeqAn C++ library implements several alignment algorithms:
http://www.sciencedirect.com/science/article/pii/S0168165617315420
* Biotite implementation of optimal alignment: https://www.biotite-python.org/apidoc/biotite.sequence.align.align_optimal.html#biotite.sequence.align.align_optimal
* Heng Li review of dynamic programming algorithms: https://figshare.com/articles/Notes_on_pairwise_alignment_with_dynamic_programming/5223973
* Parasail: SIMD-accelerated pairwise alignment algorithms https://github.com/jeffdaily/parasail
* Greedy algorithm: http://pipmaker.bx.psu.edu/dist/greedy.pdf
* [WFA](https://github.com/smarco/WFA): `O(ns)` algorithm (i.e. very fast for similar sequences)
    * WFA-inspired implementation of LV https://github.com/lh3/lv89
* Block aligner (WIP): https://github.com/Daniel-Liu-c0deb0t/block-aligner
* GPU-accelerated: https://github.com/fkallen/Accelign
#### Local

* Nucl2vec: https://github.com/prakharg24/Nucl2vec
    * local alignment only - could it be adapted to semi-global?
    * Can re-implement in SpaCy? https://spacy.io/usage/vectors-similarity
* Accelerated version of Smith-Waterman algorithm implemented in the Genomics Kernel Library: https://github.com/Intel-HLS/GKL
#### Semi-global

Semi-global alignment is global alignment with an optional gap at the end(s).

* SeqPurge algorithm: insert match algorithm that performs thresholded exhaustive
   comparison to minimize probability of incorrect alignment. Relies on the fact that
   overlapping reads share alleles and indels (i.e. no gaps are required) (in C++).
   https://github.com/imgag/ngs-bits/tree/master/src/SeqPurge.
   * Speed up sequence comparison:
      * Between adapters and overhangs: http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4080745/pdf/btu177.pdf
      * Between reads: http://bioinformatics.oxfordjournals.org/content/30/14/2000.abstract
* Quality-aware overlap alignment (in Haskell): https://hackage.haskell.org/package/bio-0.5.3/docs/Bio-Alignment-QAlign.html
* FOGSAA, modified for semi-global alignment.
   * http://www.nature.com/articles/srep01746
   * http://www.isical.ac.in/~bioinfo_miu/FOGSAA.7z
* MEM-align: https://sites.google.com/site/memalignv1/
* Atria: byte-based matching algorithm https://github.com/cihga39871/Atria
#### Global

* A*PW https://github.com/RagnarGrootKoerkamp/astar-pairwise-aligner
* https://github.com/bxskdh/TSTA
#### Barcode matching

* Pheniqs demultiplxer
   * Phred-adjusted ML for error probability: https://github.com/biosails/pheniqs/blob/master/docs/2.0/glossary.md
#### Edit Distance

* https://github.com/maxdoblas/QuickEd
* Paper on relationship between alignment score and probability models: https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btz576/5536873

## Assay Design

- MBRAnator: design of MPRA libraries https://www.genomegeek.com/

## Candidate Prioritization

- GenRank: http://bioconductor.org/packages/GenRank/
- VarFish: https://github.com/bihealth/varfish-server

## Databases

- RNA Meta Analysis (web app)
  - Perform meta-analyses of GEO microarray and RNA-Seq studies
  - Correlate resulting signatures to CMAP02/L1000 perts and 26,000+ GEO studies
  - https://www.rnama.com
- NAR catalog of databases, by subject: http://www.oxfordjournals.org/our_journals/nar/database/c/
- TransVar http://bioinformatics.mdanderson.org/transvarweb/
- dbMAE: mono-allelic expression https://mae.hms.harvard.edu/
- Disease-gene associations http://www.disgenet.org/web/DisGeNET/menu/rdf
- Human tissue-specific enhancers: http://www.enhanceratlas.org/
- Searh HLI's genome data: hli-opensearch.com
- Chromatin-state annotations + per-base functionality scores for 164 cell types: http://noble.gs.washington.edu/proj/encyclopedia/
- Feature-based classification of human transcription factors: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1349-2
- Database of epifactors (epigenetic factors): http://epifactors.autosome.ru/
- SRA metadata: http://deweylab.biostat.wisc.edu/metasra/
- CTD: drug interactions and toxicity http://ctdbase.org/
- FANTOM lncRNA catalog: http://fantom.gsc.riken.jp/cat/
- Integrated database of public ChIP-Seq datasets: http://chip-atlas.org/
- Alternative splicing: http://vastdb.crg.eu/wiki/Main_Page
- Database of cis-regulatory elements (enhancers): http://www.kostkalab.net/software.html
- Structural variation: https://www.ncbi.nlm.nih.gov/dbvar
- HIVE BioMuta: https://hive.biochemistry.gwu.edu/biomuta/readme
- CHESS: new gene annotation database
  - Contains most genes from RefSeq and Gencode, plus additional genes discovered from GTEx transcripts
  - http://ccb.jhu.edu/chess/
- Metabase: aggregates various gene annotation sources: http://metascape.org
- CIViC: https://civicdb.org/home
- Allele frequencies https://www.ncbi.nlm.nih.gov/snp/docs/gsr/alfa/
- Harmonized cancer variants https://search.cancervariants.org

## Data Formats

- Searching overlaps of large collections of interval databases: https://github.com/ryanlayer/giggle/
- FASTA/FASTQ
  - Read filtering and extraction: https://github.com/ad3002/Cookiecutter
  - fqtools: for working with FASTQ files https://github.com/alastair-droop/fqtools
  - Integrating results from multiple genome binning programs: https://github.com/songweizhi/Binning_refiner
  - Compression
    - fqzcomp: https://sourceforge.net/projects/fqzcomp/
    - Clumpify: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/
    - NAF: https://github.com/KirillKryukov/naf
    - SPRING: https://github.com/shubhamchandak94/SPRING
    - Using de Bruijn graphs: https://github.com/rongjiewang/BdBG
    - https://github.com/divonlan/genozip
- SAM/BAM/CRAM
  - Read filtering and profiling: https://github.com/jwalabroad/VariantBam
  - SQL-like querying https://github.com/maragkakislab/samql/
  - encoding Nanopore data: https://github.com/novoalab/modPhred
- VCF
  - Annotation: https://github.com/brentp/vcfanno
  - BGT: https://github.com/lh3/bgt
  - GQT
  - BCFtools: includes new tool to identify RoH http://samtools.github.io/bcftools/
  - Work with VCF in R: https://github.com/knausb/vcfR
  - More VCF tools: http://vcf-kit.readthedocs.io/
  - Validation: https://github.com/EBIvariation/vcf-validator
  - Index and subset VCF files: https://github.com/tomkurowski/tersect
  - Variant of BCF format encoding genotypes using sparse vectors https://github.com/statgen/savvy
  - Filtering
    - https://github.com/vembrane/vembrane
- BED/GFF
  - GFFTools: https://github.com/ihh/gfftools
  - Gff3Sort: https://github.com/billzt/gff3sort
  - https://github.com/alejandrogzi/gtfsort
  - Combine p-values https://github.com/brentp/combined-pvalues
  - Assessing interval overlap of multiple genomic features: https://github.com/andrew-leith/GINOM
  - Python toolkit: https://github.com/dputhier/pygtftk
  - https://github.com/lh3/bedtk
  - perturbation testing of region sets: http://bedshift.databio.org/en/latest/
- BigBED and BigWig
  - https://github.com/jackh726/bigtools
- Improved command line viewing of formatted data
  - PrettySam: http://lindenb.github.io/jvarkit/PrettySam.html
  - BioSyntax: https://biosyntax.org/install
  - DNAcol: https://github.com/koelling/dnacol
- Convert between many different formats, using a SQL database as intermediate representation: https://watson.hgen.pitt.edu/mega2/mega2r/

## Functional Enrichment/Ontology

- General:
  - R interface to DAVID: http://www.bioconductor.org/packages/release/bioc/html/RDAVIDWebService.html
  - GSEA
    - Cautions about the GSEA null model: http://bioinformatics.oxfordjournals.org/content/early/2017/01/02/bioinformatics.btw803.short
  - GiANT: uncertainty in GSEA https://cran.r-project.org/web/packages/GiANT/index.html
  - Ensemble gene set enrichment analysis: http://bioconductor.org/packages/release/bioc/html/EGSEA.html
  - Fast GSA: https://github.com/billyhw/GSALightning
  - Fast GSEA: https://github.com/ctlab/fgsea
  - Multi-dimensional GSEA: http://bioconductor.org/packages/release/bioc/html/mdgsa.html
  - GSEA with external information https://cran.r-project.org/web/packages/netgsa/index.html
  - QTest: http://statgen.snu.ac.kr/software/QTest/
  - Gene set analysis with specific alternative hypothesis: https://bioconductor.org/packages/release/bioc/html/GSAR.html
  - clusterProfiler: <https://guangchuangyu.github.io/clusterProfiler/>
  - DOSE: disease ontology <https://guangchuangyu.github.io/dose/>
  - ReactomePA: <https://guangchuangyu.github.io/reactomepa/>
  - SetRank: https://cran.r-project.org/web/packages/SetRank/index.html
  - Identify and rank significance of overlaps: https://github.com/ryanlayer/giggle
  - Generate feature vectors from ontology annotations: https://github.com/bio-ontology-research-group/onto2vec
- Gene sets
  - EGSEA: http://bioconductor.org/packages/release/data/experiment/html/EGSEAdata.html
  - ENCODE gene set hub: http://genomespot.blogspot.com.au/2017/02/introducing-encode-gene-set-hub.html
- GO term:
  - GOrilla: http://cbl-gorilla.cs.technion.ac.il/
  - LRPath: http://lrpath.ncibi.org/
  - Reduce GO term lists: 
    - Revigo: http://revigo.irb.hr/
  - clusterProfiler: http://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/
  - GO Express: https://www.bioconductor.org/packages/release/bioc/html/GOexpress.html
  - GO Extender: https://www.msu.edu/~jinchen/GOExtender/
  - DAL: http://iwera.ir/~ahmad/dal/
  - Negative GO enrichment: https://sites.google.com/site/guoxian85/neggoa
- Regional
  - LOLA: http://lola.computational-epigenetics.org
  - AnnotatR: https://github.com/rcavalcante/annotatr/
  - Compare epigenetic features in multiple samples: http://epigenome.wustl.edu/EpiCompare1/
- Trait
  - traseR: Trait-associated SNP enrichment https://www.bioconductor.org/packages/release/bioc/html/traseR.html
  - Identifying genetic heterogeneity within phenotypically defined subgroups: https://github.com/jamesliley/subtest
- Multi-omics
  - Single-sample GSA across data sets: https://www.bioconductor.org/packages/3.3/bioc/html/mogsa.html
- Network-based
  - NEAT: https://cran.r-project.org/web/packages/neat/index.html
  - EGAD: https://github.com/sarbal/EGAD
- Mutation intolerance:
  - LIMBR: https://redmine.igm.cumc.columbia.edu/projects/atav/wiki/LIMBR
  - Sub-regional: https://www.biorxiv.org/content/early/2018/09/30/431536
  
## GWAS/QTL

- Association testing
  - QTLtools: Pipeline for molecular QTL analysis https://qtltools.github.io/qtltools/
    - FastQTL: http://fastqtl.sourceforge.net/
  - RASQUAL: allele-specific QTL using phased SNPs - https://github.com/dg13/rasqual
  - Relatedness, PCA: http://zhengxwen.github.io/SNPRelate/
  - Multi-SNP, multi-trait regression https://github.com/ashlee1031/BERRRI
  - regioneR: permutation testing for association between genomic region and phenotype http://bioconductor.org/packages/release/bioc/html/regioneR.html
  - MTG2: https://sites.google.com/site/honglee0707/mtg2
  - Use local gene networks to improve trans-eQTL detection: https://github.com/PMBio/GNetLMM
  - Fast correlation testing: https://github.com/gabraham/flashpca/tree/master/flashpcaR
  - Gene and pathway association testing from summary statistics: https://cran.r-project.org/web/packages/aSPU/
  - Account for LD and functional information: https://github.com/yjingj/SFBA
  - Correcting for prediction uncertainty in TWAS: http://biorxiv.org/content/early/2017/02/14/108316
  - Using random forests https://github.com/0asa/TTree-source
  - Using k-mers: https://github.com/atifrahman/HAWK
  - BOLT-LMM: https://data.broadinstitute.org/alkesgroup/BOLT-LMM/
  - tensorQTL: https://github.com/broadinstitute/tensorqtl
  - Wavelet https://github.com/william-denault/WaveletScreening
  - Integrates total and allele-specific read counts https://github.com/hakyimlab/mixqtl
  - Apex https://corbinq.github.io/apex/
  - gVAMP https://github.com/medical-genomics-group/gVAMP
- Variance eQTL
  - veQTL mapper: https://funpopgen.github.io/veqtl-mapper/
- Multiple test correction
  - eigenMT: Efficient multiple-test correction http://montgomerylab.stanford.edu/resources/eigenMT/eigenMT.html
  - Fast multiple-test correction for LMMs: http://genetics.cs.ucla.edu/multiTrans/
  - Hierarchical eQTL MTC: http://bioinformatics.org/treeqtl/
  - Controlling bias in EWAS/TWAS using null distribution: http://bioconductor.org/packages/bacon/
- Prioritization
  - GenoWAP: Prioritization of GWAS signals using functional information http://genocanyon.med.yale.edu/GenoWAP
  - HitWalker2: https://github.com/biodev/HitWalker2
  - CRstar: https://nijingchao.github.io/CRstar/
  - Pines: http://genetics.bwh.harvard.edu/pines/
- Fine-mapping
  - Using summary statistics http://www.christianbenner.com
  - uHDset: https://sites.google.com/site/shaolongscode/home/uhdset
  - PAINTOR: fine mapping, prioritization - https://github.com/gkichaev/PAINTOR_FineMapping/
  - Genotype synthesis: https://sourceforge.net/projects/getsynth/
  - Prediction of causal variants from epigenomic annotations: https://github.mit.edu/liyue/rivieraBeta
  - DAP: Bayesian framework for QTL analysis and fine-mapping https://github.com/xqwen/dap
  - BayesFM: https://sourceforge.net/projects/bayesfm-mcmc-v1-0/
  - Determining causal genes using TADs http://biorxiv.org/content/early/2016/11/15/087718
  - Network analysis: https://github.com/YuanlongLiu/SigMod
  - Haplotype-based: http://apps.biocompute.org.uk/haprap/
  - https://github.com/liliulab/treemap
- Ancestry
  - SNPweights: https://www.hsph.harvard.edu/alkes-price/software/
  - Visualize population structure: https://github.com/kr-colab/popvae
  - https://github.com/sriramlab/SCOPE
- LD
  - LD score calculation and regression https://github.com/bulik/ldsc
  - Haplotype block detection: https://github.com/sunnyeesl/BigLD
- Find directional effect of a signed functional annotation using summary statistics: https://github.com/yakirr/sldp
- Imputation of missing phenotype information http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3513.html
- Epistasis
  - General linear-mixed model library; also includes mixed-RF method for detecting epistasis with population structure correction: https://github.com/PMBio/limix
  - GPU-accelerated detection of epistasis using Bayesian neural networks: https://github.com/beamandrew/BNN
  - MEPID: marginal epistasis test http://www.xzlab.org/software.html
- Other
  - Browser for geographical distribution of genetic variants: http://popgen.uchicago.edu/ggv/
  - Integration of GWAS with molecular QTL: https://github.com/xqwen/integrative
  - GxE interactions: https://github.com/davidaknowles/eagle
  - [DENTIST](https://github.com/Yves-CHEN/DENTIST/): QC for summary statistics 

## Microarray

- Correct for batch effects between training data and external datasets: https://cran.r-project.org/web/packages/bapred/index.html
- Impute from Affy expression arrays: http://simtk.org/home/affyimpute
- SNP
  - Call haplotypes https://cran.r-project.org/web/packages/GHap/index.html
  - Neural network for breakpoint detection for CNV calling: https://www.biorxiv.org/content/biorxiv/early/2018/06/24/354423.full.pdf
- Methylation
  - Minfi: R package for working with 450k methylation arrays
  - D3M: two-sample test of differential methylation from distribution-valued data https://cran.r-project.org/web/packages/D3M/D3M.pdf
  - Network-based approach to discovering epigenetic "modules" that can be associated with gene expression: http://bioinformatics.oxfordjournals.org/content/30/16/2360.long
  - Filtering probes using technical replicates https://cran.r-project.org/web/packages/CpGFilter/index.html
  - Imputation of genome-wide methylation: http://wanglab.ucsd.edu/star/LR450K/
  - Tutorial for analysis using bioconductor packages: http://biorxiv.org/content/biorxiv/early/2016/05/25/055087.full.pdf
  - Normalization
    - Probe design bias correction: https://www.bioconductor.org/packages/release/bioc/html/ENmix.html
  - DMR calling
    - SeqLM: https://github.com/raivokolde/seqlm
  - PyMAP: http://aminmahpour.github.io/PyMAP/
  - Interactive exploration: http://bioconductor.org/packages/release/bioc/html/shinyMethyl.html
  - eFORGE: identify cell type-specific signals in differentially methylated positions (mostly important for blood-based EWAS) http://eforge.cs.ucl.ac.uk/
  - Model for EWAS using probe signal intensities: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1347-4
  - Reference-based tissue deconvolution (three different algorithms): https://github.com/sjczheng/EpiDISH
  - Glint pipeline (qc, EWAS, population structure): http://glint-epigenetics.readthedocs.io/en/latest/
    - Bayesian extension of Refactor cell type heterogeneity correction that incorporates experimentally determined cell counts: https://github.com/cozygene/BayesCCE
  - Meffil: https://github.com/perishky/meffil/ 

## Motif/TFBS

- Fast, SNP-aware PWM matching https://www.cs.helsinki.fi/group/pssmfind/
- R package for TFBS analysis: http://bioconductor.org/packages/release/bioc/html/TFBSTools.html
- FIMO and MCAST perform best of TFBS predictors: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1298-9
- Intra-motif dependencies
  - Identifying http://www.jstacs.de/index.php/InMoDe
  - Visualizing http://bioinformaticstools.mayo.edu/circularlogo/index.html
- Dinucleotide weight tensors encode dependencies between positions in motifs: http://dwt.unibas.ch/
- Convert kernels learned by CNN to PWMs: ftp://ftp.cbi.pku.edu.cn/pub/software/CBI/k2p
- On graphs: https://github.com/pinellolab/GRAFIMO

## Network Analysis

- Identification of module with biclustering: http://web.ist.utl.pt/rmch/bicnet/temporary/index.jsp
- Functional analysis https://bioconductor.org/packages/release/bioc/html/EGAD.html
- GAGE: http://bioconductor.org/packages/release/bioc/html/gage.html
- Gene network reconstruction
  - For a set of TFs: https://sourceforge.net/projects/aracne-ap/
  - From PPI or motif sharing: https://github.com/davidvi/pypanda (is also an integrative method that can incorporate multiple sources of information)
  - Perturb-Net: Learn gene networks underlying clinical phenotypes
    - https://github.com/sssykim/Perturb-Net
    - https://www.biorxiv.org/content/early/2018/09/10/412817
    - Applied to TF binding: https://www.biorxiv.org/content/early/2018/09/10/412841
- BANFF: https://cran.r-project.org/web/packages/BANFF/index.html
- SAFE: https://bitbucket.org/abarysh/safe
- Merlin-P: https://bitbucket.org/roygroup/merlin-p
- OSS alternative to Inginuity pathway analysis: https://www.bioconductor.org/packages/release/bioc/html/QuaternaryProd.html
- Tissue-specific: https://cran.r-project.org/web/packages/GRAPE/index.html
- CLR with B-Spline for mutual information-based inference of transcriptional regulatory networks: https://bitbucket.org/Jonathan-Ish-Horowicz/fastgenemi/
- Functional module identification and enrichement through Gene Ontology terms https://www.nature.com/articles/s41598-018-23672-0

## Population genetics

- Management of allele-frequency data: https://grenaud.github.io/glactools/

## Prediction

- Gene expression
  - Predicting gene expression from methylation: http://arxiv.org/abs/1603.08386
- QTL
  - Imputation of summary statistics in multi-ethnic cohorts: http://dleelab.github.io/jepegmix/
  - Causal variant identification: http://genetics.cs.ucla.edu/caviar/
- Genetic risk
  - AnnoPred: https://github.com/yiminghu/AnnoPred
- Pathogenicity
  - Disease-specific https://github.com/samesense/pathopredictor
  - SVs https://github.com/gersteinlab/SVFX
  - https://data.broadinstitute.org/alkesgroup/LDSCORE/Kim_annotboost/
  - Frequency conservation score http://bioinfo.cnic.es/FCS
  - Rare variants: http://www.liulab.science/metarnn.html
  - SVs
    - https://github.com/gersteinlab/SVFX
    - https://github.com/andrewSharo/StrVCTVRE
  - Splice variants https://github.com/gagneurlab/MMSplice_MTSplice
  - Missense variants https://github.com/jiaying2508/LYRUS
  - Long reads: https://github.com/TheJacksonLaboratory/SvAnna
  - https://evemodel.org/
- Causal variant
  - Ensemble method: https://github.com/gifford-lab/EnsembleExpr/
  - eCAVIAR: probability that a variant is causal for both QTL and eQTL http://genetics.cs.ucla.edu/caviar/
  - qCat: https://github.com/dleelab/qcat
  - Disease-associated risk variants: https://sites.google.com/site/emorydivan/
  - Predicting gene targets from GWAS summary statistics https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4979185/
  - Orion: https://github.com/igm-team/orion-public
- Chromatin States
  - GenoSTAN: http://bioconductor.org/packages/release/bioc/html/STAN.html
  - R package for predicting chromatin states from histone marks across conditions https://github.com/ataudt/chromstaR
  - Rule-based: http://www.statehub.org/
  - Hierarchical HMM: https://github.com/gcyuan/diHMM
  - Basenji: https://github.com/calico/basenji
  - CMint: https://bitbucket.org/roygroup/cmint
- Enhancers
  - Prediction of enhancer strength from sequenced http://bioinformatics.hitsz.edu.cn/iEnhancer-2L
  - Prediction of core cell type-specific TFs from super enhancers https://bitbucket.org/young_computation/crcmapper
  - Prediction of superenhancers https://github.com/asntech/improse
  - Deep learning-based:
    - BiRen: https://github.com/wenjiegroup/BiRen
    - DeepCAPE: https://www.biorxiv.org/content/early/2018/08/24/398115
    - https://www.biorxiv.org/content/early/2018/02/14/264200
- Coding mutations
  - Predict mutation effects from sequence covariation: https://marks.hms.harvard.edu/evmutation/
  - Predicting loss-of-function from expression: https://www.nature.com/articles/s41467-017-00443-5
  - Prediction of missense variants: https://www.biorxiv.org/content/early/2018/02/02/259390
  - Disease-specific functional prediction https://sites.google.com/site/emorydivan/
  - Impact of coding SNPs: http://pantherdb.org/tools/csnpScoreForm.jsp
  - Predict disease risk from GWAS summary statistics: https://github.com/yiminghu/AnnoPred
- Regulatory variants/TF binding
  - LedPred: prediction of regulatory sequences from ChIP-seq https://github.com/aitgon/LedPred
  - GERV: prediction of regulatory variants that affect TF inding http://gerv.csail.mit.edu/
  - Score variant deleteriousness: http://cadd.gs.washington.edu/
  - BASSET: prediction of sequence activity https://github.com/davek44/Basset
  - DanQ: hybrid convolutional and recurrent neural network model for predicting the function of DNA de novo from sequence http://github.com/uci-cbcl/DanQ
  - DeepSequence: Generative model https://github.com/debbiemarkslab/DeepSequence
    - Can also be used for sequence clustering
  - [LINSIGHT](https://github.com/CshlSiepelLab/LINSIGHT)
  - Protein binding affinity: https://bitbucket.org/wenxiu/sequence-shape.git
  - Change in local frustration index: https://github.com/gersteinlab/frustration
  - TFImpute: multi-task learning from ChIP-seq data across factors and tissues to impute TF binding for an unassayed tissue/factor combination: https://bitbucket.org/feeldead/tfimpute
  - PPI https://www.ncbi.nlm.nih.gov/research/mutabind/index.fcgi/
  - Multiple instance learning of TF binding: http://www.cs.utsa.edu/~jruan/MIL/
  - Cell type-specific: https://github.com/uci-cbcl/FactorNet
  - Predict TF binding from ATAC-Seq using deep neural network:  https://github.com/hiranumn/deepatac
  - CentiSNPs: http://genome.grid.wayne.edu/centisnps/
  - Benchmark: https://github.com/Oncostat/BenchmarkNCVTools
  - PathoPredictor: https://github.com/samesense/pathopredictor
  - Cell-type agnostic regulatory activity prediction from ENCODE: http://screen.encodeproject.org/
  - Segway functional scores: https://noble.gs.washington.edu/proj/encyclopedia/
  - Predict effect of variants on epigenetic factors (chromatin accessibility, histone marks, etc): http://deepfigv.mssm.edu/downloads.html
  - HaploReg: http://www.broadinstitute.org/mammals/haploreg/haploreg.php
  - Consensus approaches:
    - PredictSNP2: https://loschmidt.chemi.muni.cz/predictsnp2/
    - PRVCS: http://jjwanglab.org/PRVCS/
  - Queryor: http://queryor.cribi.unipd.it/cgi-bin/queryor/mainpage.pl
  - Using epigenomic data https://github.com/mulin0424/cepip
  - Tissue-specific effects on expression: https://github.com/FunctionLab/ExPecto
  - https://github.com/GuanLab/Leopard
- TFBS
  - Predict TF binding affinities using open chromatin + PWMs: https://github.com/schulzlab/TEPIC
  - LR-DNAse: TFBS prediction using features derived from DNase-seq: http://biorxiv.org/content/early/2016/10/24/082594

## Security

- Secret sharing schemes for keeping patient ID secret while being able to reconstruct it from other identifiers: https://pdfs.semanticscholar.org/0307/48be9820512a1d5582351552bd0452711296.pdf
  
## Sequencing Protocols

- Single cell
  - Simultaneous RNA and methylation (and inference of CNV): http://www.nature.com/cr/journal/vaop/ncurrent/full/cr201623a.html
  - Simultaneous RNA and methylation (scM&T-seq): http://www.nature.com/nmeth/journal/v13/n3/full/nmeth.3728.html
  - Simultaneous RNA and protein measurements: http://www.sciencedirect.com/science/article/pii/S2211124715014345
- CRISPR-DS: Uses CRISPR/Cas9 excision of target sequences to improve library quality relative to hybrid capture: https://www.biorxiv.org/content/early/2017/10/23/207027

## Simulation

- Gene expression: https://www.biorxiv.org/content/10.1101/836254v1
- Human genomes: https://github.com/callum-b/H2G2/

## Variant annotation

- WGSA pipeline https://sites.google.com/site/jpopgen/wgsa/
- SnpEff: http://snpeff.sourceforge.net/
- HAIL: https://hail.is/
- HGVS
  - This is a standard notation for variants http://varnomen.hgvs.org/
  - Toolkit for working with HGVS-formatted variants: https://github.com/biocommons/hgvs
  - check/correct https://github.com/mutalyzer/mutalyzer
- Integration of multiple predictive measures of mutation deleterious effect: https://www.ncbi.nlm.nih.gov/research/snpdelscore/
- Varsome: database that aggregates variant information from many sources
  - https://varsome.com/
  - Provides a REST API that can be queried with dbSNP or HGVS IDs
- Database of human chromosomal fragile sites: http://webs.iiitd.edu.in/raghava/humcfs
- FunSeq2: http://funseq2.gersteinlab.org
- Deep learning-based: https://www.biorxiv.org/content/early/2017/12/18/235655.1
- Search engine for indexing and intersecting with large numbers of annotation BED/VCF files: https://github.com/ryanlayer/giggle/
- Gene panel/disease-specific annotation and filtering: https://www.ebi.ac.uk/gene2phenotype/g2p_vep_plugin
- Cravat: https://github.com/KarchinLab/open-cravat/tree/master/cravat
- Pharmacogenomic: https://github.com/PharmGKB/PharmCAT
- SV
  - AnnotSV: http://lbgi.fr/AnnotSV/
  - Curation/annotation by depth: https://github.com/brentp/duphold
- Sequence context: https://github.com/elolab/VarSCAT
- dna-claude-analysis: personal genome analysis toolkit for raw DNA data (23andMe/AncestryDNA) https://github.com/shmlkv/dna-claude-analysis

## Sequence Analysis

### General-purpose

- Google Genomics R API: https://followthedata.wordpress.com/2015/02/05/notes-on-genomics-apis-2-google-genomics-api/
- Fetch sequences and metadata from SRA and ENA https://wwood.github.io/kingfisher-download/
- k-mer counting
  - khmer:
    - https://github.com/ged-lab/khmer
    - https://docs.google.com/presentation/d/1biQmLkwPlCOA56mNZdUAiXa1OyGE0qyvI2nvU0qlOIE/mobilepresent?pli=1&slide=id.p58
  - KCMBT: https://github.com/abdullah009/kcmbt_mt
  - KAT https://github.com/TGAC/KAT
  - ntCard https://github.com/bcgsc/ntCard
  - KMC3 https://github.com/refresh-bio/KMC
  - Encoding counts with kmer data: https://github.com/lzhLab/kmcEx
  - Gerbil: GPU accelerated, outputs count table https://github.com/uni-halle/gerbil
  - Toolkit for working with unique kmers: https://github.com/shenwei356/unikmer
  - Across large number of datasets https://github.com/kamimrcht/REINDEER
  - denoising https://github.com/Christina-hshi/CQF-deNoise
  - differential analysis: https://github.com/tlemane/kmdiff
- Sequence dimensionality reduction 
  - Strobemers https://www.biorxiv.org/content/10.1101/2021.01.28.428549v2
- Density-based clustering: https://bitbucket.org/jerry00/densitycut_dev
- GFFutils: http://daler.github.io/gffutils/
- R package for aligned chromatin-oriented sequencing data: https://cran.r-project.org/web/packages/Pasha/
- BAMQL: query language for extracting reads from BAM files https://github.com/BoutrosLaboratory/bamql
- SAMBAMBA: samtools alternative
- BAMtools: another samtools alternative, plus some additional tools https://github.com/pezmaster31/bamtools/wiki
- DeepTools: more useful SAM/BAM operations http://deeptools.readthedocs.io/en/latest/content/list_of_tools.html
- Genomic intervals
  - bedtools http://bedtools.readthedocs.io/en/latest/
  - bedops alternative/additional BED operations http://bedops.readthedocs.io/en/latest/
  - comparison of interval sets: https://github.com/deepstanding/seqpare
  - Intersection and visualization of multiple gene/region sets: https://bitbucket.org/CBGR/intervene
  - https://github.com/noamteyssier/gia
- Demultiplexing:
  - http://gbcs.embl.de/portal/tiki-index.php?page=Je
  - https://davidsongroup.github.io/flexiplex/
  - Axe: https://axe-demultiplexer.readthedocs.io/en/latest/usage.html
  - FuzzySplit: https://github.com/Daniel-Liu-c0deb0t/Java-Fuzzy-Search
  - https://github.com/pachterlab/splitcode
- Liftover
  - SNPs: https://github.com/enormandeau/snplift
  - Variants in confident regions: https://github.com/verilylifesciences/genomewarp
  - Segments: https://github.com/baudisgroup/segment-liftover
  - Alignments
    - https://github.com/CMU-SAFARI/AirLift
    - https://github.com/milkschen/leviosam2
  - Identying unstable positions between genome versions https://github.com/cathaloruaidh/genomeBuildConversion
  - remapping reads using a chain file https://github.com/CMU-SAFARI/FastRemap
- Blacklists
  - Tools to create ENCODE blacklists, and pre-computed blacklists for model organisms: https://github.com/Boyle-Lab/Blacklist
  - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-023-03027-3
  - Cut&Run green and blacklists for endogenous normalization https://github.com/fndemello/CUT-RUN_greenlist
- Server for reference sequences and indices: https://github.com/databio/refgenie
- Coverage
  - Fast coverage estimate from BAM index: https://github.com/brentp/goleft/tree/master/indexcov
  - Quantification and normalization of coverage peaks: https://github.com/ncbi/BAMscale
  - Library for distributed coverage calculation using Spark https://github.com/ZSI-Bio/bdg-sequila
  - https://github.com/telatin/bamtocov
  - https://github.com/brentp/mosdepth
  - https://github.com/sstadick/perbase
- Additional tools built on GATK4 https://bimberlab.github.io/DISCVRSeq/toolDoc/index.html
- Sort reads and remove duplicates (for compression and to accelerate alignment): https://github.com/bioinformatics-polito/BioSeqZip
- Detect clusters of coordinates https://bioconductor.org/packages/release/bioc/html/fcScan.html
- Remap alignments from one reference to another https://github.com/CMU-SAFARI/AirLift
- Sample quality classifier https://github.com/salbrec/seqQscorer
- grabseqs: https://github.com/louiejtaylor/grabseqs
- Read depth
  - https://bioconductor.org/packages/release/bioc/html/megadepth.html
  - https://bitbucket.org/ustenzel/bathometer/src
- https://github.com/pachterlab/ffq
- Toolkit for manipulation of interleaved FASTQ files https://github.com/lskatz/fasten
- Intel accelerated versions of some key bioinformatics tools https://github.com/IntelLabs/Open-Omics-Acceleration-Framework?utm_source=pocket_reader
- Simulation
  - https://galantelab.github.io/Sandy
  - https://github.com/adamewing/bamsurgeon
  - https://github.com/WongLabPrograms/COMPASS
  - SimuG: https://github.com/yjx1217/simuG
  - Multi-modal simulator: https://www.biorxiv.org/content/10.1101/2020.02.06.936971v1
  - Reads with UMIs: https://gitlab.com/vincent-sater/umigen-master
  - https://github.com/CompGenomeLab/boquila
  - https://github.com/qasimyu/simuscop
  - https://github.com/schmeing/ReSeq
  - https://github.com/RAHenriksen/NGSNGS
  - https://github.com/ncsa/NEAT

### QC

- Qualimap2: http://qualimap.bioinfo.cipf.es/
- Determine whether two BAM files are from the same source: https://bitbucket.org/sacgf/bam-matcher
- Identify and remove UMI sequences from reads: https://github.com/CGATOxford/UMI-tools
- Integrated report from multiple tools: http://multiqc.info/
- AlmostSignificant: https://github.com/bartongroup/AlmostSignificant
- Genetic relatedness from raw reads:
  - https://github.com/kdmurray91/kwip
  - https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi
- Detecting sample swaps:
  - https://github.com/PapenfussLab/HaveYouSwappedYourSamples
  - From FASTQ https://github.com/JustinChu/ntsm
- QC Fail articles:
  - Patterned flow cells (HiSeq 3000+) have high rates of optical duplicates: https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/
- Fingerprints: http://db.systemsbiology.net/gestalt/genome_fingerprints/
- Fastq.bio: http://fastq.bio/
- Read Origin Protocol: https://github.com/smangul1/rop
- QC on Illumina run folder: https://github.com/Molmed/checkQC
- New diploid reference for benchmarking: https://www.biorxiv.org/content/early/2017/11/22/223297
- Estimate library complexity: http://smithlabresearch.org/software/preseq/
- fluff: https://fluff.readthedocs.io/en/latest/introduction.html
- GenMap: compute mappability; pre-built indices https://github.com/cpockrandt/genmap
- Identify inter-sample contamination:
  - https://github.com/Griffan/VerifyBamID
  - From VCF: https://github.com/HTGenomeAnalysisUnit/SCE-VCF
- Learn features from manual ENCODE dataset classifications and use to automate QC https://github.com/salbrec/seqQscorer
- Empirical quality filtering: https://github.com/bodeolukolu/ngsComposer
- https://quinlan-lab.github.io/bedqc
- https://github.com/pachterlab/splitcode
- https://github.com/uclahs-cds/package-PipeVal
- Quantify and predict reference bias: https://github.com/maojanlin/biastools

### Peak-based

* Handle multi-mapping reads https://github.com/seqcode/allo

#### ChIP-seq

- Pre-processing
  - Quality assessment: https://github.com/rnakato/SSP
  - Allocation of multi-mapping reads: https://github.com/keleslab/permseq
  - Mapping: https://github.com/haowenz/chromap
- Peak calling
  - GC-aware peak caller: https://bioconductor.org/packages/devel/bioc/html/gcapc.htmlR
  - GenoGAM peak caller: https://master.bioconductor.org/packages/3.3/bioc/html/GenoGAM.html
  - De-noising: https://github.com/kundajelab/TF_chipseq_pipeline
  - Peak discretizer (merging of replicates): https://github.com/nanakiksc/zerone
  - Compute error: https://github.com/tdhock/PeakError
  - Specifically for histone modifications: https://github.com/Bohdan-Khomtchouk/SUPERmerge
  - ChIPWig: https://github.com/vidarmehr/ChIPWig
  - https://github.com/biocore-ntnu/epic2
- hiddenDomains: https://sourceforge.net/projects/hiddendomains/
- Network-based identification of relationships among ChIP-seq data sets: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0925-0
- Database of labeled ChIP-seq peaks: http://cbio.ensmp.fr/thocking/chip-seq-chunk-db/ (error of peak calls computed using https://github.com/tdhock/PeakError)
- Find different modes of binding: https://narlikarlab.github.io/DIVERSITY/
- Functional analysis
  - FunChIP: http://bioconductor.org/packages/FunChIP
  - annoPeak: https://github.com/XingTang2014/annoPeak
- K-mer-based alternative to PWMs for predicting TF binding sites: http://groups.csail.mit.edu/cgs/gem/kmac/
- Recalibrate p-values from peak callers: https://github.com/theodorejperkins/RECAP
- Peak quantification (for any data type) https://github.com/ncbi/BAMscale
- Nucleotide-matched background sequences for motif enrichment https://bitbucket.org/CBGR/biasaway/src/master/
- Simulation
  - http://www.bioconductor.org/packages/release/bioc/html/ChIPsim.html
  - ChIPulate: https://github.com/vishakad/chipulate
  - Tulip: https://github.com/gymreklab/Tulip
  - Chips: https://github.com/gymreklab/chips

#### Chromatin accessibility

- HINT: http://costalab.org/publications-2/dh-hmm/ (was best out of 10 compared tools in recent NatMeth paper)
  - HINT-ATAC: method specifically for ATAC-seq footprinting
  - https://github.com/CostaLab/reg-gen
- ALTRE: https://mathelab.github.io/ALTRE/vignette.html
- LR-DNAse: TFBS prediction using features derived from DNase-seq: http://biorxiv.org/content/early/2016/10/24/082594
- Identify accessible chromatin from NOMe-seq https://sourceforge.net/projects/came/
- Predict targets of accessible regions: https://github.com/cole-trapnell-lab/cicero-release
- Correcting transposition bias in ATAC-seq: https://www.biorxiv.org/content/biorxiv/early/2019/01/
- ATAC-seq
  - MACS settings https://github.com/macs3-project/MACS/discussions/435
  - ATAC-seq peak caller: https://github.com/LiuLabUB/HMMRATAC
  - esATAC: ATAC-Seq pipeline https://bioconductor.org/packages/release/bioc/html/esATAC.html
  - ATAC-Seq QC: https://bioconductor.org/packages/release/bioc/html/ATACseqQC.html
  - Denoising of low-covereage/noisy ATAC-seq: https://github.com/clara-genomics/AtacWorks
22/525808.full.pdf
  - AIAP: https://www.biorxiv.org/content/10.1101/686808v1
    - Has associated viewer
      - https://github.com/lidaof/qATACviewer/tree/localjson
      - https://qa.targetepigenomics.org/
  - PEPATAC http://pepatac.databio.org/en/latest/
  - https://github.com/Tuteja-Lab/UnsupervisedPeakCaller
  - Cactus: ATAC-Seq and transcriptomics pipeline https://github.com/jsalignon/cactus
- Deep learning-based peak caller https://github.com/LHentges/LanceOtron
- https://github.com/Dongwon-Lee/gkmQC
- https://github.com/jvierstra/footprint-tools

### Chromatin Interactions

- Resolution enhancement: http://dna.cs.miami.edu/HiCNN/
- Binary storage format for interaction matrix: https://github.com/mirnylab/cooler
- qc: https://github.com/cerebis/qc3C
- Aligner: https://github.com/xiaochuanle/Falign
- Simulation
  - https://github.com/cerebis/sim3C
  - https://github.com/ay-lab/AveSim
  - https://github.com/yezhengSTAT/FreeHiC

### DNA

- Read Filtering
  - Removing redundant reads from deep sequencing data: https://git.informatik.uni-kiel.de/axw/Bignorm
  - Duplicate removal
    - Pardre: https://sourceforge.net/projects/pardre/
    - Nubeam-dedup: https://github.com/daihang16/nubeamdedup
    - https://github.com/yuansliu/minirmd
    - WGS and WGBS: https://github.com/huishenlab/dupsifter
- Error correction
  - BFC: https://github.com/lh3/bfc
  - RECKONER http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=reckoner&subpage=about
  - MultiRes - says it's for viral populations; not sure if applicable to humans: https://github.com/raunaq-m/MultiRes
  - http://github.com/Malfoy/BCOOL
  - Correction only of reads that contain repetative sequences: https://github.com/biointec/browniecorrector
  - https://github.com/fkallen/CARE
  - Long read
    - https://github.com/lbroseus/TALC
    - https://github.com/DecodeGenetics/Ratatosk
  - UMIs
    - https://github.com/Nik-Tavakolian/Shepherd
- Alignment
  - Recommendations:
    - https://github.com/CCDG/Pipeline-Standardization
    - These are optimized for germline variant calling, however there are a few interesting pointss
      - BWA MEM has a hidden option that enables "deterministic alignment results" (-K 100000000)
      - They recommend binning of quality scores for greater compression
        - Scores 2-6 correspond to Illumina error codes and are left as-is
        - 10, 20, 30
        - Round scores to the nearest value in probability space
      - The have some recommendations for standardizing headers
  - AGA: aligner specifically for coding sequences that takes into account the translated protein sequence: https://github.com/emweb/aga
  - Scaling alignment to massively parallel processors
    - https://www.biorxiv.org/content/early/2017/10/24/205328
    - Parallelizing hundreds of threads on Xeon Phi (64-72 core processors)
  - Alignment speed-up using PIM technologies: https://arxiv.org/abs/1711.01177
  - Whisper
    - Accelerates alignment by pre-sorting reads by sequence, and only aligning non-duplicates
    - http://sun.aei.polsl.pl/REFRESH/Whisper
    - Paper: https://www.biorxiv.org/content/early/2017/12/28/240358?rss=1&utm_source=dlvr.it&utm_medium=twitter
  - Align simultaneously against multiple reference genomes http://1001genomes.org/software/genomemapper.html
  - Compressed reference-based alignment: http://groups.csail.mit.edu/cb/cora/
  - Compression and querying of aligned haplotype data: https://github.com/richarddurbin/pbwt
  - LibGaba: https://github.com/ocxtal/libgaba
  - Minimap2: https://github.com/lh3/minimap2
  - SNP-aware: https://github.com/weiquan/SALT
  - Long read
    - https://github.com/MohammadJRS/SureMap
    - https://github.com/isovic/graphmap
    - https://github.com/ocxtal/minialign
    - for repeat regions https://github.com/marbl/Winnowmap
    - https://github.com/NWPU-903PR/smsMap
    - Vulcan - dual-mode alignment https://gitlab.com/treangenlab/vulcan
    - https://github.com/ChaissonLab/LRA
  - Graph-based:
    - https://github.com/maickrau/GraphAligner
    - Glia (mainly for local realignment): https://github.com/ekg/glia
    - Raptor: https://github.com/isovic/raptor
    - https://github.com/lh3/minigraph
    - Selecting variants for graph genome: https://www.biorxiv.org/content/early/2018/04/30/311720
    - https://github.com/maickrau/GraphAligner
    - https://github.com/ratschlab/metagraph
  - Pan-genome
    - Build pan-genome index https://github.com/SJTU-CGM/HUPAN
    - PanGenie https://bitbucket.org/jana_ebler/pangenie/src/master/
  - MECAT: https://github.com/xiaochuanle/MECAT
  - Correct sex chromosome-related errors/artifacts:
    - https://github.com/WilsonSayresLab/XYalign
    - Align samples with a Y chromosome to a YPAR-masked genome, and samples without a Y chromosome to a Y-masked genome https://www.biorxiv.org/content/10.1101/668376v1
  - Major-allele SNP references for Bowtie: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
    - Realign poorly aligning reads to population-specific major-allele references https://www.biorxiv.org/content/10.1101/2020.03.03.975219v3.full.pdf
  - Mapper with accurate MAPQ scores https://github.com/gui11aume/mmp
  - Heuristic-free aligner with guaranteed optimal alignment https://github.com/langmead-lab/vargas
  - https://github.com/ekg/wfmash
  - Other
    - Optical computing: https://www.optalysys.com/
    - FPGA: http://edicogenome.com/dragen-bioit-platform/
    - Completeness statistics for MSA: https://github.com/thomaskf/AliStat
    - Suffix array with learned data model to accelerate queries: https://github.com/mkirsche/sapling
    - Prediction of identity scores for large datasets https://github.com/BioinformaticsToolsmith/Identity
    - Puffaligner: https://github.com/COMBINE-lab/pufferfish/tree/cigar-strings
- Assembly
  - Build de Bruijn Graph from multiple genomes: https://github.com/medvedevgroup/TwoPaCo
  - Align to a de Brujn graph: https://github.com/Malfoy/BGREAT
  - Cosmo/VARI: Assembler using succinct colored de Bruijn graphs to encode population information https://github.com/cosmo-team/cosmo/tree/VARI
  - Streaming: https://github.com/Shamir-Lab/Faucet
  - BFGraph: https://github.com/pmelsted/bfgraph
  - Spark-based: https://github.com/BioHPC/SORA/wiki
  - Polishing: https://github.com/bcgsc/ntedit
  - Long-read
    - [Ra](https://github.com/lbcb-sci/ra)
    - [Raven](https://github.com/lbcb-sci/raven)
    - https://github.com/phasebook/phasebook
  - Pre-masking repeats: https://github.com/mills-lab/PALMER
  - Pan-genome graphs
    - https://github.com/lh3/minigraph
    - Extracting paths using haplotype information https://github.com/tomokveld/CHOP
    - https://github.com/LeilyR/Multi-genome-Reference
  - New measure of assembly contiguity https://lh3.github.io/2020/04/08/a-new-metric-on-assembly-contiguity
  - Linked reads: https://github.com/lauren-mak/Ariadne
  - Assemble paired reads into pseudo-long-reads: https://github.com/fkallen/CAREx
- Post-alignment
  - QC: https://github.com/achon/sramm
  - Base quality recalibration:
    - https://github.com/swainechen/lacer
    - https://github.com/broadinstitute/firepony
  - Indel realignment (also supports RNAseq): https://github.com/mozack/abra2
  - Error correction
    - Racon: https://github.com/isovic/racon
  - Mageri: https://github.com/mikessh/mageri
  - ReadStomper: https://github.com/gringer/bioinfscripts/blob/master/readstomper.pl
  - ProDuSe: https://github.com/morinlab/ProDuSe/
  - AGeNT: https://www.genomics.agilent.com/en/NGS-Data-Analysis-Software/AGeNT/?cid=AG-PT-154&tabId=prod2570007
  - Prepare for variant calling: https://github.com/ExaScience/elprep
    - replacement for sort, markdup, BQSR
  - Fast pileup of targeted regions in C/Python: https://github.com/brentp/hileup
  - Compression of quality scores:
    - Lossy https://github.com/jkbonfield/crumble
    - Lossless: https://github.com/SCUT-CCNL/LCQS
  - Denoising (error correction + quality score adjustment): https://github.com/ihwang/SAMDUDE
  - Fix reference allele bias of the Isaac aligner https://github.com/danchubb/FixVAF
  - Merge duplicate reads: https://github.com/alevar/tiebrush
  - remove cDNA contamination https://github.com/rheinbaylab/cDNA-detector
  - SQL-like read querying and filtering https://github.com/maragkakislab/samql/
  - Error stats https://github.com/google/best
  - Duplicate marking
    - https://github.com/grailbio/doppelmark
    - https://github.com/lolrenceH/TrieDedup
  - Long read
    - Realigns reads in segmental duplication regions https://gitlab.com/tprodanov/duplomap
- Variant calling
  - QC
    - https://github.com/nebfield/snpQT
  - Reference-free:
    - Kevlar: https://github.com/dib-lab/kevlar
    - RUFUS: https://github.com/jandrewrfarrell/RUFUS
  - Deep learning
    - DeepVariant: https://github.com/google/deepvariant
      - https://blog.dnanexus.com/2017-12-05-evaluating-deepvariant-googles-machine-learning-variant-caller/?utm_source=email&utm_medium=mailchimp&utm_campaign=december_newsletter
      - https://blog.dnanexus.com/2018-01-16-evaluating-the-performance-of-ngs-pipelines-on-noisy-wgs-data/
      - Currently it doesn't handle somatic calling 
      - Trio calling pipeline https://github.com/VCCRI/dv-trio/
    - Deep learning will increase accuracy in the nest GATK release: https://gatkforums.broadinstitute.org/gatk/discussion/10996/deep-learning-in-gatk4
    - Nvidia https://github.com/clara-genomics/DL4VC
    - Train DNN variant calling models (commercial): https://magnolia.sh/
    - https://github.com/kaist-ina/RUN-DVC
    - https://github.com/brendanofallon/jovian
  - Benchmarking
    - GA4GH/PrecisionFDA guidelines and tools for germline variants
    - https://www.biorxiv.org/content/early/2018/02/23/270157
    - https://github.com/ga4gh/benchmarking-tools
    - https://platform.dnanexus.com/login?scope=%7B%22full%22%3A+true%7D&redirect_uri=https%3A%2F%2Fprecision.fda.gov%2Freturn_from_login&client_id=precision_fda_gov
    - https://github.com/mazzalab-ieo/recallme
    - https://github.com/TimD1/vcfdist
  - Ococo
    - Variant calling is performed in a streaming fashion, with evidence for variation being updated from each read as it is mapped. Theoretically, new evidence can be ignored once the confidence level reaches a certain threshold. In practice, they show this is about 10% of reads. However, they only deal with SNPs. It would be interesting to investigate whether:
    - This can be extended to somatic variant calling
    - This can work with rare variants (eg cfDNA)
    - https://github.com/karel-brinda/ococo
    - Paper: https://arxiv.org/pdf/1712.01146.pdf
  - Distributed variant calling using Spark:
    - http://bdgenomics.org/projects/
    - https://github.com/atgenomix/connectedreads
    - https://halvadeforspark.readthedocs.io/en/latest/
  - De-novo variant caller: https://github.com/bgm-cwg/novoCaller
  - Create population-specific blacklists: http://lab.rockefeller.edu/casanova/BL/programs
  - From long/linked reads:
    - https://github.com/pjedge/longshot
    - https://github.com/maiziex/Aquila
    - Qscore correction prior to SNP calling from Nanopore https://github.com/rocioespci/minIONLogReg
    - ONT https://www.biorxiv.org/content/10.1101/2020.10.22.350009v1
  - Hybrid
    - https://github.com/anands-repo/hello
  - Scotch: ML approach for calling intermediate-sized indels: https://github.com/AshleyLab/scotch
  - Encoding variant information in BWT: https://www.uni-ulm.de/in/theo/research/seqana/
  - Detect variants by comparing reads between samples without genome alignment: https://github.com/akiomiyao/ped
  - Align and call simultaneously: https://github.com/hsinnan75/MapCaller
  - For reads with UMIs: https://gitlab.com/vincent-sater/umi-varcal_v3
  - Varlociraptor: call given candidate variants and filter based on FDR https://varlociraptor.github.io/docs/calling/
    - Works for somatic or an arbitrary experimental design (described in a config file)
  - Wrapper around multiple variant callers that simplifies parallelization https://github.com/lh3/unicall
  - For finding mutations in mutated sample (vs normal): https://github.com/ykat0/bivartect
  - Rare variants: https://github.com/frankMusacchia/VarGenius-HZD
  - Single-cell
    - https://github.com/prosolo/prosolo
- de novo mutations (trio variant calling)
  - https://github.com/jbelyeu/unfazed
  - https://github.com/sbg/geck
  - Prediction of de novo mutations from trio variant calling data https://github.com/james-guevara/synthdnm
  - https://github.com/yongzhuang/DNMFilter_Indel
  - Using multiple reference assemblies: https://github.com/mcveanlab/Corticall
  - filtering: https://github.com/brentp/slivar
  - DeepTrio https://github.com/google/deepvariant/tree/r1.1/deeptrio
  - https://github.com/Genome-Bioinformatics-RadboudUMC/DeNovoCNN
  - VG-DeepTrio https://github.com/vgteam/vg_wdl
  - https://github.com/TNTurnerLab/HAT
- Variant filtering/post-processing
  - Matching variant sets: https://github.com/medvedevgroup/varmatch
  - Post-processing variant calls to determine whether variants at regions with alternative loci have allele(s) from an alternate locus: https://github.com/charite/asdpex
  - Filtering low-frequency variants that likely result from DNA damage: https://github.com/eilslabs/DKFZBiasFilter
  - IonTorrent: https://www.github.com/aditya-88/intelli-ngs
  - GARFIELD: Using deep learning model to filter false-positive variants
    - https://www.biorxiv.org/content/biorxiv/early/2017/09/11/149146.full.pdf
    - Pre-trained (using GIAB) Illumina and IonTorrent models are here: https://github.com/gedoardo83/GARFIELD-NGS
  - ISOWN
    - Distinguishes somatic from germline variants
    - https://github.com/ikalatskaya/ISOWN
    - Paper: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0446-9
  - Bystro: annotation and filtering (web only)
    - https://bystro.io/
    - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1387-3
  - SevenBridges': https://github.com/sbg/smart-variant-filtering
  - Skyhawk: https://github.com/aquaskyline/Skyhawk
  - This paper describes criteria that identify 100% of FP (as determined by Sanger sequencing) in a large clinical dataset https://www.biorxiv.org/content/early/2018/05/31/335950
  - Regions subject to systematic bias: 
    - https://www.biorxiv.org/content/10.1101/679423v2
    - https://figshare.shef.ac.uk/articles/Suspect_loci_BED_files/9995945
  - Another RF classifier: https://github.com/avallonking/ForestQC
  - snpToolkit: https://snptoolkit.readthedocs.io/en/latest/index.html
  - VariFAST https://github.com/bioxsjtu/VariFAST
  - Merfin: https://github.com/arangrhie/merfin
  - Sample control variants: https://github.com/mulinlab/vSampler
  - Detect contamination https://github.com/Ancestry/VCFcontam
  - Merge SNVs into MNVs https://github.com/Sentieon/sentieon-scripts/tree/master/merge_mnp
  - Long read: https://github.com/arangrhie/merfin/
  - https://github.com/brentp/echtvar
  - Summarize, compare, visualize VCFs: https://github.com/HimanshuLab/SCI-VCF
  - https://github.com/ideateknoloji/FPDetect
- Genotyping
  - ebGenotyping: https://cran.r-project.org/web/packages/ebGenotyping/ebGenotyping.pdf
  - FastGT: http://bioinfo.ut.ee/FastGT/
  - Compression and search of genotype data: http://sun.aei.polsl.pl/REFRESH/gtc
  - Using variation graph: https://github.com/Kingsford-Group/variantstore
  - Mosaic variants: https://github.com/abyzovlab/Leucippus
  - HISAT-genotype: http://ccb.jhu.edu/hisat-genotype/index.php/Main_Page
  - Trio-based concordance: https://github.com/sbg/VBT-TrioAnalysis
  - Using a graph reference: https://github.com/bioinformatics-centre/BayesTyper
  - https://github.com/AlgoLab/malva
  - https://bioinfo.ut.ee/KATK/
  - https://github.com/ivargr/kage
  - Using haplotype-resolved pan-genome graph https://github.com/eblerjana/pangenie
  - https://github.com/kage-genotyper/kage
  - Low-pass/GBS
    - https://github.com/stevenandrewyates/GBSmode
- Imputation
  - [loimpute](https://gitlab.com/gencove/loimpute-public): impute from low-pass sequencing
  - GLIMPSE: https://odelaneau.github.io/GLIMPSE/
  - Impute parental genotypes from siblings https://github.com/AlexTISYoung/SNIPar
  - QUILT: https://github.com/rwdavies/QUILT
  - https://github.com/OpenMendel/MendelImpute.jl
  - https://github.com/shilab/STI
- SVs/CNV **calling**
  - Score SVs based on predicted functional impact https://github.com/lganel/SVScore
  - Pipelines for CNV/SV calling:
    - MARATHON: https://github.com/yuchaojiang/MARATHON
    - Anaconda: https://mcg.ustc.edu.cn/bsc/ANACONDA/
    - CODEX2: https://github.com/yuchaojiang/CODEX2
    - Ximmer: https://github.com/ssadedin/ximmer
  - Reporting
    - https://github.com/JH-Zhou/HandyCNV
    - Tool for online review of SVs: https://github.com/jbelyeu/SV-plaudit
  - Using off-target reads to call CNVs: https://ep70.eventpilotadmin.com/web/page.php?page=IntHtml&project=ASHG17&id=170120230
  - GRIDSS: detects rearrangements using du Bruijn graph assembly: https://github.com/PapenfussLab/gridss
  - ichorCNA: https://github.com/broadinstitute/ichorCNA
  - indelope: https://github.com/brentp/indelope
  - MultiGEMS: https://github.com/cui-lab/multigems
  - WHAM: CNV caller https://github.com/zeeev/wham
  - Identification of mosic events: https://github.com/asifrim/mrmosaic
    - This may be useful for calling allele-specific (i.e. heterozygous) copy number variation: https://drive.google.com/open?id=0B7E7HSjQ-SumQmlPc1Z0aUR5Sk0
  - Benchmarking of pipelines for calling complex variants: https://www.biorxiv.org/content/early/2017/11/23/218529
  - Manta and novoBreak performed best in DREME challenge for SV calling: https://www.biorxiv.org/content/early/2017/11/25/224733
  - https://github.com/Mesh89/SurVIndel
  - Allele-specific:
    - Falcon: allele-specific copy number changes https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4344483/
    - AS-GENSENG: https://sourceforge.net/projects/asgenseng/
    - ASCAT: https://github.com/Crick-CancerGenomics/ascat/tree/master/ASCAT
    - http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075350
    - TAPS (for affy arrays): https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-10-r108
  - Mappability tracks, for normalizing CNV calls http://genome-test.soe.ucsc.edu/cgi-bin/hgTrackUi?hgsid=392476445_0mwMQdJgiLvv8WtNo4YktFt41IkI&c=chr1&g=umap
  - WALDO: detection of LINE amplifications http://www.pnas.org/content/115/8/1871.short
  - Smoove: https://github.com/brentp/smoove
  - Svaba: https://github.com/walaj/svaba
  - Consensus caller: https://github.com/TheJacksonLaboratory/SVE
  - https://github.com/sandmanns/CopyDetective
  - Randomized algorithm to speed up CNV calling: https://sourceforge.net/projects/genseng/
  - CNV calling from plasma: https://www.biorxiv.org/content/early/2018/03/28/290171
  - https://github.com/kunalkathuria/SVXplorer
  - https://lh3.github.io/2021/01/11/minigraph-as-a-multi-assembly-sv-caller
  - https://github.com/hartwigmedical/hmftools/tree/master/sv-linx
  - https://github.com/ETCHING-team/ETCHING
  - https://github.com/weiyuchung/CONY
  - from haploid or diploid genome-genome alignments https://github.com/eldariont/svim-asm
  - merging, filtering, feature extraction: https://github.com/dermasugita/Viola-SV
  - https://github.com/kcleal/dysgu
  - https://github.com/tf2/CNest
  - https://github.com/alok123t/HyINDEL
  - https://github.com/yaozhong/deepIntraSV
  - https://github.com/TheJacksonLaboratory/JAX-CNV
  - https://github.com/ETCHING-team/ETCHING
  - CNVnator
    - Python extension: https://github.com/abyzovlab/CNVpytor
  - https://github.com/imgag/ClinCNV
  - https://github.com/lxwgcool/EigenDel
  - https://github.com/ITBE-Lab/MA
  - Ensemble caller:
    - https://github.com/HaoKeLab/ensembleCNV
    - Parliament2 https://github.com/dnanexus/parliament2
    - Different variant callers excel at calling different types of SVs: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1720-5
    - https://github.com/wonderful1/CNV-PG
    - https://bitbucket.org/sonnhammergroup/metacnv
    - https://github.com/SFGLab/ConsensuSV-pipeline
    - https://github.com/Mangul-Lab-USC/VISTA
  - Population-scale:
    - https://github.com/hall-lab/svtools
    - STIX: A tool for determining the frequency of SVs in population data https://github.com/ryanlayer/stix
    - https://github.com/kehrlab/PopDel
  - ML/Deep learning-based callers:
    - https://github.com/CSuperlei/DeepSV
    - Specifically for duplications and deletions: https://github.com/tomh1lll/dudeml
    - DL-CNV: https://www.aimspress.com/article/10.3934/mbe.2020011/fulltext.html
    - https://github.com/lxwgcool/EigenDel
    - https://github.com/mlinderm/npsv
    - https://github.com/luojunwei/BreakNet
    - https://github.com/SimCab-CHU/ifCNV
    - https://github.com/PopicLab/cue
    - https://github.com/xjtu-omics/SVision
    - https://github.com/mhuidong/cnnLSV
    - https://github.com/eioyuou/INSnet
  - Fusion detection: https://github.com/sgilab/JuLI
  - Genotyping
    - Benchmarks: https://academic.oup.com/gigascience/article/8/9/giz110/5565134
    - Alignment-free SV genotyping: https://github.com/Parsoa/NebulousSerendipity
    - Genotyping tandem repeats: https://github.com/mcfrith/tandem-genotypes
    - Graph-based genotyping: https://github.com/illumina/paragraph
    - Long-read: https://github.com/llecompte/SVJedi
    - Deep learning-based https://github.com/mchowdh200/samplot-ml
    - https://github.com/kage-genotyper/kage/
    - https://github.com/kehrlab/Complex-SV-Genotyping
  - Annotation
    - ACMG https://github.com/Genotek/ClassifyCNV
    - https://github.com/hartwigmedical/hmftools/tree/master/linx
    - Tandem repeats: https://github.com/TravisWheelerLab/ULTRA
  - Merging
    - [Jasmine](https://github.com/mkirsche/Jasmine)
    - https://github.com/stjude/indelPost
    - https://github.com/dermasugita/Viola-SV
    - https://github.com/maiziex/MARS
    - https://github.com/ACEnglish/truvari
  - Linked read: 
    - https://github.com/sgreer77/gemtools
    - https://github.com/BilkentCompGen/valor
  - Savvy: uses off-target reads https://github.com/rdemolgen/SavvySuite
  - svtk: utilies for consolidating, filtering, resolving, and annotating structural variants 
    - https://github.com/xjtu-omics/Makohttps://github.com/talkowski-lab/svtk
  - Alignment-free: 
    - https://github.com/smortezah/smashpp
    - https://github.com/Parsoa/Nebula (SVs)
    - https://github.com/nicolaprezza/ebwt2InDel
  - UMIs https://gitlab.com/pierrejulien.viailly/mcna/
  - Exome/targeted
    - Exon-level CNV calling in exome/targeted seqeuencing: https://github.com/SD-Genomics/DeviCNV
    - Ensemble method: https://github.com/girirajanlab/CN_Learn
    - Reference sample selection: https://www.biorxiv.org/content/early/2018/11/25/478313.1
    - Generate consensus reads (optionally using UMIs): https://github.com/OpenGene/gencore
    - Breakpoint prediction: https://github.com/SinOncology/Break
    - Atlas-CNV: https://github.com/theodorc/atlas-cnv
    - https://github.com/fredsanto/coverageMaster
    - ML model for polishing WES CNV calls https://github.com/ciceklab/DECoNT
    - https://gitlab.com/sequentiateampublic/isocnv
    - https://github.com/bihealth/clear-cnv
    - https://github.com/SimCab-CHU/ifCNV
    - Deep learning https://github.com/ciceklab/ECOLE
  - Detect non-reference sequence variants from population data https://github.com/kehrlab/PopIns2
  - https://github.com/tprodanov/parascopy
- STR
  - Pipeline to estimate mutational parameters for every STR in the human genome: https://www.nature.com/ng/journal/v49/n10/full/ng.3952.html
  - MIRMMR: MSI caller
    - https://github.com/ding-lab/MIRMMR
    - Uses methylation and mutation information
  - Targeted STR profiling: http://darwin.informatics.indiana.edu/str/
  - Krait: integrated microsatellite detection and primer design
    - https://github.com/lmdu/krait
    - https://www.biorxiv.org/content/early/2017/11/18/221754
  - MANTIS: https://github.com/OSU-SRLab/MANTIS
  - MSIsensor: https://github.com/ding-lab/msisensor
  - https://www.biorxiv.org/content/early/2018/01/10/246108
    - They generate a reference STR database using PacBio
    - They provide software to estimate STR number from short-read data using the reference database
    - This could be used to quickly extract STR-containing sequences from a BAM: https://github.com/rkmlab/perf 
  - AdVNTR: https://github.com/mehrdadbakhtiari/adVNTR
  - Tandem genotypes: https://github.com/mcfrith/tandem-genotypes
  - lobSTR: http://melissagymrek.com/lobstr-code/
  - GangSTR: Genotyping STR longer than read length: https://github.com/gymreklab/GangSTR
  - Alignment-free STR genotyping: https://github.com/jbudis/dante
  - https://github.com/rlmcl/rescan
  - https://github.com/bahlolab/superSTR
- HLA
  - arcasHLA: https://github.com/RabadanLab/arcasHLA
  - https://github.com/WansonChoi/HATK
  - https://orthanq.github.io/
- mtDNA
  - https://github.com/linzhi2013/MitoZ
  - Haplogroup classification: https://github.com/seppinho/haplogrep-cmd
  - https://github.com/KCCG/mity
- Pipelines
  - SpeedSeq: alignment/annotation pipeline - https://github.com/hall-lab/speedseq
  - GenomeVIP: cloud-based pipeline https://github.com/ding-lab/GenomeVIP/
  - superFreq: https://github.com/ChristofferFlensburg/superFreq
  - GROM: https://osf.io/6rtws/
  - Amplicon/targeted
    - Canary: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1950-z
    - pyAmpli: https://mbeyens.github.io/pyAmpli/
    - amplimap: https://github.com/koelling/amplimap
  - DNAScan
    - https://github.com/KHP-Informatics/DNAscan
    - Align with HiSAT2; realign soft-clipped and unaligned reads with BWA MEM due to poor performance of indel calling from HISAT2-aligned data
    - Variant calling with Freebayes, extract indel positions and realign with GATK HC, perform hard filtering
    - Manta and Expansion Hunter for SVs
  - PECaller/PEMapper: population caller, does not produce BAM files https://github.com/wingolab-org/pecaller
  - Imputation of personal reference https://github.com/kvaddad1/impute-first
- Ancestry and kinship analysis
  - AKT: http://illumina.github.io/akt/
  - MixFit: http://www.geenivaramu.ee/en/tools/mixfit
  - Genotype clustering and ethnicity prediction using Deep Belief Networks: https://arxiv.org/abs/1805.12218
  - Truffle: IBD and ancestry analysis https://adimitromanolakis.github.io/truffle-website/
  - ngsRelate: https://github.com/ANGSD/ngsRelate
  - Fast D-statistics: https://github.com/millanek/Dsuite
  - Compute genome fingerprint from VCF: https://github.com/gglusman/genome-fingerprints
  - https://github.com/daviddaiweizhang/fraposa
  - https://github.com/AI-sandbox/gnomix
  - https://github.com/williamslab/papi
  - https://github.com/healthdisparities/rye
  - https://github.com/Rosemeis/HaploNet
  - https://github.com/nch-igm/snvstory
  - https://github.com/GenomicRisk/aeon
- Phasing/Haplotyping
  - Eagle2: https://data.broadinstitute.org/alkesgroup/Eagle/
  - https://github.com/ikmb/eagleimp
  - Using HiC+partial haplotypes: https://github.com/YakhiniGroup/SpectraPh
  - PhaseME http://beehive.cs.princeton.edu/wiki/phaseme/
  - HapCut2 (unclear if this works with standard Illumina WGS) https://github.com/vibansal/HapCUT2
  - hap.py: Illumina tool for identifying local haplotypes (also does genotype call comparison)
    - https://github.com/Illumina/hap.py/blob/master/doc/happy.md
  - This might be useful in multiBP merging and in detecting contamination
    - For single individuals: https://github.com/jcna99/PEATH
    - https://arxiv.org/abs/1801.09864
  - Alignment-free:
    - Uses kmer frequency: https://github.com/paudano/kestrel
  - IVC
    - https://github.com/namsyvo/IVC
    - Integrates alignment and variant calling
    - Uses known variant information
    - Germline, Illumina paired-end only
  - Using decision tree: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2147-9
  - PEATH: https://github.com/jcna99/PEATH
  - SNP partitioning by halpotype block: https://bioconductor.org/packages/release/bioc/html/gpart.html
  - Using RNA-seq allele-specific expression https://github.com/0xTCG/haptreex/
  - For de novo SNVs and SVs in trios https://github.com/jbelyeu/unfazed
  - Consensus: https://github.com/ziadbkh/consHap
  - Estimation from allele frequency data https://github.com/MartaPelizzola/haploSep
  - Long reads: https://github.com/twolinin/LongPhase/
  - https://github.com/rhysf/HaplotypeTools
  - https://github.com/dmiller903/trioPhaser
  - https://github.com/twolinin/LongPhase/
  - https://ngs-smap.readthedocs.io/en/latest/
- Graph
  - PanVC: https://gitlab.com/dvalenzu/PanVC
  - Graph genome alignment and variant calling on SevenBridges: https://www.sevenbridges.com/graph/
  - pangenome graphs
    - https://github.com/lh3/minigraph
    - https://github.com/pangenome/pggb
    - https://github.com/lh3/gfatools
- Non-Illumina
  - Adapter between aligners and downstream tools for handling PacBio quality values
    - https://github.com/mchaisso/pbsamstream
  - BGI compared to Illumina: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190264
- Dimensionality reduction
  - https://github.com/yhg926/public_kssd
  - https://github.com/medvedevgroup/UST/
- Sex determination
  - https://github.com/liusihan/seGMM
  - XYalign: sex chromosome copy number estimation from sequencing https://github.com/SexChrLab/XYalign
- Simulation
  - https://github.com/jstjohn/SimSeq
  - PERsim (part of ngs-bits): enables targeting of regions of the genome using bed file
  - EAGLE: also tries to simulate different sequencers: https://github.com/sequencing/EAGLE
  - http://www.nature.com/nrg/journal/vaop/ncurrent/pdf/nrg.2016.57.pdf
    - Look at: AFG, CureSIM, MetaSim, simNGS, simhtsd, SInC, NeSSM, Grinder (amplicon sequencing), FASTQSim (can simulate STRs)
  - Bamgineer: simulate allele-specific CNVs https://github.com/pughlab/bamgineer
  - Simulate mutations in a BAM file: https://github.com/adamewing/bamsurgeon
  - SNVs
    - https://github.com/davidebolo1993/VISOR
    - Rare: https://github.com/meganmichelle/RAREsim
    - simulate multiple types of DNA variants: https://github.com/yjx1217/simuG
  - CNV
    - https://github.com/YJulyXing/SECNVs
  - SVs:
    - https://github.com/xyc0813/pysim/
    - SVGen
    - https://bitbucket.org/charade/svengine
    - https://github.com/color/clrsvsim
    - https://github.com/ndierckx/Sim-it
  - Capture sequencing
    - https://github.com/mdcao/capsim
  - From trees: https://github.com/snacktavish/TreeToReads
  - https://github.com/price0416/Simulome
  - https://github.com/adimitromanolakis/sim1000G
  - generate synthetic genomes: https://www.biorxiv.org/content/10.1101/769091v2
  - Simulator that better approximates k-mer frequencies of Illumina sequencers https://github.com/schmeing/ReSequenceR
  - Generate artificial genomes using DNN model https://www.biorxiv.org/content/10.1101/769091v2.full.pdf
  - https://github.com/umr-ds/mesa_dna_sim
- Other
  - VCF compression and data extraction: https://github.com/kedartatwawadi/GTRAC
  - Run length encoded multi-sample BWT + server: https://github.com/wtsi-svi/ReadServer
  - Generate artificial human genomes https://gitlab.inria.fr/ml_genetics/public/artificial_genomes
  - Simluate gene-disrupting variants https://github.com/wassermanlab/GeneBreaker

### Footprinting

- Centipede: http://rajanil.github.io/msCentipede/
- SeqGL: predict TF binding from DNase/ATAC-seq https://bitbucket.org/leslielab/seqgl/wiki/Home
- DeFCoM: https://bitbucket.org/bryancquach/defcom

### Gene editing

- Identify on- and off-target edits https://github.com/mpiersonsmela/seqverify

### Metagenomics

- Map reads against redundant databases: https://bitbucket.org/genomicepidemiology/kma
- Iterative bloom filter and Yara mapper for rapid updating and mapping to large metagenomic databases: https://gitlab.com/pirovc/dream_yara/
- Spaced seed hashing: https://bitbucket.org/samu661/fish/overview
- SNV calling: https://github.com/tseemann/snippy

### Methylation

- Mapping
  - https://github.com/chhylp123/BitMapperBS
  - Arioc: GPU-accelerated https://github.com/RWilton/Arioc
  - BatMeth2: https://github.com/GuoliangLi-HZAU/BatMeth2
  - Low-input: https://github.com/Dangertrip/LiBis
  - https://github.com/smithlabcode/abismal
  - https://daehwankimlab.github.io/hisat2/hisat3n
  - Handle multi-mapping reads https://github.com/lmylynn/EM-MUL
  - Abismal: mapping with a 2-letter alphabet https://github.com/smithlabcode/abismal
  - Alignment and methylation calling in single step without intermediate files https://github.com/FischerJo/FAME
- Differential methylation
  - https://github.com/tare/LuxGLM/ 
  - Metilene: Calling differential methylation http://www.bioinf.uni-leipzig.de/Software/metilene/
  - Using VAE: https://www.biorxiv.org/content/early/2018/11/07/433763
  - Differential allele-specific methylation: https://github.com/markrobinsonuzh/DAMEfinder
- Association
  - MACAU: Mixed-model association http://www.xzlab.org/software.html
  - GEM: R package for meQTL and EWAS https://bioconductor.org/packages/devel/bioc/html/GEM.html
- Pipelines
  - https://github.com/MarWoes/wg-blimp
  - https://github.com/Christensen-Lab-Dartmouth/MethylNet
  - https://github.com/jlab-code/MethylStar
  - https://github.com/NuttyLogic/BSBolt
  - https://github.com/sarahet/RLM
  - https://github.com/FischerJo/FAME
  - https://github.com/nloyfer/wgbs_tools
- Other
  - Compression tool for bedMethyl files: https://github.com/jianhao2016/METHCOMP
  - Reference-free bisulfite sequence comparison: https://github.com/thomasvangurp/epiGBS
  - Summarization of multiple bedgraph files: https://github.com/CompEpigen/methrix
  - Binary storage format https://github.com/ZhouQiangwei/DMtools
  - https://github.com/epifluidlab/FinaleMe



### Nanopore/PacBio

- Minion_qc: https://github.com/roblanf/minion_qc
- Clinical application: https://www.nature.com/articles/s41467-017-01343-4
- Scrubbing low-quality reads: https://bitbucket.org/berkeleylab/jgi-miniscrub
- Online CNA detection: https://sourceforge.net/projects/nanogladiator/
- SV calling: https://github.com/tjiangHIT/cuteSV
- Removing low-quality read segments: https://bitbucket.org/berkeleylab/jgi-miniscrub/src/master/
- Simulation
  - https://github.com/lykaust15/deep_simulator
  - https://github.com/liyu95/DeepSimulator
  - https://github.com/NWPU-903PR/NPBSS_Octave
  - http://cgm.sjtu.edu.cn/PaSS/
  - https://github.com/yukiteruono/pbsim2
  - NanoSim https://github.com/karel-brinda/NanoSim-H
  - https://github.com/aquaskyline/LRSIM

### RNA

- Reviews
  - https://arxiv.org/abs/2010.02391
- Splice junctions
  - Set of novel (i.e. missing from annotation databases) splice junctions identified from SRA datasets: https://github.com/nellore/intropolis/blob/master/README.md
  - Complete sets of splice junctions in public RNA-seq datasets: http://snaptron.cs.jhu.edu/data/
- QC
  - NOISeq - exploratory analysis of read mappings https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html
  - AuPairWise: determine replicability without replicates https://github.com/sarbal/AuPairWise
  - dupRadar: duplications https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html
  - Correcting for RNA quality: https://github.com/jengelmann/FastqPuri
  - Error correction using DBG: https://github.com/Malfoy/BCT
  - Prefilter reads to discard those that won't align to genes of interest https://github.com/AlgoLab/shark
  - Use hash of reference transcriptome to convert between transcriptome builds https://bioconductor.org/packages/tximeta
  - infer strandedness https://github.com/betsig/how_are_we_stranded_here
- Align/quantify
  - Need to re-evaluate HiSat2 + StringTie pipeline https://github.com/gpertea/stringtie
  - RapMap: stand-alone lightweight alignment library: https://github.com/COMBINE-lab/RapMap/tree/SAQuasiAlignment
  - Kallisto is a fast and accurate method for transcript quantification https://pachterlab.github.io/kallisto/
    - Sleuth is a companion R package for differential expression analysis http://pachterlab.github.io/sleuth/
    - Different models can be used in Sleuth to, for example, perform time-course experiments http://nxn.se/post/134227694720/timecourse-analysis-with-sleuth
  - tximport: R package for aggregating transcript-level quantifications for gene-level analysis: http://f1000research.com/articles/4-1521/v1
  - Wasabi: prepare Salmon/Sailfish output for Sleuth https://github.com/COMBINE-lab/wasabi
  - featureCounts:
    - http://bioinf.wehi.edu.au/featureCounts/
    - featureCounts replacement with better handling of multi-mapping reads: https://bitbucket.org/mzytnicki/multi-mapping-counter/src/master/
    - Faster version of HTSeq/featureCount: https://github.com/qinzhu/VERSE
    - https://github.com/hitaandrea/MGcount
  - D-GEX: Quantification of whole-transcriptome gene expression from landmark genes https://github.com/uci-cbcl/D-GEX
  - Quantification using both structure and abundance information: https://pypi.python.org/pypi/rsq
  - Improve transcript quantification by integrating PolII ChIP-seq data: https://github.com/pliu55/RSEM/tree/pRSEM
  - Hera: simultaneous alignment, quantification, and fusion detection https://github.com/bioturing/hera
  - Aligner calibraiton: https://bitbucket.org/irenerodriguez/fbb
  - Quantification using k-mers unique to each gene: https://github.com/informationsea/Matataki
  - Read assignment for multi-mapping reads:
    - http://gitlabscottgroup.med.usherbrooke.ca/scott-group/coco
    - https://github.com/seqcode/allo
  - TPMCalculator: https://github.com/ncbi/TPMCalculator
  - Convert between BAM and TCC: https://github.com/pachterlab/bam2tcc
  - Long read: https://github.com/hitbc/deSALT
  - XAEM: http://fafner.meb.ki.se/biostatwiki/xaem/
  - Segmentation https://github.com/HCBravoLab/yanagi
  - Group and collapse transcripts for more accurate counts https://github.com/COMBINE-lab/terminus
  - https://github.com/Kingsford-Group/subgraphquant
  - Haplotype-specific expression using pan-transcriptome graphs https://github.com/jonassibbesen/rpvg
- Correction/Normalization
  - Choosing normalization methods: https://arxiv.org/abs/1609.00959
  - TDM: cross-platform normalization https://github.com/greenelab/TDMresults
  - alpine: corrects for fragment sequence bias https://github.com/mikelove/alpine/blob/master/vignettes/alpine.Rmd
  - Filter out lowly-expressed transcripts prior to quantification decreases FP rate: http://www.genomebiology.com/2016/17/1/12
  - Partition variance between biological and technical sources: https://www.bioconductor.org/packages/3.3/bioc/vignettes/variancePartition
  - Quantify and correct for uncertainty in abundance estimates: https://github.com/PSI-Lab/BENTO-Seq
  - Simultaneous isoform discovery and quantification across multiple samples: http://cbio.ensmp.fr/flipflop
  - R package to compare normalization methods: https://github.com/Edert/NVT
  - Filtering and tissue-aware normalization: http://bioconductor.org/packages/release/bioc/html/yarn.html
  - Bias correction for transcript abundance estimation: https://www.lexogen.com/mix-square-scientific-license/
  - Replacement for htseq-counts/featureCounts that handles multi-mapping reads: https://bitbucket.org/mzytnicki/multi-mapping-counter
  - https://github.com/wgmao/DataRemix
  - Junction compatibility score for idenitifying genes whose isoform-level abundance estimates differ significantly from expectation: https://github.com/csoneson/annotation_problem_txabundance
  - Count corrector for embedded and multi-mapped genes: http://gitlabscottgroup.med.usherbrooke.ca/scott-group/coco/tree/master
  - Abundance estimates using likelihood approximation https://github.com/dcjones/polee
  - Merge paired-ends into fragments https://github.com/Shao-Group/rnabridge-align
  - Remove unwanted variation https://github.com/RMolania/TCGA_PanCancer_UnwantedVariation
- Workflows
  - Artemis (RNA-Seq workflow designed around Kallisto): https://github.com/RamsinghLab/artemis
  - Rafalab workflow: https://github.com/ririzarr/rafalib
  - Isolator: https://github.com/dcjones/isolator
  - RNACocktail: https://bioinform.github.io/rnacocktail/
  - QuickRNASeq: http://baohongz.github.io/QuickRNASeq/
  - Encode long-read: https://github.com/ENCODE-DCC/long-rna-seq-pipeline
  - Fusion detection: https://github.com/bcgsc/pavfinder
- eQTL
  - Multi-tissue
    - HT-eQTL https://github.com/reagan0323/MT-eQTL
    - MT-HESS: eQTL analysis across tissues http://www.mrc-bsu.cam.ac.uk/software/
    - eQTLBMA: cross-tissue eQTL https://github.com/timflutre/eqtlbma
    - Multi-tissue eQTL: https://cran.r-project.org/web/packages/JAGUAR/index.html
  - Multi-tissue, polygenic TWAS: https://github.com/ypark/fqtl
  - Quantile regression approach: https://xiaoyusong.shinyapps.io/QRBT/
  - Using prior knowledge: https://github.com/redsofa/LassoMP
  - CONDOR: simultaneous cis- and trans-eQTL analysis https://github.com/jplatig/condor
  - log allelic fold change (aFC) to quanitfy effect sizes of eQTL: http://biorxiv.org/content/biorxiv/early/2016/09/30/078717.full.pdf
  - Identify cis mediators of trans-eQTL: http://biorxiv.org/content/early/2016/09/30/078683
  - https://cran.r-project.org/web/packages/QRank/
  - With or w/o genotypes: https://gitlab.com/evigorito/baseqtl
- Differential expression
  - cjBitSeq: https://github.com/mqbssppe/cjBitSeq/wiki
  - Differential junction usage: https://github.com/hartleys/JunctionSeq
  - TROM: comparison of transcriptomes between species (and maybe between cell/tissue types?) https://cran.r-project.org/web/packages/TROM/index.html
  - Tissue specificity of genes, based on GTEx data: http://genetics.wustl.edu/jdlab/tsea/
  - Informative priors for Bayesian differential expression analysis using historical data: https://github.com/benliemory/IPBT
  - RNA-enrich: http://lrpath.ncibi.org/
  - Diferentially expressed region finder (also works for ChIP-seq peaks): www.bioconductor.org/packages/derfinder
  - Differentially expressed pathways using kernel MMD: https://eib.stat.ub.edu/tiki-index.php?page_ref_id=73
  - Method using dimension-reduced ANOVA: http://homepage.fudan.edu.cn/zhangh/softwares/multiDE/
  - Correct for hidden sources of variation: https://github.com/sutigit21/SVAPLSseq
  - Alternative method for estimating variances: https://github.com/mengyin/vashr
  - With few replicates: https://figshare.com/s/963e895f812d6f06468a
  - Network sub-pathways: http://bioconductor.org/packages/release/bioc/html/DEsubs.html
  - SDEAP: https://github.com/ewyang089/SDEAP/wiki
  - Local subnetworks enriched for DE genes: https://cran.rstudio.com/web/packages/LEANR/index.html
  - DiscriminantCut: https://github.com/beiyuanzhe/DiscriminantCut
  - Stage-wise DE/DT https://github.com/statOmics/stageR
  - Computing heritability of gene expression: https://cran.r-project.org/web/packages/HeritSeq/index.html
  - Power analysis: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1648-2
  - Using a beta binomial model with dynamic correction for overdispersion: https://github.com/GuoshuaiCai/BBDG
  - Kmer-based: https://github.com/Transipedia/dekupl
  - Correction for log fold change bias at the extremes: https://bioconductor.org/packages/release/bioc/html/apeglm.html
- Co-expression/networks
  - GTEx: http://biorxiv.org/content/early/2016/10/02/078741
  - BicMix: differential co-expression networks http://beehive.cs.princeton.edu/software/
- ASE
  - GeniASE: ASE without haplotypes http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758070/pdf/srep21134.pdf
  - Without phasing: http://lifecenter.sgst.cn/cisASE/
  - ASElux: https://drive.google.com/open?id=0B7E7HSjQ-SumQmlPc1Z0aUR5Sk0
  - BLMRM https://github.com/JingXieMIZZOU/BLMRM
  - https://github.com/santiagosnchez/CompMap
  - https://github.com/ncjllld/byase
  - https://github.com/anna-saukkonen/PAC
- Splicing
  - JunctionSeq: https://github.com/hartleys/JunctionSeq
  - SplicER: https://github.com/lkmklsmn/SplicER
  - Annotation of splicing types https://r-forge.r-project.org/projects/splicingtypes/
  - MAJIQ: detection of local splice variation http://majiq.biociphers.org/
  - LeafCutter: https://github.com/davidaknowles/leafcutter
  - JSplice: http://www.mhs.biol.ethz.ch/research/krek/jsplice.html
  - Splice site prediction: http://cabgrid.res.in:8080/HSplice/
  - DRIM-seq: http://bioconductor.org/packages/DRIMSeq
  - Fast quantification of differential splicing: https://github.com/comprna/SUPPA
  - Identify variant associated with splicing: https://sourceforge.net/projects/isvase/
  - Prediction of intronic splice branchpoints: https://github.com/betsig/branchpointer/
  - Proportion spliced index: https://github.com/comprna/Junckey
  - Whippet: https://github.com/timbitz/Whippet.jl
  - SparseIso: http://github.com/henryxushi/SparseIso
  - PennDiff: https://github.com/tigerhu15/PennDiff
  - ASElux: https://github.com/abl0719/ASElux
  - F1000/Bioconductor workflow on DRIMSeq and DEXSeq for DTU: https://f1000research.com/articles/7-952/v1
  - Alternative splicing https://github.com/raphaelleman/SpliceLauncher
  - splice QTL pipeline https://github.com/guigolab/sqtlseeker2-nf
- Assembly
  - CIDANE: http://ccb.jhu.edu/software/cidane/
  - transrate: evaluation of de novo assemblies http://hibberdlab.com/transrate/
  - dammit: annotator for de novo assemblies http://dammit.readthedocs.org/en/latest/
  - kma: detection of differential intron retention https://github.com/pachterlab/kma
  - RapClust: lightweight clustering of de novo transcriptomes https://pypi.python.org/pypi/rapclust/0.1
  - Shannon http://sreeramkannan.github.io/Shannon/
  - Strawberry https://github.com/ruolin/Strawberry
  - CLASS2: http://ccb.jhu.edu/people/florea/research/CLASS2/
  - BinPacker: https://sourceforge.net/projects/transcriptomeassembly/files/
  - Multi-sample transcriptome assembly: http://tacorna.github.io/
  - Clustering to decontaminate de novo assemblies: https://github.com/Lafond-LapalmeJ/MCSC_Decontamination
  - Assembly from unstranded data: http://big.hanyang.ac.kr/CAFE
  - Scallop: https://github.com/Kingsford-Group/scallop
  - Identification of transcript boundaries: https://github.com/realbigws/DeepBound
  - Corset: gene counts from a transcriptome assembly https://github.com/Oshlack/Corset/wiki
  - Consensus method: https://github.com/macmanes-lab/Oyster_River_Protocol
  - Necklace: https://github.com/Oshlack/necklace
  - Strawberry: https://github.com/ruolin/strawberry
  - BIISQ: also does differential isoform expression https://github.com/bee-hive/BIISQ
  - Cluster and annotate contigs from denovo assemblies: https://github.com/COMBINE-lab/grouper
  - https://github.com/studla/RYUTO
  - Bridging of paired-end reads prior to assembly https://github.com/Shao-Group/coral
  - Ensemble: 
    - https://sourceforge.net/projects/transcriptomeassembly/files/TransBorrow/
    - http://bioinfolab.unl.edu/emlab/consemble/
  - Graph: https://github.com/jonassibbesen/vgrna-project-paper
  - MINTIE: https://github.com/Oshlack/MINTIE
  - https://github.com/Gregor-Mendel-Institute/bookend
- Time series
  - DiceSeq: http://diceseq.sourceforge.net/
  - ctsGE: https://bioconductor.org/packages/release/bioc/html/ctsGE.html
  - Dynamic-BM: http://bioinfo.ibp.ac.cn/Dynamic-BM/
  - Differential isoform usage over time http://bioconductor.org/packages/devel/bioc/html/maSigPro.html
- Deconvolution
  - VoCAL: https://cran.r-project.org/web/packages/ComICS/index.html
  - DeconRNASeq: deconvolute expression profiles in mixed tissues http://www.bioconductor.org/packages/2.12/bioc/vignettes/DeconRNASeq/inst/doc/DeconRNASeq.pdf
- Variant calling
  - Post-processing to improve germline calling: https://github.com/inab/SmartRNASeqCaller
  - https://github.com/pmglab/RNA-SSNV
  - Structural variation
    - Squid: https://github.com/Kingsford-Group/squid
    - Identify gene expression driven by copy number alteration in samples with matched RNA-seq and CNA data: https://www.bioconductor.org/packages/release/bioc/html/iGC.html
    - Fusions
      - FusionCatcher: https://github.com/ndaniel/fusioncatcher
      - STAR: http://star-fusion.github.io
      - AF4: Alignment-free fusion detection https://github.com/grailbio/bio/tree/master/fusion
      - Using kallisto: https://github.com/pmelsted/pizzly
      - Using RapMap https://github.com/nghiavtr/FuSeq
      - https://github.com/suhrig/arriba
- Long read
  - [bambu](https://github.com/GoekeLab/bambu): reference-guided transcript discovery and quantification
  - https://github.com/vladimirsouza/lrRNAseqVariantCalling
- Simulation
  - https://omictools.com/data-simulation3-category
  - https://github.com/jaudoux/simct/
  - https://github.com/alyssafrazee/polyester
  - https://github.com/shaileshtripathi/sgnesR
  - https://github.com/benjaminjack/pinetree
  - https://github.com/CenterForStatistics-UGent/SPsimSeq (bulk and single-cell)
  - https://www.biorxiv.org/content/10.1101/836254v1
  - https://cran.r-project.org/web/packages/seqgendiff/index.html
  - https://github.com/flahertylab/scsim
  - https://github.com/CenterForStatistics-UGent/SPsimSeq
  - https://github.com/flahertylab/scsim (joint bulk and single-cell)
  - https://github.com/SindiLab/ACTIVA (single cell)
  - https://github.com/biomedbigdata/ASimulatoR
  - https://github.com/omnideconv/SimBu
- Other
  - Biclustering for gene co-expression analysis: http://bioconductor.org/packages/devel/bioc/html/QUBIC.html
  - Sample size calculation for experimental design: https://cran.r-project.org/web/packages/ssizeRNA/index.html
  - Variance estimation: http://github.com/mengyin/vashr
  - Sample expression "admixture" (can also be used for deconvolution): https://www.bioconductor.org/packages/release/bioc/html/CountClust.html
    - Essentially, this assigns samples to clusters based on similarity in expression of sets of genes determined be be most discriminating. Each sample can belong to multiple clusters (similar to an admixture analysis).
  - Identification of regulatory networks: https://sites.google.com/a/fleming.gr/rnea/
  - Predict ribosome footprint from transcripts https://sourceforge.net/projects/riboshape/
  - Phasing: https://github.com/secastel/phaser
  - Identifying the source of (almost) all RNA-seq reads: https://github.com/smangul1/rop/wiki
  - Subsampling to determine effect of read depth on downstream analyses: http://www.bio-complexity.com/samExploreR_1.0.0.tar.gz
  - Phenotype prediction: https://github.com/clabuzze/Phenotype-Prediction-Pipeline
  - Predict RNA-RNA interaction: https://github.com/satoken/ractip
  - Mitigate cell-cycle effects: http://www.nature.com/articles/srep33892
  - Fast computation of probabilities of pairwise regulation https://github.com/lingfeiwang/findr
  - Align against synthetic transcript-based reference: https://github.com/Oshlack/Lace
  - Interactive visualization http://bioconductor.org/packages/devel/bioc/html/Glimma.html
  - New factorization method for dimensionality reduction: https://github.com/brian-cleary/CS-SMAF
  - Database of public RNA-seq data sets processed using the same pipeline: https://github.com/mskcc/RNAseqDB
  - Evaluation of aligners on long reads (GMap performs best): http://biorxiv.org/content/early/2017/04/11/126656
  - Test different ML algorithms for classifying expression profiles: https://github.com/gboris/blkbox
  - Find regions of correlated expression: https://cran.r-project.org/web/packages/SegCorr/index.html
  - Tools for barcoded RNA-Seq: https://github.com/lhhunghimself/LINCS_RNAseq_cpp
    - Includes umimerge_filter, which converts barcoded RNA-Seq alignments into a 8-byte hashes (barcode + position)
  - Predict CNA from gene expression: http://wang-lab.ust.hk/software/Software.html
  - Predict age: https://github.com/LBMC/RAPToR
  - Convert a BAM from genomic to transcriptomic coordinates https://github.com/OceanGenomics/mudskipper
  - Estimate coverage: 

### Single-cell

- Comparative analysis of methods:
  - http://www.sciencedirect.com/science/article/pii/S1097276517300497
  - https://www.nature.com/nmeth/journal/v14/n4/full/nmeth.4220.html
  - Benchmarking: https://www.biorxiv.org/content/early/2018/10/08/433102
- Review of experimental design and analysis: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y
- Sean Davis' list: https://github.com/seandavi/awesome-single-cell
- Binary format and tools for storing single-cell data: https://github.com/BUStools
- Loom: http://linnarssonlab.org/loompy/index.html
- Datasets
  - scRNAseqDB: https://bioinfo.uth.edu/scrnaseqdb/
  - Analysis-ready datasets: http://imlspenticton.uzh.ch:3838/conquer/
- Platforms
  - Microwells: http://www.nature.com/articles/srep33883
- Simulation
  - SymSim: https://github.com/YosefLab/SymSim
  - Using GANs: https://github.com/imsb-uke/scGAN
  - Splatter
  - BASiCS
  - Minnow: https://github.com/COMBINE-lab/minnow/
  - https://github.com/PayamDiba/SERGIO
  - Ideas for how to simulate single-cell data http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004983
  - https://github.com/qasimyu/scssim
  - https://www.biorxiv.org/content/10.1101/2020.05.27.119594v2
- QC
  - http://www.morgridge.net/SinQC.html
  - https://github.com/YosefLab/scone
  - Parallel transcriptome/epigenome data: https://github.com/ChengchenZhao/DrSeq2
  - https://bitbucket.org/princessmaximacenter/Sharq
  - Predict coverage increase from deeper sequencing: https://github.com/smithlabcode/preseq
  - https://github.com/ayshwaryas/ddqc
- Normalization
  - ks test: http://michelebusby.tumblr.com/post/130202229486/the-ks-test-looks-pretty-good-for-single-cell
  - Accounting for technical variation: http://www.nature.com/ncomms/2015/151022/ncomms9687/full/ncomms9687.html#supplementary-information
  - ZIFA: Zero-inflated factor analysis https://github.com/epierson9/ZIFA
  - NODES
    - http://biorxiv.org/content/early/2016/04/22/049734.full.pdf+html
    - https://www.dropbox.com/s/pno78mmlj0exv7s/NODES_0.0.0.9010.tar.gz?dl=0
  - scran: http://bioconductor.org/packages/devel/bioc/html/scran.html
  - Correct for expression heterogeneity: https://github.com/PMBio/scLVM
  - Comparison of normalization methods: http://biorxiv.org/content/biorxiv/early/2016/07/17/064329.full.pdf
  - SCnorm: https://github.com/rhondabacher/SCnorm/tree/master/R
  - Factor analysis: https://github.com/UcarLab/IA-SVA/
  - qSVA corrects for RNA quality (in the SVA package)
  - ERCC https://bitbucket.org/bsblabludwig/bearscc
  - Weighting method that enables bulk RNA-seq tools to be used with single-cell: https://www.ncbi.nlm.nih.gov/pubmed/29478411
  - Benchmarking https://github.com/cellgeni/batchbench
- Cell type-specific expression
  - CellMapper: http://bioconductor.org/packages/release/bioc/html/CellMapper.html
- Clustering
  - Comparative analysis:
    - http://biorxiv.org/content/early/2016/04/07/047613
    - https://github.com/epurdom/clusterExperiment
  - SC3: consensus clustering https://github.com/hemberg-lab/sc3
  - destiny: diffusion maps for single-cell data http://bioconductor.org/packages/release/bioc/html/destiny.html
  - https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts
  - GiniClust https://github.com/lanjiangboston/GiniClust
  - pcaReduce: https://github.com/JustinaZ/pcaReduce
  - SIMLR: https://github.com/BatzoglouLabSU/SIMLR
  - CIDR: https://github.com/VCCRI/CIDR
  - Vortex: http://web.stanford.edu/~samusik/vortex/
  - Identify rare cell types: RaceID http://www.nature.com/nature/journal/v525/n7568/full/nature14966.html
  - CIDR: https://github.com/VCCRI/CIDR
  - ZIMBwave: https://github.com/drisso/zinbwave
  - Packages from the Chen lab: http://www.pitt.edu/~wec47/singlecell.html
  - Neural networks for dimensionality reduction and clustering http://sb.cs.cmu.edu/scnn/
  - DCSS: https://github.com/srmcc/dcss_single_cell
  - Cell similarity measure: https://github.com/maggiecrow/MetaNeighbor
  - Dimensionality reduction introduces distortions that can be addressed by high-dimensional PCA https://www.biorxiv.org/content/10.1101/689851v1
  - https://github.com/dnbaker/minicore
- Differential Expression
  - Monocle cole-trapnell-lab.github.io/monocle-release/ (2.0 has Census algorithm for differential transcript analysis)
  - scDD: https://github.com/kdkorthauer/scDD
  - ISOP: comparison of isoform pairs in single cells https://github.com/nghiavtr/ISOP
  - D3E: http://hemberg-lab.github.io/D3E/
  - BASiCS: https://github.com/catavallejos/BASiCS
  - Beta Poisson: https://github.com/nghiavtr/BPSC
  - Zero-inflation correct enables use of DESeq2, etc w/ single cell data: https://github.com/statOmics/zingeR
  - DESingle: https://github.com/miaozhun/DEsingle
  - VAMF: https://github.com/willtownes/vamf-paper
  - Post-clustering: https://github.com/jessemzhang/tn_test
- Transcriptome assembly
  - https://github.com/bcgsc/RNA-Bloom
- Sketching: https://github.com/bendemeo/hopper
- Splicing
  - Brie: https://github.com/huangyh09/brie    
- Time-series/ordering/lineage prediction
  - Monocle: Analysis of pseudotime uncertainty: http://cole-trapnell-lab.github.io/monocle-release/
  - ECLAIR: cell lineage prediction https://github.com/GGiecold/ECLAIR
  - Identification of ordering effects: https://github.com/lengning/OEFinder
  - Slicer: non-linear trajectories https://github.com/jw156605/SLICER
  - Wishbone: identification of bifurcations in developmental trajectories http://www.c2b2.columbia.edu/danapeerlab/html/cyt-download.html
  - SCOUP: https://github.com/hmatsu1226/SCOUP
  - Ouija: https://github.com/kieranrcampbell/ouija
  - SinCell: http://bioconductor.org/packages/release/bioc/html/sincell.html
  - SlingShot: https://github.com/kstreet13/slingshot
  - Cytoscape plugin: http://cytospade.org/
  - TSCAN: https://github.com/zji90/TSCAN
  - Scimitar: https://github.com/dimenwarper/scimitar
  - TimeSeq: https://cran.r-project.org/web/packages/timeSeq/index.html
  - cellTree: http://bioconductor.org/packages/release/bioc/html/cellTree.html
  - Diffusion pseudiotime: http://www.helmholtz-muenchen.de/icb/research/groups/machine-learning/projects/dpt/index.html
  - KBranches: https://github.com/theislab/kbranches
  - Construction co-expression networks: https://cran.r-project.org/web/packages/LEAP/index.html
  - Differentila expression between trajectories https://github.com/kieranrcampbell/switchde
  - FORKS: https://github.com/macsharma/FORKS
  - PhenoPath: https://github.com/kieranrcampbell/phenopath
- Pipelines
  - Seurat http://www.satijalab.org/seurat.html
  - SINCERA https://research.cchmc.org/pbge/sincera.html
  - MAST: https://github.com/RGLab/MAST
  - scde (differential expression + gene set over-dispersion): https://github.com/hms-dbmi/scde
  - BaSiCs: Bayesian analysis of single cell data: https://github.com/catavallejos/BASiCS
  - FastProject: https://github.com/YosefLab/FastProject/wiki
  - Citrus: http://chenmengjie.github.io/Citrus/
  - Tools from Teichmann lab (cellity, celloline, scrnatb): https://github.com/Teichlab/
  - SCell: https://github.com/diazlab/SCell
  - Scater: http://bioconductor.org/packages/scater
  - HocusPocus: https://github.com/joeburns06/hocuspocus
  - Granatum: https://gitlab.com/uhcclxgg/granatum
  - scPite: https://github.com/LuyiTian/scPipe
  - Scrat: For epigenetic data: https://zhiji.shinyapps.io/scrat/
  - scanpy: https://github.com/theislab/scanpy
  - pegasus: https://github.com/klarman-cell-observatory/pegasus
- SNVs/CNVs
  - DNA SNV calling: https://bitbucket.org/hamimzafar/monovar
  - Ginko: analysis of CNVs in single-cell data: http://qb.cshl.edu/ginkgo/?q=/XWxZEerqqY477b9i4V8F
  - CNV calling: http://genome.cshlp.org/content/early/2016/01/15/gr.198937.115.full.pdf
  - Genotyping: https://bitbucket.org/aroth85/scg/wiki/Home
  - Variant calling from RNA-seq: https://github.com/adamcornish/red_panda
  - https://github.com/HorvathLab/NGS/tree/master/SCReadCounts
  - https://github.com/KChen-lab/Monopogen
- Regulatory networks
  - Gene co-expression: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004892
  - SCODE: https://github.com/hmatsu1226/SCODE
  - SCENIC: https://github.com/aertslab/SCENIC
- Integration:
  - https://github.com/YosefLab/scVI
  - https://github.com/brianhie/scanorama
  - https://github.com/bioinfoDZ/RISC/
- Other scRNA-seq
  - Analysis of 3' tagging data: https://github.com/garber-lab/ESAT
  - StemID: Prediction of stem cells and lineage information https://github.com/dgrun/StemID
  - Phasing: https://github.com/edsgard/scphaser
  - Classification using sets of known cell type-specific genes: https://github.com/YosefLab/FastProject/wiki
  - UMI counting: https://github.com/vals/umis
  - Imputation of missing values
    - Magic: https://github.com/pkathail/magic/
    - scImpute: https://github.com/Vivianstats/scImpute
    - DrImpute: https://cran.r-project.org/web/packages/DrImpute/index.html
  - Power analysis: https://github.com/vals/umis/
  - Pooled perturbation experiments: https://github.com/asncd/MIMOSCA
  - Simulation: http://bioconductor.org/packages/splatter/
  - Comparison across experiments: https://github.com/hemberg-lab/scmap
  - Demuxlet: using natural genetic variation to demultiple 10x data https://github.com/hyunminkang/apigenome
  - Alignment of multiple single-cell data sets: https://github.com/jw156605/MATCHER
  - Topological analysis: https://github.com/RabadanLab/scTDA
  - Deconvolution of bulk RNA from scRNA: https://github.com/xuranw/MuSiC
  - Minimum metadata for reporting scRNA-seq experiments: https://arxiv.org/abs/1910.14623
- DNA-seq
  - haplotype inference: https://github.com/raphael-group/chisel
- ATAC-seq
  - Infer TF variation: https://github.com/GreenleafLab/chromVAR
  - Clustering
    - https://github.com/timydaley/scABC
    - https://github.com/QuKunLab/APEC
  - CNV https://github.com/spcdot/CopyscAT

### Somatic

- QC
  - Bias detection
    - https://www.degruyter.com/downloadpdf/j/jib.2017.14.issue-3/jib-2017-0025/jib-2017-0025.pdf
    - This method is targeted to RNA-Seq, but should work for CNV calling from Exome data as well
  - Deconvolution of signal (allele frequency, gene expression, etc) from heterogeneous tissue data: https://github.com/tedroman/WSCUnmix
  - Remove strand bias artifacts (e.g. FFPE samples): https://github.com/mikdio/SOBDetector
  - Remove DNA contamination using microhaplotypes https://www.biorxiv.org/content/10.1101/2020.12.18.423488v1
- Sentieon includes step for "co-realignment" which just seems to be merging tumor+normal and doing indel realignment on the merged sample to harmonize indel alignment for better comparability after variant calling.
- Variant Calling
  - Review of somatic variant callers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5852328/
  - FMI abstracts on variant calling from tumor-only sequencing
    - http://cancerres.aacrjournals.org/content/74/19_Supplement/1893
    - http://ascopubs.org/doi/abs/10.1200/jco.2015.33.15_suppl.11084
  - Sarek: https://github.com/SciLifeLab/Sarek
  - Comparison of 9 somatic variant callers
    - http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0151664
    - From 2016; should be updated with MuTect2 and Strelka2
  - Lancet: new somatic variant caller based on genome graphs
    - https://github.com/nygenome/lancet
    - Supports linked read data
  - Some recommended tools/pipelines for calling low-frequency variants:
    - https://github.com/tamilieberman/TB-diversity-across-organs
    - LoFreq: http://csb5.github.io/lofreq/
    - ASAP: https://github.com/TGenNorth/ASAP
  - Goby uses deep learning for somatic variant calling: http://campagnelab.org/software/goby/
  - MC3/PanCan
    - https://www.synapse.org/#!Synapse:syn7214402/wiki/406010
    - Variant call merger: https://github.com/covingto/pancanmafmerge
  - Ensemble method with random forest integration of results from multiple callers:
    - https://github.com/skandlab/SMuRF
    - Paper: https://www.biorxiv.org/content/early/2018/02/23/270413
  - http://bioinform.github.io/somaticseq/
  - TNScope: Sentieon's pipeline
    - https://www.biorxiv.org/content/early/2018/01/19/250647
    - Includes machine learning-based variant filtering
  - xATLAS
    - https://github.com/jfarek/xatlas
    - Uses a logistic regression model based on truth sets (e.g. GIAB) to assign quality scores
  - Consensus calling from multiple callers: 
    - https://www.itb.cnr.it/web/bioinformatics/isma
    - https://github.com/mingyi-wang/somatic-combiner
  - Multi-sample caller based on mpileup: https://github.com/IARCbioinfo/needlestack
  - Deep learning
    - https://github.com/jingmeng-bioinformatics/DeepSSV
    - NeuSomatic (non-commercial license)
    - Training NeuSomatic and Octopus models from SEQC2 data https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02592-9
    - https://github.com/skandlab/VarNet
    - https://github.com/JakobSkouPedersenLab/dreams
    - Validation: https://github.com/marc-vaisband/deepCNNvalid
  - Optimized for FFPE: https://www.ciscall.org/en/ciscall7.html
  - Optimized for rare variants in ultra-deep sequencing http://github.com/cibiobcg/abemus
  - Learn mutation profile to refine varinat calls: https://github.com/bmannakee/batcaver
  - Reference-free: https://github.com/izaak-coleman/GeDi
  - FFPE samples https://github.com/adamallo/ITHE
  - https://github.com/genetronhealth/uvc
  - Validating somatic variants in tumor-normal sequencing using differences in k-mer profiles between tumor and normal samples https://github.com/compbio/kmerVC
  - https://github.com/genetronhealth/uvc
  - https://github.com/LeiHaoa/RabbitVar
  - Tumor-only
    - https://github.com/pllittle/UNMASC
    - https://github.com/heiniglab/DeepSom
    - https://github.com/HKU-BAL/ClairS-TO
  - GPU-accelerated: https://github.com/wwylab/MuSE
- CNA/SV
  - Infer the tumor cell fraction of SV: https://github.com/mcmero/SVclone
  - A study of using shallow WGS on low-quality FFPE samples for CNV calling used QDNASeq for segmentation
    - https://bioconductor.org/packages/release/bioc/html/QDNAseq.html
    - Adjusts for CG content and mappability
    - Found 7M reads was optimal
    - https://www.biorxiv.org/content/early/2017/12/08/231480
  - Germline CNVs may be predictive of cancer susceptibility: https://www.biorxiv.org/content/early/2018/04/17/303339
  - Allele-specific SV calling: https://github.com/ma-compbio/Weaver
    - Phasing SVs with noisy data: http://sci-hub.tw/http://liebertpub.com/doi/pdf/10.1089/cmb.2018.0022
  - SEG: https://github.com/ZhaoS-Lab/SE
  - Multi-sample (from same patient)
    - HATCHet
    - CNValidator
  - Estimate copy number and purity: https://github.com/keyuan/ccube
  - Allele-specific CNA: https://github.com/wheelerb/hsegHMM
  - DEFOR: https://github.com/drzh/defor
  - Joint tumor-normal detection: https://github.com/yongzhuang/TumorCNV
  - Normalization: https://github.com/baudisgroup/mecan4cna
  - Allele-specific CNV calling for multi-sample tumor/normal https://github.com/imgag/ClinCNV
  - Background correction prior to CNV calling: https://github.com/mskilab/dryclean
  - https://github.com/ExpressionAnalysis/CNV_Radar
  - QC: https://caravagnalab.github.io/CNAqc/
  - Long read: https://github.com/friend1ws/nanomonsv
  - Low-depth https://github.com/crukci-bioinformatics/rascal
  - Integrate WGS and RNA-seq for gene fusion detection: https://github.com/princessmaximacenter/fusion-sq
  - Joint inference of copy number and rearrangement state in cancer whole genome sequence data https://github.com/mskilab-org/nf-jabba
  - https://github.com/songbowang125/SVision-pro
- Variant Filtering
  - TNER: https://github.com/ctDNA/TNER
  - This is a cool approach for constructing a ChIP-Seq control from integration of multiple public data sets: https://www.biorxiv.org/content/early/2018/03/08/278762. It strikes me that a similar approach could be used to generate matched controls for tumor-only samples.
  - Features used to design a RF classifier: https://www.biorxiv.org/content/biorxiv/early/2019/06/13/670687.full.pdf
  - SGZ method used by FoundationOne for filtering germline variants from tumor-only variant calls https://github.com/jsunfmi/SGZ
  - Likelihood model for high-depth sequencing https://github.com/RainyEricYe/LFMD
  - LOHGIC germline classifier and LOH prediction: https://github.com/KhiabanianLab/Genomics_Oncology_Platform/blob/master/Platform.py
  - https://github.com/cpwardell/FiNGS
  - SV simulation/filtering for FFPE data https://bioconductor.org/packages/release/bioc/html/SimFFPE.html
  - CNN: https://github.com/marc-vaisband/deepCNNvalid
  - Tumor-only SV filter https://github.com/shu-z/germline_svm
- Tumor Purity
  - IchorCNA
    - Estimates tumor cell fraction from low-pass WGS; should probably be adaptable to deep targeted sequencing
    - https://github.com/broadinstitute/ichorCNA
  - Accurity (tumor-normal): https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty043/4827681
  - Control for tumor purity differences in methylation analysis: https://www.biorxiv.org/content/early/2018/01/16/248781
  - PyLOH: https://github.com/uci-cbcl/PyLOH
  - AllFit: https://github.com/KhiabanianLab/All-FIT
  - TPES: using SNVs https://bitbucket.org/l0ka/tpes/src/master/
  - From methylation:
    - https://github.com/xjtu-omics/MEpurity
    - https://github.com/mwsill/RFpurify
  - DeepPurity https://www.biorxiv.org/content/10.1101/805135v1.full.pdf (no open-source code)
- Joint calling
  - CN, tumor purity, and LOH https://bioconductor.org/packages/release/bioc/html/PureCN.html
  - purity, ploidy, and allele-specific CN https://github.com/Crick-CancerGenomics/ascat
  - purity, ploidy, SV, and CN in tumor-normal: https://github.com/hartwigmedical/gridss-purple-linx
- Tumor Muational Burden
  - http://moat.gersteinlab.org/
  - ecTMB: https://github.com/bioinform/ecTMB
  - https://acc-bioinfo.github.io/TMBleR/
  - Normalizing mutation count by read depth https://github.com/OmnesRes/depth_norm
- Variant Annotation
  - CHASM and SNV-Box
    - CHASM Predicts functional significance of somatic missense mutations
    - SNV-Box is a database of pre-computed features of all possible amino acid substitutions in the human genome
    - http://wiki.chasmsoftware.org/index.php/Main_Page
    -  Not free for commercial use, but it was developed at JHU and Vogelstein is a co-author on at least some of the publications
  - Orchid: framework for annotation and machine learning of cancer variants: https://github.com/wittelab/orchid
  - Discover potential cancer driver elements
    - DriverPower: https://github.com/smshuai/DriverPower
    - NetSig: http://www.lagelab.org/projects/
    - CanDrA: http://bioinformatics.mdanderson.org/main/CanDrA#CanDrA
  - Predict cancer subtype using biological networks based on somatic variants: https://www.biorxiv.org/content/early/2017/12/03/228031
  - PCGR: Annotate cancer genomes with clinical information:
      - https://github.com/sigven/pcgr
    - Integrates many public databases
  - Benchmark of methods to identify pathogenic non-coding variation: https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty008/4798701?redirectedFrom=fulltext
  - FASMIC: http://bioinformatics.mdanderson.org/main/FASMIC
  - Database of TCGA alternative splicing events https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30306-4#secsectitle0080
  - Consensus driver gene identification pipeline and database of cancer driver genes in TCGA: https://www.cell.com/cell/fulltext/S0092-8674(18)30237-X
  - TCGA pathogenic signalling pathways: https://www.cell.com/cell/fulltext/S0092-8674(18)30359-3?cid=tw%26p
- Predict cancer type/signature
  - From mutations: https://www.nature.com/articles/nature12477
  - From SNV+SV calls: https://www.biorxiv.org/content/early/2018/02/18/267500
  - Machine learning method to predict cancer type from variant calls: https://www.biorxiv.org/content/early/2017/11/05/214494
  - Prediction of inherited susceptibility to 20 difference cancer types: http://www.pnas.org/content/115/6/1322.short
  - Computing polygenic risk scores and association with cancers: https://www.biorxiv.org/content/early/2017/10/19/205021
  - Predict driver mutations: https://github.com/KarchinLab/CHASMplus
  - https://bioconductor.org/packages/release/bioc/html/SparseSignatures.html
  - https://www.bioconductor.org/packages/release/bioc/html/musicatk.html
  - https://github.com/ma-compbio/MutSpace
  - https://github.com/UMCUGenetics/CUPLR
- Clonality
  - QuantumClone: https://cran.r-project.org/web/packages/QuantumClone/index.html
  - MOSEM: https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0128-6/MediaObjects/41588_2018_128_MOESM1_ESM.pdf
  - Benchmarking, plus description of CloneFinder: https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty469/5040314
  - Tracking variants across longitudinal samples: https://github.com/ChristofferFlensburg/SuperFreq
  - Reconstruction of clone- and haplotype-specific karyotypes https://github.com/aganezov/RCK
- MSI
  - Sonics: https://github.com/kzkedzierska/sonics
  - MSIsensor-pro https://github.com/xjtu-omics/msisensor-pro
  - https://github.com/YixuanWang1120/ELMSI (requires tumor-normal)
  - https://github.com/GeneCast/btMSI-CAST
  - https://github.com/mskcc/mimsi
  - https://github.com/861934367/MSIFinder
  - https://github.com/RADar-AZDelta/DeltaMSI
- Neoantigen 
  - pVACtools: https://github.com/griffithlab/pVACtools
  - NeoPredPipe: https://github.com/MathOnco/NeoPredPipe
  - https://github.com/neoanthill/neoANT-HILL
  - Using phased somatic mutations https://github.com/pdxgx/neoepiscope
  - https://github.com/koesterlab/microphaser
  - https://github.com/icbi-lab/nextNEOpi
  - https://gitlab.com/landscape-of-effective-neoantigens-software
- HRD
  - Genomic Scar Analysis - call copy number and detect HRD from targeted capture sequencing https://github.com/shaominghui/GSA
  - https://github.com/maxwell-lab/HRDex
- ctDNA/cfDNA
  - TNER: background error reduction https://github.com/ctDNA/TNER
  - SV calling https://github.com/liuhc8/Aperture
  - TruNeo pipeline https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03869-9
  - Predicting tumor content https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05478-8
  - Variant calling https://github.com/GenetronBioinfomatics/BayVarC
  - Fragment size: https://github.com/epifluidlab/FinaleToolkit
- Simulation
  - Mutation signatures: https://github.com/HillLab/SomaticSiMu
  - Simulation of cancer genomes: https://github.com/rsemeraro/XomeBlender
  - Database of simulated tumor genomes (created with BamSurgeon on NA12878 background): https://www.biorxiv.org/content/early/2018/02/07/261503
- Other
  - Bioconductor package with various methods for working with TCGA data: https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html
  - maftools: https://bioconductor.org/packages/release/bioc/vignettes/maftools/inst/doc/maftools.html
  - Sarek: end-to-end variant calling pipeline in NextFlow https://github.com/SciLifeLab/Sarek
  - Multi-regional sequencing: Mutect2 in multi-sample mode with correction step performs best https://www.biorxiv.org/content/10.1101/655605v1
  - Subclonality-aware association: https://github.com/Sun-lab/SMASH
  - Regional mutation density classifies cancer type: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006953
  - Hotspots https://github.gersteinlab.org/nimbus
  - Phase germline variants using tumor copy number https://github.com/bowentan/CNAHap
  - Somatic variant calling from singe-cell (transcriptomic, ATAC-seq) https://github.com/cortes-ciriano-lab/SComatic

### Multi-omics/integration

- Reviews:
  - http://www.biomedcentral.com/1752-0509/8/S2/I1
  - Comparison of methods http://www.biomedcentral.com/1752-0509/8/S2/S4
  - Focusing on integration of RNA-seq and ChIP-seq http://journal.frontiersin.org/article/10.3389/fcell.2014.00051/full
  - Network-based methods: http://rsif.royalsocietypublishing.org/content/12/112/20150571
- General-purpose multi-omics integration:
  - mixOmics: R package that implements several multivariate methods, including DIABLO http://mixomics.org/
  - Non-negative matrix factorization for integration of data sets https://github.com/yangzi4/iNMF
  - omicade4: Integration of multi-omics data using co-inertia analysis https://bioconductor.org/packages/release/bioc/html/omicade4.html
  - Integration of multi-omics data using random forests: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1043-4
  - Kernel-PCA: http://www.biomedcentral.com/1752-0509/8/S2/S6
  - Integrative analysis of multiple diverse omics datasets by sparse group multitask regression http://journal.frontiersin.org/article/10.3389/fcell.2014.00062/abstract
  - Joint bi-clustering of multiple data types: http://research.cs.aalto.fi/pml/software/GFAsparse/
  - Identifying covariance between sequencing data sets: http://github.com/pmb59/fCCAC/
  - PyPanda: https://github.com/davidvi/pypanda
  - SDA: integrate gene expression across multiple tissues, or multi-omics in a single tissue, for identification of trans-QTL networks: https://jmarchini.org/sda/
  - HMM for binary classification based on multivariate data: https://github.com/PetarV-/muxstep
  - OmicsIntegrator: https://github.com/fraenkel-lab/OmicsIntegrator
  - Correlation between enriched regions from different data sets: http://malone.bioquant.uni-heidelberg.de/software/mcore
  - R.Jive: https://cran.r-project.org/web/packages/r.jive/
  - SIFORM: http://bioinformatics.oxfordjournals.org/content/early/2016/07/03/bioinformatics.btw295.full
  - EpiMine: https://sourceforge.net/projects/epimine/
  - Deep learning-based framework for integrating multiple data types to predict another data type: https://github.com/ueser/FIDDLE
  - Significance-based https://www.bioconductor.org/packages/release/bioc/html/SMITE.html
  - Two-stage CCA identifies non-linear associations: https://github.com/kosyoshida/TSKCCA
  - HMM: https://link.springer.com/protocol/10.1007/978-1-4939-6753-7_10
  - PGenmi: https://github.com/KnowEnG/pgenmi
  - MOFA: https://github.com/bioFAM/MOFA
  - TensorDecomp: https://github.com/jinghan1018/tensor_decomp
  - Unsupervised: https://github.com/bioFAM/MOFA
  - Avocado: https://noble.gs.washington.edu/proj/avocado/
  - Miodi https://gitlab.com/algoromics/miodin
- Multi-tissue:
  - HDTD: http://bioconductor.org/packages/release/bioc/html/HDTD.html
  - IDEAS: https://github.com/yuzhang123/IDEAS
- Specific data types:
  - Predict gene fusions from WGS and RNA-seq: http://sourceforge.net/p/integrate-fusion/wiki/Home/
  - NuChart: layer additional omics data on Hi-C ftp://fileserver.itb.cnr.it/nuchart/
  - Predict expression from H3K27, and identify cis-regulatory elements: http://cistrome.org/MARGE/
  - GenoSkyline: predict tissue-specific functional regions from epigenomic data http://genocanyon.med.yale.edu/GenoSkyline
  - Integrate DNA and RNA: https://github.com/griffithlab/regtools
  - Association of gene expression with SVs https://github.com/chadcreighton/SVExpress
  - https://github.com/jsalignon/cactus
- Network-based:
  - Merging networks: https://github.com/maxconway/SNFtool
  - XMWAS: https://sourceforge.net/projects/xmwas/
- Causality
  - Hybrid BN/CMI approach to constructing GRNs: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005024
- Noise/bias:
  - Bias correction across different assays on the same samples: https://cran.r-project.org/web/packages/MANCIE/
  - Cell cycle heterogeneity appears to be a minor contributor to noise; instead, library size is the largerst PC by far http://www.nature.com/nbt/journal/v34/n6/full/nbt.3498.html
- Patient/disease subtyping
  - http://biorxiv.org/content/biorxiv/early/2016/09/03/073189.full.pdf (Matlab code)
- Imputation of missing data
  - MI-MFA: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1273-5
  - TensorPute: multi-dimensional imputation https://sites.google.com/site/tensortest2/
- Simulation
  - https://bitbucket.org/ConesaLab/mosim/
  - https://omicssimla.sourceforge.io/
  - InterSIM: simulate correlated multi-omics data https://cran.r-project.org/web/packages/InterSIM/index.html
- Other
  - Identify class-descriminative motifs enriched in subclasses of overlapping annotations: https://github.com/seqcode/sequnwinder
  - R class for integration algorithms: https://bioconductor.org/packages/release/bioc/html/MultiDataSet.html




# General Programming Resources

- Protocol buffers
  - Protobuf: fast cross-language/platform serialization of fixed-format messages https://developers.google.com/protocol-buffers/
  - Cap'n Proto: https://capnproto.org/
- IDEs
  - VisualStudio (now free): https://www.visualstudio.com/vs/visual-studio-mac/
- Diff/patch/merge for data tables
  - js/python: http://paulfitz.github.io/daff/
  - R: https://github.com/edwindj/daff
- JSON Diff: http://www.jsondiff.com/
- Google C++ and python libraries: https://github.com/abseil/abseil-py
- Share a terminal as a web page: https://github.com/yudai/gotty
- Spec for CSV files with metadata: http://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html

## C/C++

- kmers
  - Streaming kmer counting https://github.com/bcgsc/ntCard
  - kmer bloom filters: https://github.com/Kingsford-Group/kbf
  - index: https://github.com/Malfoy/BLight
- High-performance concurrent hash table (C++11): https://github.com/efficient/libcuckoo
- BWT that incorporates genetic variants: https://github.com/iqbal-lab/gramtools
- C++ interface to htslib, BWA-MEM, and Fermi (local assembly) (would be useful to build python bindings for this): https://github.com/walaj/SeqLib
- Minimal perfect hash function: fast, large data sets https://github.com/rizkg/BBHash
- rr debugger: record and playback executions http://rr-project.org/
- clip: command line parser https://github.com/muellan/clipp
- klib: generic library with fast implementations of striped Smith-Waterman and other NGS-related algorithms https://github.com/attractivechaos/klib
- Sketch data structures (includes Python bindings): https://github.com/dnbaker/sketch

## R

- Tidy data cheatsheet: http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
- Multiple variable assignment: https://github.com/nteetor/zeallot
- Summary statistics: https://github.com/ropenscilabs/skimr

### Find packages

- Awesome R: https://github.com/qinwf/awesome-R#integrated-development-environment
- Cranberries: http://dirk.eddelbuettel.com/cranberries/
- METACRAN: identify R packages http://www.r-pkg.org/

### Database

- MonetDB - embeddable column-store DB with R integration (MonetDB.R)

### Data Cleaning

- dplyr
  - chunked processing of large files: https://github.com/edwindj/chunked
- Join tables on inexact matching: https://github.com/dgrtwo/fuzzyjoin
- Tidy text: http://juliasilge.com/blog/Life-Changing-Magic/

### Reporting

- MarkDeep for R documentation: http://casual-effects.com/markdeep/

### Misc

- Advanced table formatting in knitr: https://github.com/renkun-ken/formattable
- Access data frames using SQL: sqldf package
- Work with PDF files: https://cran.r-project.org/web/packages/pdftools/index.html
- Language-agnostic data frame format: https://github.com/wesm/feather
- Work with genomic ranges using dplyr-like syntax. 

## Python

- Data structures/formats
  - Chunked, compressed, disk-based arrays: https://github.com/alimanfoo/zarr
  - Tabular data
    - Working with tabular data: http://docs.python-tablib.org/en/latest/
    - Apache Arrow
    - Pandas
    - Vaesx:https://github.com/vaexio/vaex
  - GFA: https://github.com/ggonnella/gfapy/tree/master/gfapy
- API for interacting with databases: https://github.com/kennethreitz/records
- Debugging
  - boltons.debugutils: The entire boltons package has lots of useful stuff, but debugutils is particularly cool - you can add one line of code to enable you to drop into a debugger on signal (e.g. Ctrl-C): https://boltons.readthedocs.io/en/latest/debugutils.html
  - Web-based debugger: https://github.com/alexmojaki/birdseye
- Stats
  - Non-negative matrix factorization: https://github.com/ccshao/nimfa
  - Statsmodels: http://www.statsmodels.org/stable/index.html
  - R formulas in python: https://github.com/pydata/patsy
- pyrasite: code injection into running applications
- Dexy: documentation
- Event loops for asynchronous programming
  - curio
  - gevent
- dill: alternative serialization
- arrow: alternative to datetime
- Disabling python garbage collection speeds up programs: https://engineering.instagram.com/dismissing-python-garbage-collection-at-instagram-4dca40b29172#.ri55nyjdu (only safe when the lifecycle is straight-forward for all objects, an thus reference counting is sufficient for memory management)
- Easily implementing function proxies/wrappers: http://wrapt.readthedocs.io/en/latest/
- Cache system: https://bitbucket.org/zzzeek/dogpile.cache
- Parse TOML (an enhanced config-file spec): https://github.com/uiri/toml
- Web scraping
  - Goose: https://github.com/grangier/python-goose
  - ScraPy: https://scrapy.org/
  - Requestium: https://github.com/tryolabs/requestium
- Visualize python code execution time as a heatmap in a Jupyter notebook: https://github.com/csurfer/pyheatmagic
- Graph genomes
  - High-level API for working with GFA (i.e. graph genomes): https://github.com/ggonnella/gfapy
  - https://github.com/jambler24/GenGraph
- Fast string matching: https://bergvca.github.io/2017/10/14/super-fast-string-matching.html
- Trace system calls in multiprocessing context: https://github.com/pinterest/ptracer
- Extract keywords from text: https://github.com/vi3k6i5/flashtext
- Probability distributions and other related functions: https://pomegranate.readthedocs.io/en/latest/
- CuPy: Numpy-like implementation of GPU-accelerated array operations https://github.com/cupy/cupy
- Command-line interface to download and manage reference genomes: https://github.com/simonvh/genomepy
- API for working with HGVS-formatted variants: https://github.com/biocommons/hgvs
- Chunked, compressed ndarrays http://zarr.readthedocs.io/en/latest/index.html
- Easy colored/styled printing in command line apps https://github.com/UltimateHackers/hue
- Launching a subprocess in a pseudo-terminal (e.g. for accepting passwords) https://github.com/pexpect/ptyprocess
- VCF:
  - ReadVCF: http://alimanfoo.github.io/2017/06/14/read-vcf.html
  - Faster alternative to pyvcf: https://github.com/brentp/cyvcf2
- BAM:
  - Pure python BAM parser https://github.com/betteridiot/bamnostic
- Alternative to pandas for larger-than-memory datasets: https://vaex.readthedocs.io/en/latest/
- Load PLINK data into Dask arrays https://github.com/limix/pandas-plink

## HPC

- Spark:
  - https://spark.apache.org/downloads.html
  - https://amplab-extras.github.io/SparkR-pkg/
- Ibis: http://blog.cloudera.com/blog/2015/07/ibis-on-impala-python-at-scale-for-data-science/
- Flink: https://flink.apache.org/
- Dask: http://dask.pydata.org/en/latest/

## Command Line (OSX/Linux)

- Diff tables: https://github.com/paulfitz/daff
- Miller - work with tables http://johnkerl.org/miller/doc/build.html

## Databases

- CockroachDB: based on Google's distributed database https://github.com/cockroachdb/cockroach

## Reproducibility/Containerization

- Packrat: http://rstudio.github.io/packrat/walkthrough.html
- Docker
  - http://arxiv.org/pdf/1410.0846v1.pdf
  - http://bioboxes.org/available-bioboxes/
  - http://ivory.idyll.org/blog//2015-docker-and-replicating-papers.html
  - GUI for running Docker images locally: https://kitematic.com/
  - Convert Docker images to AWS Lambda container images: https://github.com/awslabs/aws-lambda-container-image-converter
  - Programmatic access to Docker from Python: https://github.com/docker/docker-py
- Jupyter notebooks
  - http://jupyter.org/
  - http://mybinder.org/
  - RISE: presentations from Jupyter notebooks https://github.com/damianavila/RISE
  - Convert to markdown: https://github.com/aaren/notedown
  - Jupyter notebooks on Azure: https://notebooks.azure.com/
  - Diff/merge: https://github.com/jupyter/nbdime
  - Web-based IDE built on Jupyter: https://github.com/jupyterlab/jupyterlab
  - Convert Jupyter notebook to Sphinx documentation: https://matthew-brett.github.io/nb2plots
  - Diagram editor: https://blog.jupyter.org/a-diagram-editor-for-jupyterlab-a254121ff919
- Jupyter alternatives
  - https://nteract.io/
  - http://www.unnotebook.com/
- Bioinformatics software containers
  - BioContainers: https://github.com/BioContainers
- Popper: http://falsifiable.us/

### Building Pipelines

- Luigi https://github.com/spotify/luigi
- Flo http://flo.readthedocs.org/en/latest/index.html
- Nextflow and Nextflow Workbench: https://www.nextflow.io/ and http://campagnelab.org/software/nextflow-workbench/
  - https://github.com/assemblerflow/assemblerflow
- SUSHI: https://github.com/uzh/sushi
- WorkflowR: https://github.com/jdblischak/workflowr
- Data manager for R: https://cran.r-project.org/web/packages/repo/index.html
- dgsh: bash variant that has primitives for parallelizing tasks https://www2.dmst.aueb.gr/dds/sw/dgsh
- Rabix: visual editor for CWL pipelines http://rabix.io/ (deprecated?)
- List of workflow systems: https://github.com/common-workflow-language/common-workflow-language/wiki/Existing-Workflow-systems
- Experimental directory structure: alternative to HDF5 https://github.com/CINPLA/exdir/
- Structured projects that can be run on multiple workflow engines: http://reana.readthedocs.io/
- Invoke: http://docs.pyinvoke.org/en/latest/
- Snakemake
  - Pipelines for several genomic/epigenomic analyses https://github.com/maxplanck-ie/snakepipes
- drake: https://github.com/ropensci/drake



## Other

- Automated system for sequencing core facilities using microservices architecture: https://www.biorxiv.org/content/early/2017/11/06/214858

# Statistics/Machine Learning

- Search for papers: https://www.semanticscholar.org/
- Common probability distributions http://blog.cloudera.com/blog/2015/12/common-probability-distributions-the-data-scientists-crib-sheet/
- How to share data with a statistician: https://github.com/jtleek/datasharing
- Precision-recall curves:
  - https://cran.r-project.org/web/packages/precrec/index.html
  - Alternative: http://link.springer.com/article/10.1007/s10994-009-5119-5
- R package for subsetting data into training/testing/validation sets: https://cran.r-project.org/web/packages/STPGA/index.html
- PMML: interchange format for trained models https://en.wikipedia.org/wiki/Predictive_Model_Markup_Language

## Methods/algorithms

- Lists
  - https://github.com/rushter/MLAlgorithms
- Decision tree methods
  - R randomForest package
  - FuzzyForests are an extension of random forests for classification in which subsets of variables/features are highly correlated https://github.com/OHDSI/FuzzyForest
  - CatBoost: https://github.com/catboost
  - 2D Boosting: https://github.com/kundajelab/boosting2D/
  - New performance metric: https://cran.r-project.org/web/packages/IPMRF/IPMRF.pdf
- Modrian forests https://scikit-garden.github.io/examples/MondrianTreeRegressor/
- Sequence embedding
  - https://github.com/krishnanlab/pecanpy
- Clustering
  - Help selecting biclustering algorithm: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1487-1
  - DBScan https://cran.r-project.org/web/packages/dbscan/dbscan.pdf
  - KODAMA: https://cran.r-project.org/web/packages/KODAMA/index.html
  - t-SNE: alternative to PCA and MDS: http://lvdmaaten.github.io/tsne/
    - http://distill.pub/2016/misread-tsne/
  - UMAP: https://github.com/lmcinnes/umap
  - Compressive k-means: https://arxiv.org/pdf/1610.08738.pdf
  - Sparse convex: https://arxiv.org/pdf/1601.04586.pdf
  - Estimate `k`: https://github.com/RoheLab/gdim
  - https://github.com/sinc-lab/clustermatch
- Multivariate
  - CRAN packages: http://cran.r-project.org/web/views/Multivariate.html
  - Multivariate analysis of covariance (MANCOVA): http://en.wikipedia.org/wiki/MANCOVA
  - Correlation https://www.biorxiv.org/content/biorxiv/early/2018/05/25/330498.full.pdf
- Nonnegative Matrix Factorization: https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
- Tensor factorization: https://cran.r-project.org/package=tensorBF
- Bayesian regression modeling:
  - brms and rstanarm are R packages based on stan
  - JAGS http://jeromyanglim.blogspot.com/2012/04/getting-started-with-jags-rjags-and.html
  - MCMC http://www.stat.umn.edu/geyer/mcmc/library/mcmc/doc/demo.pdf
- Multiple test correction
  - FDR for multi-dimensional pairwise comparisons (e.g. RNA-seq): http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0937-5
  - Local false signal rate: an alternative to FDR that operates on standard error estimates rather than p-values: https://github.com/stephens999/ashr
  - MTC weighted by variant: effect http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3507.html
  - Parallelizable FDR correction: http://bioinformatics.oxfordjournals.org/content/early/2016/02/25/bioinformatics.btw029.short
  - IHW: http://bioconductor.org/packages/release/bioc/html/IHW.html
  - BonEV: Bonferonni FDR correction https://cran.r-project.org/web/packages/BonEV/index.html
  - Myriads: http://myriads.webs.uvigo.es/
- Fast Bayesian alternative to lasso and ElasticNet for feature selection and effect estimation: https://cran.r-project.org/web/packages/EBglmnet/index.html
- Genome-wide generalized addative models: https://master.bioconductor.org/packages/3.3/bioc/html/GenoGAM.html
- Causal inference test: https://cran.r-project.org/web/packages/cit/index.html
- Reed-Sololmon error correction
  - https://github.com/tomerfiliba/reedsolomon
  - https://github.com/brownan/Reed-Solomon
  - https://github.com/catid/wirehair
- Automated generation of ML pipelines: https://github.com/rhiever/tpot/tree/tpot-mdr
- Topological analysis:
  - https://github.com/MLWave/kepler-mapper
  - https://github.com/RabadanLab/sakmapper
  - http://danifold.net/mapper/
  - https://github.com/paultpearson/TDAmapper
- Feature selection
  - Workflows: https://github.com/enriquea/feseR
  - FeatureSelect: https://github.com/LBBSoft/FeatureSelect
- Compares method for outlier detection: https://cran.r-project.org/web/packages/OutliersO3/index.html
- Comparison of inference methods: http://jenniferhill7.wixsite.com/acic-2016/competition
- Block random forest method for multi-omic data: https://github.com/bips-hb/blockForest



## Deep Learning

- Reading
  - Papers
    - https://github.com/songrotek/Deep-Learning-Papers-Reading-Roadmap
  - Books
    - https://hackerlists.com/free-machine-learning-books/
    - http://www.deeplearningbook.org/contents/intro.html
- Libraries
  - List: http://www.teglor.com/b/deep-learning-libraries-language-cm569/
  - Keras: https://github.com/fchollet/keras
  - Chainer: https://www.oreilly.com/learning/complex-neural-networks-made-easy-by-chainer
  - biologicaly-focused neural networks https://github.com/kundajelab/dragonn/tree/master/dragonn
  - analysis of features in deep neural networks https://github.com/kundajelab/deeplift
  - API to add fuzzy logic: https://fuzzy.ai/docs
  - Edward: probabalistic modeling, inference, and criticism; build on TensorFlow https://github.com/blei-lab/edward
  - VectorFlow: specifically designed for sparse data https://github.com/Netflix/vectorflow
  - secml: https://gitlab.com/secml/secml
  - https://github.com/BIMSBbioinfo/janggu
- Architectures:
  - http://www.asimovinstitute.org/neural-network-zoo/
  - Deep residual:
    - http://image-net.org/challenges/talks/ilsvrc2015_deep_residual_learning_kaiminghe.pdf
    - https://arxiv.org/abs/1512.03385
  - Wide residual: https://arxiv.org/abs/1605.07146
  - Time-delay http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=809100&tag=1
     - Semi-supervised classification of nodes in a graph: https://github.com/tkipf/gcn
  - Adversarial: https://arxiv.org/abs/1406.2661
  - Recurrent scalable deep kernels: https://arxiv.org/abs/1610.08936
  - Multiplicative LSTM: https://arxiv.org/abs/1609.07959
  - Graph CNN: https://arxiv.org/abs/1609.02907
  - Encoding variable-length sequences: 
    - https://arxiv.org/abs/1505.01504
    - dna2vec: https://pnpnpn.github.io/dna2vec/
  - Interactive sequence generation from RNNs: https://arxiv.org/abs/1612.04687
  - DNNs with external memory: http://www.nature.com/nature/journal/v538/n7626/full/nature20101.html
  - The Predictron: https://arxiv.org/abs/1612.08810
  - Associative LSTM: https://arxiv.org/abs/1602.03032
  - Conditional variational autoencoders: http://ijdykeman.github.io/ml/2016/12/21/cvae.html
  - Group equivariant CNN: http://jmlr.org/proceedings/papers/v48/cohenc16.pdf
  - QuickNet: https://arxiv.org/abs/1701.02291
  - Collection of pre-trained tensorflow models: https://github.com/tensorflow/models/tree/master/research
  - Generalized framework for genomic signal recognition: https://zenodo.org/record/1117159
  - Using autoencoders to impute missing values: https://github.com/gevaertlab/DAPL
  - Explainable: https://github.com/wassermanlab/ExplaiNN
- Tools
  - Generate synthetic training data: https://hazyresearch.github.io/snorkel/
  - Tensorflow playground: http://playground.tensorflow.org/
  - DeepVis: http://yosinski.com/deepvis
  - https://medium.com/@shivon/the-current-state-of-machine-intelligence-3-0-e4d305da032e#.gw0ywcpkv
  - Visualize hyperparameters: https://github.com/ContextLab/hypertools
  - Convolutional layers for DNA sequence: https://github.com/koonimaru/DeepGMAP
  - Interpret neuron importance: https://arxiv.org/abs/1807.09946
- Non-bio networks that might be applied
  - srez: https://github.com/david-gpu/srez
  - Deep learning with text: https://explosion.ai/blog/deep-learning-formula-nlp

## Web APIs

- AWS ML: https://aws.amazon.com/machine-learning



## Text classification

- Automatic text summarization: https://pypi.python.org/pypi/sumy
- Word2vec: models for predicting missing words https://en.m.wikipedia.org/wiki/Word2vec

# Visualization

- Types of plots: 
  - http://www.informationisbeautifulawards.com/showcase/611-the-graphic-continuum
  - Line graph + heat map: http://www.fastcodesign.com/3052450/an-easy-intuitive-tool-for-making-sense-of-your-data
  - http://datavizproject.com/
- Examples: http://www.visualcomplexity.com/vc/
- VDA Lab: https://bitbucket.org/vda-lab/
- Visualization of GO results: http://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html
- Fluff: publication-quality genomics plots http://fluff.readthedocs.org
- Comparison of different types of heatmaps: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1442-6
- Interactively build simulations and animated interaction diagrams: http://ncase.me/loopy/
- Plotly: Charting library with multiple language bindings https://github.com/plotly/plotly.py
- Open-access visualization research: http://oavis.steveharoz.com/
- Carto Colors: https://carto.com/carto-colors/
- Generate flow charts, etc using markdown-like syntax: https://github.com/knsv/mermaid
- SamPlot: https://github.com/ryanlayer/samplot
- Visual inspection:
  - Viper: https://github.com/MarWoes/viper
  - https://www.biorxiv.org/content/early/2018/02/21/266262
- Venn diagrams: https://github.com/vqf/nVenn
- BBC Visual and Data Journalism guide
  - Cookbook: https://bbc.github.io/rcookbook/
  - bbplot package: https://github.com/bbc/bbplot
- https://genomespy.app/

## Genome browsers

* Graph
  * https://github.com/MoMI-G/MoMI-G
  * https://github.com/docpaa/mumdex/

## Networks

- Gephi: http://gephi.org/

## Phylogenetic trees

- web-based http://microreact.org/showcase/
- JavaScript libraries: https://github.com/tntvis

## R

- Comparison of libraries:
  - http://lisacharlotterost.github.io/2016/05/17/one-chart-code/
  - http://ouzor.github.io/blog/2014/11/21/interactive-visualizations.html
- D3 from R: http://christophergandrud.github.io/networkD3/
- SVG device: https://github.com/hadley/svglite/blob/master/README.md
- Multilayer data plotted on a Hilbert curve: http://www.bioconductor.org/packages/devel/bioc/html/HilbertCurve.html
- Visualize local epigenetic neighborhood of a SNP: http://bioconductor.org/packages/release/bioc/html/SNPhood.html
- CIRCOS plots: https://cggl.horticulture.wisc.edu/software/
- Make multiplanel figures from a combination of plot types: https://cran.r-project.org/web/packages/multipanelfigure/
- KaryotypeR: http://bioconductor.org/packages/release/bioc/html/karyoploteR.html
- PCA+clustering: http://r-forge.r-project.org/projects/thresher/
- clustree: Plot custer trees, for determining number of clusters to use: https://cran.r-project.org/web/packages/clustree/

### ggplot2

- Themes
  - Cowplot - improve default ggplot: http://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
  - ggplot2 theme for publication-quality figures: https://github.com/robertwilson190/ggplot2-theme
  - HBR Themes: https://github.com/hrbrmstr/hrbrthemes
  - xkcd-style plots: http://xkcd.r-forge.r-project.org/
- ggtree: phylogenetic trees https://bioconductor.org/packages/release/bioc/html/ggtree.html
- geomnet: network visualization
- ggforce: many extensions to ggplot
- ggalt: many extensions to ggplot
- ggraph: plotting graphs/networks
- ggedit: interactive plot editor (Shiny gadget)
- ggRandomForests https://cran.r-project.org/web/packages/ggRandomForests/index.html
- biplots https://github.com/vqv/ggbiplot
- logos: https://github.com/omarwagih/ggseqlogo
- Make any ggplot interactive: 
  - http://rpubs.com/turnersd/ggplotly-demo
  - ggiraph: https://github.com/davidgohel/ggiraph
- ggtern: Ternary diagrams https://bitbucket.org/nicholasehamilton/ggtern
- ggpubr: Publication-ready plots, including
  - barplot alternatives: http://www.sthda.com/english/rpkgs/ggpubr/
  - ggarrange, for flexible multi-panel figures
- Marginal plots: https://twitter.com/ClausWilke/status/900776341494276096
- gggenes: https://github.com/wilkox/gggenes
- egg: combine multiple plots https://cran.rstudio.com/web/packages/egg/
- ggSeqLogo: https://omarwagih.github.io/ggseqlogo/
- Plotting background data: https://drsimonj.svbtle.com/plotting-background-data-for-groups-with-ggplot2
- Patchwork: https://github.com/thomasp85/patchwork
- ggupset: https://github.com/const-ae/ggupset

### Plot Types

- Correlograms: corrgram package
- DiagrammR http://rich-iannone.github.io/DiagrammeR/
- Gene word clouds: http://genomespot.blogspot.co.uk/2014/10/geneclouds-unconventional-genetics-data.html?m=1
- Upset plots: https://cran.r-project.org/web/packages/UpSetR/index.html
- Genomic data: https://bioconductor.org/packages/GenVisR
- hextri: multiclass hexagonal bins https://cran.r-project.org/web/packages/hextri/vignettes/hexbin-classes.html
- trellis: https://www.bioconductor.org/packages/release/bioc/html/gtrellis.html
- Complex heat maps: http://www.bioconductor.org/packages/devel/bioc/html/ComplexHeatmap.html
  - Includes ability to generate upset plots
- Scatterplot Matrix: http://bl.ocks.org/mbostock/4063663
- Beeswarm: https://flowingdata.com/2016/09/08/beeswarm-plot-in-r-to-show-distributions/
- http://flowingdata.com/2016/10/25/r-graph-gallery/
- Tilegrams: http://flowingdata.com/2016/10/13/tilegrams-in-r/
- Karyotype plots: http://bioconductor.org/packages/devel/bioc/html/karyoploteR.html
- Ridgeplots (for multi-sample/category time-series data): 
  - https://github.com/halhen/viz-pub/blob/master/sports-time-of-day/2_gen_chart.R (TODO: make a dedicated geom for this)
  - In python: https://github.com/mwaskom/seaborn/pull/1238
- Half-violin plot: https://gist.github.com/jbburant/b3bd4961f3f5b03aeb542ed33a8fe062

### Data Types

- EpiViz: visualization of epigenomic data sets in R, http://epiviz.github.io/
- Differential expression from RNA-seq: http://bioconductor.org/packages/devel/bioc/html/Glimma.html

### Interactive

- https://gist.github.com/jcheng5/cbcc3b439a949deb544b
- Interactive charts: https://benjaminlmoore.wordpress.com/2015/05/19/interactive-charts-in-r/
- http://www.htmlwidgets.org/showcase_leaflet.html
- Interactive ROC plots: https://github.com/sachsmc/plotROC
- Dull: create interactive web applications - https://github.com/nteetor/dull

## Python

- Graphic of python libraries for visualization: https://pbs.twimg.com/media/C9QBxNsU0AEokR_.jpg
- Seaborn: http://web.stanford.edu/~mwaskom/software/seaborn/index.html
  - Tutorial: http://www.jesshamrick.com/2016/04/13/reproducible-plots/
- Pyxley (Shiny for python): http://multithreaded.stitchfix.com/blog/2015/07/16/pyxley/
- XKCD-style plots: http://jakevdp.github.io/blog/2012/10/07/xkcd-style-plots-in-matplotlib/?utm_content=buffera9a76&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer
- Phylogenetic trees: http://etetoolkit.org/
- Altair: https://altair-viz.github.io/
- Visualization dashboard: https://github.com/facebookresearch/visdom
- Holoviews: http://holoviews.org/
- GGplot clone: https://github.com/has2k1/plotnine
- https://python-graph-gallery.com/
- Plot genome tracks: https://github.com/deeptools/pyGenomeTracks
- PyGal: http://pygal.org/en/stable/
- Web based/interactive
  - Gleam: https://github.com/dgrtwo/gleam
  - Spyre: https://github.com/adamhajari/spyre
  - Bokeh: http://bokeh.pydata.org/docs/user_guide/charts.html
  - Dash: for building data-centric web apps https://github.com/plotly/dash
- Dashboards: https://panel.pyviz.org/
- https://github.com/harbourlab/SparK

## Javascript

- vis.js
- [Javascript libraries](http://blog.nextgenetics.net/?e=62)
- Scraping/JS rendering: https://splash.readthedocs.org/en/latest/
- VegaLite: ggplot-like framework built on top of Vega, which is built on D3.js: https://vega.github.io/vega-lite/
- GPU rendering: https://stardustjs.github.io
- Semiotic: https://emeeks.github.io/semiotic/#/

## Examples

- Cool interactive visualization of differential data: http://graphics.wsj.com/gender-pay-gap/

# Publication/Archiving

- Licenses: http://choosealicense.com/licenses/
- APIs for literature search: http://libguides.mit.edu/apis
- Assessing credit for bioinformatics software authorship: http://depsy.org/
- Icons for presentations: http://cameronneylon.net/blog/some-slides-for-granting-permissions-or-not-in-presentations/
- A nice template for bootstrapping your own academic website using GitHub: https://academicpages.github.io/
- Word choices: [](https://pbs.twimg.com/media/DB4-B__XkAA3Nrw.jpg)
- Prepare papers for any journal format: https://typeset.io/
- Generate manuscripts on GitHub: https://github.com/greenelab/manubot-rootstock
- General-purpose archival format: https://github.com/fair-research/bdbag/

## Code/Data sharing

- DOI for code
  - http://zenodo.org/
  - https://guides.github.com/activities/citable-code/
  - https://mozillascience.github.io/code-research-object/
  - Ruby library for fetching metadata for DOI: https://rubygems.org/gems/terrier
- CodeOcean: web platform to run algorithms https://codeocean.com
- Amazon CodeCommit: http://aws.amazon.com/codecommit/
- Data sharing
  - Dat: https://datproject.org/
  - Academic Torrents: http://academictorrents.com/
  - FigShare: https://figshare.com/
  - Patterns for data sharing: http://project-if.github.io/data-permissions-catalogue/
  - Globus is an open source toolkit for transferring large data files; it implements the GridFTP protocol https://www.globus.org/
  - bbcp is multi-stream scp (for point-to-point large file transfer) https://www.olcf.ornl.gov/kb_articles/transferring-data-with-bbcp/
  - Quilt: package manager for data https://quiltdata.com/
- Data project management: https://www.datazar.com
- Execute code from GitHub project with Jupyter notebooks: http://mybinder.org/
- DataRetriever: http://www.data-retriever.org/
- ResearchObject: http://www.researchobject.org/
- DataLAD: http://datalad.org/
- Publishing/sharing Jupyter notebooks: https://kyso.io/
- Git for data (tables): https://github.com/liquidata-inc/dolt
- Data version control https://dvc.org/

## Journals

- Distill: Online, ML-focused journal http://distill.pub/journal/
- JORS: https://openresearchsoftware.metajnl.com/
- JOSS: http://joss.theoj.org/
- Open source publishing platform: https://libero.pub/

## Writing

- Scripts to identify "bad smells" in science writing (would want to convert this to python): http://matt.might.net/articles/shell-scripts-for-passive-voice-weasel-words-duplicates/
- Collaborative writing
  - DraftIn: https://draftin.com/documents
  - Quip: http://quip.com
  - Canvas: https://usecanvas.com/
- Templates:
  - InDesign template for preprint: https://github.com/cleterrier/ManuscriptTools/blob/master/biorxiv_template_CC2015.indd
  - Rmarkdown templates for journal articles https://github.com/rstudio/rticles
  - GitHub template for authoring papers: https://github.com/peerj/paper-now
  - Two-column rmarkdown template: http://dirk.eddelbuettel.com/code/pinp.html
- Publication-quality tables for markdown: https://cran.r-project.org/web/packages/flextable/index.html
- Convert DOI to Bibtex entry
- Online equation editor: https://www.mathcha.io/
- Simplified markup language for equations: http://asciimath.org/
- Texture: https://github.com/substance/texture
- Scientific writing template for R Markdown: https://rstudio.github.io/radix/
- Extension of MarkDown that handles diagraMs, equations, etc https://casual-effects.com/markdeep/







# Promising methods without software implementation

- SDR gene set analysis http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0928-6
- EpiTensor (3D genomes from 1D data): http://www.nature.com.ezproxy.nihlibrary.nih.gov/ncomms/2016/160310/ncomms10812/full/ncomms10812.html
- MOCHA: identifying modulators of transcriptional regulation from gene expression http://www.nature.com.ezproxy.nihlibrary.nih.gov/articles/srep22656
- HiC deconvolution: http://www.pnas.org.ezproxy.nihlibrary.nih.gov/content/early/2016/03/04/1512577113.full
- Method for multi-omics integration: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1122-6
- Dynamic Bayesian network for predicting TF binding from DNase-seq: http://bioinformatics.oxfordjournals.org/content/26/12/i334.long
- De-noising ChIP-seq using DNNs: http://biorxiv.org/content/biorxiv/early/2016/05/07/052118.full.pdf
- DNA clustering: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-321
- Bayesian identification of bifurcations in single-cell data: http://biorxiv.org/content/biorxiv/early/2016/09/21/076547.full.pdf
- Prediction of promoter-enhancer interactions: http://biorxiv.org/content/biorxiv/early/2016/11/02/085241.full.pdf
- Mocap: TFBS prediction from ATAC-seq http://biorxiv.org/content/biorxiv/early/2016/10/27/083998.full.pdf
- Statistical method to determine deviation from time linearity in epigenetic aging: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005183
- Multi-tissue eQTL meta-analysis: http://biorxiv.org/content/early/2017/01/16/100701
- Missing value imputation evaluator: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1429-3
- Simulated annealing for biological network alignment: https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btx090/2996219/SANA-Simulated-Annealing-far-outperforms-many
- TF binding site prediction using memory-matching networks: https://arxiv.org/pdf/1702.06760.pdf
- V-ALIGN: alignmnet on genome graphs http://www.biorxiv.org/content/biorxiv/early/2017/04/06/124941.full.pdf
- Self-organizing maps for single-cell http://www.biorxiv.org/content/biorxiv/early/2017/04/05/124693.full.pdf
- New locality sensitive hashing method: http://www.biorxiv.org/content/biorxiv/early/2017/08/25/180471.full.pdf
- Extracting knowledge from RFs: https://www.biorxiv.org/content/early/2017/11/10/217695
- Bloom filters without false-positives across a constrained set of values: http://lendulet.tmit.bme.hu/lendulet_website/wp-content/papercite-data/pdf/kiss2018bloom.pdf
- Adaptive pruning of boosted trees: https://arxiv.org/pdf/1805.07592.pdf
- Calculations to predict VAF for a mutation in a somatic sample: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6148761/
- Zip Trees: https://arxiv.org/pdf/1806.06726.pdf
- Order Min Hash: https://www.biorxiv.org/content/10.1101/534446v1
- Phasing by VAF in tumor-normal samples: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2753-1
- Variance stabilization of genomic signals https://www.biorxiv.org/content/10.1101/2020.01.31.929174v1
- Fast computation of edit distance https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-3184-8
- Local leveling for bias reduction in RNA-seq https://www.biorxiv.org/content/10.1101/2023.10.31.564992v1
- https://www.biorxiv.org/content/10.1101/2024.02.20.581311v1
