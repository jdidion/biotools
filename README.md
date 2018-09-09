This is a collection of bioinformatics tools I have sourced from recent literature, organized by topic. I have not used most of these tools.

**Table of Contents**

<!-- toc -->

- [Discovery](#discovery)
- [Data Sets](#data-sets)
- [Genomics](#genomics)
  * [General Information](#general-information)
  * [Algorithms](#algorithms)
  * [Assay Design](#assay-design)
  * [Candidate Prioritization](#candidate-prioritization)
  * [Databases](#databases)
  * [Data Formats](#data-formats)
  * [Footprinting](#footprinting)
  * [Functional Enrichment](#functional-enrichment)
  * [GWAS/QTL](#gwasqtl)
  * [Methylation Array](#methylation-array)
  * [Network Analysis](#network-analysis)
  * [Prediction](#prediction)
  * [Sequencing Protocols](#sequencing-protocols)
  * [Simulation](#simulation)
  * [SNP Array](#snp-array)
  * [Variant annotation](#variant-annotation)
  * [Sequence Analysis](#sequence-analysis)
    + [General-purpose](#general-purpose)
    + [QC](#qc)
    + [ATAC-seq/DNase-seq](#atac-seqdnase-seq)
    + [ChIP-seq](#chip-seq)
    + [Chromatin Interactions](#chromatin-interactions)
    + [DNA](#dna)
    + [Methylation](#methylation)
    + [RNA](#rna)
    + [Single-cell](#single-cell)
    + [Integrated Methods](#integrated-methods)
  * [Upcomming methods to watch for](#upcomming-methods-to-watch-for)
- [General Programming Resources](#general-programming-resources)
  * [C++](#c)
  * [R](#r)
    + [Find packages](#find-packages)
    + [Database](#database)
    + [Data Cleaning](#data-cleaning)
    + [Reporting](#reporting)
    + [Misc](#misc)
  * [Python](#python)
  * [HPC](#hpc)
  * [Command Line (OSX/Linux)](#command-line-osxlinux)
  * [Reproducibility/Containerization](#reproducibilitycontainerization)
    + [Building Pipelines](#building-pipelines)
- [Statistics/Machine Learning](#statisticsmachine-learning)
  * [Methods/algorithms](#methodsalgorithms)
  * [Python](#python-1)
  * [Deep Learning](#deep-learning)
  * [Web APIs](#web-apis)
  * [Data Sets](#data-sets-1)
  * [Text classification](#text-classification)
  * [Misc](#misc-1)
- [Visualization](#visualization)
  * [Networks](#networks)
  * [Phylogenetic trees](#phylogenetic-trees)
  * [R](#r-1)
    + [ggplot2](#ggplot2)
    + [Plot Types](#plot-types)
    + [Data Types](#data-types)
    + [Interactive](#interactive)
  * [Python](#python-2)
  * [Javascript](#javascript)
  * [Examples](#examples)
- [Publication/Archiving](#publicationarchiving)
  * [Writing](#writing)
- [Promising methods without software implementation](#promising-methods-without-software-implementation)

<!-- tocstop -->

# Clinical data

## EHR

* BlueButton related tools https://github.com/amida-tech

# Data Management

* Send stdout to a Google Sheet: https://github.com/kren1/tosheets

# Data Sets

* 1000 Genomes (RNA-Seq, ChIP-Seq): http://archive.gersteinlab.org/docs/2015/06.04/1kg_fun_studies.htm
* Reference panel of ~250 Dutch families http://biorxiv.org/content/early/2016/01/18/036897
* FANTOM consortium has CAGE (5' single molecule RNA counting) data from ~1000 human cell/tissue/cell-line samples from ~300 different cell/tissue types
* Full text of all PMC papers from 2008-present: ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/
* Phased genome sequences
    * >100 fully-phased: http://gigascience.biomedcentral.com/articles/10.1186/s13742-016-0148-z
    * Statistically phased: http://www.haplotype-reference-consortium.org/
    * Phased variants: http://genome.cshlp.org/content/early/2016/11/25/gr.210500.116.full.pdf+html
* Clinical trials: http://vivli.org/
* Targets for drug discovery: https://www.targetvalidation.org/
* Large medical datasets for ML: https://github.com/beamandrew/medical-data
* Species images: http://phylopic.org/image/browse/
* Public RNA-Seq data: 
    * https://jhubiostatistics.shinyapps.io/recount/ 
    * with phenotype predictions https://bioconductor.org/packages/release/bioc/html/recount.html
    * human and mouse: http://amp.pharm.mssm.edu/archs4
* Whole genomes of 150 Danish individuals: http://www.nature.com/nature/journal/v548/n7665/full/nature23264.html
* Simons Genome Diversity Project: https://www.simonsfoundation.org/simons-genome-diversity-project/
* BloodPAC: https://doi.bloodpac.org/BLOODPAC.0001/
* Reference haplotypes for STR phasing: https://www.biorxiv.org/content/early/2018/03/06/277673
* http://www.docm.info/
* Normalized, batch-corrected RNA-seq for both normal (GTEx) and cancer (TCGA) data: https://www.nature.com/articles/sdata201861
* Hartwig WGS (requires data request): https://www.hartwigmedicalfoundation.nl/en/wgs-database/
* Orthogonally validated TP and FP CNVs: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5473400/
* Microarray transcriptome profiles prepared for ML applications: https://www.biorxiv.org/content/biorxiv/early/2018/06/22/353698.full.pdf
* Cloud resource for filtering 1kg data
    * http://bamsi.research.it.uu.se/
    * Python API: https://github.com/NGDSG/BAMSI-API
* Paired mRNA and protein abundance in 26 human tissues: https://www.biorxiv.org/content/biorxiv/early/2018/06/27/357137.full.pdf

# Discovery

When looking for a bioinformatics tool for a specific application:

* http://omictools.com/
* http://www.gitxiv.com/?cat%5B0%5D=bioinformatics
* https://bio.tools/
* https://biosharing.org

# Genomics

## General Information

* Migrating to GRCh38: https://software.broadinstitute.org/gatk/blog?id=8180
* Mappings between contig names in different assemblies: https://github.com/dpryan79/ChromosomeMappings
* Command-line interface to download and manage reference genomes: https://github.com/simonvh/genomepy

## Algorithms

* Suffix arrays: http://almob.biomedcentral.com/articles/10.1186/s13015-016-0068-6
* gapped k-mer SVM: classifiers for DNA and protein sequences http://www.beerlab.org/gkmsvm/
    * Version for large-scale data: https://github.com/Dongwon-Lee/lsgkm/
* Fast BWT creation: https://github.com/hitbc/deBWT

## Assay Design

* Choosing assays based on complementarity to existing data: https://github.com/melodi-lab/Submodular-Selection-of-Assays
* MBRAnator: design of MPRA libraries https://www.genomegeek.com/

## Candidate Prioritization

* http://bioconductor.org/packages/GenRank/
* Bayesian prioritizaiton of rare functional variants using RNA-seq data: https://github.com/ipw012/RIVER
* http://snp-nexus.org/IW-Scoring/

## Databases

* RNA Meta Analysis (web app)
	* Perform meta-analyses of GEO microarray and RNA-Seq studies
	* Correlate resulting signatures to CMAP02/L1000 perts and 26,000+ GEO studies
	* https://www.rnama.com
* NAR catalog of databases, by subject: http://www.oxfordjournals.org/our_journals/nar/database/c/
* Super-Enhancer Archive: http://www.bio-bigdata.com/SEA/
* GWAS database: http://jjwanglab.org/gwasdb
* rVarBase: regulatory features of human genetic variants http://rv.psych.ac.cn/
* TransVar http://bioinformatics.mdanderson.org/transvarweb/
* dbMAE: mono-allelic expression https://mae.hms.harvard.edu/
* Disease-gene associations http://www.disgenet.org/web/DisGeNET/menu/rdf
* BISQUE: convert between database identifiers http://bisque.yulab.org/
* Human tissue-specific enhancers: http://www.enhanceratlas.org/
* Searh HLI's genome data: hli-opensearch.com
* Chromatin-state annotations + per-base functionality scores for 164 cell types: http://noble.gs.washington.edu/proj/encyclopedia/
* Feature-based classification of human transcription factors: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1349-2
* Database of epifactors (epigenetic factors): http://epifactors.autosome.ru/
* Database of disease-associated methylation: http://202.97.205.78/diseasemeth/
* SRA metadata: http://deweylab.biostat.wisc.edu/metasra/
* Human histone modifications: http://www.tongjidmb.com/human/index.html
* iRegNet3D: SNP-focused catalog of TF-TF, TF-DNA, and DNA-DNA interactions http://iregnet3d.yulab.org/index/
* CTD: drug interactions and toxicity http://ctdbase.org/
* FANTOM lncRNA catalog: http://fantom.gsc.riken.jp/cat/
* Integrated database of public ChIP-Seq datasets: http://chip-atlas.org/
* Alternative splicing: http://vastdb.crg.eu/wiki/Main_Page
* Interactive multi-omics tissue assay database: https://ccb-web.cs.uni-saarland.de/imota/
* Database of cis-regulatory elements (enhancers): http://www.kostkalab.net/software.html
* Database of genetic variant effects on gene expression: https://xhaubem01.u.hpc.mssm.edu/gwas2genes/
* Structural variation: https://www.ncbi.nlm.nih.gov/dbvar
* PharmacoDB
    * Search multiple cancer pharmacogenomic databases with a single query
    * Software is GPL licensed; target databases have various licenses
    * http://pharmacodb.pmgenomics.ca
* https://hive.biochemistry.gwu.edu/biomuta/readme
* Database of druggable variant information (web only): http://depo-dinglab.ddns.net/
* CHESS: new gene annotation database
    * Contains most genes from RefSeq and Gencode, plus additional genes discovered from GTEx transcripts
    * http://ccb.jhu.edu/chess/

## Data Formats

* Fast ordered sampling of records: https://github.com/tonymugen/sampleSNPs/
* Searching overlaps of large collections of interval databases: https://github.com/ryanlayer/giggle/
* FASTA/FASTQ
    * Read filtering and extraction: https://github.com/ad3002/Cookiecutter
    * fqtools: for working with FASTQ files https://github.com/alastair-droop/fqtools
    * Integrating results from multiple genome binning programs: https://github.com/songweizhi/Binning_refiner
* SAM/BAM/CRAM
    * Read filtering and profiling: https://github.com/jwalabroad/VariantBam
* VCF
    * Annotation: https://github.com/brentp/vcfanno
    * https://github.com/lh3/bgt
    * GQT
    * BCFtools: includes new tool to identify RoH http://samtools.github.io/bcftools/
    * Work with VCF in R: https://github.com/knausb/vcfR
    * More VCF tools: http://vcf-kit.readthedocs.io/
    * Validation: https://github.com/EBIvariation/vcf-validator
* BED/GFF
    * https://github.com/ihh/gfftools
    * https://github.com/billzt/gff3sort
    * Combine p-values https://github.com/brentp/combined-pvalues
    * Assessing interval overlap of multiple genomic features: https://github.com/andrew-leith/GINOM
* Improved command line viewing of formatted data
    * http://lindenb.github.io/jvarkit/PrettySam.html
    * https://biosyntax.org/install
    * https://github.com/koelling/dnacol
* Generate NCBI submission data: https://github.com/genomeannotation/GAG
* Convert between many different formats, using a SQL database as intermediate representation: https://watson.hgen.pitt.edu/mega2/mega2r/

## Footprinting

* http://rajanil.github.io/msCentipede/
* NucID: nucleosome positioning from DNase-seq https://jianlingzhong.github.io/NucID/
* SeqGL: predict TF binding from DNase/ATAC-seq https://bitbucket.org/leslielab/seqgl/wiki/Home
* DeFCoM: https://bitbucket.org/bryancquach/defcom

## Functional Enrichment/Ontology

* General:
    * R interface to DAVID: http://www.bioconductor.org/packages/release/bioc/html/RDAVIDWebService.html
    * GSEA
        * Cautions about the GSEA null model: http://bioinformatics.oxfordjournals.org/content/early/2017/01/02/bioinformatics.btw803.short
    * GiANT: uncertainty in GSEA https://cran.r-project.org/web/packages/GiANT/index.html
    * Ensemble gene set enrichment analysis: http://bioconductor.org/packages/release/bioc/html/EGSEA.html
    * Fast GSA: https://github.com/billyhw/GSALightning
    * Fast GSEA: https://github.com/ctlab/fgsea
    * Multi-dimensional GSEA: http://bioconductor.org/packages/release/bioc/html/mdgsa.html
    * GSEA with external information https://cran.r-project.org/web/packages/netgsa/index.html
    * QTest: http://statgen.snu.ac.kr/software/QTest/
    * Gene set analysis with specific alternative hypothesis: https://bioconductor.org/packages/release/bioc/html/GSAR.html
    * clusterProfiler: <https://guangchuangyu.github.io/clusterProfiler/>
    * DOSE: disease ontology <https://guangchuangyu.github.io/dose/>
    * ReactomePA: <https://guangchuangyu.github.io/reactomepa/>
    * https://cran.r-project.org/web/packages/SetRank/index.html
    * Identify and rank significance of overlaps: https://github.com/ryanlayer/giggle
    * Generate feature vectors from ontology annotations: https://github.com/bio-ontology-research-group/onto2vec
* Gene sets
    * http://bioconductor.org/packages/release/data/experiment/html/EGSEAdata.html
    * http://genomespot.blogspot.com.au/2017/02/introducing-encode-gene-set-hub.html
* GO term:
    * http://cbl-gorilla.cs.technion.ac.il/
    * http://lrpath.ncibi.org/
    * Reduce GO term lists: 
    	+ http://revigo.irb.hr/
	+ <http://guangchuangyu.github.io/2015/10/use-simplify-to-remove-redundancy-of-enriched-go-terms/>
    * GO Express: https://www.bioconductor.org/packages/release/bioc/html/GOexpress.html
    * GO Extender: https://www.msu.edu/~jinchen/GOExtender/
    * http://iwera.ir/~ahmad/dal/
    * Negative GO enrichment: https://sites.google.com/site/guoxian85/neggoa
* Variant Set
    * https://cran.r-project.org/web/packages/VSE/vignettes/my-vignette.html
    * Functional enrichment with LD correction
* MESH
    * MeSH over-representation: http://www.bioconductor.org/packages/release/bioc/vignettes/meshr/inst/doc/MeSH.pdf
    * meshes: <https://guangchuangyu.github.io/meshes/>
* Regional
    * LOLA: http://lola.computational-epigenetics.org
    * AnnotatR: https://github.com/rcavalcante/annotatr/
    * Compare epigenetic features in multiple samples: http://epigenome.wustl.edu/EpiCompare1/
* Trait
    * traseR: Trait-associated SNP enrichment https://www.bioconductor.org/packages/release/bioc/html/traseR.html
    * Identifying genetic heterogeneity within phenotypically defined subgroups: https://github.com/jamesliley/subtest
* Multi-omics
    * Single-sample GSA across data sets: https://www.bioconductor.org/packages/3.3/bioc/html/mogsa.html
* Network-based
    * https://cran.r-project.org/web/packages/neat/index.html
    * https://github.com/sarbal/EGAD

## GWAS/QTL

* Association testing
    * QTLtools: Pipeline for molecular QTL analysis https://qtltools.github.io/qtltools/
        * FastQTL: http://fastqtl.sourceforge.net/
    * RASQUAL: allele-specific QTL using phased SNPs - https://github.com/dg13/rasqual
    * Relatedness, PCA: http://zhengxwen.github.io/SNPRelate/
    * Multi-SNP, multi-trait regression https://github.com/ashlee1031/BERRRI
    * regioneR: permutation testing for association between genomic region and phenotype http://bioconductor.org/packages/release/bioc/html/regioneR.html
    * https://sites.google.com/site/honglee0707/mtg2
    * Use local gene networks to improve trans-eQTL detection: https://github.com/PMBio/GNetLMM
    * Fast correlation testing: https://github.com/gabraham/flashpca/tree/master/flashpcaR
    * Gene and pathway association testing from summary statistics: https://cran.r-project.org/web/packages/aSPU/
    * Account for LD and functional information: https://github.com/yjingj/SFBA
    * Correcting for prediction uncertainty in TWAS: http://biorxiv.org/content/early/2017/02/14/108316
    * Using random forests https://github.com/0asa/TTree-source
    * Using k-mers: https://github.com/atifrahman/HAWK
    * BOLT-LMM: https://data.broadinstitute.org/alkesgroup/BOLT-LMM/
* Variance eQTL
    * https://funpopgen.github.io/veqtl-mapper/
* Multiple test correction
    * eigenMT: Efficient multiple-test correction http://montgomerylab.stanford.edu/resources/eigenMT/eigenMT.html
    * Fast multiple-test correction for LMMs: http://genetics.cs.ucla.edu/multiTrans/
    * Hierarchical eQTL MTC: http://bioinformatics.org/treeqtl/
    * Controlling bias in EWAS/TWAS using null distribution: http://bioconductor.org/packages/bacon/
* Prioritization
    * GenoWAP: Prioritization of GWAS signals using functional information http://genocanyon.med.yale.edu/GenoWAP
    * HitWalker2: https://github.com/biodev/HitWalker2
    * https://nijingchao.github.io/CRstar/
    * http://genetics.bwh.harvard.edu/pines/
* Fine-mapping
    * Using summary statistics http://www.christianbenner.com
    * http://bioinformatics.oxfordjournals.org/content/32/3/330.full
    * PAINTOR: fine mapping, prioritization - https://github.com/gkichaev/PAINTOR_FineMapping/
    * Genotype synthesis: https://sourceforge.net/projects/getsynth/
    * Prediction of causal variants from epigenomic annotations: https://github.mit.edu/liyue/rivieraBeta
    * DAP: Bayesian framework for QTL analysis and fine-mapping https://github.com/xqwen/dap
    * BayesFM: https://sourceforge.net/projects/bayesfm-mcmc-v1-0/
    * http://apps.biocompute.org.uk/haprap/
    * Determining causal genes using TADs http://biorxiv.org/content/early/2016/11/15/087718
    * Network analysis: https://github.com/YuanlongLiu/SigMod
    * Haplotype-based: http://apps.biocompute.org.uk/haprap/
* LD score calculation and regression https://github.com/bulik/ldsc
* Find directional effect of a signed functional annotation using summary statistics: https://github.com/yakirr/sldp
* Imputation of missing phenotype information http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3513.html
* Epistasis
    * General linear-mixed model library; also includes mixed-RF method for detecting epistasis with population structure correction: https://github.com/PMBio/limix
    * GPU-accelerated detection of epistasis using Bayesian neural networks: https://github.com/beamandrew/BNN
    * MEPID: marginal epistasis test http://www.xzlab.org/software.html
* Other
    * Browser for geographical distribution of genetic variants: http://popgen.uchicago.edu/ggv/
    * Integration of GWAS with molecular QTL: https://github.com/xqwen/integrative
    * GxE interactions: https://github.com/davidaknowles/eagle

## Microarray

* Correct for batch effects between training data and external datasets: https://cran.r-project.org/web/packages/bapred/index.html
* Impute from Affy expression arrays: http://simtk.org/home/affyimpute
* SNP
    * Call haplotypes https://cran.r-project.org/web/packages/GHap/index.html
    * Neural network for breakpoint detection for CNV calling: https://www.biorxiv.org/content/biorxiv/early/2018/06/24/354423.full.pdf
* Methylation
    * Minfi: R package for working with 450k methylation arrays
    * D3M: two-sample test of differential methylation from distribution-valued data https://cran.r-project.org/web/packages/D3M/D3M.pdf
    * Network-based approach to discovering epigenetic "modules" that can be associated with gene expression: http://bioinformatics.oxfordjournals.org/content/30/16/2360.long
    * Filtering probes using technical replicates https://cran.r-project.org/web/packages/CpGFilter/index.html
    * Imputation of genome-wide methylation: http://wanglab.ucsd.edu/star/LR450K/
    * Tutorial for analysis using bioconductor packages: http://biorxiv.org/content/biorxiv/early/2016/05/25/055087.full.pdf
    * Normalization
        * Probe design bias correction: https://www.bioconductor.org/packages/release/bioc/html/ENmix.html
    * DMR calling
        * https://github.com/raivokolde/seqlm
    * http://aminmahpour.github.io/PyMAP/
    * Interactive exploration: http://bioconductor.org/packages/release/bioc/html/shinyMethyl.html
    * eFORGE: identify cell type-specific signals in differentially methylated positions (mostly important for blood-based EWAS) http://eforge.cs.ucl.ac.uk/
    * Model for EWAS using probe signal intensities: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1347-4
    * Reference-based tissue deconvolution (three different algorithms): https://github.com/sjczheng/EpiDISH
    * Glint pipeline (qc, EWAS, population structure): http://glint-epigenetics.readthedocs.io/en/latest/
    	* Bayesian extension of Refactor cell type heterogeneity correction that incorporates experimentally determined cell counts: https://github.com/cozygene/BayesCCE
    * https://github.com/perishky/meffil/ 

## Motif/TFBS

* Fast, SNP-aware PWM matching https://www.cs.helsinki.fi/group/pssmfind/
* Cell type-specific TFBS analysis (focuses analysis on TFs expressed in cell type of interest): https://github.com/Danko-Lab/rtfbs_db
* Bayesian motif discovery: https://github.com/soedinglab/BaMMmotif
* R package for TFBS analysis: http://bioconductor.org/packages/release/bioc/html/TFBSTools.html
* FIMO and MCAST perform best of TFBS predictors: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1298-9
* Intra-motif dependencies
    * Identifying http://www.jstacs.de/index.php/InMoDe
    * Visualizing http://bioinformaticstools.mayo.edu/circularlogo/index.html
* Dinucleotide weight tensors encode dependencies between positions in motifs: http://dwt.unibas.ch/
* Circular logos: http://bioinformaticstools.mayo.edu/circularlogo/index.html
* Convert kernels learned by CNN to PWMs: ftp://ftp.cbi.pku.edu.cn/pub/software/CBI/k2p
* https://github.com/schulter/crbm
* https://github.com/daquang/YAMDA

## Network Analysis

* Identification of module with biclustering: http://web.ist.utl.pt/rmch/bicnet/temporary/index.jsp
* Functional analysis https://bioconductor.org/packages/release/bioc/html/EGAD.html
* GAGE: http://bioconductor.org/packages/release/bioc/html/gage.html
* PAXToolsR: http://bioconductor.org/packages/release/bioc/html/paxtoolsr.html
* ancGWAS: post GWAS association with protein-protein interaction networks http://www.cbio.uct.ac.za/~emile/software.html
* Gene network reconstruction
    * For a set of TFs: https://sourceforge.net/projects/aracne-ap/
    * From PPI or motif sharing: https://github.com/davidvi/pypanda (is also an integrative method that can incorporate multiple sources of information)
* BANFF: https://cran.r-project.org/web/packages/BANFF/index.html
* https://bitbucket.org/abarysh/safe
* https://bitbucket.org/roygroup/merlin-p
* OSS alternative to Inginuity pathway analysis: https://www.bioconductor.org/packages/release/bioc/html/QuaternaryProd.html
* SAFE: spatial analysis of functional enrichment https://bitbucket.org/abarysh/safe
* Similarity search: https://github.com/zhangjiaobxy/nssrfPackage
* Tissue-specific: https://cran.r-project.org/web/packages/GRAPE/index.html
* CLR with B-Spline for mutual information-based inference of transcriptional regulatory networks: https://bitbucket.org/Jonathan-Ish-Horowicz/fastgenemi/
* Functional module identification and enrichement through Gene Ontology terms https://www.nature.com/articles/s41598-018-23672-0

## Population genetics

* Population history from unphased whole-genomes: https://github.com/popgenmethods/smcpp
* Management of allele-frequency data: https://grenaud.github.io/glactools/

## Prediction

* QTL
    * Imputation of summary statistics in multi-ethnic cohorts: http://dleelab.github.io/jepegmix/
    * Causal variant identification: http://genetics.cs.ucla.edu/caviar/
* eQTL
    * Imputation of gene expression from genotype data : https://github.com/hriordan/PrediXcan
* Genetic risk
    * https://github.com/yiminghu/AnnoPred
* Causal variant
    * Ensemble method: https://github.com/gifford-lab/EnsembleExpr/
    * eCAVIAR: probability that a variant is causal for both QTL and eQTL http://genetics.cs.ucla.edu/caviar/
    * https://github.com/dleelab/qcat
    * Disease-associated risk variants: https://sites.google.com/site/emorydivan/
    * Predicting gene targets from GWAS summary statistics https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4979185/
    * https://github.com/igm-team/orion-public
* Chromatin States
    * GenoSTAN: http://bioconductor.org/packages/release/bioc/html/STAN.html
    * R package for predicting chromatin states from histone marks across conditions https://github.com/ataudt/chromstaR
    * Rule-based: http://www.statehub.org/
    * Hierarchical HMM: https://github.com/gcyuan/diHMM
    * https://github.com/calico/basenji
    * https://bitbucket.org/roygroup/cmint
* Enhancers
    * Prediction of enhancer strength from sequenced http://bioinformatics.hitsz.edu.cn/iEnhancer-2L
    * Prediction of core cell type-specific TFs from super enhancers https://bitbucket.org/young_computation/crcmapper
    * Prediction of superenhancers https://github.com/asntech/improse
    * Deep learning-based:
        * https://github.com/wenjiegroup/BiRen
        * https://www.biorxiv.org/content/early/2018/08/24/398115
        * https://www.biorxiv.org/content/early/2018/02/14/264200
* Coding mutations
    * Predict mutation effects from sequence covariation: https://marks.hms.harvard.edu/evmutation/
    * Predicting loss-of-function from expression: https://www.nature.com/articles/s41467-017-00443-5
    * Prediction of missense variants: https://www.biorxiv.org/content/early/2018/02/02/259390
    * Disease-specific functional prediction https://sites.google.com/site/emorydivan/
    * Impact of coding SNPs: http://pantherdb.org/tools/csnpScoreForm.jsp
    * Predict disease risk from GWAS summary statistics: https://github.com/yiminghu/AnnoPred
* Regulatory variants/TF binding
    * LedPred: prediction of regulatory sequences from ChIP-seq https://github.com/aitgon/LedPred
    * GERV: prediction of regulatory variants that affect TF inding http://gerv.csail.mit.edu/
    * Score variant deleteriousness: http://cadd.gs.washington.edu/
    * BASSET: prediction of sequence activity https://github.com/davek44/Basset
    * DanQ: hybrid convolutional and recurrent neural network model for predicting the function of DNA de novo from sequence http://github.com/uci-cbcl/DanQ
    * [LINSIGHT](https://github.com/CshlSiepelLab/LINSIGHT)
    * Protein binding affinity: https://bitbucket.org/wenxiu/sequence-shape.git
    * Change in local frustration index: https://github.com/gersteinlab/frustration
    * TFImpute: multi-task learning from ChIP-seq data across factors and tissues to impute TF binding for an unassayed tissue/factor combination: https://bitbucket.org/feeldead/tfimpute
    * PPI https://www.ncbi.nlm.nih.gov/research/mutabind/index.fcgi/
    * Multiple instance learning of TF binding: http://www.cs.utsa.edu/~jruan/MIL/
    * Cell type-specific: https://github.com/uci-cbcl/FactorNet
    * Predict TF binding from ATAC-Seq using deep neural network:  https://github.com/hiranumn/deepatac
    * http://genome.grid.wayne.edu/centisnps/
    * Benchmark: https://github.com/Oncostat/BenchmarkNCVTools
    * https://github.com/samesense/pathopredictor
    * Cell-type agnostic regulatory activity prediction from ENCODE: http://screen.encodeproject.org/
    * Segway functional scores: https://noble.gs.washington.edu/proj/encyclopedia/
    * Predict effect of variants on epigenetic factors (chromatin accessibility, histone marks, etc): http://deepfigv.mssm.edu/downloads.html
    * HaploReg: http://www.broadinstitute.org/mammals/haploreg/haploreg.php
    * Consensus approaches:
        * http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004962
        * http://jjwanglab.org/PRVCS/
    * http://queryor.cribi.unipd.it/cgi-bin/queryor/mainpage.pl
    * Using epigenomic data https://github.com/mulin0424/cepip
* Methylation
    * CpGenie: predicts methylation from sequence, predicts impact of variants on nearby methylation https://github.com/gifford-lab/CpGenie
* Chromatin accessibility
    * http://scm.csail.mit.edu/
* TFBS
    * Predict TF binding affinities using open chromatin + PWMs: https://github.com/schulzlab/TEPIC
    * LR-DNAse: TFBS prediction using features derived from DNase-seq: http://biorxiv.org/content/early/2016/10/24/082594
* Classification of cis-regulatory modules: https://github.com/weiyangedward/IMMBoost
* Imputation: https://github.com/tdurham86/PREDICTD
* Variable selection for random forest: https://github.com/jomayer/SMuRF

## Security
    * Secret sharing schemes for keeping patient ID secret while being able to reconstruct it from other identifiers: https://pdfs.semanticscholar.org/0307/48be9820512a1d5582351552bd0452711296.pdf
    
## Sequencing Protocols

* Single cell
    * Simultaneous RNA and methylation (and inference of CNV): http://www.nature.com/cr/journal/vaop/ncurrent/full/cr201623a.html
    * Simultaneous RNA and methylation (scM&T-seq): http://www.nature.com/nmeth/journal/v13/n3/full/nmeth.3728.html
    * Simultaneous RNA and protein measurements: http://www.sciencedirect.com/science/article/pii/S2211124715014345
* CRISPR-DS: Uses CRISPR/Cas9 excision of target sequences to improve library quality relative to hybrid capture: https://www.biorxiv.org/content/early/2017/10/23/207027

## Simulation

* InterSIM: simulate correlated multi-omics data https://cran.r-project.org/web/packages/InterSIM/index.html
* Simulate capture sequencing: https://github.com/mdcao/capsim
* Database of simulated tumor genomes (created with BamSurgeon on NA12878 background): https://www.biorxiv.org/content/early/2018/02/07/261503
* Simulation of cancer genomes: https://github.com/rsemeraro/XomeBlender
* Simulate CNVs: https://github.com/pughlab/bamgineer
* Simulate mutations in a BAM file: https://github.com/adamewing/bamsurgeon

## Variant annotation

* WGSA pipeline https://sites.google.com/site/jpopgen/wgsa/
* http://snpeff.sourceforge.net/
* Normalization of SNP ID's from literature: https://github.com/rockt/SETH
* https://hail.is/
* Tissue-specific https://github.com/kevinVervier/TiSAn
* VCF visualization with Circos plot: http://legolas.ariel.ac.il/~tools/CircosVCF/
* HGVS
    * This is a standard notation for variants http://varnomen.hgvs.org/
    * Toolkit for working with HGVS-formatted variants: https://github.com/biocommons/hgvs
* Integration of multiple predictive measures of mutation deleterious effect: https://www.ncbi.nlm.nih.gov/research/snpdelscore/
* Varsome: database that aggregates variant information from many sources
    * https://varsome.com/
    * Provides a REST API that can be queried with dbSNP or HGVS IDs
* Constrained coding regions: https://github.com/quinlan-lab/ccrhtml
* Database of human chromosomal fragile sites: http://webs.iiitd.edu.in/raghava/humcfs
* FunSeq2: http://funseq2.gersteinlab.org
* Deep learning-based: https://www.biorxiv.org/content/early/2017/12/18/235655.1
* Search engine for indexing and intersecting with large numbers of annotation BED/VCF files: https://github.com/ryanlayer/giggle/
* SV: http://lbgi.fr/AnnotSV/
* Clinical annotation: https://github.com/ding-lab/CharGer

## Sequence Analysis

### General-purpose

* Google Genomics R API: https://followthedata.wordpress.com/2015/02/05/notes-on-genomics-apis-2-google-genomics-api/
* k-mer counting
    * khmer:
        * https://github.com/ged-lab/khmer
        * https://docs.google.com/presentation/d/1biQmLkwPlCOA56mNZdUAiXa1OyGE0qyvI2nvU0qlOIE/mobilepresent?pli=1&slide=id.p58
    * https://github.com/abdullah009/kcmbt_mt
    * KAT https://github.com/TGAC/KAT
    * ntCard https://github.com/bcgsc/ntCard
    * KMC3 https://github.com/refresh-bio/KMC
    * https://github.com/uni-halle/gerbil
    * Toolkit for working with unique kmers: https://github.com/shenwei356/unikmer
* Density-based clustering: https://bitbucket.org/jerry00/densitycut_dev
* chopBAI: segment BAM indexes by region for faster access https://github.com/DecodeGenetics/chopBAI
* GFFutils: http://daler.github.io/gffutils/
* R package for aligned chromatin-oriented sequencing data: https://cran.r-project.org/web/packages/Pasha/
* MMR: resolve multi-mapping reads https://github.com/ratschlab/mmr
* BAMQL: query language for extracting reads from BAM files https://github.com/BoutrosLaboratory/bamql
* SAMBAMBA: samtools alternative
* BAMtools: another samtools alternative, plus some additional tools https://github.com/pezmaster31/bamtools/wiki
* DeepTools: more useful SAM/BAM operations http://deeptools.readthedocs.io/en/latest/content/list_of_tools.html
* bedtools http://bedtools.readthedocs.io/en/latest/
* bedops alternative/additional BED operations http://bedops.readthedocs.io/en/latest/
* Normalization:
    * https://github.com/allenxhcao/glscale
    * https://github.com/stephaniehicks/qsmooth
    * https://github.com/SchulzLab/ORNA
    * https://github.com/jpvert/suquan
* Demultiplexing/deduping barcoded reads w/ UMIs: http://gbcs.embl.de/portal/tiki-index.php?page=Je
* Hardware acceleration of alignment (requires $5k FPGA module): https://github.com/BilkentCompGen/GateKeeper
* Detection and removement of barcode swapping (issue on Illumina sequencers that used patterned flow cells: https://github.com/MarioniLab/BarcodeSwapping2017
* Data processing pipelines for many types of omics data, built using NextFlow and Singularity: https://github.com/c-guzman/cipher-workflow-platform
* Index and fetch data from BGZF-compressed files.
* https://github.com/baudisgroup/segment-liftover
* Recover unaligned reads: https://github.com/VCCRI/Scavenger
* Intersection and visualization of multiple gene/region sets: https://bitbucket.org/CBGR/intervene
* Mantis: index of raw-read datasets for efficient and exacty queries https://github.com/splatlab/mantis
* Machine learning method for determining sequence identity: https://github.com/TulsaBioinformaticsToolsmith/FASTCAR

### Demultiplexing

* https://axe-demultiplexer.readthedocs.io/en/latest/usage.html

### QC

* Qualimap2: http://qualimap.bioinfo.cipf.es/
* Determine whether two BAM files are from the same source: https://bitbucket.org/sacgf/bam-matcher
* DOGMA: Measure completeness of a transcriptome or proteome assembly https://ebbgit.uni-muenster.de/domainWorld/DOGMA/
* Identify and remove UMI sequences from reads: https://github.com/CGATOxford/UMI-tools
* Integrated report from multiple tools: http://multiqc.info/
* Batch effects: 
    * https://github.com/mani2012/BatchQC
    * Correct batch effects using residual neural net: https://github.com/ushaham/BatchEffectRemoval
* AlmostSignificant: https://github.com/bartongroup/AlmostSignificant
* Genetic relatedness from raw reads:
    * https://github.com/kdmurray91/kwip
    * https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi
* Fast coverage estimate from BAM index: https://github.com/brentp/goleft/tree/master/indexcov
* Detecting sample swaps: https://github.com/PapenfussLab/HaveYouSwappedYourSamples
* QC Fail articles:
    * Patterned flow cells (HiSeq 3000+) have high rates of optical duplicates: https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/
* http://samstat.sourceforge.net/
* Fingerprints: http://db.systemsbiology.net/gestalt/genome_fingerprints/
* http://fastq.bio/
* https://github.com/nilesh-tawari/ChronQC
* Read Origin Protocol: https://github.com/smangul1/rop
* QC on Illumina run folder: https://github.com/Molmed/checkQC
* New diploid reference for benchmarking: https://www.biorxiv.org/content/early/2017/11/22/223297
* Historical tracking: https://github.com/nilesh-tawari/ChronQC
* Estimate library complexity: http://smithlabresearch.org/software/preseq/
* Protocol-specific pre-processing: https://github.com/lh3/pre-pe

### ChIP-seq

* Pre-processing
    * Quality assessment: https://github.com/rnakato/SSP
    * Allocation of multi-mapping reads: https://github.com/keleslab/permseq
* Peak calling
    * GC-aware peak caller: https://bioconductor.org/packages/devel/bioc/html/gcapc.html
    * GenoGAM peak caller: https://master.bioconductor.org/packages/3.3/bioc/html/GenoGAM.html
    * De-noising: https://github.com/kundajelab/TF_chipseq_pipeline
    * Peak discretizer (merging of replicates): https://github.com/nanakiksc/zerone
    * Compute error: https://github.com/tdhock/PeakError
    * Specifically for histone modifications: https://github.com/Bohdan-Khomtchouk/SUPERmerge
    * https://github.com/vidarmehr/ChIPWig
* hiddenDomains: https://sourceforge.net/projects/hiddendomains/
* Network-based identification of relationships among ChIP-seq data sets: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0925-0
* Motif assessment: http://www.bioinf.ict.ru.ac.za/
* https://github.com/cfce/chilin
* Web-based tool to compute enrichment at a variety of genomic features: http://liulab.dfci.harvard.edu/CEAS/
* Database of labeled ChIP-seq peaks: http://cbio.ensmp.fr/thocking/chip-seq-chunk-db/ (error of peak calls computed using https://github.com/tdhock/PeakError)
* EM algorithm for cooperatively bound TFs: https://github.com/vishakad/cpi-em
* Find different modes of binding: https://narlikarlab.github.io/DIVERSITY/
* Functional analysis
    * http://bioconductor.org/packages/FunChIP
    * https://github.com/XingTang2014/annoPeak
* Pipelines
    * https://github.com/milospjanic/ChIPSeqFPro
* Shape motifs: https://github.com/h-samee/shape-motif
* K-mer-based alternative to PWMs for predicting TF binding sites: http://groups.csail.mit.edu/cgs/gem/kmac/

### Chromatin accessibility

* DNase footprinting: https://github.com/ajank/Romulus
* HINT: http://costalab.org/publications-2/dh-hmm/ (was best out of 10 compared tools in recent NatMeth paper)
* ALTRE: https://mathelab.github.io/ALTRE/vignette.html
* Predict TF binding affinities using open chromatin + PWMs: https://github.com/schulzlab/TEPIC
* LR-DNAse: TFBS prediction using features derived from DNase-seq: http://biorxiv.org/content/early/2016/10/24/082594
* Identify accessible chromatin from NOMe-seq https://sourceforge.net/projects/came/
* Nucleotide-specific bias adjustment: https://github.com/txje/sequence-bias-adjustment
* ATAC-seq peak caller: https://github.com/LiuLabUB/HMMRATAC
* esATAC: ATAC-Seq pipeline https://bioconductor.org/packages/release/bioc/html/esATAC.html
* Predict targets of accessible regions: https://github.com/cole-trapnell-lab/cicero-release
* HINT-ATAC: TF footprinting from ATAC-Seq data http://www.regulatory-genomics.org/hint/introduction/
* ATAC-Seq QC: https://bioconductor.org/packages/release/bioc/html/ATACseqQC.html

### Chromatin Interactions

* Model 3D chromosome structure from Hi-C contact maps + optional FISH constraints: https://github.com/yjzhang/FISH_MDS.jl, https://github.com/yjzhang/3DC-Browser
* Predict enhancer targets: https://github.com/shwhalen/targetfinder
* Predicting TADs from histone modifications: https://cb.utdallas.edu/CITD/index.htm#ajax=home

### DNA

* Filtering
    * Removing redundant reads from deep sequencing data: https://git.informatik.uni-kiel.de/axw/Bignorm
* Error correction
    * https://github.com/lh3/bfc
    * RECKONER http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=reckoner&subpage=about
    * MultiRes - says it's for viral populations; not sure if applicable to humans: https://github.com/raunaq-m/MultiRes
* Duplicate removal
    * https://sourceforge.net/projects/pardre/
* Alignment
    * Recommendations:
	* https://github.com/CCDG/Pipeline-Standardization
	* These are optimized for germline variant calling, however there are a few interesting points
	* BWA MEM has a hidden option that enables "deterministic alignment results" (-K 100000000)
	* They recommend binning of quality scores for greater compression
    	    * Scores 2-6 correspond to Illumina error codes and are left as-is
            * 10, 20, 30
            * Round scores to the nearest value in probability space
        * The have some recommendations for standardizing headers
    * AGA: aligner specifically for coding sequences that takes into account the translated protein sequence: https://github.com/emweb/aga
    * Scaling alignment to massively parallel processors
        * https://www.biorxiv.org/content/early/2017/10/24/205328
        * Parallelizing hundreds of threads on Xeon Phi (64-72 core processors)
    * Alignment speed-up using PIM technologies: https://arxiv.org/abs/1711.01177
    * Whisper
        * Accelerates alignment by pre-sorting reads by sequence, and only aligning non-duplicates
        * http://sun.aei.polsl.pl/REFRESH/Whisper
        * Paper: https://www.biorxiv.org/content/early/2017/12/28/240358?rss=1&utm_source=dlvr.it&utm_medium=twitter
    * Other
        * Optical computing: https://www.optalysys.com/
        * FPGA: http://edicogenome.com/dragen-bioit-platform/
        * Selecting variants for graph genome: https://www.biorxiv.org/content/early/2018/04/30/311720
    * Align simultaneously against multiple reference genomes http://1001genomes.org/software/genomemapper.html
    * Compressed reference-based alignment: http://groups.csail.mit.edu/cb/cora/
    * Compression and querying of aligned haplotype data: https://github.com/richarddurbin/pbwt
    * Graph-based (mainly for local realignment): https://github.com/ekg/glia
    * https://github.com/ocxtal/libgaba
    * Long read
        * https://github.com/MohammadJRS/SureMap
	* https://github.com/lh3/minimap2
	    * Python binding: https://pypi.python.org/pypi/mappy
	* https://github.com/ocxtal/minialign
	* Graph-based: https://github.com/isovic/graphmap
	* https://github.com/xiaochuanle/MECAT
    * Correct sex chromosome-related errors/artifacts: https://github.com/WilsonSayresLab/XYalign
    * Lossy compression of quality scores: https://github.com/jkbonfield/crumble
* Assembly
    * Build de Bruijn Graph from multiple genomes: https://github.com/medvedevgroup/TwoPaCo
    * Align to a de Brujn graph: https://github.com/Malfoy/BGREAT
    * Cosmo/VARI: Assembler using succinct colored de Bruijn graphs to encode population information https://github.com/cosmo-team/cosmo/tree/VARI
    * Streaming: https://github.com/Shamir-Lab/Faucet
    * https://github.com/pmelsted/bfgraph
* Post-alignment
    * Base quality recalibration: https://github.com/swainechen/lacer
    * Indel realignment: https://github.com/mozack/abra
    * Error correction
	* https://github.com/isovic/racon
	* https://github.com/mikessh/mageri
	* https://github.com/gringer/bioinfscripts/blob/master/readstomper.pl
	* https://github.com/morinlab/ProDuSe/
	* https://www.genomics.agilent.com/en/NGS-Data-Analysis-Software/AGeNT/?cid=AG-PT-154&tabId=prod2570007 
* Variant calling
    * Matching variant sets: https://github.com/medvedevgroup/varmatch
    * Post-processing variant calls to determine whether variants at regions with alternative loci have allele(s) from an alternate locus: https://github.com/charite/asdpex
    * Filtering low-frequency variants that likely result from DNA damage: https://github.com/eilslabs/DKFZBiasFilter
    * Reference-free:
        * https://github.com/dib-lab/kevlar
        * https://github.com/jandrewrfarrell/RUFUS
    * DeepVariant is now open source
	* https://github.com/google/deepvariant
	* https://blog.dnanexus.com/2017-12-05-evaluating-deepvariant-googles-machine-learning-variant-caller/?utm_source=email&utm_medium=mailchimp&utm_campaign=december_newsletter
	* https://blog.dnanexus.com/2018-01-16-evaluating-the-performance-of-ngs-pipelines-on-noisy-wgs-data/
	* Currently it doesn't handle somatic calling 
    * Deep learning will increase accuracy in the nest GATK release: https://gatkforums.broadinstitute.org/gatk/discussion/10996/deep-learning-in-gatk4
    * Benchmarking
	* GA4GH/PrecisionFDA guidelines and tools for germline variants
	* https://www.biorxiv.org/content/early/2018/02/23/270157
	* https://github.com/ga4gh/benchmarking-tools
	* https://platform.dnanexus.com/login?scope=%7B%22full%22%3A+true%7D&redirect_uri=https%3A%2F%2Fprecision.fda.gov%2Freturn_from_login&client_id=precision_fda_gov
	* Using trios: https://github.com/sbg/geck
    * Ococo
    	* Variant calling is performed in a streaming fashion, with evidence for variation being updated from each read as it is mapped. Theoretically, new evidence can be ignored once the confidence level reaches a certain threshold. In practice, they show this is about 10% of reads. However, they only deal with SNPs. It would be interesting to investigate whether:
	* This can be extended to somatic variant calling
	* This can work with rare variants (eg cfDNA)
	* https://github.com/karel-brinda/ococo
	* Paper: https://arxiv.org/pdf/1712.01146.pdf
    * Distributed variant calling using Spark: http://bdgenomics.org/projects/
    * De-novo variant caller: https://github.com/bgm-cwg/novoCaller
* Genotyping
    * https://cran.r-project.org/web/packages/ebGenotyping/ebGenotyping.pdf
    * http://bioinfo.ut.ee/FastGT/
    * Compression of genotype data: http://sun.aei.polsl.pl/REFRESH/gtc
    * Base quality recalibration: https://github.com/swainechen/lacer
    * Mosaic variants: https://github.com/abyzovlab/Leucippus
    * HISAT-genotype: http://ccb.jhu.edu/hisat-genotype/index.php/Main_Page
    * Trio-based concordance: https://github.com/sbg/VBT-TrioAnalysis
    * Using a graph reference: https://github.com/bioinformatics-centre/BayesTyper
* SVs/CNV calling
    * Score SVs based on predicted functional impact https://github.com/lganel/SVScore
    * Pipelines for CNV calling:
	* https://www.biorxiv.org/content/early/2017/09/28/195230
	* https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1833-3
	* https://github.com/yuchaojiang/CODEX2
    * STIX: A tool for determining the frequency of SVs in population data: https://github.com/ryanlayer/stix
    * Tool for online review of SVs: https://github.com/jbelyeu/SV-plaudit
    * Using off-target reads to call CNVs: https://ep70.eventpilotadmin.com/web/page.php?page=IntHtml&project=ASHG17&id=170120230
    * GRIDSS: detects rearrangements using du Bruijn graph assembly
    * https://github.com/PapenfussLab/gridss
    * https://github.com/broadinstitute/ichorCNA
    * https://github.com/brentp/indelope
    * https://github.com/cui-lab/multigems
    * WHAM: CNV caller https://github.com/zeeev/wham
    * Identification of mosic events: https://github.com/asifrim/mrmosaic
        * This may be useful for calling allele-specific (i.e. heterozygous) copy number variation: https://drive.google.com/open?id=0B7E7HSjQ-SumQmlPc1Z0aUR5Sk0
    * Benchmarking of pipelines for calling complex variants: https://www.biorxiv.org/content/early/2017/11/23/218529
    * Manta and novoBreak performed best in DREME challenge for SV calling: https://www.biorxiv.org/content/early/2017/11/25/224733
    * https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4344483/
    * https://academic.oup.com/nar/article/43/14/e90/2902627
    * http://www.pnas.org/content/107/39/16910.full
    * http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075350
    * https://genomebiology.biomedcentral.com/articles/10.1186/gb-2011-12-10-r108
    * https://www.biorxiv.org/content/early/2017/10/10/119636
    * Mappability tracks, for normalizing CNV calls http://genome-test.soe.ucsc.edu/cgi-bin/hgTrackUi?hgsid=392476445_0mwMQdJgiLvv8WtNo4YktFt41IkI&c=chr1&g=umap
    * http://www.pnas.org/content/115/8/1871.short
    * https://github.com/brentp/smoove
    * https://github.com/walaj/svaba
    * Consensus caller: https://github.com/TheJacksonLaboratory/SVE
    * Randomized algorithm to speed up CNV calling: https://sourceforge.net/projects/genseng/
    * CNV calling from plasma: https://www.biorxiv.org/content/early/2018/03/28/290171
    * CODEX2: https://github.com/yuchaojiang/CODEX2
    * Ensemble caller: https://github.com/HaoKeLab/ensembleCNV
* STR
    * Pipeline to estimate mutational parameters for every STR in the human genome: https://www.nature.com/ng/journal/v49/n10/full/ng.3952.html
    * MIRMMR: MSI caller
	* https://github.com/ding-lab/MIRMMR
	* Uses methylation and mutation information
    * Targeted STR profiling: http://darwin.informatics.indiana.edu/str/
    * Krait: integrated microsatellite detection and primer design
	* https://github.com/lmdu/krait
	* https://www.biorxiv.org/content/early/2017/11/18/221754
    * MANTIS: https://github.com/OSU-SRLab/MANTIS
    * MSIsensor: https://github.com/ding-lab/msisensor
	* https://www.biorxiv.org/content/early/2018/01/10/246108
	* They generate a reference STR database using PacBio
	* They provide software to estimate STR number from short-read data using the reference database
    * This could be used to quickly extract STR-containing sequences from a BAM: https://github.com/rkmlab/perf 
    * https://github.com/mehrdadbakhtiari/adVNTR
    * https://github.com/mcfrith/tandem-genotypes
    * http://melissagymrek.com/lobstr-code/
    * Genotyping STR longer than read length: https://github.com/gymreklab/GangSTR
* Variant filtering
    * Using deep learning model to filter false-positive variants
	* https://www.biorxiv.org/content/biorxiv/early/2017/09/11/149146.full.pdf
	* Pre-trained (using GIAB) Illumina and IonTorrent models are here: https://github.com/gedoardo83/GARFIELD-NGS
    * ISOWN
	* Distinguishes somatic from germline variants
	* https://github.com/ikalatskaya/ISOWN
	* Paper: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-017-0446-9
    * Bystro: annotation and filtering (web only)
	* https://bystro.io/
	* https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1387-3
    * SevenBridges': https://github.com/sbg/smart-variant-filtering
    * https://github.com/aquaskyline/Skyhawk
    * This paper describes criteria that identify 100% of FP (as determined by Sanger sequencing) in a large clinical dataset https://www.biorxiv.org/content/early/2018/05/31/335950?rss=1
* Repeat calling
    * REPdenovo: https://github.com/Reedwarbler/REPdenovo
* Pipelines
    * SpeedSeq: alignment/annotation pipeline - https://github.com/hall-lab/speedseq
    * GenomeVIP: cloud-based pipeline https://github.com/ding-lab/GenomeVIP/
* Ancestry and kinship analysis
    * http://illumina.github.io/akt/
    * http://www.geenivaramu.ee/en/tools/mixfit
    * Genotype clustering and ethnicity prediction using Deep Belief Networks: https://arxiv.org/abs/1805.12218
* Phasing/Haplotyping
    * Eagle2: https://data.broadinstitute.org/alkesgroup/Eagle/
    * Using HiC+partial haplotypes: https://github.com/YakhiniGroup/SpectraPh
    * PhaseME http://beehive.cs.princeton.edu/wiki/phaseme/
    * HapCut2 (unclear if this works with standard Illumina WGS) https://github.com/vibansal/HapCUT2
    * hap.py: Illumina tool for identifying local haplotypes
        * https://github.com/Illumina/hap.py/blob/master/doc/happy.md
    * This might be useful in multiBP merging and in detecting contamination
        * For single individuals: https://github.com/jcna99/PEATH
        * https://arxiv.org/abs/1801.09864
    * Alignment-free:
        * Uses kmer frequency: https://github.com/paudano/kestrel
    * IVC
        * https://github.com/namsyvo/IVC
        * Integrates alignment and variant calling
        * Uses known variant information
        * Germline, Illumina paired-end only
    * Using decision tree: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2147-9
    * https://github.com/jcna99/PEATH
* Graph
    * https://gitlab.com/dvalenzu/PanVC
    * Graph genome alignment and variant calling on SevenBridges: https://www.sevenbridges.com/graph/
* Full Pipelines
    * superFreq: https://github.com/ChristofferFlensburg/superFreq
    * GROM: https://osf.io/6rtws/
    * Amplicon
        * Canary: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1950-z
        * https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1985-1
    * DNAScan
        * https://github.com/KHP-Informatics/DNAscan
        * Align with HiSAT2; realign soft-clipped and unaligned reads with BWA MEM due to poor performance of indel calling from HISAT2-aligned data
        * Variant calling with Freebayes, extract indel positions and realign with GATK HC, perform hard filtering
        * Manta and Expansion Hunter for SVs
* Non-Illumina
    * Adapter between aligners and downstream tools for handling PacBio quality values
        * https://github.com/mchaisso/pbsamstream
    * BGI compared to Illumina: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0190264
    * Nanopore
        * https://github.com/nanoporetech/wub
	* Detect analytes from raw event data: https://github.com/nanoporetech/mako
* Other
    * VCF compression and data extraction: https://github.com/kedartatwawadi/GTRAC
    * Run length encoded multi-sample BWT + server: https://github.com/wtsi-svi/ReadServer
    * Amplicon
	* Pisces: https://github.com/Illumina/Pisces/

### Metagenomics

* Alignment-free functional binning and abundance estimation: https://github.com/snz20/carnelian
* Map reads against redundant databases: https://bitbucket.org/genomicepidemiology/kma
* Taxonomic classificatino using pseudoalignment: https://github.com/mreppell/Karp

### Methylation

* BS-SNPer: fast SNP calling from bisulfite-converted sequencing reads https://github.com/hellbelly/BS-Snper
* MACAU: Mixed-model association http://www.xzlab.org/software.html
* ME-plot: Error detection and correction in bisulfite-converted sequencing reads https://github.com/joshuabhk/methylsuite
* Metilene: Calling differential methylation http://www.bioinf.uni-leipzig.de/Software/metilene/
* Reference-free bisulfite sequence comparison: https://github.com/thomasvangurp/epiGBS
* Correction for cell-type composition: http://www.cs.tau.ac.il/~heran/cozygene/software/refactor.html
* Predicting gene expression from methylation: http://arxiv.org/abs/1603.08386
* GEM: R package for meQTL and EWAS https://bioconductor.org/packages/devel/bioc/html/GEM.html
* Call CNVs from methylation array data: https://github.com/mknoll/cnAnalysis450k
* Compression tool for bedMethyl files: https://github.com/jianhao2016/METHCOMP
* https://github.com/RWilton/Arioc

### MNase-seq

* Alternative/differential nuleosome positioning: https://github.com/airoldilab/cplate

### Nanopore

* https://github.com/roblanf/minion_qc
* Clinical application: https://www.nature.com/articles/s41467-017-01343-4
* Signal-level barcode demux: https://github.com/rrwick/Deepbinner

### RNA

* Annotations
    * Reassembly and annotation of public data for multi-tissue transcriptome map: http://big.hanyang.ac.kr/CAFE
* Splice junctions
    * Set of novel (i.e. missing from annotation databases) splice junctions identified from SRA datasets: https://github.com/nellore/intropolis/blob/master/README.md
    * Complete sets of splice junctions in public RNA-seq datasets: http://snaptron.cs.jhu.edu/data/	
* QC
    * NOISeq - exploratory analysis of read mappings https://www.bioconductor.org/packages/release/bioc/html/NOISeq.html
    * AuPairWise: determine replicability without replicates https://github.com/sarbal/AuPairWise
    * dupRadar: duplications https://www.bioconductor.org/packages/release/bioc/html/dupRadar.html
    * Correcting for RNA quality
* Align/quantify:
    * Need to re-evaluate HiSat2 + StringTie pipeline https://github.com/gpertea/stringtie
    * RapMap: stand-alone lightweight alignment library: https://github.com/COMBINE-lab/RapMap/tree/SAQuasiAlignment
    * Kallisto is a fast and accurate method for transcript quantification https://pachterlab.github.io/kallisto/
        * Sleuth is a companion R package for differential expression analysis http://pachterlab.github.io/sleuth/
        * Different models can be used in Sleuth to, for example, perform time-course experiments http://nxn.se/post/134227694720/timecourse-analysis-with-sleuth
    * tximport: R package for aggregating transcript-level quantifications for gene-level analysis: http://f1000research.com/articles/4-1521/v1
    * Wasabi: prepare Salmon/Sailfish output for Sleuth https://github.com/COMBINE-lab/wasabi
    * featureCounts: read summarization http://bioinf.wehi.edu.au/featureCounts/
    * D-GEX: Quantification of whole-transcriptome gene expression from landmark genes https://github.com/uci-cbcl/D-GEX
    * Faster version of HTSeq/featureCount: https://github.com/qinzhu/VERSE
    * Quantification using both structure and abundance information: https://pypi.python.org/pypi/rsq
    * Improve transcript quantification by integrating PolII ChIP-seq data: https://github.com/pliu55/RSEM/tree/pRSEM
    * Hera: simultaneous alignment, quantification, and fusion detection https://github.com/bioturing/hera
    * Aligner calibraiton: https://bitbucket.org/irenerodriguez/fbb
    * Quantification using k-mers unique to each gene: https://github.com/informationsea/Matataki
* Correction/Normalization
    * Choosing normalization methods: https://arxiv.org/abs/1609.00959
    * TDM: cross-platform normalization https://github.com/greenelab/TDMresults
    * alpine: corrects for fragment sequence bias https://github.com/mikelove/alpine/blob/master/vignettes/alpine.Rmd
    * Filter out lowly-expressed transcripts prior to quantification decreases FP rate: http://www.genomebiology.com/2016/17/1/12
    * Partition variance between biological and technical sources: https://www.bioconductor.org/packages/3.3/bioc/vignettes/variancePartition
    * Quantify and correct for uncertainty in abundance estimates: https://github.com/PSI-Lab/BENTO-Seq
    * Simultaneous isoform discovery and quantification across multiple samples: http://cbio.ensmp.fr/flipflop
    * R package to compare normalization methods: https://github.com/Edert/NVT
    * Filtering and tissue-aware normalization: http://bioconductor.org/packages/release/bioc/html/yarn.html
    * Bias correction for transcript abundance estimation: https://www.lexogen.com/mix-square-scientific-license/
    * Replacement for htseq-counts/featureCounts that handles multi-mapping reads: https://bitbucket.org/mzytnicki/multi-mapping-counter
    * https://github.com/wgmao/DataRemix
    * Junction compatibility score for idenitifying genes whose isoform-level abundance estimates differ significantly from expectation: https://github.com/csoneson/annotation_problem_txabundance
* Workflows:
    * Artemis (RNA-Seq workflow designed around Kallisto): https://github.com/RamsinghLab/artemis
    * https://github.com/ririzarr/rafalib
    * Isolator: https://github.com/dcjones/isolator
    * https://bioinform.github.io/rnacocktail/
* eQTL
    * Multi-tissue:
        * HT-eQTL https://github.com/reagan0323/MT-eQTL
        * MT-HESS: eQTL analysis across tissues http://www.mrc-bsu.cam.ac.uk/software/
        * eQTLBMA: cross-tissue eQTL https://github.com/timflutre/eqtlbma
        * Multi-tissue eQTL: https://cran.r-project.org/web/packages/JAGUAR/index.html
	* Multi-tissue, polygenic TWAS: https://github.com/ypark/fqtl
    * Quantile regression approach: https://xiaoyusong.shinyapps.io/QRBT/
    * Using prior knowledge: https://github.com/redsofa/LassoMP
    * CONDOR: simultaneous cis- and trans-eQTL analysis https://github.com/jplatig/condor
    * log allelic fold change (aFC) to quanitfy effect sizes of eQTL: http://biorxiv.org/content/biorxiv/early/2016/09/30/078717.full.pdf
    * Identify cis mediators of trans-eQTL: http://biorxiv.org/content/early/2016/09/30/078683
    * https://cran.r-project.org/web/packages/QRank/
* Differential expression:
    * cjBitSeq: https://github.com/mqbssppe/cjBitSeq/wiki
    * Differential junction usage: https://github.com/hartleys/JunctionSeq
    * TROM: comparison of transcriptomes between species (and maybe between cell/tissue types?) https://cran.r-project.org/web/packages/TROM/index.html
    * Tissue specificity of genes, based on GTEx data: http://genetics.wustl.edu/jdlab/tsea/
    * Informative priors for Bayesian differential expression analysis using historical data: https://github.com/benliemory/IPBT
    * RNA-enrich: http://lrpath.ncibi.org/
    * Diferentially expressed region finder (also works for ChIP-seq peaks): www.bioconductor.org/packages/derfinder
    * Differentially expressed pathways using kernel MMD: https://eib.stat.ub.edu/tiki-index.php?page_ref_id=73
    * Method using dimension-reduced ANOVA: http://homepage.fudan.edu.cn/zhangh/softwares/multiDE/
    * Correct for hidden sources of variation: https://github.com/sutigit21/SVAPLSseq
    * Alternative method for estimating variances: https://github.com/mengyin/vashr
    * With few replicates: https://figshare.com/s/963e895f812d6f06468a
    * Network sub-pathways: http://bioconductor.org/packages/release/bioc/html/DEsubs.html
    * https://github.com/ewyang089/SDEAP/wiki
    * Local subnetworks enriched for DE genes: https://cran.rstudio.com/web/packages/LEANR/index.html
    * https://github.com/beiyuanzhe/DiscriminantCut
    * Stage-wise DE/DT https://github.com/statOmics/stageR
    * Computing heritability of gene expression: https://cran.r-project.org/web/packages/HeritSeq/index.html
    * https://github.com/bee-hive/BIISQ
    * Power analysis: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1648-2
    * Using a beta binomial model with dynamic correction for overdispersion: https://github.com/GuoshuaiCai/BBDG
    * Kmer-based: https://github.com/Transipedia/dekupl
* Differential transcript usage
    * https://github.com/bartongroup/Rats
    * Bayesian extenstion to BitSeq for differential transcript usage: https://github.com/mqbssppe/cjBitSeq
* Co-expression/networks
    * http://biorxiv.org/content/early/2016/10/02/078741
    * BicMix: differential co-expression networks http://beehive.cs.princeton.edu/software/
* ASE
    * GeniASE: ASE without haplotypes http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758070/pdf/srep21134.pdf
    * Without phasing: http://lifecenter.sgst.cn/cisASE/
    * ASElux: https://drive.google.com/open?id=0B7E7HSjQ-SumQmlPc1Z0aUR5Sk0
* Splicing
    * https://github.com/hartleys/JunctionSeq
    * https://github.com/lkmklsmn/SplicER
    * Annotation of splicing types https://r-forge.r-project.org/projects/splicingtypes/
    * MAJIQ: detection of local splice variation http://majiq.biociphers.org/
    * https://github.com/davidaknowles/leafcutter
    * http://www.mhs.biol.ethz.ch/research/krek/jsplice.html
    * Splice site prediction: http://cabgrid.res.in:8080/HSplice/
    * DRIM-seq: http://bioconductor.org/packages/DRIMSeq
    * Fast quantification of differential splicing: https://github.com/comprna/SUPPA
    * Identify variant associated with splicing: https://sourceforge.net/projects/isvase/
    * Prediction of intronic splice branchpoints: https://github.com/betsig/branchpointer/
    * Proportion spliced index: https://github.com/comprna/Junckey
    * https://github.com/timbitz/Whippet.jl
    * http://github.com/henryxushi/SparseIso
    * https://github.com/tigerhu15/PennDiff
    * https://github.com/abl0719/ASElux
    * F1000/Bioconductor workflow on DRIMSeq and DEXSeq for DTU: https://f1000research.com/articles/7-952/v1
* Assembly:
    * CIDANE: http://ccb.jhu.edu/software/cidane/
    * transrate: evaluation of de novo assemblies http://hibberdlab.com/transrate/
    * dammit: annotator for de novo assemblies http://dammit.readthedocs.org/en/latest/
    * kma: detection of differential intron retention https://github.com/pachterlab/kma
    * RapClust: lightweight clustering of de novo transcriptomes https://pypi.python.org/pypi/rapclust/0.1
    * Shannon http://sreeramkannan.github.io/Shannon/
    * Strawberry https://github.com/ruolin/Strawberry
    * CLASS2: http://ccb.jhu.edu/people/florea/research/CLASS2/
    * https://sourceforge.net/projects/transcriptomeassembly/files/
    * Multi-sample transcriptome assembly: http://tacorna.github.io/
    * Clustering to decontaminate de novo assemblies: https://github.com/Lafond-LapalmeJ/MCSC_Decontamination
    * Assembly from unstranded data: http://big.hanyang.ac.kr/CAFE
    * https://github.com/Kingsford-Group/scallop
    * Identification of transcript boundaries: https://github.com/realbigws/DeepBound
    * Corset: gene counts from a transcriptome assembly https://github.com/Oshlack/Corset/wiki
    * Consensus method: https://github.com/macmanes-lab/Oyster_River_Protocol
    * https://github.com/Oshlack/necklace
    * https://github.com/ruolin/strawberry
    * https://github.com/bee-hive/BIISQ
* Time series
    * http://diceseq.sourceforge.net/
    * https://bioconductor.org/packages/release/bioc/html/ctsGE.html
    * http://bioinfo.ibp.ac.cn/Dynamic-BM/
    * Differential isoform usage over time http://bioconductor.org/packages/devel/bioc/html/maSigPro.html
* Deconvolution:
    * VoCAL: https://cran.r-project.org/web/packages/ComICS/index.html
    * DeconRNASeq: deconvolute expression profiles in mixed tissues http://www.bioconductor.org/packages/2.12/bioc/vignettes/DeconRNASeq/inst/doc/DeconRNASeq.pdf
* Search
    * http://www.cs.cmu.edu/∼ckingsf/software/bloomtree/
    * SBTs https://github.com/medvedevgroup/bloomtree-allsome
* Structural variation
    * Squid: https://github.com/Kingsford-Group/squid
    * Identify gene fusions: https://github.com/ndaniel/fusioncatcher
    * Fusion genes: http://star-fusion.github.io
    * Identify gene expression driven by copy number alteration in samples with matched RNA-seq and CNA data: https://www.bioconductor.org/packages/release/bioc/html/iGC.html
    * Fusion identification using kallisto: https://github.com/pmelsted/pizzly
* Other
    * Biclustering for gene co-expression analysis: http://bioconductor.org/packages/devel/bioc/html/QUBIC.html
    * Sample size calculation for experimental design: https://cran.r-project.org/web/packages/ssizeRNA/index.html
    * Variance estimation: http://github.com/mengyin/vashr
    * Sample expression "admixture" (can also be used for deconvolution): https://www.bioconductor.org/packages/release/bioc/html/CountClust.html
        * Essentially, this assigns samples to clusters based on similarity in expression of sets of genes determined be be most discriminating. Each sample can belong to multiple clusters (similar to an admixture analysis).
    * Identification of regulatory networks: https://sites.google.com/a/fleming.gr/rnea/
    * Predict ribosome footprint from transcripts https://sourceforge.net/projects/riboshape/
    * Phasing: https://github.com/secastel/phaser
    * Identifying the source of (almost) all RNA-seq reads: https://github.com/smangul1/rop/wiki
    * Subsampling to determine effect of read depth on downstream analyses: http://www.bio-complexity.com/samExploreR_1.0.0.tar.gz
    * Phenotype prediction: https://github.com/clabuzze/Phenotype-Prediction-Pipeline
    * Predict RNA-RNA interaction: https://github.com/satoken/ractip
    * Mitigate cell-cycle effects: http://www.nature.com/articles/srep33892
    * Fast computation of probabilities of pairwise regulation https://github.com/lingfeiwang/findr
    * Align against synthetic transcript-based reference: https://github.com/Oshlack/Lace
    * Interactive visualization http://bioconductor.org/packages/devel/bioc/html/Glimma.html
    * New factorization method for dimensionality reduction: https://github.com/brian-cleary/CS-SMAF
    * Database of public RNA-seq data sets processed using the same pipeline: https://github.com/mskcc/RNAseqDB
    * Evaluation of aligners on long reads (GMap performs best): http://biorxiv.org/content/early/2017/04/11/126656
    * Test different ML algorithms for classifying expression profiles: https://github.com/gboris/blkbox
    * Find regions of correlated expression: https://cran.r-project.org/web/packages/SegCorr/index.html

### Single-cell

* Comparative analysis of methods: 
    * http://www.sciencedirect.com/science/article/pii/S1097276517300497
    * https://www.nature.com/nmeth/journal/v14/n4/full/nmeth.4220.html
* Review of experimental design and analysis: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y
* Sean Davis' list: https://github.com/seandavi/awesome-single-cell
* Datasets:
    * scRNAseqDB: https://bioinfo.uth.edu/scrnaseqdb/
    * Analysis-ready datasets: http://imlspenticton.uzh.ch:3838/conquer/
* Platforms
    * Microwells: http://www.nature.com/articles/srep33883
* Simulation
    * https://github.com/YosefLab/SymSim
* QC:
    * http://www.morgridge.net/SinQC.html
    * https://github.com/YosefLab/scone
    * Parallel transcriptome/epigenome data: https://github.com/ChengchenZhao/DrSeq2
    * https://bitbucket.org/princessmaximacenter/Sharq
* Normalization
    * http://michelebusby.tumblr.com/post/130202229486/the-ks-test-looks-pretty-good-for-single-cell
    * Accounting for technical variation: http://www.nature.com/ncomms/2015/151022/ncomms9687/full/ncomms9687.html#supplementary-information
    * ZIFA: Zero-inflated factor analysis https://github.com/epierson9/ZIFA
    * http://biorxiv.org/content/early/2016/04/22/049734.full.pdf+html
        * https://www.dropbox.com/s/pno78mmlj0exv7s/NODES_0.0.0.9010.tar.gz?dl=0
    * scran: http://bioconductor.org/packages/devel/bioc/html/scran.html
    * Correct for expression heterogeneity: https://github.com/PMBio/scLVM
    * Comparison of normalization methods: http://biorxiv.org/content/biorxiv/early/2016/07/17/064329.full.pdf
    * https://github.com/rhondabacher/SCnorm/tree/master/R
    * Factor analysis: https://github.com/UcarLab/IA-SVA/
    * qSVA corrects for RNA quality (in the SVA package)
    * ERCC https://bitbucket.org/bsblabludwig/bearscc
    * Weighting method that enables bulk RNA-seq tools to be used with single-cell: https://www.ncbi.nlm.nih.gov/pubmed/29478411
* Gene/Transcript counting
    * Modified version of Kallisto: https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts
    * DISCO: https://pbtech-vc.med.cornell.edu/git/mason-lab/disco/tree/master
    * http://garberlab.umassmed.edu/software/esat/
* Cell type-specific expression
    * http://bioconductor.org/packages/release/bioc/html/CellMapper.html
* Clustering
    * Comparative analysis:
        * http://biorxiv.org/content/early/2016/04/07/047613
        * https://github.com/epurdom/clusterExperiment
    * SC3: consensus clustering https://github.com/hemberg-lab/sc3
    * destiny: diffusion maps for single-cell data http://bioconductor.org/packages/release/bioc/html/destiny.html
    * https://github.com/govinda-kamath/clustering_on_transcript_compatibility_counts
    * GiniClust https://github.com/lanjiangboston/GiniClust
    * pcaReduce: https://github.com/JustinaZ/pcaReduce
    * https://github.com/BatzoglouLabSU/SIMLR
    * CIDR: https://github.com/VCCRI/CIDR
    * Vortex: http://web.stanford.edu/~samusik/vortex/
    * Identify rare cell types: RaceID http://www.nature.com/nature/journal/v525/n7568/full/nature14966.html
    * https://github.com/VCCRI/CIDR
    * https://github.com/drisso/zinbwave
    * http://www.pitt.edu/~wec47/singlecell.html
    * Neural networks for dimensionality reduction and clustering http://sb.cs.cmu.edu/scnn/
    * https://github.com/srmcc/dcss_single_cell
    * Cell similarity measure: https://github.com/maggiecrow/MetaNeighbor
* Differential Expression
    * Monocle cole-trapnell-lab.github.io/monocle-release/ (2.0 has Census algorithm for differential transcript analysis)
    * scDD: https://github.com/kdkorthauer/scDD
    * ISOP: comparison of isoform pairs in single cells https://github.com/nghiavtr/ISOP
    * D3E: http://hemberg-lab.github.io/D3E/
    * BASiCS: https://github.com/catavallejos/BASiCS
    * Beta Poisson: https://github.com/nghiavtr/BPSC
    * Zero-inflation correct enables use of DESeq2, etc w/ single cell data: https://github.com/statOmics/zingeR
    * DESingle: https://github.com/miaozhun/DEsingle
    * https://github.com/willtownes/vamf-paper
* Allele-specific expression
    * SCALE accounts for "burstiness" of transcription: https://github.com/yuchaojiang/SCALE
    * https://github.com/churchill-lab/scBASE
* Splicing
    * https://github.com/huangyh09/brie    
* Time-series/ordering/lineage prediction
    * Monocle
        * Analysis of pseudotime uncertainty: http://biorxiv.org/content/biorxiv/early/2016/04/05/047365.full.pdf
    * ECLAIR: cell lineage prediction https://github.com/GGiecold/ECLAIR
    * Identification of ordering effects: https://github.com/lengning/OEFinder
    * Slicer: non-linear trajectories https://github.com/jw156605/SLICER
    * Wishbone: identification of bifurcations in developmental trajectories http://www.c2b2.columbia.edu/danapeerlab/html/cyt-download.html
    * SCOUP: https://github.com/hmatsu1226/SCOUP
    * Ouija: https://github.com/kieranrcampbell/ouija
    * http://bioconductor.org/packages/release/bioc/html/sincell.html
    * https://github.com/kstreet13/slingshot
    * Cytoscape plugin: http://cytospade.org/
    * https://github.com/zji90/TSCAN
    * https://github.com/dimenwarper/scimitar
    * https://cran.r-project.org/web/packages/timeSeq/index.html
    * http://bioconductor.org/packages/release/bioc/html/cellTree.html
    * Diffusion pseudiotime: http://www.helmholtz-muenchen.de/icb/research/groups/machine-learning/projects/dpt/index.html
    * https://github.com/theislab/kbranches
    * Construction co-expression networks: https://cran.r-project.org/web/packages/LEAP/index.html
    * Differentila expression between trajectories https://github.com/kieranrcampbell/switchde
    * FORKS: https://github.com/macsharma/FORKS
    * https://github.com/kieranrcampbell/phenopath
* Pipelines
    * Seurat http://www.satijalab.org/seurat.html
    * SINCERA https://research.cchmc.org/pbge/sincera.html
    * MAST: https://github.com/RGLab/MAST
    * scde (differential expression + gene set over-dispersion): https://github.com/hms-dbmi/scde
    * BaSiCs: Bayesian analysis of single cell data: https://github.com/catavallejos/BASiCS
    * FastProject: https://github.com/YosefLab/FastProject/wiki
    * Citrus: http://chenmengjie.github.io/Citrus/
    * Tools from Teichmann lab (cellity, celloline, scrnatb): https://github.com/Teichlab/
    * SCell: https://github.com/diazlab/SCell
    * http://bioconductor.org/packages/scater
    * https://github.com/joeburns06/hocuspocus
    * https://gitlab.com/uhcclxgg/granatum
    * https://github.com/LuyiTian/scPipe
    * For epigenetic data: https://zhiji.shinyapps.io/scrat/
    * scanpy: https://github.com/theislab/scanpy
* SNVs/CNVs
    * DNA SNV calling: https://bitbucket.org/hamimzafar/monovar
    * Ginko: analysis of CNVs in single-cell data: http://qb.cshl.edu/ginkgo/?q=/XWxZEerqqY477b9i4V8F
    * CNV calling: http://genome.cshlp.org/content/early/2016/01/15/gr.198937.115.full.pdf
    * Genotyping: https://bitbucket.org/aroth85/scg/wiki/Home
* Regulatory networks:
    * Gene co-expression: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004892
    * https://github.com/hmatsu1226/SCODE
    * https://github.com/aertslab/SCENIC
* Other scRNA-seq
    * Analysis of 3' tagging data: https://github.com/garber-lab/ESAT
    * StemID: Prediction of stem cells and lineage information https://github.com/dgrun/StemID
    * Phasing: https://github.com/edsgard/scphaser
    * Classification using sets of known cell type-specific genes: https://github.com/YosefLab/FastProject/wiki
    * UMI counting: https://github.com/vals/umis
    * Imputation of missing values
        * https://github.com/pkathail/magic/
	* https://github.com/Vivianstats/scImpute
	* https://cran.r-project.org/web/packages/DrImpute/index.html
    * Power analysis: https://github.com/vals/umis/
    * Pooled perturbation experiments: https://github.com/asncd/MIMOSCA
    * Simulation: http://bioconductor.org/packages/splatter/
    * Comparison across experiments: https://github.com/hemberg-lab/scmap
    * Demuxlet: using natural genetic variation to demultiple 10x data https://github.com/hyunminkang/apigenome
    * Alignment of multiple single-cell data sets: https://github.com/jw156605/MATCHER
    * Topological analysis: https://github.com/RabadanLab/scTDA
    * Deconvolution of bulk RNA from scRNA: https://github.com/xuranw/MuSiC
* Methylation
    * Prediction of missing information: https://github.com/cangermueller/deepcpg
* ATAC-seq
    * Infer TF variation: https://github.com/GreenleafLab/chromVAR
    * Clustering: https://github.com/timydaley/scABC

### Somatic

* QC
    * Bias detection
        * https://www.degruyter.com/downloadpdf/j/jib.2017.14.issue-3/jib-2017-0025/jib-2017-0025.pdf
        * This method is targeted to RNA-Seq, but should work for CNV calling from Exome data as well
    * Deconvolution of signal (allele frequency, gene expression, etc) from heterogeneous tissue data: https://github.com/tedroman/WSCUnmix
    * Remove strand bias artifacts (e.g. FFPE samples): https://github.com/mikdio/SOBDetector
* Variant Calling
    * Review of somatic variant callers: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5852328/
    * FMI abstracts on variant calling from tumor-only sequencing
    * http://cancerres.aacrjournals.org/content/74/19_Supplement/1893
    * http://ascopubs.org/doi/abs/10.1200/jco.2015.33.15_suppl.11084
    * https://github.com/SciLifeLab/Sarek
    * Comparison of 9 somatic variant callers
        * http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0151664
        * From 2016; should be updated with MuTect2 and Strelka2
    * Lancet: new somatic variant caller based on genome graphs
        * https://github.com/nygenome/lancet
        * Non-commercial license
    * Some recommended tools/pipelines for calling low-frequency variants:
	* https://github.com/tamilieberman/TB-diversity-across-organs
	* http://csb5.github.io/lofreq/
	* https://github.com/TGenNorth/ASAP
	* https://www.ncbi.nlm.nih.gov/pubmed/25201145
	* Goby uses deep learning for somatic variant calling: http://campagnelab.org/software/goby/
    * MC3/PanCan
        * https://www.synapse.org/#!Synapse:syn7214402/wiki/406010
        * Variant call merger: https://github.com/covingto/pancanmafmerge
    * Ensemble method with random forest integration of results from multiple callers:
	* https://github.com/skandlab/SMuRF
	* Paper: https://www.biorxiv.org/content/early/2018/02/23/270413
    * TNScope: Sentieon's pipeline
	* https://www.biorxiv.org/content/early/2018/01/19/250647
	* Includes machine learning-based variant filtering
    * xATLAS
	* https://github.com/jfarek/xatlas
	* Uses a logistic regression model based on truth sets (e.g. GIAB) to assign quality scores
* CNA/SV
    * Infer the tumor cell fraction of SV: https://github.com/mcmero/SVclone
    * A study of using shallow WGS on low-quality FFPE samples for CNV calling used QDNASeq for segmentation
	* https://bioconductor.org/packages/release/bioc/html/QDNAseq.html
	* Adjusts for CG content and mappability
	* Found 7M reads was optimal
	* https://www.biorxiv.org/content/early/2017/12/08/231480
    * Germline CNVs may be predictive of cancer susceptibility: https://www.biorxiv.org/content/early/2018/04/17/303339
    * https://github.com/yuchaojiang/MARATHON
    * Allele-specific SV calling: https://github.com/ma-compbio/Weaver
        * Phasing SVs with noisy data: http://sci-hub.tw/http://liebertpub.com/doi/pdf/10.1089/cmb.2018.0022
* Variant Filtering
    * https://github.com/ctDNA/TNER
    * This is a cool approach for constructing a ChIP-Seq control from integration of multiple public data sets: https://www.biorxiv.org/content/early/2018/03/08/278762. It strikes me that a similar approach could be used to generate matched controls for tumor-only samples.
* Tumor Purity
    * IchorCNA
	* Estimates tumor cell fraction from low-pass WGS; should probably be adaptable to deep targeted sequencing
	* https://github.com/broadinstitute/ichorCNA
    * Accurity (tumor-normal)
	* https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty043/4827681
	* Control for tumor purity differences in methylation analysis: https://www.biorxiv.org/content/early/2018/01/16/248781
    * PyLOH: https://github.com/uci-cbcl/PyLOH
    * Tumor Muational Burden
	* http://moat.gersteinlab.org/
* Variant Annotation
    * CHASM and SNV-Box
	* CHASM Predicts functional significance of somatic missense mutations
	* SNV-Box is a database of pre-computed features of all possible amino acid substitutions in the human genome
	* http://wiki.chasmsoftware.org/index.php/Main_Page
	*  Not free for commercial use, but it was developed at JHU and Vogelstein is a co-author on at least some of the publications
    * Orchid: framework for annotation and machine learning of cancer variants: https://github.com/wittelab/orchid
    * Discover potential cancer driver elements
	* https://github.com/smshuai/DriverPower
	* http://www.lagelab.org/projects/
	* http://bioinformatics.mdanderson.org/main/CanDrA#CanDrA
    * Predict cancer subtype using biological networks based on somatic variants: https://www.biorxiv.org/content/early/2017/12/03/228031
    * PCGR: Annotate cancer genomes with clinical information:
	* https://github.com/sigven/pcgr
	* Integrates many public databases
    * Benchmark of methods to identify pathogenic non-coding variation: https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty008/4798701?redirectedFrom=fulltext
    * http://bioinformatics.mdanderson.org/main/FASMIC
    * Database of TCGA alternative splicing events https://www.cell.com/cancer-cell/fulltext/S1535-6108(18)30306-4#secsectitle0080
    * Consensus driver gene identification pipeline and database of cancer driver genes in TCGA: https://www.cell.com/cell/fulltext/S0092-8674(18)30237-X
    * TCGA pathogenic signalling pathways: https://www.cell.com/cell/fulltext/S0092-8674(18)30359-3?cid=tw%26p
* Predict cancer type/signature
    * From mutations: https://www.nature.com/articles/nature12477
    * From SNV+SV calls: https://www.biorxiv.org/content/early/2018/02/18/267500
    * Machine learning method to predict cancer type from variant calls: https://www.biorxiv.org/content/early/2017/11/05/214494
    * Prediction of inherited susceptibility to 20 difference cancer types: http://www.pnas.org/content/115/6/1322.short
    * Computing polygenic risk scores and association with cancers: https://www.biorxiv.org/content/early/2017/10/19/205021
    * Predict driver mutations: https://github.com/KarchinLab/CHASMplus
* Clonality
    * https://cran.r-project.org/web/packages/QuantumClone/index.html
    * https://static-content.springer.com/esm/art%3A10.1038%2Fs41588-018-0128-6/MediaObjects/41588_2018_128_MOESM1_ESM.pdf
    * Benchmarking, plus description of CloneFinder: https://academic.oup.com/bioinformatics/advance-article-abstract/doi/10.1093/bioinformatics/bty469/5040314
* MSI
    * https://github.com/kzkedzierska/sonics
* Other
    * Bioconductor package with various methods for working with TCGA data: https://bioconductor.org/packages/release/bioc/html/TCGAbiolinks.html

### Integrated Methods

* Reviews:
    * http://www.biomedcentral.com/1752-0509/8/S2/I1
    * Comparison of methods http://www.biomedcentral.com/1752-0509/8/S2/S4
    * Focusing on integration of RNA-seq and ChIP-seq http://journal.frontiersin.org/article/10.3389/fcell.2014.00051/full
    * Network-based methods: http://rsif.royalsocietypublishing.org/content/12/112/20150571
* General-purpose multi-omics integration:
    * mixOmics: R package that implements several multivariate methods, including DIABLO http://mixomics.org/
    * Non-negative matrix factorization for integration of data sets https://github.com/yangzi4/iNMF
    * omicade4: Integration of multi-omics data using co-inertia analysis https://bioconductor.org/packages/release/bioc/html/omicade4.html
    * Integration of multi-omics data using random forests: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1043-4
    * Kernel-PCA: http://www.biomedcentral.com/1752-0509/8/S2/S6
    * Integrative analysis of multiple diverse omics datasets by sparse group multitask regression http://journal.frontiersin.org/article/10.3389/fcell.2014.00062/abstract
    * Joint bi-clustering of multiple data types: http://research.cs.aalto.fi/pml/software/GFAsparse/
    * Identifying covariance between sequencing data sets: http://github.com/pmb59/fCCAC/
    * https://github.com/davidvi/pypanda
    * SDA: integrate gene expression across multiple tissues, or multi-omics in a single tissue, for identification of trans-QTL networks: https://jmarchini.org/sda/
    * HMM for binary classification based on multivariate data: https://github.com/PetarV-/muxstep
    * https://github.com/fraenkel-lab/OmicsIntegrator
    * Correlation between enriched regions from different data sets: http://malone.bioquant.uni-heidelberg.de/software/mcore
    * https://cran.r-project.org/web/packages/r.jive/
    * SIFORM: http://bioinformatics.oxfordjournals.org/content/early/2016/07/03/bioinformatics.btw295.full
    * https://sourceforge.net/projects/epimine/
    * Deep learning-based framework for integrating multiple data types to predict another data type: https://github.com/ueser/FIDDLE
    * Significance-based https://www.bioconductor.org/packages/release/bioc/html/SMITE.html
    * Two-stage CCA identifies non-linear associations: https://github.com/kosyoshida/TSKCCA
    * HMM: https://link.springer.com/protocol/10.1007/978-1-4939-6753-7_10
    * https://github.com/KnowEnG/pgenmi
    * https://github.com/bioFAM/MOFA
    * https://github.com/jinghan1018/tensor_decomp
    * Unsupervised: https://github.com/bioFAM/MOFA
    * https://noble.gs.washington.edu/proj/avocado/
* Multi-tissue:
    * http://bioconductor.org/packages/release/bioc/html/HDTD.html
    * https://github.com/yuzhang123/IDEAS
* Specific data types:
    * Predict gene fusions from WGS and RNA-seq: http://sourceforge.net/p/integrate-fusion/wiki/Home/
    * NuChart: layer additional omics data on Hi-C ftp://fileserver.itb.cnr.it/nuchart/
    * Predict expression from H3K27, and identify cis-regulatory elements: http://cistrome.org/MARGE/
    * GenoSkyline: predict tissue-specific functional regions from epigenomic data http://genocanyon.med.yale.edu/GenoSkyline
* Network-based:
    * Merging networks: https://github.com/maxconway/SNFtool
    * https://sourceforge.net/projects/xmwas/
* Causality
    * Hybrid BN/CMI approach to constructing GRNs: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005024
* Noise/bias:
    * Bias correction across different assays on the same samples: https://cran.r-project.org/web/packages/MANCIE/
    * Cell cycle heterogeneity appears to be a minor contributor to noise; instead, library size is the largerst PC by far http://www.nature.com/nbt/journal/v34/n6/full/nbt.3498.html
* Patient/disease subtyping
    * http://biorxiv.org/content/biorxiv/early/2016/09/03/073189.full.pdf (Matlab code)
* Imputation of missing data
    * http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1273-5
    * TensorPute: multi-dimensional imputation https://sites.google.com/site/tensortest2/
* Other
    * Identify class-descriminative motifs enriched in subclasses of overlapping annotations: https://github.com/seqcode/sequnwinder
    * R class for integration algorithms: https://bioconductor.org/packages/release/bioc/html/MultiDataSet.html

# General Programming Resources

* Generate data type-specific compression formats: http://algorithms.cnag.cat/cargo/
* Protocol buffers
    * Protobuf: fast cross-language/platform serialization of fixed-format messages https://developers.google.com/protocol-buffers/
    * https://capnproto.org/
* IDEs
    * VisualStudio (now free): https://www.visualstudio.com/vs/visual-studio-mac/
* Parameter optimization
    * https://github.com/beamandrew/Spearmint
* Diff/patch/merge for data tables
    * js/python: http://paulfitz.github.io/daff/
    * R: https://github.com/edwindj/daff
* Pipe output of a shell command to a website (unfortunately can't be used in NIH HPC since nodes do not allow network connections): https://seashells.io/
* Debugging
    * Sandsifter: Fuzzer https://github.com/xoreaxeaxeax/sandsifter
    * https://github.com/agermanidis/livepython
* JSON Diff: http://www.jsondiff.com/
* Google C++ and python libraries: https://github.com/abseil/abseil-py
* Share a terminal as a web page: https://github.com/yudai/gotty
* Spec for CSV files with metadata: http://docs.astropy.org/en/stable/api/astropy.io.ascii.Ecsv.html

## C/C++

* kmers
    * Streaming kmer counting https://github.com/bcgsc/ntCard
    * kmer bloom filters: https://github.com/Kingsford-Group/kbf
* High-performance concurrent hash table (C++11): https://github.com/efficient/libcuckoo
* BWT that incorporates genetic variants: https://github.com/iqbal-lab/gramtools
* Fast bitwise operations on nucleotide sequences: https://github.com/kloetzl/biotwiddle
* C++ interface to htslib, BWA-MEM, and Fermi (local assembly) (would be useful to build python bindings for this): https://github.com/walaj/SeqLib
* Minimal perfect hash function: fast, large data sets https://github.com/rizkg/BBHash
* Counting quotient filter: https://github.com/splatlab/cqf
* Succinct de Bruijn Graphs: http://alexbowe.com/succinct-debruijn-graphs/
* Blazing signature filter: fast pairwise comparison of e.g. gene expression matrices https://github.com/PNNL-Comp-Mass-Spec/bsf-py
* rr debugger: record and playback executions http://rr-project.org/
* clip: command line parser https://github.com/muellan/clipp
* klib: generic library with fast implementations of striped Smith-Waterman and other NGS-related algorithms https://github.com/attractivechaos/klib

## R

* Tidy data cheatsheet: http://www.rstudio.com/wp-content/uploads/2015/02/data-wrangling-cheatsheet.pdf
* Multiple variable assignment: https://github.com/nteetor/zeallot
* Summary statistics: https://github.com/ropenscilabs/skimr

### Find packages

* https://github.com/qinwf/awesome-R#integrated-development-environment
* http://www.computerworld.com/article/2497464/business-intelligence/business-intelligence-60-r-resources-to-improve-your-data-skills.html
* http://dirk.eddelbuettel.com/cranberries/
* METACRAN: identify R packages http://www.r-pkg.org/

### Database

* MonetDB - embeddable column-store DB with R integration (MonetDB.R)

### Data Cleaning

* daff: diff/merge for data frames - https://github.com/edwindj/daff
* dplyr
    * chunked processing of large files: https://github.com/edwindj/chunked
* magrittr/pipes
    * debugging: https://github.com/gaborcsardi/tamper
* Join tables on inexact matching: https://github.com/dgrtwo/fuzzyjoin
* Tidy text: http://juliasilge.com/blog/Life-Changing-Magic/

### Reporting

* MarkDeep for R documentation: http://casual-effects.com/markdeep/
* https://confluence.broadinstitute.org/display/GDAC/Nozzle

### Misc

* Access to Google spreadsheets from R: https://github.com/jennybc/googlesheets
* Advanced table formatting in knitr: https://github.com/renkun-ken/formattable
* Access data frames using SQL: sqldf package
* Developing R packages: https://github.com/jtleek/rpackages
* Work with PDF files: https://cran.r-project.org/web/packages/pdftools/index.html
* Language-agnostic data frame format: https://github.com/wesm/feather
* Make for R: https://github.com/richfitz/maker
* Find root of current package: https://krlmlr.github.io/here/
* Work with genomic ranges using dplyr-like syntax. 

## Python

* Data structures/formats
    * Chunked, compressed, disk-based arrays: https://github.com/alimanfoo/zarr
    * Tabular data
        * Working with tabular data: http://docs.python-tablib.org/en/latest/
    	* Watch for Apache Arrow
    	* Pandas
    * GFA: https://github.com/ggonnella/gfapy/tree/master/gfapy
* A regular expression scanner: https://github.com/mitsuhiko/python-regex-scanner
* API for interacting with databases: https://github.com/kennethreitz/records
* RStudio for python: https://www.yhat.com/products/rodeo
* Debugging
    * boltons.debugutils: The entire boltons package has lots of useful stuff, but debugutils is particularly cool - you can add one line of code to enable you to drop into a debugger on signal (e.g. Ctrl-C): https://boltons.readthedocs.io/en/latest/debugutils.html
    * Web-based debugger: https://github.com/alexmojaki/birdseye
* Stats
    * Non-negative matrix factorization: https://github.com/ccshao/nimfa
    * http://www.statsmodels.org/stable/index.html
    * R formulas in python: https://github.com/pydata/patsy
* pyrasite: code injection into running applications
* Dexy: documentation
* Event loops for asynchronous programming
    * curio
    * gevent
* Fast microservices: https://github.com/squeaky-pl/japronto
* dill: alternative serialization
* arrow: alternative to datetime
* Template for scientific projects: https://github.com/uwescience/shablona
* FFI
    * GO transplier: https://github.com/google/grumpy
    * Calling Rust libraries from python: https://medium.com/@caulagi/complementing-python-with-rust-657a8cb3d066#.6in8v0bte
    * pyjamas: javascript bridge
* Disabling python garbage collection speeds up programs: https://engineering.instagram.com/dismissing-python-garbage-collection-at-instagram-4dca40b29172#.ri55nyjdu (only safe when the lifecycle is straight-forward for all objects, an thus reference counting is sufficient for memory management)
* Easily implementing function proxies/wrappers: http://wrapt.readthedocs.io/en/latest/
* Cache system: https://bitbucket.org/zzzeek/dogpile.cache
* Parse TOML (an enhanced config-file spec): https://github.com/uiri/toml
* Web scraping
    * https://github.com/grangier/python-goose
    * https://scrapy.org/
    * https://github.com/tryolabs/requestium
* Visualize python code execution time as a heatmap in a Jupyter notebook: https://github.com/csurfer/pyheatmagic
* High-level API for working with GFA (i.e. graph genomes): https://github.com/ggonnella/gfapy
* Fast and memory-efficient object counting: https://github.com/RaRe-Technologies/bounter
* Fast string matching: https://bergvca.github.io/2017/10/14/super-fast-string-matching.html
* Trace system calls in multiprocessing context: https://github.com/pinterest/ptracer
* skbio: http://scikit-bio.org/docs/0.4.1/index.html
* Extract keywords from text: https://github.com/vi3k6i5/flashtext
* Working with time-series data: https://github.com/RJT1990/pyflux
* Probability distributions and other related functions: https://pomegranate.readthedocs.io/en/latest/
* Global optimization: https://github.com/chrisstroemel/Simple
* HPAT: framework for automatically parallizing numpy and pandas code using MPI https://github.com/IntelLabs/hpat
* CuPy: Numpy-like implementation of GPU-accelerated array operations https://github.com/cupy/cupy
* Command-line interface to download and manage reference genomes: https://github.com/simonvh/genomepy
* API for working with HGVS-formatted variants: https://github.com/biocommons/hgvs
* B+ Tree: Hash table backed by on-disk storage https://github.com/NicolasLM/bplustree
* Chunked, compressed ndarrays http://zarr.readthedocs.io/en/latest/index.html
* Easy colored/styled printing in command line apps https://github.com/UltimateHackers/hue
* http://alimanfoo.github.io/2017/06/14/read-vcf.html
* Launching a subprocess in a pseudo-terminal (e.g. for accepting passwords) https://github.com/pexpect/ptyprocess
* Launch an editor from python: https://github.com/fmoo/python-editor
* Faster alternative to pyvcf: https://github.com/brentp/cyvcf2

## HPC

* Spark:
    * https://spark.apache.org/downloads.html
    * https://amplab-extras.github.io/SparkR-pkg/
* Ibis: http://blog.cloudera.com/blog/2015/07/ibis-on-impala-python-at-scale-for-data-science/
* Petuum: http://petuum.github.io/
* Flink: https://flink.apache.org/
* Dask: http://dask.pydata.org/en/latest/
* Efficient tabular storage: http://matthewrocklin.com/blog/work/2015/08/28/Storage/
* Common runtime for various libraries (e.g. Pandas, TensorFlow) that speeds up execution when interfacing the libraries with each other: https://weld-project.github.io/

## Command Line (OSX/Linux)

* Diff tables: https://github.com/paulfitz/daff
* Miller - work with tables http://johnkerl.org/miller/doc/build.html

## Databases

* CockroachDB: based on Google's distributed database https://github.com/cockroachdb/cockroach
* In-memory key-value db in python: https://github.com/paxos-bankchain/subconscious

## Reproducibility/Containerization

* http://rstudio.github.io/packrat/walkthrough.html
* Docker:
	* http://arxiv.org/pdf/1410.0846v1.pdf
	* http://bioboxes.org/available-bioboxes/
	* http://ivory.idyll.org/blog//2015-docker-and-replicating-papers.html
	* GUI for running Docker images locally: https://kitematic.com/
* Jupyter notebooks:
    * http://jupyter.org/
    * http://mybinder.org/
    * RISE: presentations from Jupyter notebooks https://github.com/damianavila/RISE
    * Convert to markdown: https://github.com/aaren/notedown
    * Jupyter notebooks on Azure: https://notebooks.azure.com/
    * Diff/merge: https://github.com/jupyter/nbdime
    * Web-based IDE built on Jupyter: https://github.com/jupyterlab/jupyterlab
    * Convert Jupyter notebook to Sphinx documentation: https://matthew-brett.github.io/nb2plots
    * Diagram editor: https://blog.jupyter.org/a-diagram-editor-for-jupyterlab-a254121ff919
* Jupyter alternatives:
    * https://nteract.io/
    * http://www.unnotebook.com/
* GUI frontend to IPython: http://nwhitehead.github.io/pineapple/
* Stencila: interesting alternative to Jupyter notebooks and R markdown https://stenci.la/
* Bioinformatics software containers
    * https://github.com/BioContainers
* Continuous analysis: https://github.com/greenelab/continuous_analysis
* Creating reproducible workflows with R Markdown documents: https://jdblischak.github.io/workflowr/

### Building Pipelines

* Luigi https://github.com/spotify/luigi
* Flo http://flo.readthedocs.org/en/latest/index.html
* Qsubsec: template language for defining SGE workflows https://github.com/alastair-droop/qsubsec
* Nextflow and Nextflow Workbench: https://www.nextflow.io/ and http://campagnelab.org/software/nextflow-workbench/
    * https://github.com/assemblerflow/assemblerflow
* SUSHI: https://github.com/uzh/sushi
* https://github.com/jdblischak/workflowr
* Data manager for R: https://cran.r-project.org/web/packages/repo/index.html
* SnakeMake
    * SnakeChunks: components for SnakeMake https://github.com/SnakeChunks/SnakeChunks
    * Create command line programs from SnakeMake workflows: https://github.com/nh13/snakeparse
    * GUI: http://sequana.readthedocs.io/en/master/sequanix.html
* dgsh: bash variant that has primitives for parallelizing tasks https://www2.dmst.aueb.gr/dds/sw/dgsh
* Rabix: visual editor for CWL pipelines http://rabix.io/
* List of workflow systems: https://github.com/common-workflow-language/common-workflow-language/wiki/Existing-Workflow-systems
* Experimental directory structure: alternative to HDF5 https://github.com/CINPLA/exdir/
* Structured projects that can be run on multiple workflow engines: http://reana.readthedocs.io/
* Invoke: http://docs.pyinvoke.org/en/latest/
* Toil: http://toil.readthedocs.io/en/latest/installation.html
* Snakemake
    * Pipelines for several genomic/epigenomic analyses https://github.com/maxplanck-ie/snakepipes
* Reflow: workflow management system developed by Grail https://github.com/grailbio/reflow
* drake: https://github.com/ropensci/drake

### Software Distribution

* Executable archives (xar): https://github.com/facebookincubator/xar


## Other

* Automated system for sequencing core facilities using microservices architecture: https://www.biorxiv.org/content/early/2017/11/06/214858
* Aether manages bidding for AWS and Azure credits: http://aether.kosticlab.org/

# Statistics/Machine Learning

* Search for papers: https://www.semanticscholar.org/
* Common probability distributions http://blog.cloudera.com/blog/2015/12/common-probability-distributions-the-data-scientists-crib-sheet/
* How to share data with a statistician: https://github.com/jtleek/datasharing
* Precision-recall curves:
    * https://cran.r-project.org/web/packages/precrec/index.html
    * Alternative: http://link.springer.com/article/10.1007/s10994-009-5119-5
* R package for subsetting data into training/testing/validation sets: https://cran.r-project.org/web/packages/STPGA/index.html
* Second-generation p-values: https://github.com/lucymcgowan/sgpvalue
* PMML: interchange format for trained models https://en.wikipedia.org/wiki/Predictive_Model_Markup_Language

## Methods/algorithms

* Lists
	* https://github.com/rushter/MLAlgorithms
* Decision tree methods
    * R randomForest package
    * FuzzyForests are an extension of random forests for classification in which subsets of variables/features are highly correlated https://github.com/OHDSI/FuzzyForest
    * https://github.com/catboost
    * https://github.com/kundajelab/boosting2D/
    * New performance metric: https://cran.r-project.org/web/packages/IPMRF/IPMRF.pdf
* Modrian forests https://scikit-garden.github.io/examples/MondrianTreeRegressor/
* Clustering
    * Help selecting biclustering algorithm: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1487-1
    * DBScan https://cran.r-project.org/web/packages/dbscan/dbscan.pdf
    * https://cran.r-project.org/web/packages/KODAMA/index.html
    * t-SNE: alternative to PCA and MDS: http://lvdmaaten.github.io/tsne/
	* http://distill.pub/2016/misread-tsne/
    * UMAP: https://github.com/lmcinnes/umap
    * Compressive k-means: https://arxiv.org/pdf/1610.08738.pdf
    * Sparse convex: https://arxiv.org/pdf/1601.04586.pdf
* Multivariate
    * http://cran.r-project.org/web/views/Multivariate.html
    * Multivariate analysis of covariance (MANCOVA): http://en.wikipedia.org/wiki/MANCOVA
    * Correlation https://www.biorxiv.org/content/biorxiv/early/2018/05/25/330498.full.pdf
* Nonnegative Matrix Factorization: https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
* Tensor factorization: https://cran.r-project.org/package=tensorBF
* Identification of correlated features within or between datasets: https://github.com/siskac/discordant
* Bayesian alternatives to standard R functions: https://github.com/rasmusab/bayesian_first_aid
* Bayesian regression modeling:
    * brms and rstanarm are R packages based on stan
    * JAGS http://jeromyanglim.blogspot.com/2012/04/getting-started-with-jags-rjags-and.html
    * MCMC http://www.stat.umn.edu/geyer/mcmc/library/mcmc/doc/demo.pdf
* Multiple test correction
    * FDR for multi-dimensional pairwise comparisons (e.g. RNA-seq): http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0937-5
    * Local false signal rate: an alternative to FDR that operates on standard error estimates rather than p-values: https://github.com/stephens999/ashr
    * MTC weighted by variant: effect http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3507.html
    * Parallelizable FDR correction: http://bioinformatics.oxfordjournals.org/content/early/2016/02/25/bioinformatics.btw029.short
    * IHW: http://bioconductor.org/packages/release/bioc/html/IHW.html
    * BonEV: Bonferonni FDR correction https://cran.r-project.org/web/packages/BonEV/index.html
    * Myriads: http://myriads.webs.uvigo.es/
* Analysis of mutual information: https://github.com/jkleinj/arMI
* Fast Bayesian alternative to lasso and ElasticNet for feature selection and effect estimation: https://cran.r-project.org/web/packages/EBglmnet/index.html
* Genome-wide generalized addative models: https://master.bioconductor.org/packages/3.3/bioc/html/GenoGAM.html
* Causal inference test: https://cran.r-project.org/web/packages/cit/index.html
* Iterative denoising tree: https://github.com/youngser/behaviotypes/blob/master/doidt.r
* Reed-Sololmon error correction
    * https://github.com/tomerfiliba/reedsolomon
    * https://github.com/brownan/Reed-Solomon
    * https://github.com/catid/wirehair
* Fast exact calculation of p-values in Friedman rank sum test: http://www.ru.nl/publish/pages/726696/friedmanrsd.zip
    * The Friedman test is for testing whether any columns are consistently different from other columns in a matrix
* Automated generation of ML pipelines: https://github.com/rhiever/tpot/tree/tpot-mdr
* Forecasting from time-series data: https://facebookincubator.github.io/prophet/ (this is a retail-centric model from Facebook, but could be adapted to biological data)
* Visualizing ML features: https://github.com/pair-code/facets
* Topological analysis:
    * https://github.com/MLWave/kepler-mapper
    * https://github.com/RabadanLab/sakmapper
    * http://danifold.net/mapper/
    * https://github.com/paultpearson/TDAmapper
* Feature selection workflows: https://github.com/enriquea/feseR
* Iterative random forest: https://github.com/sumbose/iRF
* Compares method for outlier detection: https://cran.r-project.org/web/packages/OutliersO3/index.html
* Comparison of inference methods: http://jenniferhill7.wixsite.com/acic-2016/competition
* Polynomial regression: https://github.com/matloff/polyreg

## Python

* http://www.clips.ua.ac.be/pages/pattern
* Linear mixed-model solver https://github.com/nickFurlotte/pylmm
* Exact inference in HMMs with large number of hidden states: https://github.com/regevs/factorial_hmm

## Deep Learning

* Reading
	* Papers
		* https://github.com/songrotek/Deep-Learning-Papers-Reading-Roadmap
	* Books
		* https://hackerlists.com/free-machine-learning-books/
		* http://www.deeplearningbook.org/contents/intro.html
* Platforms
	* http://aetros.com/
* Libraries
	* http://www.teglor.com/b/deep-learning-libraries-language-cm569/
	* https://github.com/fchollet/keras
	* Chainer: https://www.oreilly.com/learning/complex-neural-networks-made-easy-by-chainer
	* biologicaly-focused neural networks https://github.com/kundajelab/dragonn/tree/master/dragonn
	* analysis of features in deep neural networks https://github.com/kundajelab/deeplift
	* API to add fuzzy logic: https://fuzzy.ai/docs
	* Edward: probabalistic modeling, inference, and criticism; build on TensorFlow https://github.com/blei-lab/edward
	* VectorFlow: specifically designed for sparse data https://github.com/Netflix/vectorflow
* Architectures:
	* http://www.asimovinstitute.org/neural-network-zoo/
	* Deep residual:
		* http://image-net.org/challenges/talks/ilsvrc2015_deep_residual_learning_kaiminghe.pdf
		* https://arxiv.org/abs/1512.03385
	* Wide residual: https://arxiv.org/abs/1605.07146
	* Time-delay http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=809100&tag=1
   	* Semi-supervised classification of nodes in a graph: https://github.com/tkipf/gcn
	* Adversarial: https://arxiv.org/abs/1406.2661
	* Recurrent scalable deep kernels: https://arxiv.org/abs/1610.08936
	* Multiplicative LSTM: https://arxiv.org/abs/1609.07959
	* Graph CNN: https://arxiv.org/abs/1609.02907
	* Encoding variable-length sequences: 
	    * https://arxiv.org/abs/1505.01504
	    * https://pnpnpn.github.io/dna2vec/
	* Interactive sequence generation from RNNs: https://arxiv.org/abs/1612.04687
	* DNNs with external memory: http://www.nature.com/nature/journal/v538/n7626/full/nature20101.html
	* The Predictron: https://arxiv.org/abs/1612.08810
	* Associative LSTM: https://arxiv.org/abs/1602.03032
	* Conditional variational autoencoders: http://ijdykeman.github.io/ml/2016/12/21/cvae.html
	* Group equivariant CNN: http://jmlr.org/proceedings/papers/v48/cohenc16.pdf
	* QuickNet: https://arxiv.org/abs/1701.02291
	* Collection of pre-trained tensorflow models: https://github.com/tensorflow/models/tree/master/research
	* Generalized framework for genomic signal recognition: https://zenodo.org/record/1117159
	* Using autoencoders to impute missing values: https://github.com/gevaertlab/DAPL
* Tools
	* Tensorflow playground: http://playground.tensorflow.org/
	* http://yosinski.com/deepvis
	* https://medium.com/@shivon/the-current-state-of-machine-intelligence-3-0-e4d305da032e#.gw0ywcpkv
	* Visualize hyperparameters: https://github.com/ContextLab/hypertools
	* Convolutional layers for DNA sequence: https://github.com/koonimaru/DeepGMAP
	* Interpret neuron importance: https://arxiv.org/abs/1807.09946
* Non-bio networks that might be applied
	* https://github.com/david-gpu/srez
	* Deep learning with text: https://explosion.ai/blog/deep-learning-formula-nlp

## Web APIs

* https://cloud.google.com/prediction
* https://aws.amazon.com/machine-learning

## Data Sets

* https://medium.com/@olivercameron/20-weird-wonderful-datasets-for-machine-learning-c70fc89b73d5#.9e5byk1mo

## Text classification

* Automatic text summarization: https://pypi.python.org/pypi/sumy
* https://github.com/facebookresearch/fastText
* Word2vec: models for predicting missing words https://en.m.wikipedia.org/wiki/Word2vec

# Visualization

* Breve is a mac application that displays large tables in a way that makes it easy to identify patterns and missing data http://breve.designhumanities.org/
* Types of plots: 
    * http://www.informationisbeautifulawards.com/showcase/611-the-graphic-continuum
    * Line graph + heat map: http://www.fastcodesign.com/3052450/an-easy-intuitive-tool-for-making-sense-of-your-data
    * http://datavizproject.com/
* Making colorblind-friendly figures: http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/
* http://www.informationisbeautifulawards.com/showcase?acategory=free-tool&action=index&award=2015&controller=showcase&page=1&pcategory=long-list&type=awards
* Examples: http://www.visualcomplexity.com/vc/
* Feedback: http://helpmeviz.com/
* https://bitbucket.org/vda-lab/
* Visualization of GO results: http://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html
* Fluff: publication-quality genomics plots http://fluff.readthedocs.org
* Visualization of feature density along the genome: https://github.com/sguizard/DensityMap
* Circos-like visualization of chromosome structure with support for multiple data types https://rondo.ws
* Comparison of different types of heatmaps: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1442-6
* Interactively build simulations and animated interaction diagrams: http://ncase.me/loopy/
* Plotly: Charting library with multiple language bindings https://github.com/plotly/plotly.py
* Open-access visualization research: http://oavis.steveharoz.com/
* https://carto.com/carto-colors/
* libvips: fast imgae processing library (C++ and python bindings) https://github.com/jcupitt/libvips
* Generate flow charts, etc using markdown-like syntax: https://github.com/knsv/mermaid
* CircosVCF
    * Visualization of variant calls from a VCF on a CIRCOS plot
    * http://legolas.ariel.ac.il/~tools/CircosVCF
* Interactive reports: https://github.com/Genomon-Project/paplot
    * Make plots from regions of BAM files
    * Useful for very deep sequencing
* https://github.com/ryanlayer/samplot
* Visual inspection:
    * https://github.com/MarWoes/viper
    * https://www.biorxiv.org/content/early/2018/02/21/266262
* Venn diagrams: https://github.com/vqf/nVenn

## Networks

* Human-like Orthogonal Layout: https://github.com/skieffer/hola
* http://gephi.org/

## Phylogenetic trees

* scripts from Holt lab: https://github.com/katholt/plotTree
* web-based http://microreact.org/showcase/
* https://github.com/allendecid/TreeLink
* JavaScript libraries: https://github.com/tntvis

## R

* Comparison of libraries:
    * http://lisacharlotterost.github.io/2016/05/17/one-chart-code/
    * http://ouzor.github.io/blog/2014/11/21/interactive-visualizations.html
* D3 from R: http://christophergandrud.github.io/networkD3/
* SVG device: https://github.com/hadley/svglite/blob/master/README.md
* Multilayer data plotted on a Hilbert curve: http://www.bioconductor.org/packages/devel/bioc/html/HilbertCurve.html
* Visualize local epigenetic neighborhood of a SNP: http://bioconductor.org/packages/release/bioc/html/SNPhood.html
* CIRCOS plots: https://cggl.horticulture.wisc.edu/software/
* Nice looking boxplots: https://github.com/mw55309/perceptions
* Make multiplanel figures from a combination of plot types: https://cran.r-project.org/web/packages/multipanelfigure/
* http://bioconductor.org/packages/release/bioc/html/karyoploteR.html
* https://github.com/jkruppa/kmerPyramid
* https://github.com/venyao/shinyCircos
* PCA+clustering: http://r-forge.r-project.org/projects/thresher/
* Create pallete from colors in an image: https://github.com/AndreaCirilloAC/paletter

### ggplot2

* Gallery of extensions: http://www.ggplot2-exts.org/gallery/
* Themes
    * Cowplot - improve default ggplot: http://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
    * ggplot2 theme for publication-quality figures: https://github.com/robertwilson190/ggplot2-theme
    * https://github.com/hrbrmstr/hrbrthemes
    * xkcd-style plots: http://xkcd.r-forge.r-project.org/
* Color palattes
    * https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html
    * Color scales with clustering (would want to adapt this to ggplot): https://github.com/schne
    * https://ggsci.net/
* ggtree: phylogenetic trees https://bioconductor.org/packages/release/bioc/html/ggtree.html
* geomnet: network visualization
* ggrepel: displaying text labels with minimal overlapping https://github.com/slowkow/ggrepelrd/d3-scale-cluster
* ggforce: many extensions to ggplot
* ggalt: many extensions to ggplot
* ggraph: plotting graphs/networks
* ggedit: interactive plot editor (Shiny gadget)
* ggRandomForests https://cran.r-project.org/web/packages/ggRandomForests/index.html
* superheat: pretty heatmaps https://github.com/rlbarter/superheat
* biplots https://github.com/vqv/ggbiplot
* logos: https://github.com/omarwagih/ggseqlogo
* Make any ggplot interactive: 
    * http://rpubs.com/turnersd/ggplotly-demo
    * https://github.com/davidgohel/ggiraph
* https://github.com/eclarke/ggbeeswarm: plot overlapping points without jitter
* ggtern: Ternary diagrams https://bitbucket.org/nicholasehamilton/ggtern
* ggpubr: Publication-ready plots, including
    * barplot alternatives: http://www.sthda.com/english/rpkgs/ggpubr/
    * ggarrange, for flexible multi-panel figures
* joyplots: https://github.com/clauswilke/ggjoy
* colorbindr: test effect of colorblindness on readability of plots https://github.com/clauswilke/colorblindr
* Marginal plots: https://twitter.com/ClausWilke/status/900776341494276096
* gggenes: https://github.com/wilkox/gggenes
* egg: combine multiple plots https://cran.rstudio.com/web/packages/egg/
* https://omarwagih.github.io/ggseqlogo/
* Plotting background data: https://drsimonj.svbtle.com/plotting-background-data-for-groups-with-ggplot2
* Theme using IBM Plex fonts: https://github.com/juliasilge/silgelib/blob/master/R/graphing.R
* https://github.com/thomasp85/patchwork

### Plot Types

* Correlograms: corrgram package
* DiagrammR http://rich-iannone.github.io/DiagrammeR/
* Gene word clouds: http://genomespot.blogspot.co.uk/2014/10/geneclouds-unconventional-genetics-data.html?m=1
* Upset plots: https://cran.r-project.org/web/packages/UpSetR/index.html
* Genomic data: https://bioconductor.org/packages/GenVisR
* hextri: multiclass hexagonal bins https://cran.r-project.org/web/packages/hextri/vignettes/hexbin-classes.html
* trellis: https://www.bioconductor.org/packages/release/bioc/html/gtrellis.html
* Complex heat maps: http://www.bioconductor.org/packages/devel/bioc/html/ComplexHeatmap.html
* Scatterplot Matrix: http://bl.ocks.org/mbostock/4063663
* Beeswarm: https://flowingdata.com/2016/09/08/beeswarm-plot-in-r-to-show-distributions/
* http://flowingdata.com/2016/10/25/r-graph-gallery/
* Tilegrams: http://flowingdata.com/2016/10/13/tilegrams-in-r/
* Karyotype plots: http://bioconductor.org/packages/devel/bioc/html/karyoploteR.html
* Joyplots (for multi-sample/category time-series data): 
    * https://github.com/halhen/viz-pub/blob/master/sports-time-of-day/2_gen_chart.R (TODO: make a dedicated geom for this)
    * In python: https://github.com/mwaskom/seaborn/pull/1238
* Half-violin plot: https://gist.github.com/jbburant/b3bd4961f3f5b03aeb542ed33a8fe062

### Data Types

* Shushi: publication-quality figures from multiple data types https://github.com/dphansti/Sushi
* EpiViz: visualization of epigenomic data sets in R, http://epiviz.github.io/
* Interaction data: https://github.com/kcakdemir/HiCPlotter
* Network visualization from R using vis.js: http://dataknowledge.github.io/visNetwork/
* Differential expression from RNA-seq: http://bioconductor.org/packages/devel/bioc/html/Glimma.html

### Interactive

* https://gist.github.com/jcheng5/cbcc3b439a949deb544b
* Interactive charts: https://benjaminlmoore.wordpress.com/2015/05/19/interactive-charts-in-r/
* http://www.htmlwidgets.org/showcase_leaflet.html
* Interactive ROC plots: https://github.com/sachsmc/plotROC
* Dull: create interactive web applications - https://github.com/nteetor/dull

## Python

* Graphic of python libraries for visualization: https://pbs.twimg.com/media/C9QBxNsU0AEokR_.jpg
* Seaborn: http://web.stanford.edu/~mwaskom/software/seaborn/index.html
    * Tutorial: http://www.jesshamrick.com/2016/04/13/reproducible-plots/
* Vincent: https://github.com/wrobstory/vincent
* Pyxley (Shiny for python): http://multithreaded.stitchfix.com/blog/2015/07/16/pyxley/
* https://github.com/svaksha/pythonidae/blob/master/Computer-Graphics.md
* XKCD-style plots: http://jakevdp.github.io/blog/2012/10/07/xkcd-style-plots-in-matplotlib/?utm_content=buffera9a76&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer
* Phylogenetic trees: http://etetoolkit.org/
* Altair: https://altair-viz.github.io/
* Visualization dashboard: https://github.com/facebookresearch/visdom
* http://holoviews.org/
* GGplot clone: https://github.com/has2k1/plotnine
* https://python-graph-gallery.com/
* Plot genome tracks: https://github.com/deeptools/pyGenomeTracks
* http://pygal.org/en/stable/
* Web based/interactive
    * https://github.com/dgrtwo/gleam
    * https://github.com/adamhajari/spyre
    * Bokeh: http://bokeh.pydata.org/docs/user_guide/charts.html
    * Dash: for building data-centric web apps https://github.com/plotly/dash

## Javascript

* vis.js
* [Javascript libraries](http://blog.nextgenetics.net/?e=62)
* http://blog.webkid.io/javascript-chart-libraries/
* D3.js: https://pub.beakernotebook.com/#/publications/560c9f9b-14e6-4d95-8e78-cc0a60bf4e5a?fullscreen=false
* Scraping/JS rendering: https://splash.readthedocs.org/en/latest/
* Circos for Javascript: http://bioinfo.ibp.ac.cn/biocircos/
* VegaLite: ggplot-like framework built on top of Vega, which is built on D3.js: https://vega.github.io/vega-lite/
* GPU rendering: https://stardustjs.github.io
* Dataviz components built on top of D3: http://nivo.rocks/?ref=producthunt#/components
* https://emeeks.github.io/semiotic/#/
* Genome browser components: https://github.com/Zhong-Lab-UCSD/Genomic-Interactive-Visualization-Engine

## Examples

* Cool interactive visualization of differential data: http://graphics.wsj.com/gender-pay-gap/

# Publication/Archiving

* Licenses: http://choosealicense.com/licenses/
* APIs for literature search: http://libguides.mit.edu/apis
* Assessing credit for bioinformatics software authorship: http://depsy.org/
* Icons for presentations: http://cameronneylon.net/blog/some-slides-for-granting-permissions-or-not-in-presentations/
* Continuous analysis: https://github.com/greenelab/continuous_analysis
* A nice template for bootstrapping your own academic website using GitHub: https://academicpages.github.io/
* Desktop app for searching/managing papers: https://github.com/codeforscience/sciencefair
* Recommend papers to cite based on your bibliography: http://labs.semanticscholar.org/citeomatic/
* Word choices: [](https://pbs.twimg.com/media/DB4-B__XkAA3Nrw.jpg)
* Prepare papers for any journal format: https://typeset.io/
* Slideboards: Mashup of slides and FAQ to explain a publication http://slideboard.herokuapp.com/
* Generate manuscripts on GitHub: https://github.com/greenelab/manubot-rootstock
* General-purpose archival format: https://github.com/fair-research/bdbag/

## Code/Data sharing

* DOI for code
    * http://zenodo.org/
    * https://guides.github.com/activities/citable-code/
    * https://mozillascience.github.io/code-research-object/
    * Ruby library for fetching metadata for DOI: https://rubygems.org/gems/terrier
* CodeOcean: web platform to run algorithms https://codeocean.com
* GitHub: https://github.com/blog/1986-announcing-git-large-file-storage-lfs
* Amazon CodeCommit: http://aws.amazon.com/codecommit/
* Data sharing
    * Dat: https://datproject.org/
    * http://academictorrents.com/
    * https://figshare.com/
    * Patterns for data sharing: http://project-if.github.io/data-permissions-catalogue/
    * Globus is an open source toolkit for transferring large data files; it implements the GridFTP protocol https://www.globus.org/
    * bbcp is multi-stream scp (for point-to-point large file transfer) https://www.olcf.ornl.gov/kb_articles/transferring-data-with-bbcp/
    * Quilt: package manager for data https://quiltdata.com/
* Git plugin for version-control of data files: https://github.com/ctjacobs/git-rdm
* OSF API: https://test-api.osf.io/v2/docs/
* Data project management: https://www.datazar.com
* Execute code from GitHub project with Jupyter notebooks: http://mybinder.org/
* http://www.data-retriever.org/
* https://thinklab.com/
* http://www.researchobject.org/
* http://datalad.org/

## Journals

* Distill: Online, ML-focused journal http://distill.pub/journal/
* JORS: https://openresearchsoftware.metajnl.com/
* JOSS: http://joss.theoj.org/

## Writing

* Scripts to identify "bad smells" in science writing (would want to convert this to python): http://matt.might.net/articles/shell-scripts-for-passive-voice-weasel-words-duplicates/
* Collaborative writing
	* https://draftin.com/documents
	* http://quip.com
	* https://usecanvas.com/
* Templates:
    * InDesign template for preprint: https://github.com/cleterrier/ManuscriptTools/blob/master/biorxiv_template_CC2015.indd
    * Rmarkdown templates for journal articles https://github.com/rstudio/rticles
    * GitHub template for authoring papers: https://github.com/peerj/paper-now
    * Two-column rmarkdown template: http://dirk.eddelbuettel.com/code/pinp.html
* Convert between (R)markdown and iPython notebooks: https://github.com/aaren/notedown
* Publication-quality tables for markdown: https://cran.r-project.org/web/packages/flextable/index.html
* Pandoc scholar: https://github.com/pandoc-scholar/pandoc-scholar
    * Nice paper showing example of generating manuscripts for two different journals: https://peerj.com/preprints/2648.pdf
* Convert DOI to Bibtex entry
* Online equation editor: https://www.mathcha.io/
* Simplified markup language for equations: http://asciimath.org/
* https://editoria.pub/
* https://github.com/substance/texture

# Promising methods without software implementation

* SDR gene set analysis http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-0928-6
* ABBA http://abba.systems-genetics.net
* EpiTensor (3D genomes from 1D data): http://www.nature.com.ezproxy.nihlibrary.nih.gov/ncomms/2016/160310/ncomms10812/full/ncomms10812.html
* MOCHA: identifying modulators of transcriptional regulation from gene expression http://www.nature.com.ezproxy.nihlibrary.nih.gov/articles/srep22656
* HiC deconvolution: http://www.pnas.org.ezproxy.nihlibrary.nih.gov/content/early/2016/03/04/1512577113.full
* MR_eQTL: https://github.com/PrincetonUniversity/MR_eQTL
* Method for multi-omics integration: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1122-6
* Dynamic Bayesian network for predicting TF binding from DNase-seq: http://bioinformatics.oxfordjournals.org/content/26/12/i334.long
* De-noising ChIP-seq using DNNs: http://biorxiv.org/content/biorxiv/early/2016/05/07/052118.full.pdf
* DNA clustering: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-15-321
* Bayesian identification of bifurcations in single-cell data: http://biorxiv.org/content/biorxiv/early/2016/09/21/076547.full.pdf
* Prediction of promoter-enhancer interactions: http://biorxiv.org/content/biorxiv/early/2016/11/02/085241.full.pdf
* Mocap: TFBS prediction from ATAC-seq http://biorxiv.org/content/biorxiv/early/2016/10/27/083998.full.pdf
* Statistical method to determine deviation from time linearity in epigenetic aging: http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005183
* Multi-tissue eQTL meta-analysis: http://biorxiv.org/content/early/2017/01/16/100701
* Missing value imputation evaluator: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1429-3
* Simulated annealing for biological network alignment: https://academic.oup.com/bioinformatics/article/doi/10.1093/bioinformatics/btx090/2996219/SANA-Simulated-Annealing-far-outperforms-many
* TF binding site prediction using memory-matching networks: https://arxiv.org/pdf/1702.06760.pdf
* V-ALIGN: alignmnet on genome graphs http://www.biorxiv.org/content/biorxiv/early/2017/04/06/124941.full.pdf
* Self-organizing maps for single-cell http://www.biorxiv.org/content/biorxiv/early/2017/04/05/124693.full.pdf
* New locality sensitive hashing method: http://www.biorxiv.org/content/biorxiv/early/2017/08/25/180471.full.pdf
* Extracting knowledge from RFs: https://www.biorxiv.org/content/early/2017/11/10/217695
* Bloom filters without false-positives across a constrained set of values: http://lendulet.tmit.bme.hu/lendulet_website/wp-content/papercite-data/pdf/kiss2018bloom.pdf
* Adaptive pruning of boosted trees: https://arxiv.org/pdf/1805.07592.pdf
