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

# Discovery

When looking for a bioinformatics tool for a specific application:

* http://omictools.com/
* http://gitxiv.com/category/bioinformatics
* https://bio.tools/
* https://biosharing.org

# Data Sets

* 1000 Genomes (RNA-Seq, ChIP-Seq): http://archive.gersteinlab.org/docs/2015/06.04/1kg_fun_studies.htm
* Reference panel of ~250 Dutch families http://biorxiv.org/content/early/2016/01/18/036897
* FANTOM consortium has CAGE (5' single molecule RNA counting) data from ~1000 human cell/tissue/cell-line samples from ~300 different cell/tissue types
* Full text of all PMC papers from 2008-present: ftp://ftp.ncbi.nlm.nih.gov/pub/pmc/manuscript/
* Human tissue-specific enhancers: http://www.enhanceratlas.org/

# Genomics

## General Information

* Migrating to GRCh38: https://software.broadinstitute.org/gatk/blog?id=8180
* Mappings between contig names in different assemblies: https://github.com/dpryan79/ChromosomeMappings

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

## Databases

* NAR catalog of databases, by subject: http://www.oxfordjournals.org/our_journals/nar/database/c/
* Super-Enhancer Archive: http://www.bio-bigdata.com/SEA/
* GWAS database: http://jjwanglab.org/gwasdb
* rVarBase: regulatory features of human genetic variants http://rv.psych.ac.cn/
* TransVar http://bioinformatics.mdanderson.org/transvarweb/
* dbMAE: mono-allelic expression https://mae.hms.harvard.edu/
* Disease-gene associations http://www.disgenet.org/web/DisGeNET/menu/rdf
* BISQUE: convert between database identifiers http://bisque.yulab.org/

## Data Formats

* FASTA/FASTQ
    * Read filtering and extraction: https://github.com/ad3002/Cookiecutter
    * fqtools: for working with FASTQ files https://github.com/alastair-droop/fqtools
* SAM/BAM/CRAM
    * Read filtering and profiling: https://github.com/jwalabroad/VariantBam
* VCF
    * Annotation: https://github.com/brentp/vcfanno
    * https://github.com/lh3/bgt
    * GQT
    * BCFtools: includes new tool to identify RoH http://samtools.github.io/bcftools/
    * Work with VCF in R: https://github.com/knausb/vcfR

## Footprinting

* http://rajanil.github.io/msCentipede/
* NucID: nucleosome positioning from DNase-seq https://jianlingzhong.github.io/NucID/
* SeqGL: predict TF binding from DNase/ATAC-seq https://bitbucket.org/leslielab/seqgl/wiki/Home
* DeFCoM: https://bitbucket.org/bryancquach/defcom

## Functional Enrichment

* General:
    * R interface to DAVID: http://www.bioconductor.org/packages/release/bioc/html/RDAVIDWebService.html
    * GSEA
    * GiANT: uncertainty in GSEA https://cran.r-project.org/web/packages/GiANT/index.html
    * Ensemble gene set enrichment analysis: https://bitbucket.org/malhamdoosh/egsea
    * Fast GSA: https://github.com/billyhw/GSALightning
    * Fast GSEA: https://github.com/ctlab/fgsea
    * Multi-dimensional GSEA: http://bioconductor.org/packages/release/bioc/html/mdgsa.html
    * GSEA with external information https://cran.r-project.org/web/packages/netgsa/index.html
    * QTest: http://statgen.snu.ac.kr/software/QTest/
* GO term:
    * http://cbl-gorilla.cs.technion.ac.il/
    * http://lrpath.ncibi.org/
    * Reduce GO term lists: http://revigo.irb.hr/
    * GO Express: https://www.bioconductor.org/packages/release/bioc/html/GOexpress.html
    * GO Extender: https://www.msu.edu/~jinchen/GOExtender/
    * http://iwera.ir/~ahmad/dal/
    * Negative GO enrichment: https://sites.google.com/site/guoxian85/neggoa
* Variant Set
    * https://cran.r-project.org/web/packages/VSE/vignettes/my-vignette.html
    * Functional enrichment with LD correction
* MESH
    * MeSH over-representation: http://www.bioconductor.org/packages/release/bioc/vignettes/meshr/inst/doc/MeSH.pdf
* Regional
    * LOLA: http://lola.computational-epigenetics.org
    * AnnotatR: https://github.com/rcavalcante/annotatr/
* Trait
    * traseR: Trait-associated SNP enrichment https://www.bioconductor.org/packages/release/bioc/html/traseR.html
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
* LD score calculation and regression https://github.com/bulik/ldsc
* Imputation of missing phenotype information http://www.nature.com/ng/journal/vaop/ncurrent/full/ng.3513.html
* Epistasis
    * General linear-mixed model library; also includes mixed-RF method for detecting epistasis with population structure correction: https://github.com/PMBio/limix
    * GPU-accelerated detection of epistasis using Bayesian neural networks: https://github.com/beamandrew/BNN
    * MEPID: marginal epistasis test http://www.xzlab.org/software.html
* Other
    * Browser for geographical distribution of genetic variants: http://popgen.uchicago.edu/ggv/

## Microarray

* Correct for batch effects between training data and external datasets: https://cran.r-project.org/web/packages/bapred/index.html
* Impute from Affy expression arrays: http://simtk.org/home/affyimpute
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

## Network Analysis

* Identification of module with biclustering: http://web.ist.utl.pt/rmch/bicnet/temporary/index.jsp
* Functional analysis https://bioconductor.org/packages/release/bioc/html/EGAD.html
* GAGE: http://bioconductor.org/packages/release/bioc/html/gage.html
* PAXToolsR: http://bioconductor.org/packages/release/bioc/html/paxtoolsr.html
* ancGWAS: post GWAS association with protein-protein interaction networks http://www.cbio.uct.ac.za/~emile/software.html
* Gene network reconstruction
    * For a set of TFs: https://sourceforge.net/projects/aracne-ap/
    * From PPI or motif sharing: https://github.com/davidvi/pypanda (is also an integrative method that can incorporate multiple sources of information)
* BicMix: differential co-expression networks http://beehive.cs.princeton.edu/software/
* BANFF: https://cran.r-project.org/web/packages/BANFF/index.html
* https://bitbucket.org/abarysh/safe
* https://bitbucket.org/roygroup/merlin-p

## Prediction

* QTL
    * Imputation of summary statistics in multi-ethnic cohorts: http://dleelab.github.io/jepegmix/
    * Causal variant identification: http://genetics.cs.ucla.edu/caviar/
* eQTL
    * Imputation of gene expression from genotype data : https://github.com/hriordan/PrediXcan
* Causal variant
    * Ensemble method: https://github.com/gifford-lab/EnsembleExpr/
    * eCAVIAR: probability that a variant is causal for both QTL and eQTL http://genetics.cs.ucla.edu/caviar/
    * https://github.com/dleelab/qcat
    * Disease-associated risk variants: https://sites.google.com/site/emorydivan/
* Chromatin States
    * GenoSTAN: http://bioconductor.org/packages/release/bioc/html/STAN.html
* Enhancers
    * Prediction of enhancer strength from sequenced http://bioinformatics.hitsz.edu.cn/iEnhancer-2L
    * Prediction of core cell type-specific TFs from super enhancers https://bitbucket.org/young_computation/crcmapper
* Regulatory variants/TF binding
    * LedPred: prediction of regulatory sequences from ChIP-seq https://github.com/aitgon/LedPred
    * GERV: prediction of regulatory variants that affect TF inding http://gerv.csail.mit.edu/
    * Score variant deleteriousness: http://cadd.gs.washington.edu/
    * BASSET: prediction of sequence activity https://github.com/davek44/Basset
    * DanQ: hybrid convolutional and recurrent neural network model for predicting the function of DNA de novo from sequence http://github.com/uci-cbcl/DanQ
    * [LINSIGHT](http://genome-mirror.cshl.edu/cgi-bin/hgTables?hgsid=110380_RYHnZaJsP5sXJRVSjw7Z6nhgc0oE&clade=mammal&org=Human&db=hg19&hgta_group=allTracks&hgta_track=LINSIGHT&hgta_table=0&hgta_regionType=genome&position=chr21%3A33031597-33041570&hgta_outputType=primaryTable&hgta_outFileName=)
    * Protein binding affinity: https://bitbucket.org/wenxiu/sequence-shape.git
    * Change in local frustration index: https://github.com/gersteinlab/frustration

## Sequencing Protocols

* Single cell
    * Simultaneous RNA and methylation (and inference of CNV): http://www.nature.com/cr/journal/vaop/ncurrent/full/cr201623a.html
    * Simultaneous RNA and methylation (scM&T-seq): http://www.nature.com/nmeth/journal/v13/n3/full/nmeth.3728.html
    * Simultaneous RNA and protein measurements: http://www.sciencedirect.com/science/article/pii/S2211124715014345

## Simulation

* InterSIM: simulate correlated multi-omics data https://cran.r-project.org/web/packages/InterSIM/index.html

## SNP Array

* Call haplotypes https://cran.r-project.org/web/packages/GHap/index.html

## Variant annotation

* WGSA pipeline https://sites.google.com/site/jpopgen/wgsa/
* http://snpeff.sourceforge.net/
* Normalization of SNP ID's from literature: https://github.com/rockt/SETH
* https://hail.is/
* Prediction of functional impact
    * HaploReg: http://www.broadinstitute.org/mammals/haploreg/haploreg.php
    * Several tools/score sets: CADD, DANN, etc
    * Consensus approaches:
        * http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004962
        * http://jjwanglab.org/PRVCS/
    * Impact of coding SNPs: http://pantherdb.org/tools/csnpScoreForm.jsp
    * Predict disease risk from GWAS summary statistics: https://github.com/yiminghu/AnnoPred

## Sequence Analysis

### General-purpose

* Google Genomics R API: https://followthedata.wordpress.com/2015/02/05/notes-on-genomics-apis-2-google-genomics-api/
* k-mer counting
    * khmer:
        * https://github.com/ged-lab/khmer
        * https://docs.google.com/presentation/d/1biQmLkwPlCOA56mNZdUAiXa1OyGE0qyvI2nvU0qlOIE/mobilepresent?pli=1&slide=id.p58
    * https://github.com/abdullah009/kcmbt_mt
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
* Demultiplexing/deduping barcoded reads w/ UMIs: http://gbcs.embl.de/portal/tiki-index.php?page=Je

### QC

* Qualimap2: http://qualimap.bioinfo.cipf.es/
* Determine whether two BAM files are from the same source: https://bitbucket.org/sacgf/bam-matcher
* DOGMA: Measure completeness of a transcriptome or proteome assembly https://ebbgit.uni-muenster.de/domainWorld/DOGMA/
* Identify and remove UMI sequences from reads: https://github.com/CGATOxford/UMI-tools
* Integrated report from multiple tools: http://multiqc.info/
* Identify batch effects: https://github.com/mani2012/BatchQC
* AlmostSignificant: https://github.com/bartongroup/AlmostSignificant
* Genetic relatedness from raw reads: 
    * https://github.com/kdmurray91/kwip
    * https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/Software.cgi
* Fast coverage estimate from BAM index: https://github.com/brentp/goleft/tree/master/indexcov

### ATAC-seq/DNase-seq

* DNase footprinting: https://github.com/ajank/Romulus
* HINT: http://costalab.org/publications-2/dh-hmm/ (was best out of 10 compared tools in recent NatMeth paper)
* ALTRE: https://mathelab.github.io/ALTRE/vignette.html

### ChIP-seq

* Allocation of multi-mapping reads: https://github.com/keleslab/permseq
* R package for TFBS analysis: http://bioconductor.org/packages/release/bioc/html/TFBSTools.html
* R package for predicting chromatin states from histone marks across conditions https://github.com/ataudt/chromstaR
* hiddenDomains: https://sourceforge.net/projects/hiddendomains/
* GenoGAM peak caller: https://master.bioconductor.org/packages/3.3/bioc/html/GenoGAM.html
* Bayesian motif discovery: https://github.com/soedinglab/BaMMmotif
* Network-based identification of relationships among ChIP-seq data sets: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0925-0
* De-noising: https://github.com/kundajelab/TF_chipseq_pipeline
* Peak discretizer (merging of replicates): https://github.com/nanakiksc/zerone
* Cell type-specific TFBS analysis (focuses analysis on TFs expressed in cell type of interest): https://github.com/Danko-Lab/rtfbs_db
* FIMO and MCAST perform best of TFBS predictors: http://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-016-1298-9
* Motif assessment: http://www.bioinf.ict.ru.ac.za/
* https://github.com/cfce/chilin

### Chromatin Interactions

* Model 3D chromosome structure from Hi-C contact maps + optional FISH constraints: https://github.com/yjzhang/FISH_MDS.jl, https://github.com/yjzhang/3DC-Browser
* Predict enhancer targets: https://github.com/shwhalen/targetfinder
* Predicting TADs from histone modifications: https://cb.utdallas.edu/CITD/index.htm#ajax=home

### DNA

* Error correction
    * https://github.com/lh3/bfc
    * RECKONER http://sun.aei.polsl.pl/REFRESH/index.php?page=projects&project=reckoner&subpage=about
    * MultiRes - says it's for viral populations; not sure if applicable to humans: https://github.com/raunaq-m/MultiRes
* Duplicate removal
    * https://sourceforge.net/projects/pardre/
* Alignment
    * Align simultaneously against multiple reference genomes http://1001genomes.org/software/genomemapper.html
    * Compressed reference-based alignment: http://groups.csail.mit.edu/cb/cora/
    * Compression and querying of aligned haplotype data: https://github.com/richarddurbin/pbwt
* Assembly
    * Build de Bruijn Graph from multiple genomes: https://github.com/medvedevgroup/TwoPaCo
    * Align to a de Brujn graph: https://github.com/Malfoy/BGREAT
* Genotyping
    * https://cran.r-project.org/web/packages/ebGenotyping/ebGenotyping.pdf
    * http://bioinfo.ut.ee/FastGT/
    * STR genotyping from NGS: http://melissagymrek.com/lobstr-code/
* CNV calling
    * https://github.com/cui-lab/multigems
    * WHAM: CNV caller https://github.com/zeeev/wham
    * Identification of mosic events: https://github.com/asifrim/mrmosaic
* Repeat calling
    * REPdenovo: https://github.com/Reedwarbler/REPdenovo
* Pipelines
    * SpeedSeq: alignment/annotation pipeline - https://github.com/hall-lab/speedseq
* Ancestry and kinship analysis: http://illumina.github.io/akt/
* Phasing
    * Eagle2: https://data.broadinstitute.org/alkesgroup/Eagle/
    * Using HiC+partial haplotypes: https://github.com/YakhiniGroup/SpectraPh
* Other
    * VCF compression and data extraction: https://github.com/kedartatwawadi/GTRAC
    * Run length encoded multi-sample BWT + server: https://github.com/wtsi-svi/ReadServer

### Methylation

* BS-SNPer: fast SNP calling from bisulfite-converted sequencing reads https://github.com/hellbelly/BS-Snper
* MACAU: Mixed-model association http://www.xzlab.org/software.html
* ME-plot: Error detection and correction in bisulfite-converted sequencing reads https://github.com/joshuabhk/methylsuite
* Metilene: Calling differential methylation http://www.bioinf.uni-leipzig.de/Software/metilene/
* Reference-free bisulfite sequence comparison: https://github.com/thomasvangurp/epiGBS
* Correction for cell-type composition: http://www.cs.tau.ac.il/~heran/cozygene/software/refactor.html
* Predicting gene expression from methylation: http://arxiv.org/abs/1603.08386
* GEM: R package for meQTL and EWAS https://bioconductor.org/packages/devel/bioc/html/GEM.html

### RNA

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
    * Wasabi: prepare Salmon/Sailfish output for Kallisto https://github.com/COMBINE-lab/wasabi
    * featureCounts: read summarization http://bioinf.wehi.edu.au/featureCounts/
    * D-GEX: Quantification of whole-transcriptome gene expression from landmark genes https://github.com/uci-cbcl/D-GEX
    * Faster version of HTSeq/featureCount: https://github.com/qinzhu/VERSE
    * Quantification using both structure and abundance information: https://pypi.python.org/pypi/rsq
    * Improve transcript quantification by integrating PolII ChIP-seq data: https://github.com/pliu55/RSEM/tree/pRSEM
* Correction/Normalization
    * TDM: cross-platform normalization https://github.com/greenelab/TDMresults
    * alpine: corrects for fragment sequence bias https://github.com/mikelove/alpine/blob/master/vignettes/alpine.Rmd
    * Filter out lowly-expressed transcripts prior to quantification decreases FP rate: http://www.genomebiology.com/2016/17/1/12
    * Partition variance between biological and technical sources: https://www.bioconductor.org/packages/3.3/bioc/vignettes/variancePartition
    * Quantify and correct for uncertainty in abundance estimates: https://github.com/PSI-Lab/BENTO-Seq
    * Simultaneous isoform discovery and quantification across multiple samples: http://cbio.ensmp.fr/flipflop
    * R package to compare normalization methods: https://github.com/Edert/NVT
* Workflows:
    * Artemis (RNA-Seq workflow designed around Kallisto): https://github.com/RamsinghLab/artemis
    * https://github.com/ririzarr/rafalib
    * Isolator: https://github.com/dcjones/isolator
* eQTL
    * MT-HESS: eQTL analysis across tissues http://www.mrc-bsu.cam.ac.uk/software/
    * eQTLBMA: cross-tissue eQTL https://github.com/timflutre/eqtlbma
    * Multi-tissue eQTL: https://cran.r-project.org/web/packages/JAGUAR/index.html
    * Quantile regression approach: https://xiaoyusong.shinyapps.io/QRBT/
    * Using prior knowledge: https://github.com/redsofa/LassoMP
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
* ASE
    * GeniASE: ASE without haplotypes http://www.ncbi.nlm.nih.gov/pmc/articles/PMC4758070/pdf/srep21134.pdf
    * Without phasing: http://lifecenter.sgst.cn/cisASE/
* Splicing
    * https://github.com/hartleys/JunctionSeq
    * https://github.com/lkmklsmn/SplicER
    * Annotation of splicing types https://r-forge.r-project.org/projects/splicingtypes/
    * MAJIQ: detection of local splice variation http://majiq.biociphers.org/
    * https://github.com/davidaknowles/leafcutter
    * http://www.mhs.biol.ethz.ch/research/krek/jsplice.html
    * Splice site prediction: http://cabgrid.res.in:8080/HSplice/
    * DRIM-seq: http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html
    * Fast quantification of differential splicing: https://github.com/comprna/SUPPA
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
* Time series
    * http://diceseq.sourceforge.net/
* Deconvolution:
    * VoCAL: https://cran.r-project.org/web/packages/ComICS/index.html
    * DeconRNASeq: deconvolute expression profiles in mixed tissues http://www.bioconductor.org/packages/2.12/bioc/vignettes/DeconRNASeq/inst/doc/DeconRNASeq.pdf
* Search
    * http://www.cs.cmu.edu/∼ckingsf/software/bloomtree/
    * SBTs https://github.com/medvedevgroup/bloomtree-allsome
* Other:
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
    * Identify gene fusions: https://github.com/ndaniel/fusioncatcher
    * Prediction of intronic splice branchpoints: https://github.com/betsig/branchpointer/

### Single-cell

* Comparative analysis of methods: http://biorxiv.org/content/early/2016/01/05/035758
* Review of experimental design and analysis: http://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0927-y
* QC:
    * http://www.morgridge.net/SinQC.html
    * https://github.com/YosefLab/scone
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
* Differential Expression
    * Monocle cole-trapnell-lab.github.io/monocle-release/
    * scDD: https://github.com/kdkorthauer/scDD
    * ISOP: comparison of isoform pairs in single cells https://github.com/nghiavtr/ISOP
    * D3E: http://hemberg-lab.github.io/D3E/
    * BASiCS: https://github.com/catavallejos/BASiCS
    * Beta Poisson: https://github.com/nghiavtr/BPSC
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
* SNVs/CNVs
    * DNA SNV calling: https://bitbucket.org/hamimzafar/monovar
    * Ginko: analysis of CNVs in single-cell data: http://qb.cshl.edu/ginkgo/?q=/XWxZEerqqY477b9i4V8F
    * CNV calling: http://genome.cshlp.org/content/early/2016/01/15/gr.198937.115.full.pdf
* Regulatory networks:
    * Gene co-expression: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1004892
    * https://github.com/hmatsu1226/SCODE
* Other
    * Analysis of 3' tagging data: https://github.com/garber-lab/ESAT
    * StemID: Prediction of stem cells and lineage information https://github.com/dgrun/StemID
    * Phasing: https://github.com/edsgard/scphaser
    * Classification using sets of known cell type-specific genes: https://github.com/YosefLab/FastProject/wiki
    * UMI counting: https://github.com/vals/umis
* Methylation
    * Prediction of missing information: https://github.com/cangermueller/deepcpg

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
    * Multi-tissue: http://bioconductor.org/packages/release/bioc/html/HDTD.html
    * SDA: integrate gene expression across multiple tissues, or multi-omics in a single tissue, for identification of trans-QTL networks: https://jmarchini.org/sda/
    * HMM for binary classification based on multivariate data: https://github.com/PetarV-/muxstep
    * https://github.com/fraenkel-lab/OmicsIntegrator
    * Correlation between enriched regions from different data sets: http://malone.bioquant.uni-heidelberg.de/software/mcore
    * https://cran.r-project.org/web/packages/r.jive/
    * SIFORM: http://bioinformatics.oxfordjournals.org/content/early/2016/07/03/bioinformatics.btw295.full
    * https://sourceforge.net/projects/epimine/
    * http://bioconductor.org/packages/release/bioc/html/fCCAC.html
* Specific data types:
    * Predict gene fusions from WGS and RNA-seq: http://sourceforge.net/p/integrate-fusion/wiki/Home/
    * NuChart: layer additional omics data on Hi-C ftp://fileserver.itb.cnr.it/nuchart/
    * Predict expression from H3K27, and identify cis-regulatory elements: http://cistrome.org/MARGE/
    * GenoSkyline: predict tissue-specific functional regions from epigenomic data http://genocanyon.med.yale.edu/GenoSkyline
* Network-based:
    * Merging networks: https://github.com/maxconway/SNFtool
* Causality
    * Hybrid BN/CMI approach to constructing GRNs: http://journals.plos.org/ploscompbiol/article?id=10.1371%2Fjournal.pcbi.1005024
* Noise/bias:
    * Bias correction across different assays on the same samples: https://cran.r-project.org/web/packages/MANCIE/
    * Cell cycle heterogeneity appears to be a minor contributor to noise; instead, library size is the largerst PC by far http://www.nature.com/nbt/journal/v34/n6/full/nbt.3498.html

## Upcomming methods to watch for

* From ASHG:
    * Wang - SplineAdjust DE correction for gene length
    * FIRE prediction of cis-eqtls
    * Hiccups Hi-C processing
* phaseME http://beehive.cs.princeton.edu/wiki/phaseme/
* ccRemover: remove cell-cycle effect from single-cell data https://arxiv.org/ftp/arxiv/papers/1605/1605.04492.pdf
* NETAM: network-based GWAS http://bioinformatics.oxfordjournals.org/content/32/12/i164.short

# General Programming Resources

* Generate data type-specific compression formats: http://algorithms.cnag.cat/cargo/

## C++

* High-performance concurrent hash table (C++11): https://github.com/efficient/libcuckoo
* kmer bloom filters: https://github.com/Kingsford-Group/kbf
* BWT that incorporates genetic variants: https://github.com/iqbal-lab/gramtools
* Fast bitwise operations on nucleotide sequences: https://github.com/kloetzl/biotwiddle

## R

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

## Python

* Data structures/formats
    * Chunked, compressed, disk-based arrays: https://github.com/alimanfoo/zarr
    * Tabular data
        * Working with tabular data: http://docs.python-tablib.org/en/latest/
    	* Watch for Apache Arrow
    	* Pandas
* Pipelines
    * Invoke: http://docs.pyinvoke.org/en/latest/
    * Toil: http://toil.readthedocs.io/en/latest/installation.html
    * Snakemake
* A regular expression scanner: https://github.com/mitsuhiko/python-regex-scanner
* API for interacting with databases: https://github.com/kennethreitz/records
* R formulas in python: https://github.com/pydata/patsy
* RStudio for python: https://www.yhat.com/products/rodeo
* boltons.debugutils: The entire boltons package has lots of useful stuff, but debugutils is particularly cool - you can add one line of code to enable you to drop into a debugger on signal (e.g. Ctrl-C): https://boltons.readthedocs.io/en/latest/debugutils.html
* Non-negative matrix factorization: https://github.com/ccshao/nimfa
* pyjamas: javascript bridge
* pyrasite: code injection into running applications
* Dexy: documentation
* Event loops for asynchronous programming
    * curio
    * gevent
* dill: alternative serialization
* arrow: alternative to datetime
* Template for scientific projects: https://github.com/uwescience/shablona

## HPC

* Spark:
    * https://spark.apache.org/downloads.html
    * https://amplab-extras.github.io/SparkR-pkg/
* Ibis: http://blog.cloudera.com/blog/2015/07/ibis-on-impala-python-at-scale-for-data-science/
* Petuum: http://petuum.github.io/
* Flink: https://flink.apache.org/
* Dask: http://dask.pydata.org/en/latest/
* Efficient tabular storage: http://matthewrocklin.com/blog/work/2015/08/28/Storage/

## Command Line (OSX/Linux)

* Diff tables: https://github.com/paulfitz/daff
* Miller - work with tables http://johnkerl.org/miller/doc/build.html

## Reproducibility/Containerization

* http://rstudio.github.io/packrat/walkthrough.html
* Docker:
	* http://arxiv.org/pdf/1410.0846v1.pdf
	* http://bioboxes.org/available-bioboxes/
	* http://ivory.idyll.org/blog//2015-docker-and-replicating-papers.html
	* GUI for running Docker images locally: https://kitematic.com/
* NextFlow: http://www.nextflow.io/
* Jupyter notebooks:
    * http://jupyter.org/
    * http://mybinder.org/
    * http://nwhitehead.github.io/pineapple/
    * RISE: presentations from Jupyter notebooks https://github.com/damianavila/RISE
* Stencila: interesting alternative to Jupyter notebooks and R markdown https://stenci.la/

### Building Pipelines

* Luigi https://github.com/spotify/luigi
* Flo http://flo.readthedocs.org/en/latest/index.html
* Qsubsec: template language for defining SGE workflows https://github.com/alastair-droop/qsubsec
* Nextflow and Nextflow Workbench: http://campagnelab.org/software/nextflow-workbench/
* SUSHI: https://github.com/uzh/sushi

# Statistics/Machine Learning

* Search for papers: https://www.semanticscholar.org/
* Common probability distributions http://blog.cloudera.com/blog/2015/12/common-probability-distributions-the-data-scientists-crib-sheet/
* How to share data with a statistician: https://github.com/jtleek/datasharing
* Precision-recall curves: https://cran.r-project.org/web/packages/precrec/index.html

## Methods/algorithms

* Lists
	* https://github.com/rushter/MLAlgorithms
* Random Forests in R: randomForest package
    * FuzzyForests are an extension of random forests for classification in which subsets of variables/features are highly correlated https://github.com/OHDSI/FuzzyForest
* Clustering
    * DBScan https://cran.r-project.org/web/packages/dbscan/dbscan.pdf
    * https://cran.r-project.org/web/packages/KODAMA/index.html
    * t-SNE: alternative to PCA and MDS: http://lvdmaaten.github.io/tsne/
	* http://distill.pub/2016/misread-tsne/
* Multivariate analysis http://cran.r-project.org/web/views/Multivariate.html
* Multivariate analysis of covariance (MANCOVA): http://en.wikipedia.org/wiki/MANCOVA
* Nonnegative Matrix Factorization: https://cran.r-project.org/web/packages/NMF/vignettes/NMF-vignette.pdf
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
* Analysis of mutual information: https://github.com/jkleinj/arMI
* Fast Bayesian alternative to lasso and ElasticNet for feature selection and effect estimation: https://cran.r-project.org/web/packages/EBglmnet/index.html
* Genome-wide generalized addative models: https://master.bioconductor.org/packages/3.3/bioc/html/GenoGAM.html
* Causal inference test: https://cran.r-project.org/web/packages/cit/index.html
* Iterative denoising tree: https://github.com/youngser/behaviotypes/blob/master/doidt.r
* Reed-Sololmon error correction
    * https://github.com/tomerfiliba/reedsolomon
    * https://github.com/brownan/Reed-Solomon
    * https://github.com/catid/wirehair

## Python

* http://www.clips.ua.ac.be/pages/pattern
* Linear mixed-model solver https://github.com/nickFurlotte/pylmm

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
	* Encoding variable-length sequences: https://arxiv.org/abs/1505.01504
	* Interactive sequence generation from RNNs: https://arxiv.org/abs/1612.04687
* Tools
	* Tensorflow playground: http://playground.tensorflow.org/
	* http://yosinski.com/deepvis
	* https://medium.com/@shivon/the-current-state-of-machine-intelligence-3-0-e4d305da032e#.gw0ywcpkv
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

# Visualization

* Breve is a mac application that displays large tables in a way that makes it easy to identify patterns and missing data http://breve.designhumanities.org/
* Types of plots: http://www.informationisbeautifulawards.com/showcase/611-the-graphic-continuum
    * Line graph + heat map: http://www.fastcodesign.com/3052450/an-easy-intuitive-tool-for-making-sense-of-your-data
* Making colorblind-friendly figures: http://bconnelly.net/2013/10/creating-colorblind-friendly-figures/
* http://www.informationisbeautifulawards.com/showcase?acategory=free-tool&action=index&award=2015&controller=showcase&page=1&pcategory=long-list&type=awards
* Examples: http://www.visualcomplexity.com/vc/
* Feedback: http://helpmeviz.com/
* https://bitbucket.org/vda-lab/
* Visualization of GO results: http://cran.r-project.org/web/packages/GOplot/vignettes/GOplot_vignette.html
* Fluff: publication-quality genomics plots http://fluff.readthedocs.org
* Visualization of feature density along the genome: https://github.com/sguizard/DensityMap
* Circos-like visualization of chromosome structure with support for multiple data types https://rondo.ws

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

### ggplot2

* Cowplot - improve default ggplot: http://cran.r-project.org/web/packages/cowplot/vignettes/introduction.html
* ggplot2 theme for publication-quality figures: https://github.com/robertwilson190/ggplot2-theme
* ggtree: phylogenetic trees https://bioconductor.org/packages/release/bioc/html/ggtree.html
* geomnet: network visualization
* ggrepel: displaying text labels with minimal overlapping https://github.com/slowkow/ggrepel
* Color scales with clustering (would want to adapt this to ggplot): https://github.com/schnerd/d3-scale-cluster
* ggforce: many extensions to ggplot
* ggalt: many extensions to ggplot
* ggraph: plotting graphs/networks
* ggedit: interactive plot editor (Shiny gadget)

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

* Seaborn: http://web.stanford.edu/~mwaskom/software/seaborn/index.html
    * Tutorial: http://www.jesshamrick.com/2016/04/13/reproducible-plots/
* Bokeh: http://bokeh.pydata.org/docs/user_guide/charts.html
* Vincent: https://github.com/wrobstory/vincent
* Pyxley (Shiny for python): http://multithreaded.stitchfix.com/blog/2015/07/16/pyxley/
* https://github.com/svaksha/pythonidae/blob/master/Computer-Graphics.md
* XKCD-style plots: http://jakevdp.github.io/blog/2012/10/07/xkcd-style-plots-in-matplotlib/?utm_content=buffera9a76&utm_medium=social&utm_source=twitter.com&utm_campaign=buffer
* Phylogenetic trees: http://etetoolkit.org/

## Javascript

* vis.js
* [Javascript libraries](http://blog.nextgenetics.net/?e=62)
* http://blog.webkid.io/javascript-chart-libraries/
* D3.js: https://pub.beakernotebook.com/#/publications/560c9f9b-14e6-4d95-8e78-cc0a60bf4e5a?fullscreen=false
* Scraping/JS rendering: https://splash.readthedocs.org/en/latest/
* Circos for Javascript: http://bioinfo.ibp.ac.cn/biocircos/
* VegaLite: ggplot-like framework built on top of Vega, which is built on D3.js: https://vega.github.io/vega-lite/

## Examples

* Cool interactive visualization of differential data: http://graphics.wsj.com/gender-pay-gap/

# Publication/Archiving

* DOI for code
	* https://guides.github.com/activities/citable-code/
	* https://mozillascience.github.io/code-research-object/
* Licenses: http://choosealicense.com/licenses/
* Large file storage:
	* GitHub: https://github.com/blog/1986-announcing-git-large-file-storage-lfs
	* Amazon CodeCommit: http://aws.amazon.com/codecommit/
* Templates:
	* InDesign template for preprint: https://github.com/cleterrier/ManuscriptTools/blob/master/biorxiv_template_CC2015.indd
	* Rmarkdown templates for journal articles https://github.com/rstudio/rticles
	* GitHub template for authoring papers: https://github.com/peerj/paper-now
* OSF API: https://test-api.osf.io/v2/docs/
* APIs for literature search: http://libguides.mit.edu/apis
* Places to archive research output:
	* http://zenodo.org/
	* https://figshare.com/
	* DAT: distributed data sharing http://dat-data.com/blog/2016-02-01-dat-1
	* Git plugin for version-control of data files: https://github.com/ctjacobs/git-rdm
* Assessing credit for bioinformatics software authorship: http://depsy.org/
* Icons for presentations: http://cameronneylon.net/blog/some-slides-for-granting-permissions-or-not-in-presentations/
* Patterns for data sharing: http://project-if.github.io/data-permissions-catalogue/
* Ruby library for fetching metadata for DOI: https://rubygems.org/gems/terrier
* Data project management: https://www.datazar.com
* Convert between (R)markdown and iPython notebooks: https://github.com/aaren/notedown

## Writing

* Scripts to identify "bad smells" in science writing (would want to convert this to python): http://matt.might.net/articles/shell-scripts-for-passive-voice-weasel-words-duplicates/
* Collaborative writing
	* https://draftin.com/documents
	* http://quip.com

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
