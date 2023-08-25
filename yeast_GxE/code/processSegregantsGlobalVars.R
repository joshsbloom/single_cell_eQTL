# for parallelizing the HMM

#base code directory for additional necessary functions
code.dir='/data/single_cell_eQTL/yeast/code/'
source(paste0(code.dir, 'rqtl.R'))
source(paste0(code.dir, 'getGenoInformativeCounts.R'))
source(paste0(code.dir, 'additional_fxs.R'))
source(paste0(code.dir, 'estimateEmissionProbs.R'))
source(paste0(code.dir, 'HMM_fxs.R'))
source(paste0(code.dir, 'runHMM_custom.R'))

source(paste0(code.dir, 'ASE_fxs.R'))


source(paste0(code.dir, 'processSegregantsSetup.R'))


ncores=16
registerDoMC(cores=ncores)

# set some variables---------------------------------------------------------
chroms=paste0('chr', as.roman(1:16)) 
#c('I','II', 'III', 'IV', 'V', 'X')

#some genome annotation
sacCer3=BSgenome.Scerevisiae.UCSC.sacCer3
reference.dir='/data/single_cell_eQTL/yeast/reference/'
sgd.granges=import.gff(paste0(reference.dir, 'saccharomyces_cerevisiae.gff'))
sgd=as.data.frame(sgd.granges)
gcoord.key= build.gcoord.key(paste0(reference.dir, 'sacCer3.fasta'))
sgd.genes=getSGD_GeneIntervals(reference.dir)
sgd.genes$chromosome_name=seqnames(sgd.genes)
sgd.genes$start_position=start(sgd.genes)

# for now load all crosses--------------------------------------------------
load(paste0(reference.dir, 'cross.list.RData'))
load(paste0(reference.dir, 'parents.list.RData'))
#reorganize
crossL=list('A'=cross.list[['A']], 
            'B'=cross.list[['B']], 
            '3004'=cross.list[['3004']],
            '393'=cross.list[['393']]

)
parentsL=list('A'=parents.list[['A']], 
              'B'=parents.list[['B']], 
              '3004'=parents.list[['3004']],
              '393'=parents.list[['393']]
)
rm(cross.list)
rm(parents.list)

# need genetic map
genetic.maps=lapply(crossL, extractGeneticMapDf)
#----------------------------------------------------------------------------
# 


base.dir='/data/single_cell_eQTL/yeast_GxE/'
cranger.dir='filtered_feature_bc_matrix/'

#GT.YJM145x : num [1:69385] 0 1 0 1 1 0 0 1 1 2 ...
#  .. ..$ GT.YPS163a

#.$ GT.CBS2888a : num [1:78041] 0 1 1 2 0 0 0 0 0 0 ...
#  .. ..$ GT.YJM981x  : num [1:78041] 1 0 0 0 1 1 1 1 1 1 ...

#cList=list(
#           '01_2444_44_1-2'='B',
#           '02_2444_44-'='B',
#           '03_ByxRM_51-_1-2'='A',
#           '04_ByxRM_51-'='A',
#           '00_BYxRM_480MatA_1'='A',
#           '00_BYxRM_480MatA_2'='A',
#           '08_2444_cross_10k_Feb_21'='B',
#           '09_2444_cross_5k_Feb_21'='B',
#           '10_3004_cross_10k_Feb_21'='3004',
#           '11_3004_cross_5k_Feb_21'='3004',
#           '07_2444-cross-1'='B',
#           '07_2444-cross-2'='B',
#           '19_393_10k_May10'='393',
#           '20_393_20k_May10'='393',
#           '21_3004_10k_May10'='3004',
#           '22_3004_20k_May10'='3004'
#)
#
#hList=list(
#           '01_2444_44_1-2'=2^6.1,
#           '02_2444_44-'= 2^7.5,
#           '03_ByxRM_51-_1-2'=2^7,
#           '04_ByxRM_51-'=2^7.4,
#           '00_BYxRM_480MatA_1'=2^6.5,
#           '00_BYxRM_480MatA_2'=2^6.5,
#           '08_2444_cross_10k_Feb_21'=2^5.2,
#           '09_2444_cross_5k_Feb_21'= 2^4.9,
#           '10_3004_cross_10k_Feb_21'= 2^5.6,
#           '11_3004_cross_5k_Feb_21'=  2^4.9,
#           '07_2444-cross-1'  = 2^6,
#           '07_2444-cross-2'  = 2^5.5,
#           '19_393_10k_May10' = 2^5.2,
#           '20_393_20k_May10' = 2^5.1,
#           '21_3004_10k_May10'= 2^6,
#           '22_3004_20k_May10'= 2^6.7
#)

cList=list(
           'NaCl_0.7M_t0_3051_rep1'='A',
           'NaCl_0.7M_t30_3051_rep1'='A')
hList=list(
           'NaCl_0.7M_t0_3051_rep1'=2^5,
           'NaCl_0.7M_t30_3051_rep1'=2^5 )



experiments=names(cList)
cdesig=as.vector(sapply(cList, function(x) x))
het.thresholds=as.vector(sapply(hList, function(x) x))
data.dirs=paste0(base.dir, 'processed/', experiments, '/')
het.across.cells.threshold=.1

cl <- makeCluster(36)
clusterEvalQ(cl, {
       library(fastglm)
       library(Matrix)
       library(MASS) 
       library(nebula)
       NULL  })

# changed from 10k
nUMI_thresh=16384 #000
minInformativeCellsPerTranscript=128

nperm=5

#note sgd gene def is updated to include 3' utr               
#sgd.genes=sgd.granges[sgd.granges$type=='gene',]
#sgd.genes$gcoord=gcoord.key[as.character(seqnames(sgd.genes))]+start(sgd.genes)

# RUN HMM 

# combine experiments

#good YPSxYJM was 11 and 12
#7,8,

#each separate
sets=list('A_NaCl_p7M_t0'=c(1),
          'A_NaCl_p7M_t30'=c(2))

#together
#sets=list('A_NaCl_p7M'=c(1,2))        
#sets=list(
#    '3004'=c(9,10),
#    'A'=c(3,4),
#    'B'=c(1,2,11,12)
#    'Ap'=c(5,6)
#)
#sets=list(
#          'Ap'=c(5,6)
#          )
cycle.cats=c('G1', 'G1:S', 'S', 'G2:M', 'M:G1')

