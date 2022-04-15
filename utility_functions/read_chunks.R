Read_chunks_csv <- function(directory){
  library(tidyverse)
  print(directory)
	fileNames <- Sys.glob(file.path(directory,"*.csv"))
	combined <- read.csv(fileNames[1])
	for (fileName in fileNames[2:length(fileNames)]){
		combined <- rbind(combined, read.csv(fileName))
	}
	return(combined)
}

#setwd('/gpfs/gsfs5/users/TCR/10X_Genomics/scRNAseq_P1_HUMAN_GEX_V2')
#df <- Read_chunks_csv('./data/interim/regression/5.0.6_ME_Int_OneGrp_logNorm_AllStats_scripts_saveChunks_cluster')
