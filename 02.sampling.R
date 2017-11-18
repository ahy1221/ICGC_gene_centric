#!/usr/bin/env Rscript
rm(list = ls())
suppressPackageStartupMessages(library(optparse))
parser <- OptionParser(usage = "%prog [options] histology_name output_dir")
parser <- add_option(parser, c("-c", "--CPU"), action = "store", type = "integer", default = 4, 
                     help = "CPU cores to use for parallel [default %default]")
parser <- add_option(parser, c("-r", "--replicates"), action = "store", type = "integer", default = 10000,
                     help = "How many replicates for sampling [default %default]")

args.list <- parse_args(parser, positional_arguments = 2)


M.CORES = args.list$options$CPU
REPLICATES = args.list$options$replicates
CANCER = args.list$args[1]
OUTDIR = args.list$args[2]

cat("Running arguments: \n")
print(args.list)

#======= Library and settring
source("./R/functions.R") # Analysis functions
#source("./gene_sets.R") # Analysis used gene sets
INFO <- sprintf("[INFO] Loading packages ...\n")
cat(INFO)
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Matrix))
suppressPackageStartupMessages(library(doSNOW))


INFO <- sprintf("[INFO] Loading data ...\n")
cat(INFO)
load("./data/gene.centric.rds")


#===== Get alterations 
INFO <- sprintf("[INFO] Get alteration hits matrics...\n")
cat(INFO)
alterations <- get_all_alterations(gt, alterations = names(gt)[6:14], alteration.names = names(gt)[6:14])

#==== Summary donor hits
INFO <- sprintf("[INFO] Summarising total donor hits...\n")
cat(INFO)

count.aggregate.donors<- lapply(alterations, function(x) {setNames(s<-rowSums(x), rownames(x))}) %>% 
  do.call("cbind", .) %>% 
  as.data.frame()



#===== Prepare for sampling
ndonors <- donor2cancer.use %>% 
  filter(histology_abbreviation == CANCER) %>% 
  nrow()
  
INFO <- sprintf("[INFO] Start sampling %s with %i donors for %i replicates with %i cores ...\n", CANCER, ndonors, REPLICATES, M.CORES)
cat(INFO)
cat("[TimeStamp] ", "Start sampling at ", as.character(Sys.time()), " \n")

for(idx in seq.int(alterations)) {
  alt.name <- names(alterations)[idx]
  INFO <- sprintf("[INFO] Start sampling alterations %s ...\n", alt.name )
  cat(INFO)
  sampling.hits <- sampling_hits_matrix(alterations[[idx]], n = ndonors, n.repeat = REPLICATES, 
                        group.table = donor2cancer.use, seed.use = 319, is.balance = T, n.cores = M.CORES)
  obj.name = paste0(CANCER, ".", alt.name, ".random.", REPLICATES, ".hits.mat")
  assign(obj.name, sampling.hits)
  outfile <- paste0(OUTDIR, "/",obj.name, ".rds")
  cat("[TimeStamp] ", alt.name, "is done at ", as.character(Sys.time()), "\n")
  INFO <- sprintf("[INFO] Writing output to %s ...\n", outfile)
  cat(INFO)
  
  save(list = obj.name, file = outfile)
  rm(obj.name)
}
  
  
#=================
cat("Task done at ", as.character(Sys.time()), "\n")



q(status = 0)