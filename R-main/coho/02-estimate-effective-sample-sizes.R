library(readr)
library(dplyr)
library(stringr)

source("R/ccc-sonc-functions.R")

# In this script, we read in the estimated pedigree structure of each sample
# and drop genes down it to see what the allele frequency variance is in the 
# sample, and use that to compute an effective sample size.  





# compute effective sample sizes for Colony-Run-1
pops <- dir("slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections")
paths <- dir("slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections", full.names = TRUE)
names(paths) <- pops
tmp <- sapply(paths, function(x) {
  ped <- read_best_config(file.path(x, "Colony-Run-1", "output.BestConfig"))
  sibgroup_eff_sample_size(ped, reps = 1000)
})
eff_sizes <- data_frame(collection = names(tmp), eff_num_gc = tmp)



# now we are going to get the allele freqs
AFreqs <- get_alle_freqs_from_slg_pipe(path = "slg_pipe/arena/COHO_FIRST_RUN/alle_freqs.txt")



# now, we 

# from this long format it should not be too hard to attach
# the effective sample size, and then write out the CoNe files.
# note that the first sample can be reals and the second sample will
# have to be rounded.