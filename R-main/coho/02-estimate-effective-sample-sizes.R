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
eff_sizes <- data_frame(Collection = names(tmp), eff_num_gc = tmp)



# now we are going to get the allele freqs, and add to them the effective sample sizes
# and also the effective counts.  Note that for the first sample these can be real numbers
# but the for the second sample they have to be rounded and they must sum to an even number
# (that is just how CoNe works...)
AFreqs <- get_alle_freqs_from_slg_pipe(path = "slg_pipe/arena/COHO_FIRST_RUN/alle_freqs.txt") %>%
  left_join(., eff_sizes) %>%
  mutate(eff_counts = ifelse(Year == "a", freq * eff_num_gc, round(freq * eff_num_gc))) %>%
  group_by(Pop, Year, Locus) %>%
  mutate(eff_sum = sum(eff_counts),
         eff_counts_even = ifelse(Year == "b" & eff_sum %% 2 == 1, 
                                  eff_counts + rmultinom(1, 1, prob = eff_counts)[,1],  # we always add one to make it even.  Otherwise we get loci in which the sample size goes to zero.
                                  eff_counts)) %>%
  mutate(eff_counts_str = ifelse(Year == "a", sprintf("%02f", eff_counts_even), sprintf("%d", as.integer(eff_counts_even))))


# we need to toss those populations that don't have two samples
toss_pops <- AFreqs %>%
  group_by(Pop, Year) %>%
  tally() %>%
  group_by(Pop) %>%
  filter(n() == 1) %>%
  select(Pop)

AFreqs2 <- anti_join(AFreqs, toss_pops)

# at the end of that, eff_counts_str is what we are going to want to print out
# into the CoNe file.
write_cone_eff_count_files(AFreqs2)


# now, we 

# from this long format it should not be too hard to attach
# the effective sample size, and then write out the CoNe files.
# note that the first sample can be reals and the second sample will
# have to be rounded.