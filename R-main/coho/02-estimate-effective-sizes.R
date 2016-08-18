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
write_cone_eff_count_files(AFreqs2, pathprefix = "CoNe_area/arena")


# now, we run all those pops through CoNe...
pops <- sort(unique(AFreqs2$Pop))
names(pops) <- pops
CoNe_results_list <- lapply(pops, function(x) runCoNe(x))

# now, find those that we drop because their sample sizes were too small:
drop_pops <- names(CoNe_results_list)[sapply(CoNe_results_list, is.null)]

results <- CoNe_results_list[!(names(CoNe_results_list) %in% drop_pops)]

# now, bung them into a few data frames
Ne_Logls <- lapply(results, function(x) x$logl) %>%
  bind_rows(.id = "parent_n_factor") %>%
  mutate(pop = pop) %>%
  select(pop, everything())

Ne_mles <- lapply(results, function(x) x$max_etc) %>%
  bind_rows(.id = "parent_n_factor") %>%
  mutate(pop = pop) %>%
  mutate(UpperSuppLim = ifelse(UpperSuppLim == -999.9990, Inf, UpperSuppLim)) %>%
  select(pop, everything())

