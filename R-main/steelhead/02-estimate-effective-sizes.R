library(tidyverse)
library(afblue)

source("R/ccc-sonc-functions.R")

# In this script, we read in the estimated pedigree structure of each sample
# and drop genes down it to see what the allele frequency variance is in the 
# sample, and use that to compute an effective sample size.  




#### DO some stuff to assess the BLUE ####
# compute effective sample sizes for Colony-Run-1
pops <- dir("slg_pipe/arena/STEELHEAD_FIRST_RUN/ColonyArea/Collections")
paths <- dir("slg_pipe/arena/STEELHEAD_FIRST_RUN/ColonyArea/Collections", full.names = TRUE)
names(paths) <- pops

#### First do the simulations as if the colony-inferred pedigree is the truth ####
tmp <- lapply(paths, function(x) {
  # do the straightforward and the BLUE estimator
  ped <- read_best_config(file.path(x, "Colony-Run-1", "output.BestConfig"))
  
  bcped <- read_colony_best_config(path = file.path(x, "Colony-Run-1", "output.BestConfig"))
  samples <- bcped$id[!is.na(bcped$dad)]
  if (length(samples) > 1) {
    L <- matrix_L_from_pedigree(bcped, samples)
    weights <- weights_from_matrix_L(L)
  } else {
    weights <- NULL
  }
  
  sibgroup_eff_sample_size(ped, reps = 10000, wts = weights)
})
eff_sizes <- bind_rows(tmp, .id = "Collection")

#### Then do the simulations assuming that every individual is totally unrelated, and colony inferred the pedigree from permuted data ####
tmp <- lapply(paths, function(x) {
  # do the straightforward and the BLUE estimator
  ped <- read_best_config(file.path(x, "Permed-Run-1", "output.BestConfig"))
  
  bcped <- read_colony_best_config(path = file.path(x, "Permed-Run-1", "output.BestConfig"))
  samples <- bcped$id[!is.na(bcped$dad)]
  if (length(samples) > 1) {
    L <- matrix_L_from_pedigree(bcped, samples)
    weights <- weights_from_matrix_L(L)
  } else {
    weights <- NULL
  }
  
  sibgroup_eff_sample_size(ped, reps = 10000, wts = weights, force_unrelated = TRUE)
})
eff_sizes_unrel <- bind_rows(tmp, .id = "Collection")

names(eff_sizes_unrel)[-1] <- paste("perm-unrel-", names(eff_sizes_unrel)[-1], sep = "")

full_blue_results <- inner_join(eff_sizes, eff_sizes_unrel)

### we add an eff_num_gc (effective number of gene copies) column in eff_sizes. We use the Naive estimator for now...
eff_sizes <- eff_sizes %>%
  mutate(eff_num_gc = 2 * EffNumKidsNaive)



# now we are going to get the allele freqs, and add to them the effective sample sizes
# and also the effective counts.  Note that for the first sample these can be real numbers
# but the for the second sample they have to be rounded and they must sum to an even number
# (that is just how CoNe works...)
AFreqs <- get_steelhead_alle_freqs_from_slg_pipe(path = "slg_pipe/arena/STEELHEAD_FIRST_RUN/alle_freqs.txt") %>%
  left_join(., eff_sizes) %>%
  mutate(eff_counts = ifelse(Year == "A", freq * eff_num_gc, round(freq * eff_num_gc))) %>%
  group_by(Pop, Year, Locus) %>%
  mutate(eff_sum = sum(eff_counts),
         eff_counts_even = ifelse(Year == "B" & eff_sum %% 2 == 1, 
                                  eff_counts + rmultinom(1, 1, prob = eff_counts)[,1],  # we always add one to make it even.  Otherwise we get loci in which the sample size goes to zero.
                                  eff_counts)) %>%
  mutate(eff_counts_str = ifelse(Year == "a", sprintf("%02f", eff_counts_even), sprintf("%d", as.integer(eff_counts_even)))) %>%
  ungroup()


# we need to toss those populations that don't have two samples
toss_pops <- AFreqs %>%
  group_by(Pop, Year) %>%
  tally() %>%
  group_by(Pop) %>%
  filter(n() == 1) %>%
  select(Pop)

AFreqs2 <- anti_join(AFreqs, toss_pops)

# now, we need to drop loci that have 0 (less than or equal to 1) gene copies  sampled in any time
# period in a population
toss_pop_loci <- AFreqs2 %>%
  filter(eff_sum < 1) %>%
  group_by(Pop, Locus) %>% 
  tally() %>%
  select(Pop, Locus)

AFreqs3 <- anti_join(AFreqs2, toss_pop_loci)



# at the end of that, eff_counts_str is what we are going to want to print out
# into the CoNe file.

# now remove any existing CoNe input files and write new ones
# then for reproducibility, set a seed (not really necessary cuz they are biallelic markers....)
system("cd CoNe_area/steelhead_arena; rm *.txt *.out;  echo 12345 678910 > cone_seeds")

dump <- write_cone_eff_count_files(AFreqs3, pathprefix = "CoNe_area/steelhead_arena")


# now, we run all those pops through CoNe...
pops <- sort(unique(AFreqs2$Pop))
names(pops) <- pops
CoNe_results_list <- lapply(pops, function(x) runCoNe(x, arena_path = "CoNe_area/steelhead_arena/"))

# now, find those that we drop because their sample sizes were too small:
drop_pops <- names(CoNe_results_list)[sapply(CoNe_results_list, is.null)]

results <- CoNe_results_list[!(names(CoNe_results_list) %in% drop_pops)]

# now, bung them into a few data frames
Ne_Logls <- lapply(results, function(x) x$logl) %>%
  bind_rows(.id = "pop") %>%
  select(pop, everything())

Ne_mles <- lapply(results, function(x) x$max_etc) %>%
  bind_rows(.id = "pop") %>%
  mutate(LowerSuppLim = ifelse(LowerSuppLim == -999.9990, Inf, LowerSuppLim),
         MLE = ifelse(MLE == -999.9990, Inf, MLE),
         UpperSuppLim = ifelse(UpperSuppLim == -999.9990, Inf, UpperSuppLim)
         ) %>%
  select(pop, LowerSuppLim, MLE, UpperSuppLim)




#### Now, we also want to record the sample size and the inferred number of parents contributing to it ####
inferred_num_parents <- lapply(paths, function(x) {
  # do the straightforward and the BLUE estimator
  ped <- read_best_config(file.path(x, "Colony-Run-1", "output.BestConfig"))
  numpar <- length(unique(c(ped$pa, ped$ma)))
  numoffs <- nrow(ped)
  data_frame(inferred_num_parents = numpar, sample_size_n = numoffs)
}) %>%
 bind_rows(.id = "Collection")



#### Write out the results files to send to Libby and Carlos ####
dir.create("outputs")
write_csv(eff_sizes, path = "outputs/steelhead-effective_sample_sizes.csv")
write_csv(Ne_mles, path = "outputs/steelhead-temporal_method_ne_estimates.csv")
write_csv(inferred_num_parents, path = "outputs/steelhead-inferred_num_parents.csv")


# now, while we are at it, let's make the inferred_num_parents easier
# to interpret
inf_num_wide <- inferred_num_parents %>%
  filter(sample_size_n >= 20) %>%   # toss any collection with fewer than 20 indivs
  mutate(year_letter = str_sub(Collection, nchar(Collection), nchar(Collection))) %>%
  mutate(pop = str_sub(Collection, 1, nchar(Collection) - 1)) %>%
  mutate(parents_per_sample = inferred_num_parents / sample_size_n) %>%
  select(pop, year_letter, parents_per_sample) %>%
  tidyr::spread(data = ., key = year_letter, value = parents_per_sample)

g <- ggplot(inf_num_wide, aes(x = A)) +
  geom_point(aes(y = B), colour = "red") +
  geom_point(aes(y = C), colour = "blue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  xlab("First Episode Inferred Num Parents / Sample Size") + 
  ylab("Subsequent Episode Inferred Num Parents / Sample Size")

ggsave(g, filename = "outputs/steelhead-num-parents-scatter.pdf")

tmpA <- tmp %>% 
  filter(year_letter == "A")
tmpO <- tmp %>%
  filter(year_letter != "A")