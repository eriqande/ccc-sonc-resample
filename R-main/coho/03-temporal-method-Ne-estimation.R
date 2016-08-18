library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(ggplot2)

#### some useful functions ####

# write out some CoNe_input
# note that this only works for biallelic markers.
#' @param par_n_factor how much do you want to decrease the observed number of parents?
CoNe_input <- function(afreq1 = "slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections/ALBa/Colony-Run-1/output.AlleleFreq", 
                       afreq2 = "slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections/ALBb/Colony-Run-1/output.AlleleFreq",
                       bestcluster1 = "slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections/ALBa/Colony-Run-1/output.BestConfig",
                       bestcluster2 = "slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections/ALBb/Colony-Run-1/output.BestConfig",
                       outf = "CoNe_infile.txt",
                       par_n_factor = 1) {
 
  af1 <- read_table(afreq1)
  af2 <- read_table(afreq2)
  
  bc1 <- read_table(bestcluster1)
  bc2 <- read_table(bestcluster2)
  
  c1 <- af1 %>% 
    mutate(N = length(unique(c(bc1$FatherID, bc1$MotherID))),
           pcount = round(2 * N * UpdatedFreq * par_n_factor))  %>%
    group_by(MarkerID) %>%
    mutate(pcount_even = pcount) # sum of allele counts at each locus need not be even in the first sample
  
  
  c2 <- af2 %>% 
    mutate(N = length(unique(c(bc2$FatherID, bc2$MotherID))),
           pcount = round(2 * N * UpdatedFreq * par_n_factor))  %>%
    group_by(MarkerID) %>%
    mutate(pcount_even = pcount - ((sum(pcount) %% 2) == 1) * as.vector(rmultinom(1, 1, UpdatedFreq))) # this ensures the count at each locus is even
  
  both <- bind_rows(c1, c2, .id = "sample")
  
  # now deal with cases where one allele did not show up at all in one of the years
  filled <- both %>%
    select(MarkerID, sample, AlleleID, pcount_even) %>%
    tidyr::complete(sample, nesting(AlleleID, MarkerID)) %>%
    arrange(MarkerID, sample, AlleleID)
  
  filled$pcount_even[is.na(filled$pcount_even)] <- 0
  
  # now, prepare strings to write out the CoNe file
  lines <- filled %>%
    group_by(MarkerID) %>%
    filter(n() == 4)  %>%   # toss out monomorhic loci
    arrange(sample, AlleleID) %>%
    summarise(str = paste(2, "\n", pcount_even[1], pcount_even[2], "\n", pcount_even[3], pcount_even[4]))
  
  # then write that dude out
  cat("0\n2",nrow(lines), sep = "\n", file = outf)
  cat(lines$str, sep = "\n", file = outf, append = TRUE)
}


#' Do a CoNe run and store the stdout in a file, then pull out what we want and return it
#' 
#' It will do this run in the CoNe_area/arena and should
#' be called from the top directory of the repository. It also assumes that
#' the colony results 
#' @param pop  The name of the population without the a or b on the end that
#' indicates the first and second samples
#' @run_name  the name of the directory in which the colony run results you want to get are
#' @parn vector of values of par_n_factor that you want to use
runCoNe <- function(pop, run_name = "Colony-Run-1", parn = c(0.5, 1)) {
  apath <- paste("slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections/", pop, "a/", run_name, sep = "")
  bpath <- paste("slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections/", pop, "b/", run_name, sep = "")
  
  afile <- file.path(apath, "output.AlleleFreq")
  bfile <- file.path(bpath, "output.AlleleFreq")
  bestcluster_a <- file.path(apath, "output.BestConfig")
  bestcluster_b <- file.path(bpath, "output.BestConfig")
  
  
  names(parn) <- parn
  results <- lapply(parn, function(p) {
    CoNe_input(afile, bfile, bestcluster_a, bestcluster_b, 
               outf = paste("CoNe_area/arena/", pop, ".txt", p,  sep = ""), 
               par_n_factor = p)
    outf <- paste(pop, "_cone.out", p, sep = "")
    
    system(paste("cd CoNe_area/arena; ../bin/CoNe -f ", pop, ".txt", p, " -p ../probs/ -T 4 -m 10 -n 2 5000 1 > ", outf, sep = ""))
    
    # now slurp those data in
    x <- readLines(file.path("CoNe_area/arena", outf))
    
    if(length(x) >  12) {
      tmp <- x[str_detect(x, "^NE_LOGLIKE")] %>%
        str_split_fixed(., "  *", 9) 
      tmp <- tmp[, 3:5]
      
      header <- tmp[1,]
      
      tmp <- tmp[-1,]
      mode(tmp) <- "numeric"
      
      # get the logl curve
      logl <- as.data.frame(tmp) %>%
        setNames(header) %>%
        tbl_df()
      
      # now pick out the max and the support limits
      tmp <- x[str_detect(x, "MaxByParabolicInterp|LowerSupportLimit|UpperSupportLimit")] %>%
        str_split_fixed(., "  *", 7)
      
      max_etc <- data_frame(MLE = as.numeric(tmp[1,3]),
                            LowerSuppLim = as.numeric(tmp[2,3]),
                            UpperSuppLim = as.numeric(tmp[3,3]))
      
      ret <- list(logl = logl, max_etc = max_etc)
    } else {
      ret <- NULL
    }
    ret
  })
  
  if(all(sapply(results, is.null))) {
    ret <- NULL
  } else {
    # now, bung them into a few data frames
    Logls <- lapply(results, function(x) x$logl) %>%
      bind_rows(.id = "parent_n_factor") %>%
      mutate(pop = pop) %>%
      select(pop, everything())
    
    mles <- lapply(results, function(x) x$max_etc) %>%
      bind_rows(.id = "parent_n_factor") %>%
      mutate(pop = pop) %>%
      mutate(UpperSuppLim = ifelse(UpperSuppLim == -999.9990, Inf, UpperSuppLim)) %>%
      select(pop, everything())
    
    ret <- list(logl = Logls, mle = mles)
  }
  ret
  
}



####  Now we are ready to run them all ####
# first pick out the names of the pops that have both an a and a b
pops <- dir("slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections/") %>%
  str_sub(.,start = 1, end = 3)
two_pops <- names(table(pops))[table(pops) == 2]


all_results <- lapply(two_pops, function(x) {print(x); runCoNe(x)})

# drop the ones that were full NULL
all_results <- all_results[!sapply(all_results, is.null)]

# get tidy data frames
logls <- lapply(all_results, function(x) x$logl) %>%
  bind_rows()

mles <- lapply(all_results, function(x) x$mle) %>%
  bind_rows()

mles[mles == -999.99900] <- Inf


mles %>% 
  arrange(parent_n_factor, MLE) %>% 
  write.csv(file = "ccc-ne-estimates.csv") 

