
#### SLG PIPE DOWNLOADING ####
#' download the slg_pipe repo and binaries and install it
#' 
#' Just a convenience function that wraps up a few different commands.  If the directory
#' slg_pipe already exists, this function does nothing.
#' @param dir  the directory in which to stick slg_pipe.  Default is .
#' @param commit  The SHA-1 hash for the commit of slg_pipe that you want to
#' get.
#' @param binary_pack  Web address of the binary pack to download.
get_slg_pipe <- function(DIR = ".", 
             commit = "e2140a5876901db29328dab94ab9f7e0cf5cb160", 
             binary_pack = "https://dl.dropboxusercontent.com/u/19274778/slg_pipe_binaries-2016-04-21.tar.gz") {
  
  if(file.exists(file.path(DIR, "slg_pipe"))) {
    message(paste("Directory slg_pipe already exists at", DIR, "   Leaving it untouched..."))
    return(NULL)
  }
  
  curdir <- getwd()
  setwd(DIR)
  
  
 
# this stuff below didn't work because all the permissions were hosed.  None of the
# scripts were executable, for goodness sake!
#  SLG <- paste("https://github.com/eriqande/slg_pipe/archive/", commit, ".zip", sep = "")
#  # get slg_pipe
#  message("Downloading slg_pipe from GitHub")
#  download.file(SLG, destfile = "tmp.zip")
#  unzip("tmp.zip")
#  file.rename(paste("slg_pipe-", commit, sep = ""), "slg_pipe") 
  
  # so, instead of the above we are going to just clone the thing...
  system("git clone https://github.com/eriqande/slg_pipe.git")
  system(paste("cd slg_pipe; git checkout ", commit, "; git checkout -b coho-working-branch;"))
  
  # get the binary pack
  message("Downloading binaries from Dropbox and rsyncing them into place")
  download.file(binary_pack, destfile = "slg_pipe_binaries.tar.gz")
  system("gunzip slg_pipe_binaries.tar.gz;
          tar -xvf slg_pipe_binaries.tar;
          rsync -avh slg_pipe_binaries/* slg_pipe")
  
  message("Removing temporary download files")
  file.remove("slg_pipe_binaries.tar")
  unlink("slg_pipe_binaries", recursive = TRUE)

  
  # change back to original working directory
  setwd(curdir)
}




#### EFFECTIVE SAMPLE SIZES ####
#' drop genes down a pedigree into sibling groups and return an effective sample size
#' 
#' This just assumes an allele freq of p and gives the parents genotypes
#' then segregates genes to their offspring according to a pedigree and then
#' it estimates the allele freq amongst the offspring.  Doing this multiple
#' times over and computing the variance of the allele freq amongst the offspring
#' we can derive an effective number of gene copies (effective sample size) which
#' we then can use while estimating Ne
#' @param ped  a data frame that gives the pedigree
#' @param reps the number of iterations to do for the Monte Carlo estimate
#' @param p the initial freq (default = 0.5)
sibgroup_eff_sample_size <- function(ped, reps, p = 0.5) {
  
  if(nrow(ped) == 1) {
    return(NA)
  }
  # get a vector of pa's and ma's
  Dads <- unique(ped$pa)
  Moms <- unique(ped$ma)
  
  # now simulate reps iterations of genotypes for all of those
  q <- 1 - p
  gDad <- matrix(sample(x = 0:2, size = length(Dads) * reps, prob = c(q^2, 2 * p * q, p^2), replace = TRUE),
                 ncol = reps)
  gMom <- matrix(sample(x = 0:2, size = length(Moms) * reps, prob = c(q^2, 2 * p * q, p^2), replace = TRUE),
                 ncol = reps)
  rownames(gDad) <- Dads
  rownames(gMom) <- Moms
  
  # now make a matrix of those where the rows are replicated as necessary for each child
  # the ma or pa had
  gD <- gDad[ped$pa,]
  gM <- gMom[ped$ma,]
  
  # now we segregate gametes from them
  hD <- gD/2
  hD[hD > 0.25 & hD < 0.75]  <- sample(0:1, size = sum(hD > 0.25 & hD < 0.75), replace = TRUE)
  
  # now we segregate gametes from them
  hM <- gM/2
  hM[hM > 0.25 & hM < 0.75]  <- sample(0:1, size = sum(hM > 0.25 & hM < 0.75), replace = TRUE)
  
  # now we make the kids by adding the haplotypes together
  gKids <- hD + hM
  
  TheVar <- mean(((colMeans(gKids) / 2) - p)^2)
  
  # this is the effective number of gene copies
  eff_num_gc <- p * (1 - p) / TheVar 
  
  eff_num_gc
}


#' read a Colony BestConfig file into a pedigree with columns pa, ma, kid
#' @param path the path to the BestConfig file you want to read.
read_best_config <- function(path) {
  read.table(path,
             header = TRUE,
             comment = "",
             stringsAsFactors = FALSE) %>%
    tbl_df %>%
    mutate(pa = str_replace(FatherID, "\\*", "pa_"),
           ma = str_replace(MotherID, "\\#", "ma_")
    ) %>%
    rename(kid = OffspringID) %>%
    select(pa, ma, kid)
}





#### CREATING CONE INPUT FILES ####

#' Read the alle_freqs from an slg_pipe run and return them in long format
#' @param path the path that the alle_freq file is found at.  For example
#' "slg_pipe/arena/COHO_FIRST_RUN/alle_freqs.txt"
get_alle_freqs_from_slg_pipe <- function(path) {
  af <- readLines(path)
  af_long_mat <- paste(af[c(T,F)], af[c(F,T)], sep = "   ") %>%
    str_split_fixed(., pattern = "[\\t ]+", n = 6)
  
  
  af_long_mat[, c(2,3,4,5)] %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    setNames(c("Locus", "Collection", "alle1", "alle2")) %>%
    tbl_df() %>%
    mutate(alle1 = as.numeric(alle1),
           alle2 = as.numeric(alle2)) %>%
    tidyr::gather(., key = "allele", value = "freq", alle1, alle2) %>%
    mutate(Pop = str_sub(Collection, 1, 3),
           Year = str_sub(Collection, 4, 4)) %>%
    select(Pop, Year, Collection, Locus, allele, freq) %>%
    arrange(Pop, Locus, Year, Collection, allele)

}



#' given a data frame D with the eff_counts_str (and freqs), group it by pop and for each one, write
#' out a CoNe input file.
write_cone_eff_count_files <- function(D, pathprefix = "tmp") {
  # first, toss loci that are monomorphic within populations
  tossers <- D %>% 
    group_by(Pop, Locus, allele) %>%
    summarise(allesum = sum(freq)) %>%
    group_by(Pop, Locus) %>%
    summarise(tossit = any(allesum < 0.000001)) %>% # this is testing for any freqs equal to zero
    filter(tossit == TRUE) %>%
    select(Pop, Locus)
  
  keepers <- ungroup(D) %>%
    anti_join(., tossers) %>%
    arrange(Pop, Locus, Year) %>%
    select(Pop, Year, Locus, allele, eff_counts_str) %>%
    tidyr::spread(., allele, value = eff_counts_str) %>%
    arrange(Pop, Locus, Year)
  
  keep_split <- split(keepers, keepers$Pop)
  
  # then just lapply over these to write the CoNe files:
  lapply(keep_split, function(x) {
    file <- file.path(pathprefix, paste(x$Pop[1], ".txt", sep = ""))
    cat(c("0", "2", nrow(x) / 2), sep = "\n", file = file)
    cat("\n\n", file = file, append = TRUE)
    tmp <- cbind(c(2, ""), x$alle1, x$alle2)
    write.table(tmp, row.names = FALSE, col.names = FALSE, file = file, append = TRUE, quote = FALSE)
    NULL
  })
}


#### FOR RUNNING CONE AND SLURPING BACK THE OUTPUT ####


#' Do a CoNe run and store the stdout in a file, then pull out what we want and return it
#' 
#' It will do this run in the CoNe_area/arena and should
#' be called from the top directory of the repository.
#' @param pop  The name of the CoNe file (without the .txt extension).  
#' @param path The directory in which to the pop.txt file exists and in which you want to run CoNe.
runCoNe <- function(pop, path = "CoNe_area/arena") {
  
    inf <- paste(pop, ".txt", sep = "")
    outf <- paste(pop, "_cone.out", sep = "")
    
    system(paste("cd ", path, ";"   " ../bin/CoNe -f ", pop, ".txt", p, " -p ../probs/ -T 4 -m 10 -n 2 5000 1 > ", outf, sep = ""))
    
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
