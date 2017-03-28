
# This script executes the initial pipeline analyses on the data, and then it 
# sets up the ColonyRuns so the user can launch them.

# This script should be run at the top level of the ccc-sonc-resample respository (the 
# directory in which ccc-sonc-resample.Rproj lives.)

# NOTE! YOU CAN'T REALLY RUN THIS ALL STRAIGHT THROUGH IN R BECAUSE A LOT OF THE SYSTEM
# CALLS ARE RUN AS BACKGROUND JOBS.  I JUST PUT THOSE IN AS SYSTEM CALLS BUT I MOSTLY
# JUST PASTE THEM INTO THE COMMAND LINE AND LET THEM RUN.  OTHERWISE THINGS WON'T BE DONE
# WHEN R TRIES TO MOVE ON TO THE NEXT ITEM.

library(dplyr)
library(stringr)
source("R/ccc-sonc-functions.R")


#### Setup slg_pipe  ####
# set up the slg_pipe working area (must be connected to the internet)
get_slg_pipe()

# make a directory to store the input and some log files
INP <- "slg_pipe/arena/steelhead-inputs"
dir.create(INP)



#### Get the data prepared ####
metageno <- read.csv("data/steelhead-meta-and-genos.csv.gz", stringsAsFactors = FALSE) %>%
  tbl_df() %>%
  filter(!(str_detect(Duplicated_Sample_Flag, "omit")))
  

## Retain only the YOYs from 2001 ##
# this involves filtering using stream-specific length criteria
sz_grps <- suppressWarnings(prepare_sth_sizegrp_data())

# the sz_grps file that brian gave me has some sort of redundant
# columns.  We can get everything that we need from yoy_upper_mm and
# one_plus_lower
size_cats <- metageno %>%
  select(STREAM_NUMBER_SPENCE, Pop_ID, Year, AGE, LENGTH) %>%
  left_join(sz_grps %>% select(STREAM_NUMBER_SPENCE, yoy_upper_mm, one_plus_lower)) %>%
  mutate(age_category = ifelse(Year != "A",
                               ifelse(AGE == 0, "yoy", "one_plus"),
                               ifelse(LENGTH <= yoy_upper_mm,
                                      "yoy",
                                      ifelse(LENGTH < one_plus_lower, "uncertain", "one_plus")
                               )
  )) %>%
  mutate(age_category = factor(age_category, levels = c("yoy", "uncertain", "one_plus"))) # get em in the right order

         
# now, it turns out that CarRC, ElCrB, and NaInB all have nothing but yoy's (Brian checked this)
# but they do not have an entry in the AGE column.  These are, in fact, the only fish from
# years B or C that don't have age listed.  So, let's just change that.
size_cats$age_category[size_cats$Year %in% c("B", "C") & is.na(size_cats$AGE)] <- "yoy"
   
size_cats %>% 
  group_by(Year, age_category) %>%
  tally()
  
# this lets us make a little thing that will show we have gotten the right number 
# of different age categories
compare_it <- size_cats %>%
  filter(Year == "A") %>%
  group_by(STREAM_NUMBER_SPENCE, Year, age_category) %>%
  tally() %>%
  ungroup() %>% 
  tidyr::spread(data = ., key = age_category, value = n) %>%
  left_join(sz_grps, .) %>%
  mutate(yoy_inconsistency = YOY_cnt != yoy)

write.csv(compare_it, "outputs/compare-age-category-nums.csv")

### HAVE GOTTA PICK OUT THE YOY's and PROCEED NOW!!

#####################################
# be warned that there is a column at the end of all the genos giving the amount of 
# missing data there.  That is easily removed below.

# get just the genotypes we want and the IDs
genos <- metageno %>%
  select(Pop_ID:SH112876.45.1) %>%
  select(-Sample_ID, -SexID, -SexID.1)

# turns out that we are going to want to toss out some individuals that
# have a lot of missing data.  Because we are going to be tossing out 
# some apparent siblings.  It is important to have fairly complete genotypes.

# it looks like three of the loci are untyped in everyone, too.  So lets quickly
# take care of that:
geno_mat <- as.matrix(genos[, -1])
rownames(geno_mat) <- genos$Pop_ID

# find untyped loci
loc_type_counts <- colSums(geno_mat != 0)
toss_these_locus_columns <- names(loc_type_counts)[loc_type_counts == 0]

# now, let's also toss anyone that has more than 10 loci missing 
# (in addition to the two that are missing in everyone).  So, that
# translates to more than 24 0's in a row
twice_num_missing_loci <- rowSums(geno_mat == 0)
keep_these_individuals <- names(twice_num_missing_loci)[twice_num_missing_loci < 24]

# that will toss about 500 indivs

# and now we filter the individuals
genos <- genos %>%
  filter(Pop_ID %in% keep_these_individuals)

# and then toss out the loci that we don't want
genos <- genos[, !(names(genos) %in% toss_these_locus_columns)]

names(genos)[1] <- ""  # that first column has be to empty


# now make sure that the locus name is the same in each column
slg_genos <- genos
names(slg_genos)[-1][c(F,T)] <- names(slg_genos)[-1][c(T,F)]

write.table(slg_genos, file = file.path(INP, "steelhead-first-run-genos.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
            

# grab the populations in the order that they appear in the summary file:
pops <- read.csv("data/steelhead-pops-summary-S-to-N.csv", stringsAsFactors = FALSE, skip = 1) %>%
  tbl_df() %>%
  select(Population.Code, Basin, Sampling.Location)

write.table(pops, file = file.path(INP, "steelhead-first-run-pops.txt"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# make the locus file, using all the loci
loci <- names(genos)[-1][c(T, F)]
cat(loci, file = file.path(INP, "steelhead-first-run-loci.txt"),
    sep = "\n")



#### Launch the initial slg_pipe run ####

# note that the coho-first-pipe-run.sh settings file is fine to use 
# for steelhead
system(paste("cd slg_pipe/arena;",
             "../script/Do_standard_analyses.sh steelhead-inputs/steelhead-first-run-genos.txt",
             "steelhead-inputs/steelhead-first-run-pops.txt",
             "steelhead-inputs/steelhead-first-run-loci.txt",
             "STEELHEAD_FIRST_RUN",
             "../../inputs/coho-first-pipe-run.sh | tee steelhead-inputs/STEELHEAD_FIRST_RUN-log.txt")
)



# after that, we can launch the colony runs.  Note that I typically just do this on the 
# unix command line...
message("Launching colony runs")
system("cd slg_pipe/arena/STEELHEAD_FIRST_RUN/ColonyArea; ./script/RunAllColony.sh  Colony-Run-1   0  20  &")

# I didn't do the following yet
#system("cd slg_pipe/arena/STEELHEAD_FIRST_RUN/ColonyArea; ./script/RunAllColony.sh  Permed-Run-1   1  20  &")



#### Then we can do the sibyanking procedure
source("slg_pipe/R/slg_pipe_r_funcs.R")
sibyankout <- "slg_pipe/arena/STEELHEAD_FIRST_RUN/ColonyArea/Colony-Run-1-SibYankedDataSets"
dir.create(sibyankout)
set.seed(444)  # set this here for reproducibility
yanked_sibs_list <- yank_sibs(genos = file.path(INP, "steelhead-first-run-genos.txt"),
          CollDir = "slg_pipe/arena/STEELHEAD_FIRST_RUN/ColonyArea/Collections",
          Run = "Colony-Run-1",
          the_pops = file.path(INP, "steelhead-first-run-pops.txt"),
          Cutoff = 3,
          num_kept = 2,
          Num = 4,
          OutDir = sibyankout
          )

saveRDS(yanked_sibs_list, file = "slg_pipe/arena/STEELHEAD_FIRST_RUN/ColonyArea/Colony-Run-1-yanked-sibs-list.rds", compress = "xz")

#### And here we prep the 4 sib-yanked data sets for 2 runs of structure, each, at K = 2,...,10
system("cd slg_pipe/arena/STEELHEAD_FIRST_RUN/; export SLG_PATH=../..; ../../script/Prepare_StructureArea.sh ../../../inputs/coho-structure-setup-input.sh ColonyArea/Colony-Run-1-SibYankedDataSets/*.txt ")


#### Then launch those structure runs on 20 processors
system("cd slg_pipe/arena/STEELHEAD_FIRST_RUN/StructureArea/arena; nohup ../script/ExecuteStructureRuns.sh  20  > BIG_LOG.txt  2>&1 &")


#### Once they are done, we clump and distruct them
system("cd slg_pipe/arena/STEELHEAD_FIRST_RUN/StructureArea/clump_and_distruct; ./script/ClumpAndDistructAll.sh 6")

#### Then Latex that stuff.  This creates a file: slg_pipe/arena/STEELHEAD_FIRST_RUN/StructureArea/clump_and_distruct/steelhead_struct.pdf
# Must have a TeX installation.
system("cd slg_pipe/arena/STEELHEAD_FIRST_RUN/StructureArea/clump_and_distruct; ./script/LaTeXify.sh -b ./final_pdf \"2 3 4 5 6 7 8 9 10\" > steelhead_struct.tex; pdflatex steelhead_struct.tex; pdflatex steelhead_struct.tex;")
