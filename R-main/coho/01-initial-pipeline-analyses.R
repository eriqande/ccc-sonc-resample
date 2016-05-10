
# This script executes the initial pipeline analyses on the data, and then it 
# sets up the ColonyRuns so the user can launch them.

# This script should be run at the top level of the ccc-sonc-resample respository (the 
# directory in which ccc-sonc-resample.Rproj lives.)

# NOTE! YOU CAN'T REALLY RUN THIS ALL STRAIGHT THROUGH IN R BECAUSE A LOT OF THE SYSTEM
# CALLS ARE RUN AS BACKGROUND JOBS.  I JUST PUT THOSE IN AS SYSTEM CALLS BUT I MOSTLY
# JUST PASTE THEM INTO THE COMMAND LINE AND LET THEM RUN.  OTHERWISE THINGS WON'T BE DONE
# WHEN R TRIES TO MOVE ON TO THE NEXT ITEM.


library(readxl)
library(dplyr)
source("R/ccc-sonc-functions.R")


#### Setup slg_pipe  ####
# set up the slg_pipe working area (must be connected to the internet)
get_slg_pipe()

# make a directory to store the input and some log files
INP <- "slg_pipe/arena/coho-inputs"
dir.create(INP)



#### Get the data prepared ####
# slurp the data out of Libby's excel files and prepare the slg_pipe run.
# this assumes that the pipe-input file is in sheet #2.
genos <- read_excel("data/coho_coastw_2003_2015_genotypes.xlsx", sheet = 2) %>%
  tbl_df
names(genos)[1] <- ""  # that first column has be to empty
genos <- genos[!is.na(genos[[1]]), ]  # the excel file ends with three rows of NA's.  This removes them.
# write that to a text file for input
write.table(genos, file = file.path(INP, "coho-first-run-genos.txt"), 
            sep = "\t", row.names = FALSE, quote = FALSE)
            

# grab the populations in the order that libby has put them
pops <- read_excel("data/coho_coastw_pops_N_S.xlsx", sheet = 1)
pops <- pops[, c(2,4,5)]
write.table(pops, file = file.path(INP, "coho-first-run-pops.txt"), 
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# make the locus file, using all the loci
loci <- names(genos)[c(F,T)]
cat(loci, file = "slg_pipe/arena/coho-inputs/coho-first-run-loci.txt",
    sep = "\n")



#### Launch the initial slg_pipe run ####
system(paste("cd slg_pipe/arena;",
             "../script/Do_standard_analyses.sh coho-inputs/coho-first-run-genos.txt",
             "coho-inputs/coho-first-run-pops.txt",
             "coho-inputs/coho-first-run-loci.txt",
             "COHO_FIRST_RUN",
             "../../inputs/coho-first-pipe-run.sh | tee coho-inputs/COHO_FIRST_RUN-log.txt")
)



# after that, we can launch the colony runs.  Note that I typically just do this on the 
# unix command line...
message("Launching colony runs")
system("cd slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/; ./script/RunAllColony.sh  Colony-Run-1   0  20  &")
system("cd slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/; ./script/RunAllColony.sh  Permed-Run-1   1  20  &")



#### Then we can do the sibyanking procedure
source("slg_pipe/R/slg_pipe_r_funcs.R")
sibyankout <- "slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Colony-Run-1-SibYankedDataSets"
dir.create(sibyankout)
set.seed(555)  # set this here for reproducibility
yanked_sibs_list <- yank_sibs(genos = file.path(INP, "coho-first-run-genos.txt"),
          CollDir = "slg_pipe/arena/COHO_FIRST_RUN/ColonyArea/Collections",
          Run = "Colony-Run-1",
          the_pops = file.path(INP, "coho-first-run-pops.txt"),
          Cutoff = 3,
          Num = 4,
          OutDir = sibyankout
          )


#### And here we prep the 4 sib-yanked data sets for 2 runs of structure, each, at K = 2,...,10
system("cd slg_pipe/arena/COHO_FIRST_RUN/; export SLG_PATH=../..; ../../script/Prepare_StructureArea.sh ../../../inputs/coho-structure-setup-input.sh ColonyArea/Colony-Run-1-SibYankedDataSets/*.txt ")


#### Then launch those structure runs on 20 processors
system("cd slg_pipe/arena/COHO_FIRST_RUN/StructureArea/arena; nohup ../script/ExecuteStructureRuns.sh  20  > BIG_LOG.txt  2>&1 &")
