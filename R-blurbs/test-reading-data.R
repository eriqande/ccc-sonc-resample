

# just a thing to make sure we can read the data libby put
# here in XL files.

library(readxl)
library(dplyr)


# note that column names (like for loci) are not unique, so be careful there.

pops <- read_excel("data/coho_coastw_pops_TABLE.xlsx", sheet = 1) %>%
  tbl_df


geno_repo <- read_excel("data/coho_coastw_2003_2015_genotypes.xlsx", sheet = 1, col_names = TRUE) %>%
  tbl_df

geno_infile <- read_excel("data/coho_coastw_2003_2015_genotypes.xlsx", sheet = 2, col_names = TRUE) %>%
  tbl_df 

geno_discards <- read_excel("data/coho_coastw_2003_2015_genotypes.xlsx", sheet = 3, col_names = TRUE) %>%
  tbl_df 