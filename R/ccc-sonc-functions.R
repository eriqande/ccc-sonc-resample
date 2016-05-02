

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