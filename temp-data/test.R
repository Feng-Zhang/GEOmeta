# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
project_dir <- rprojroot::find_rstudio_root_file()
temp_dir <- file.path(project_dir,"temp-data")
if(!file.exists(temp_dir)) dir.create(temp_dir)

# library and set default parameter  ----------------
library(tidyverse)
library(GEOquery)
library(SummarizedExperiment)

# test function ---------------
gse <- getGEO("GSE35640", destdir = temp_dir, getGPL = T)
eSet <- gse[[1]]
gpl <- Table(getGEO("GPL25929", destdir = temp_dir))
update_gpl_list("GPL25929",probeID="ID",symbolID="ORF",overwrite = T)
load(file.path(project_dir,"data","GPL_list.rda"))


# test package install ---------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
library(devtools)
load_all()
devtools::run_examples()

library(roxygen2)
roxygen2::roxygenize()
check()
build()


#再对函数进行单元测试
usethis::use_testthat() #创建tests文件夹
usethis::use_test("DEseqObj.R") #在testthat下创建test文件
devtools::test()


