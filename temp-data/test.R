# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
project_dir <- rprojroot::find_rstudio_root_file()
temp_dir <- file.path(project_dir,"temp-data")
if(!file.exists(temp_dir)) dir.create(temp_dir)

# library and set default parameter  ----------------
library(tidyverse)
library(GEOquery)
library(SummarizedExperiment)
library(devtools)

# test function ---------------
gse_id <- "GSE218847"
gse <- getGEO(gse_id, destdir = temp_dir, getGPL = T)
eSet <- gse[[1]]
gpl_id <- "GPL29371"
gpl <- Table(getGEO(gpl_id, destdir = temp_dir)); head(gpl)
data( "GPL_list" )
update_gpl_list(gpl_id,probeID="ID",symbolID="GeneSymbol",overwrite = T)
load(file.path(project_dir,"data","GPL_list.rda"))
GPL_list

# test package install ---------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
load_all(project_dir)
devtools::run_examples()
library(roxygen2)
roxygen2::roxygenize()
check()
build()


#再对函数进行单元测试
usethis::use_testthat() #创建tests文件夹
usethis::use_test("DEseqObj.R") #在testthat下创建test文件
devtools::test()


