# set path ------------------------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
project_dir <- rprojroot::find_rstudio_root_file()
temp_dir <- file.path(project_dir,"temp-data")
if(!file.exists(temp_dir)) dir.create(temp_dir)

# library and set default parameter  ----------------
library(tidyverse)
library(GEOquery)
library(SummarizedExperiment)
gse <- getGEO("GSE35640", destdir = temp_dir, getGPL = T)
eSet <- gse[[1]]


# test package install ---------------------------
rm(list=ls());options(stringsAsFactors=FALSE)
library(devtools)
load_all()
library(roxygen2)
roxygen2::roxygenize()
check()
devtools::run_examples()
build()


#再对函数进行单元测试
usethis::use_testthat() #创建tests文件夹
usethis::use_test("DEseqObj.R") #在testthat下创建test文件
devtools::test()


