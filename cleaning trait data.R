### Read trait data 
###
### Authors: kevin wilcox (kevin.wilcox@uwyo.edu)
### created: April 22, 2019;  last updated: April 22, 2019; 

file.choose()
### Set up workspace
rm(list=ls())
setwd("C:\\Users\\wilco\\Desktop\\Working groups\\Grazing\\trait data\\")
library(tidyverse)
library(data.table)

### REad in data
try_data <- fread("try data_raw.txt", header = T, sep = "\t", dec = ".", quote = "", data.table = T)












