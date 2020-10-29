## packages I want loaded for all pages 
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
})


## knitr options I want set as default for all ('global') code chunks
knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE)