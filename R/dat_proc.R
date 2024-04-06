rm(list=ls(all=TRUE))

library(googledrive)
library(here)

drive_deauth()

# Read EPC water quality data (xlsx) from google drive and export as RData

  # google drive file url, make sure permissions are editor
  fl1url <- 'https://docs.google.com/spreadsheets/d/1KCJKYHA24Lc3a0ajYQn1kqVV7df-vlsN/edit?usp=sharing&ouid=111030724967489561225&rtpof=true&sd=true'

  # download to temp file
  fl1 <- tempfile(fileext = '.xlsx')
  drive_download(fl1url, path = fl1)
  
  # save as RData object and remove temp file
  epcwq <- readxl::read_excel(fl1)
  save(epcwq, file = here('data/epcwq.RData'))
  file.remove(fl1)

  
# Read WIN water quality data (txt) from google drive and export as RData
  
  # google drive file url, make sure permissions are editor
  fl2url <- 'https://drive.google.com/file/d/1yQQ9Y4qhneyt4m35s6w4Qd_27EhAUQPQ/view?usp=sharing'

  # download to temp file
  fl2 <- tempfile(fileext = '.txt')
  drive_download(fl2url, path = fl2)
  
  # save as RData object and remove temp file
  winwq <- read.table( fl2, sep = "|", header = TRUE, fill = TRUE, skip = 8 )
  save(winwq, file = here('data/winwq.RData'))
  file.remove(fl2)
  