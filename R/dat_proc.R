library(googledrive)
library(here)

drive_deauth()

# example file download from google drive -----------------------------------------------------

# google drive file url, make sure permissions are editor
flurl <- 'https://docs.google.com/spreadsheets/d/1KCJKYHA24Lc3a0ajYQn1kqVV7df-vlsN/edit?usp=sharing&ouid=111030724967489561225&rtpof=true&sd=true'

# download to temp file
fl <- tempfile(fileext = '.xlsx')
drive_download(flurl, path = fl)

# save as RData object and remove temp file
epcwq <- readxl::read_excel(fl)
save(epcwq, file = here('data/epcwq.RData'))
file.remove(fl)