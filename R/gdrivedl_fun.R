# download and load RData files from Google Drive ---------------------------------------------

gdrivedl_fun <- function(url){
  
  # temp file
  tmpfl <- tempfile(fileext = ".RData")
  
  # get downloadable URL
  id <- sub(".*d/([^/]+).*", "\\1", url)
  direct_url <- paste0("https://drive.google.com/uc?export=download&id=", id)
  
  # download and load
  download.file(direct_url, tmpfl, mode = "wb")
  load(tmpfl)
  
  # find the object name and remove temp file
  lenv <- ls()
  obj <- lenv[!lenv %in% c('direct_url', 'id', 'tmpfl', 'url')]
  unlink(tmpfl)
  
  # out
  out <- get(obj)
  
  return(out)
  
}

# example file
url <- 'https://drive.google.com/file/d/1-hSMlq9NwhNZGuVcte22Og5qppH2ouIn/view?usp=sharing'
dat <- gdrivedl_fun(url)
