
######################################################################
##########                 Helper Functions                 ##########
######################################################################

google.red <- '#EA4335'
google.yellow <- '#FBBC05'
google.green <- '#34A853'
google.blue <- '#4285F4'


# ============== Download and Preprocess TCGA-PRAD Data ============ #

downloadMethy <- function(manifest=manifest,directory) {
  
  ### download gdc-client
  if (! file.exists('gdc-client') & !file.exists('gdc-client.exe')) {
    downloadClientFun(Sys.info()[1])
  }
  
  Sys.chmod('gdc-client')
  
  manifestDa <- read.table(manifest, sep='\t', header=TRUE, 
                           stringsAsFactors = FALSE)
  ex <- manifestDa$filename %in% dir(paste(directory, 
                                           dir(directory), sep='/'))
  nonex <- ! ex
  numFiles <- sum(ex)
  
  if(numFiles > 0) {
    message (paste('Already exists', numFiles, 'files !', sep=' '))
    
    if (sum(nonex) > 0 ) {
      message (paste('Download the other', 
                     sum(nonex), 'files !', sep=' '))
      
      manifestDa <- manifestDa[nonex,]
      manifest <- paste(manifestDa$id, collapse =' ')
      system(paste('./gdc-client download ', manifest, sep=''))
    } else {
      return(invisible())
    }
    
    
  } else {
    system(paste('./gdc-client download -m ', manifest, sep=''))
  }
  
  
  #### move to the directory
  files <- manifestDa$id
  if (directory == 'Data') {
    if (! dir.exists('Data')) {
      dir.create('Data')
    }
  } else {
    if (! dir.exists(directory)) {
      dir.create(directory, recursive = TRUE)
    }
  }
  
  file.move(files, directory)
}


#################################################

file.move <- function(files, directory) {
  file.copy(from=files, to=directory, recursive = TRUE)
  unlink(files, recursive=TRUE)
}

#################################################

capitalize <- function(x, tow.lower=TRUE) {
  if (to.lower==TRUE) {
    x <- tolower(x)
  }
  substr(x,1,1) <- toupper(substr(x,1,1))
  return(x)
}


#################################################

generateCV <- function(n, nfold=10, repeats=10) {
  foldidID<-lapply(1:repeats,function(i){
    sample(rep(1:nfold,ceiling(n/nfold))[1:n])
  })
  return (foldidID)
}

### 10-fold CV
# foldidList <- generateCV(n=285, nfold=10, repeats=10)

### Traning-Test
# foldidList <- generateCV(n=285, nfold=2, repeats=10)

### LOOCV
# foldidList <- generateCV(n=285, nfold=285, repeats=10)

