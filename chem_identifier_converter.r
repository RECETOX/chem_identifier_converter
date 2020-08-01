#!/usr/bin/env Rscript

# Script for translating InChIKey/SMILES/ChemicalName
# 
# INPUT: InChIKey/SMILES/ChemicalName, CSV file, mixed input is allowed
# OUTPUT: InChIKey/SMILES/ChemicalName, CSV file with all identifiers for each line in input
#
# RUN: Rscript chem_identificator_convertor.r --inFile input_mix.csv --outFile test_output.csv
options(warn = -1)

library("optparse")
library("webchem")
library("ChemmineOB")
library("stringr")
library("xml2")
library("XML")
library("RJSONIO")

option_list = list(
  make_option(c("-f", "--inFile"), action="store", type="character", default=NULL, 
              help="CSV file with input SMILES/InChIKey/Chemical Name. Each molecule identifier should be on separate line, mixed identifier formats are allowed."),
  make_option(c("-t", "--outFile"), action="store", type="character", default=NULL, 
              help="CSV file with output SMILES/InChIKey/Chemical Name.")
)

opt = parse_args(OptionParser(option_list=option_list))

##############################################################################
# HELPER FUNCTIONS -> translation between different identifiers

# INCHIKEY -> COMMON NAME
# Extract common/preferred name from given table (result of cts_compinfo())
getName <- function(table, my_inchikey){
  tryCatch({
    my_name <- table[[my_inchikey]][["synonyms"]][table[[my_inchikey]][["synonyms"]]$type=="Synonym",]$name
    return(my_name)
  }, error=function(e)
    return("NA")
  )
}

# Extract common/preferred name from NIH
getName2 <- function(inchikey){
  tryCatch({
    qurl <- paste0("https://chem.nlm.nih.gov/chemidplus/inchikey/",inchikey,"?DT_START_ROW=0&DT_ROWS_PER_PAGE=50",sep="")
    Sys.sleep(runif(1, 3.0, 6.0))
    ttt <- read_html(qurl)
    txt <- xml_text(ttt)
    name <- as.vector(str_match(txt, "Substance Name:(.*?)RN:"))[2]
    name <- sub("\\[.*", "", name)
    return(substring(name,2))
  }, error=function(e){
    print("Chemical name for this compound could not be retrieved due to connection error -> database could not be reached...")
    return("NA")
  })
}

# Extract common/preferred name from NIST
getName3 <- function(table, inchikey){
  cas <- table[[inchikey]]$externalIds[table[[inchikey]]$externalIds$name=="CAS",]$value
  qurl <- paste0("https://webbook.nist.gov/cgi/cbook.cgi?ID=",cas[1],sep="")
  Sys.sleep(runif(1, 3.0, 6.0))
  ttt <- read_html(qurl)
  txt <- xml_text(ttt)
  name <- sub("\\\n.*", "", txt)
  return(name)
}

# INCHIKEY -> IUPAC NAME
# Extract IUPAC name from given table (result of cts_compinfo())
getIUPACName <- function(table, my_inchikey){
  tryCatch({
    my_iupac <- table[[my_inchikey]][["synonyms"]][table[[my_inchikey]][["synonyms"]]$type=="IUPAC Name (Preferred)",]$name
    return(my_iupac)
  }, error=function(e)
    return("NA")
  )
}

# INCHI -> SMILES
# Obtain canonical SMILES from open babel
getCanonicalSmiles <- function(inchi){
  return(str_sub(convertFormat("inchi","can",source = inchi), end = -3))
}

# Obtain SMILES from open babel
getSmiles <- function(inchi){
  return(str_sub(convertFormat("inchi","smiles",source = inchi), end = -3))
}

# NAME -> INCHIKEY
# Translate chemical name to inchikey with NIH
nameToInchikey1 <- function(name){
  qurl <- paste0("https://chem.nlm.nih.gov/chemidplus/name/",name,sep="")
  qurl <- URLencode(qurl)
  Sys.sleep(runif(1, 3.0, 6.0))
  ttt <- try(read_html(qurl), silent=T)
  if(class(ttt) == "try-error"){
    return(nameToInchikey2(name))
  } else {
  txt <- xml_text(ttt)
  inchikey <- as.vector(str_match(txt, "inchikey='(.*?)',"))[2]
  if(is.na(inchikey) | inchikey=="null" | identical(inchikey, character(0))){
    message(paste("1"))
    return(nameToInchikey2(name))
  } else {
  message(paste("2"))
  return(inchikey)
  }
  }
}

# Translate chemical name to inchikey with CTS
nameToInchikey2 <- function(name){
  cnv <- cts_convert(name,'Chemical Name','inchikey', verbose = F)
  cnv <- as.data.frame(unlist(cnv))
  inchikey <- as.character(cnv[1,1])
  return(inchikey)
}

#########################################################################################
# MAIN FUNCTIONS -> based on user input, get all the other identifiers

# INPUT IS SMILES -> get all other identifiers
inpSmiles <- function(item){
  
  inchi <- c()
  inchi <- cs_smiles_inchi(item, verbose=F)
  inchikey <- c()
  inchikey <- cs_inchi_inchikey(inchi, verbose = F)

  tmp <- cts_compinfo(inchikey, verbose = F)
  
  if(is.na(tmp[1])){
    print(paste0("For SMILES ",item,", chemical name was not found -> possibly the entry is not SMILES but chemical name, trying running entry as chemical name..."))
    inpChemName(item) }
  else  {
    print("Extracting name from NIH...")
    name <- getName2(inchikey)
    if(is.na(name)){
      print("Extracting name from NIH failed, trying CTS...")
      name <- getName(tmp,inchikey)[1]
    }
    if(is.na(name)){
      print("Extracting name from CTS failed, trying NIST...")
      name <- getName3(tmp,inchikey)
    }
    
    iupac <- getIUPACName(tmp,inchikey)
    
    namesFull <<- c(namesFull, name[1])
    iupacFull <<- c(iupacFull, iupac[1])
    smilesCanFull <<- c(smilesCanFull,getCanonicalSmiles(inchi))
    smilesFull <<- c(smilesFull,getSmiles(inchi))
    inchikeyFull <<- c(inchikeyFull,inchikey)
  }
}

# INPUT IS INCHIKEY -> get all other identifiers
inpInchikey <- function(item){
  inchi <- c()
  inchi <- cs_inchikey_inchi(item, verbose = F)

  tmp <- cts_compinfo(item, verbose = F)
  
  if (is.na(tmp[[item]])){
    inchikeyFull <<- c(inchikeyFull, item)
    smilesCanFull <<- c(smilesCanFull,"NA")
    smilesFull <<- c(smilesFull,"NA")
    namesFull <<- c(namesFull, "NA")
    iupacFull <<- c(iupacFull, "NA")
  } else {
  print("Extracting name from NIH...")
  name <- getName2(item)

  if(is.na(name)){
    print("Extracting name from NIH failed, trying CTS...")
    name <- getName(tmp,item)[1]
  }
  if(is.na(name)){
    print("Extracting name from CTS failed, trying NIST...")
    name <- getName3(tmp,item)
  }
  
  iupac <- getIUPACName(tmp,item)
  
  inchikeyFull <<- c(inchikeyFull, item)
  smilesCanFull <<- c(smilesCanFull,getCanonicalSmiles(inchi))
  smilesFull <<- c(smilesFull,getSmiles(inchi))
  namesFull <<- c(namesFull, name[1])
  iupacFull <<- c(iupacFull, iupac[1])
  }
}

# INPUT IS CHEMICAL NAME -> get all other identifiers
inpChemName <- function(item){

  if(grepl('\"', item, fixed = TRUE)){
    print("There is banned character present in the string so I am removing it...")
    item <- gsub('\"',"",item)
  }
  inchikey <- c()
  inchikey <- nameToInchikey1(item)

  if(identical(inchikey,character(0))){
    inchikey<-NA
  }

  if(is.na(inchikey)){
    print(paste0("Chemical name ",item," is not correct OR was not found in any database.",sep=""))
    inchikeyFull <<- c(inchikeyFull, "NA")
    namesFull <<- c(namesFull, item)
    iupacFull <<- c(iupacFull, "NA")
    smilesCanFull <<- c(smilesCanFull, "NA")
    smilesFull <<- c(smilesFull,"NA")
  } else {
    tmp <- cts_compinfo(inchikey = inchikey, verbose = F)
    
    print("Extracting name from NIH...")
    name <- getName2(inchikey)
    if(is.na(name)){
      print("Extracting name from NIH failed, trying CTS...")
      name <- getName(tmp,inchikey)[1]
      if(is.na(name)){
        name <- item
      }
    }
    
    iupac <- getIUPACName(tmp,inchikey)
    
    inchikeyFull <<- c(inchikeyFull, inchikey)
    namesFull <<- c(namesFull, name[1])
    iupacFull <<- c(iupacFull, iupac[1])
    smilesCanFull <<- c(smilesCanFull, getCanonicalSmiles(cs_inchikey_inchi(inchikey = inchikey, verbose=F)))
    smilesFull <<- c(smilesFull,getSmiles(cs_inchikey_inchi(inchikey = inchikey, verbose=F)))
    
  }
}

intab <- readLines(opt$inFile)

smilesCanFull <- c()
smilesFull <- c()
inchikeyFull <- c()
namesFull <- c()
iupacFull <- c()

for(item in intab){
  if(item==""){
    "Input line is empty."
  }
  else if(suppressWarnings(rcdk::parse.smiles(item))[1] != "NULL"){
    print(paste0("Input ",item," is SMILES.", sep=""))
    inpSmiles(item)
  } else if(is.inchikey(item, verbose = F)){
    print(paste0("Input ",item," is InChIKey.", sep=""))
    inpInchikey(item)
  } else #if(!is.na(cir_query(item, representation = 'smiles', verbose = F))){
    {
    print(paste0("Input ",item," is chemical name.",sep=""))
    inpChemName(item)
  }
}

out <- data.frame(smiles_canonical=smilesCanFull,smiles=smilesFull,inchikey=inchikeyFull,common_name=namesFull,iupac_name=iupacFull)

write.csv(out, opt$outFile, row.names = F)
#########################################################################################

# #intab <- readLines("/home/kaja/metacentrum/chem_identification_converter/input_mix.csv")
# intab <- readLines(opt$inFile)
# #intab
# 
# smilesCanFull <- c()
# smilesFull <- c()
# inchikeyFull <- c()
# namesFull <- c()
# iupacFull <- c()
# item <- "NICOTINIC ACID ADENINE DINUCLEOTIDE PHOSPHATE"
# for(item in intab){
#   if(suppressWarnings(rcdk::parse.smiles(item))[1] != "NULL"){
#     print(paste0("Input ",item," is SMILES.", sep=""))
# 
#     inchi <- c()
#     inchi <- cs_smiles_inchi(item, verbose=F)
#     
#     smilesCanFull <- c(smilesCanFull,getCanonicalSmiles(inchi))
#     smilesFull <- c(smilesFull,getSmiles(inchi))
#     
#     inchikey <- c()
#     inchikey <- cs_inchi_inchikey(inchi, verbose = F)
#     inchikeyFull <- c(inchikeyFull,inchikey)
#     
#     tmp <- cts_compinfo(inchikey, verbose = F)
#     
#     print("Extracting name from NIH...")
#     name <- getName2(inchikey)
#     if(is.na(name)){
#       print("Extracting name from NIH failed, trying CTS...")
#       name <- getName(tmp,inchikey)[1]
#     }
#     if(is.na(name)){
#       print("Extracting name from CTS failed, trying NIST...")
#       name <- getName3(tmp,inchikey)
#     }
#     
#     iupac <- getIUPACName(tmp,inchikey)
#     namesFull <- c(namesFull, name[1])
#     iupacFull <- c(iupacFull, iupac[1])
#   }
#   else if(is.inchikey(item, verbose = F)){
#     print(paste0("Input ",item," is InChIKey.", sep=""))
#     #item <- "MGSRCZKZVOBKFT-UHFFFAOYSA-N"
#     inchikeyFull <- c(inchikeyFull, item)
# 
#     inchi <- c()
#     inchi <- cs_inchikey_inchi(item, verbose = F)
# 
#     
#     smilesCanFull <- c(smilesCanFull,getCanonicalSmiles(inchi))
#     smilesFull <- c(smilesFull,getSmiles(inchi))
#     
#     tmp <- cts_compinfo(item, verbose = F)
#   
#     print("Extracting name from NIH...")
#     name <- getName2(item)
#     if(is.na(name)){
#       print("Extracting name from NIH failed, trying CTS...")
#       name <- getName(tmp,item)[1]
#     }
#     if(is.na(name)){
#       print("Extracting name from CTS failed, trying NIST...")
#       name <- getName3(tmp,item)
#     }
# 
#     iupac <- getIUPACName(tmp,item)
# 
#     namesFull <- c(namesFull, name[1])
#     iupacFull <- c(iupacFull, iupac[1])
# 
#   }
#   # if user puts in incorrect chemical name (does not matter if common or IUPAC), this will fail
#   else #if(!is.na(cir_query(item, representation = 'smiles', verbose = F))){
#   {
#     print(paste0("Input ",item," is chemical name.",sep=""))
#     if(grepl('\"', item, fixed = TRUE)){
#       print("There is banned character present in the string so I am removing it...")
#       item <- gsub('\"',"",item)
#     }
# 
#     inchikey <- c()
# 
#     inchikey <- nameToInchikey1(item)
# 
#     if(is.na(inchikey)){
#       print(paste0("Chemical name ",item," is not correct OR was not found in any database.",sep=""))
#       inchikeyFull <- c(inchikeyFull, "NA")
#       namesFull <- c(namesFull, item)
#       iupacFull <- c(iupacFull, "NA")
#       smilesCanFull <- c(smilesCanFull, "NA")
#       smilesFull <- c(smilesFull,"NA")
#     } else {
#     
#     inchikeyFull <- c(inchikeyFull, inchikey)
#     
#     tmp <- cts_compinfo(inchikey = inchikey, verbose = F)
#   
#     print("Extracting name from NIH...")
#     name <- getName2(inchikey)
# 
#     if(is.na(name)){
#       print("Extracting name from NIH failed, trying CTS...")
#       name <- getName(tmp,inchikey)[1]
#     }
# 
#     # if(is.na(name)){
#     #   name3 <- getName3(tmp,inchikey)
#     # }
#     
#     iupac <- getIUPACName(tmp,inchikey)
#     namesFull <- c(namesFull, name[1])
#     iupacFull <- c(iupacFull, iupac[1])
#     
#     smilesCanFull <- c(smilesCanFull, getCanonicalSmiles(cs_inchikey_inchi(inchikey = inchikey, verbose=F)))
#     smilesFull <- c(smilesFull,getSmiles(cs_inchikey_inchi(inchikey = inchikey, verbose=F)))
#     }
#   } #else {
#     #stop("Supplied identifier is neither of SMILES, InChIKey or common/IUPAC name.")
#   #}
# }
# 
# out <- data.frame(smiles_canonical=smilesCanFull,smiles=smilesFull,inchikey=inchikeyFull,common_name=namesFull,iupac_name=iupacFull)
# 
# write.csv(out, opt$outFile, row.names = F)

# # Deprecated tests
# 
# inchikey <- c("RZVAJINKPMORJF-UHFFFAOYSA-N")
# cts_compinfo("RZVAJINKPMORJF-UHFFFAOYSA-N")
# View(d[[inchikey]][["synonyms"]])
# 
# #get Chemical Name
# d[[inchikey]][["synonyms"]][d[[inchikey]][["synonyms"]]$type=="IUPAC Name (Traditional)",]$name
# 
# #get IUPAC 
# d[[inchikey]][["synonyms"]][d[[inchikey]][["synonyms"]]$type=="IUPAC Name (Preferred)",]$name

# print(smilesCanFull)
# print("--")
# print(smilesFull)
# print("--")
# print(inchikeyFull)
# print("--")
# print(namesFull)
# print("--")
# print(iupacFull)
