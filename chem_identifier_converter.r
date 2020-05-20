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

# Extract common/preferred name from given table (result of cts_compinfo())
getName <- function(table, my_inchikey){
  my_name <- table[[my_inchikey]][["synonyms"]][table[[my_inchikey]][["synonyms"]]$type=="Synonym",]$name
  return(my_name)
}

# Extract common/preferred name from NIH
getName2 <- function(inchikey){
  qurl <- paste0("https://chem.sis.nlm.nih.gov/chemidplus/inchikey/",inchikey,"?DT_START_ROW=0&DT_ROWS_PER_PAGE=50",sep="")
  ttt <- read_html(qurl)
  txt <- xml_text(ttt)
  name <- as.vector(str_match(txt, "Substance Name:(.*?)RN:"))[2]
  name <- sub("\\[.*", "", name)
  return(substring(name,2))
}

# Extract common/preferred name from NIST
getName3 <- function(table, inchikey){
  cas <- table[[inchikey]]$externalIds[table[[inchikey]]$externalIds$name=="CAS",]$value
  qurl <- paste0("https://webbook.nist.gov/cgi/cbook.cgi?ID=",cas[1],sep="")
  ttt <- read_html(qurl)
  txt <- xml_text(ttt)
  name <- sub("\\\n.*", "", txt)
  return(name)
}

# Extract IUPAC name from given table (result of cts_compinfo())
getIUPACName <- function(table, my_inchikey){
  my_iupac <- table[[my_inchikey]][["synonyms"]][table[[my_inchikey]][["synonyms"]]$type=="IUPAC Name (Preferred)",]$name
  return(my_iupac)
}

# Obtain canonical SMILES from open babel
getCanonicalSmiles <- function(inchi){
  return(str_sub(convertFormat("inchi","can",source = inchi), end = -3))
}

# Obtain SMILES from open babel
getSmiles <- function(inchi){
  return(str_sub(convertFormat("inchi","smiles",source = inchi), end = -3))
}

# Translate chemical name to inchikey
nameToInchikey <- function(name){
  qurl <- paste0("https://chem.nlm.nih.gov/chemidplus/name/",name,sep="")
  qurl <- URLencode(qurl)
  ttt <- try(read_html(qurl), silent=T)
  if(class(ttt) == "try-error"){
    cnv <- cts_convert(name,'Chemical Name','inchikey', verbose = F)
    cnv <- as.data.frame(unlist(cnv))
    inchikey <- as.character(cnv[1,1])
    return(inchikey)
  } else {
  txt <- xml_text(ttt)
  inchikey <- as.vector(str_match(txt, "inchikey='(.*?)',"))[2]
  return(inchikey)
  }
}

#intab <- readLines("/home/kaja/metacentrum/chem_identification_converter/input_mix.csv")
intab <- readLines(opt$inFile)
#intab

smilesCanFull <- c()
smilesFull <- c()
inchikeyFull <- c()
namesFull <- c()
iupacFull <- c()

for(item in intab){
  print(item)
  if(suppressWarnings(rcdk::parse.smiles(item))[1] != "NULL"){
    print(paste0("Input ",item," is SMILES.", sep=""))
    
    inchi <- c()
    inchi <- cs_smiles_inchi(item, verbose=F)
    
    smilesCanFull <- c(smilesCanFull,getCanonicalSmiles(inchi))
    smilesFull <- c(smilesFull,getSmiles(inchi))
    
    inchikey <- c()
    inchikey <- cs_inchi_inchikey(inchi, verbose = F)
    inchikeyFull <- c(inchikeyFull,inchikey)
    
    tmp <- cts_compinfo(inchikey, verbose = F)
    name <- getName2(inchikey)
    if(is.na(name)){
      name <- getName(tmp,inchikey)[1]
    }
    if(is.na(name)){
      name <- getName3(tmp,inchikey)
    }
    
    iupac <- getIUPACName(tmp,inchikey)
    namesFull <- c(namesFull, name[1])
    iupacFull <- c(iupacFull, iupac[1])
  }
  else if(is.inchikey(item, verbose = F)){
    print(paste0("Input ",item," is InChIKey.", sep=""))
    
    inchikeyFull <- c(inchikeyFull, item)

    inchi <- c()
    inchi <- cs_inchikey_inchi(item, verbose = F)

    smilesCanFull <- c(smilesCanFull,getCanonicalSmiles(inchi))
    smilesFull <- c(smilesFull,getSmiles(inchi))
    
    tmp <- cts_compinfo(item, verbose = F)

    name <- getName2(item)
    if(is.na(name)){
      name <- getName(tmp,item)[1]
    }
    if(is.na(name)){
      name <- getName3(tmp,item)
    }

    iupac <- getIUPACName(tmp,item)

    namesFull <- c(namesFull, name[1])
    iupacFull <- c(iupacFull, iupac[1])

  }
  # if user puts in incorrect chemical name (does not matter if common or IUPAC), this will fail
  else if(!is.na(cir_query(URLencode(item), representation = 'smiles', verbose = F))){
    print(paste0("Input ",item," is chemical name.",sep=""))
    
    inchikey <- c()
    inchikey <- nameToInchikey(item)
    if(is.na(inchikey)){
      print(paste0("Chemical name ",item," is not correct OR was not found in any database.",sep=""))
      inchikeyFull <- c(inchikeyFull, "NA")
      namesFull <- c(namesFull, item)
      iupacFull <- c(iupacFull, "NA")
      smilesCanFull <- c(smilesCanFull, "NA")
      smilesFull <- c(smilesFull,"NA")
    } else {
    inchikeyFull <- c(inchikeyFull, inchikey)
    
    tmp <- cts_compinfo(inchikey = inchikey, verbose = F)
    name <- getName2(inchikey)
    if(is.na(name)){
      name <- getName(tmp,inchikey)[1]
    }
    # if(is.na(name)){
    #   name3 <- getName3(tmp,inchikey)
    # }
    iupac <- getIUPACName(tmp,inchikey)
    namesFull <- c(namesFull, name)
    iupacFull <- c(iupacFull, iupac[1])
    smilesCanFull <- c(smilesCanFull, getCanonicalSmiles(cs_inchikey_inchi(inchikey = inchikey, verbose=F)))
    smilesFull <- c(smilesFull,getSmiles(cs_inchikey_inchi(inchikey = inchikey, verbose=F)))
    }
  } else {
    stop("Supplied identifier is neither of SMILES, InChIKey or common/IUPAC name.")
  }
}

out <- data.frame(smiles_canonical=smilesCanFull,smiles=smilesFull,inchikey=inchikeyFull,common_name=namesFull,iupac_name=iupacFull)

write.csv(out, opt$outFile, row.names = F)

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
