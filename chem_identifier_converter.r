#!/usr/bin/env Rscript

# Script for translating InChIKey/SMILES/ChemicalName
# 
# INPUT: InChIKey/SMILES/ChemicalName, CSV file, mixed input is allowed
# OUTPUT: InChIKey/SMILES/ChemicalName, CSV file with all identifiers for each line in input
#
# RUN: Rscript chem_identificator_convertor.r --inFile input_mix.csv --outFile test_output.csv
library("optparse")
library("webchem")
library("ChemmineOB")
library("stringr")

option_list = list(
  make_option(c("-f", "--inFile"), action="store", type="character", default=NULL, 
              help="CSV file with input SMILES/InChIKey/Chemical Name. Each molecule identifier should be on separate line, mixed identifier formats are allowed."),
  make_option(c("-t", "--outFile"), action="store", type="character", default=NULL, 
              help="CSV file with output SMILES/InChIKey/Chemical Name.")
)

opt = parse_args(OptionParser(option_list=option_list))

# Extract common/preferred name from given table (result of cts_compinfo())
getName <- function(table, my_inchikey){
  my_name <- table[[my_inchikey]][["synonyms"]][table[[my_inchikey]][["synonyms"]]$type=="IUPAC Name (Traditional)",]$name
  return(my_name)
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

#intab <- readLines("/home/kaja/metacentrum/chem_identification_converter/input_mix.csv")
intab <- readLines(opt$inFile)
#intab

smilesFull <- c()
inchikeyFull <- c()
namesFull <- c()
iupacFull <- c()

for(item in intab){
  if(rcdk::parse.smiles(item)[1] != "NULL"){
    print(paste0("Input ",item," is SMILES.", sep=""))
    
    inchi <- c()
    inchi <- cs_smiles_inchi(item, verbose=F)
    
    smilesFull <- c(smilesFull,getCanonicalSmiles(inchi))
    
    inchikey <- c()
    inchikey <- cs_inchi_inchikey(inchi, verbose = F)
    inchikeyFull <- c(inchikeyFull,inchikey)
    
    tmp <- cts_compinfo(inchikey, verbose = F)
    name <- getName(tmp,inchikey)
    iupac <- getIUPACName(tmp,inchikey)
    namesFull <- c(namesFull, name)
    iupacFull <- c(iupacFull, iupac)
  }
  else if(is.inchikey(item, verbose = F)){
    print(paste0("Input ",item," is InChIKey.", sep=""))
    
    inchikeyFull <- c(inchikeyFull, item)
    
    inchi <- c()
    inchi <- cs_inchikey_inchi(item, verbose = F)

    smilesFull <- c(smilesFull,getCanonicalSmiles(inchi))

    tmp <- cts_compinfo(item, verbose = F)
    tmp <- cts_compinfo("MGSRCZKZVOBKFT-UHFFFAOYSA-N", verbose = F)
    name <- getName(tmp,item)
    iupac <- getIUPACName(tmp,item)
    namesFull <- c(namesFull, name)
    iupacFull <- c(iupacFull, iupac)
    
  }
  # if user puts in incorrect chemical name (does not matter if common or IUPAC), this will fail
  else if(!is.na(cir_query(item, representation = 'smiles', verbose = F))){
    print(paste0("Input ",item," is chemical name.",sep=""))
    
    inchikey <- c()

    cnv <- cts_convert(item,'Chemical Name','inchikey', verbose = F)
    cnv <- as.data.frame(unlist(cnv))
    inchikey <- as.character(cnv[1,1])
    
    inchikeyFull <- c(inchikeyFull, inchikey)
    
    tmp <- cts_compinfo(inchikey = inchikey, verbose = F)
    
    name <- getName(tmp,inchikey)
    iupac <- getIUPACName(tmp,inchikey)
    namesFull <- c(namesFull, name)
    iupacFull <- c(iupacFull, iupac)
    
    smilesFull <- c(smilesFull, getCanonicalSmiles(cs_inchikey_inchi(inchikey = inchikey, verbose=F)))
  } else {
    stop("Supplied identifier is neither of SMILES, InChIKey or common/IUPAC name.")
  }
}

out <- data.frame(smiles=smilesFull,inchikey=inchikeyFull,common_name=namesFull,iupac_name=iupacFull)

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


