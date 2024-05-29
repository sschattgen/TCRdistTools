require(tidyverse)
require(Biostrings)

#stuck at 1477

white_list_species <- c(
  "Homo sapiens", 
  "Mus musculus", 
  "Macaca mulatta", 
  "Bos taurus"
  )

#download imgt refs ====
tmp <- tempdir()
AAseq_file <- paste(tmp, 'IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP', sep = "")

download.file("https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP", 
              AAseq_file, mode = "wb")

nucseq_file <-  paste(tmp, 'IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP', sep = "")
download.file("https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP", 
              nucseq_file, mode = "wb")

#read in and filter ====
in_AA <- readAAStringSet(AAseq_file)
in_NT <- readDNAStringSet(nucseq_file)

keep_list <- which(str_split(names(in_AA), pattern = '[|]',simplify = T)[,3] %in% white_list_species)

in_AA <-in_AA[keep_list]
in_NT <-in_NT[keep_list]

col_vals <- c('id','organism','chain',
              'region','nucseq','frame',
              'aligned_protseq',
              'cdr_columns','cdrs')

rough_stats_df <- as.data.frame(
  str_split(
    names(in_AA), 
    pattern = '[|]',
    simplify = T)
  )

# functions to parse and annotate entries ====
species_assigner <- function(gene, species){
  
  chain2 <- case_when(
    grepl('IG',gene) ~ '_ig',
    grepl('TR[DG]',gene) ~ '_gd',
    .default = '')
    

  species2 <- case_when(
    species =='Homo sapiens' ~'human',
    species =='Mus musculus' ~'mouse',
    species =='Macaca mulatta' ~'rhesus',
    species =='Bos taurus' ~ 'bovine',
    )
  new_species <- paste0(species2, chain2)
  

  return(new_species)
}

V_gene_parser <- function(seq, species, len){ 
  
  cdr1 <- substr(seq, 27, 38) #IMGT gap settings
  cdr2 <-substr(seq, 56, 65) #IMGT gap settings
  cdr2.5 <- substr(seq, 81,86) #IMGT gap settings
  cdr3_start <- tail(str_locate_all(seq, 'C')[[1]][,1], n=1)
  cdr3 <- substr(seq, cdr3_start, len)
  #gaps_to_add <- (cdr3_stop-cdr3_start)-nchar(cdr3)
  
  #cdr3 <- paste0(cdr3, paste(rep('.',gaps_to_add),collapse = ''))
  
  CDR_columns <- paste(
    paste(27, 38, sep = '-'),
    paste(56, 65, sep = '-'),
    paste(81, 86, sep = '-'),
    paste(cdr3_start, len, sep = '-'),
    sep = ';'
  )
  if (grepl('ig',species)){
    CDRs <- paste(cdr1,cdr2,cdr2.5,cdr3,sep = ';')
  } else {
    CDRs <- paste(cdr1,cdr2,cdr3,sep = ';')
  }
  
  
  return(list(CDR_columns, CDRs))
} #something off with gaps here

J_gene_parser <- function(seq, species, chain){
  
  ending <- case_when(
    species =='human' & chain == 'A' ~ 12,
    species =='human' & chain == 'B' ~ 8,
    species =='mouse' & chain == 'A' ~ 13,
    species =='mouse' & chain == 'B' ~ 8,
    species =='rhesus' & chain == 'A' ~ 11,
    species =='rhesus' & chain == 'B' ~ 8,
    species =='bovine' & chain == 'A' ~ 9,
    species =='bovine' & chain == 'B' ~ 7,
    species =='mouse_gd' & chain == 'A' ~ 9,
    species =='mouse_gd' & chain == 'B' ~ 9,
    species =='human_gd' & chain == 'A' ~ 11,
    species =='human_gd' & chain == 'B' ~ 9,
    species =='rhesus_gd' & chain == 'A' ~ 6,
    species =='rhesus_gd' & chain == 'B' ~ 7,
    species =='bovine_gd' & chain == 'B' ~ 8,
    species =='human_ig' & chain == 'A' ~ 11,
    species =='human_ig' & chain == 'B' ~ 9,
    species =='mouse_ig' & chain == 'A' ~ 3,
    species =='mouse_ig' & chain == 'B' ~ 5,
    species =='bovine_ig' & chain == 'A' ~ 3,
    species =='bovine_ig' & chain == 'B' ~ 5
  )
  

  seq_out <- subseq(seq, 1, ending)
  CDR_columns <- paste0('1-', ending)

  
  return(list(CDR_columns, seq_out))
}


#run through the list and parse entries
paired_NT_AA <- list()
for (i in seq_along(in_AA)){

  name_in <- str_split(names(in_AA)[i], pattern = '[|]',simplify = T)
  
  if( name_in[[1,4]] == "F" & 
      grepl("partial",name_in[[1,14]]) == F &
      grepl('[VJ]-REGION', name_in[[1,5]]) ) {
    
    gene <- name_in[[1,2]]
    speciesx <- name_in[[1,3]]
    frame <- name_in[[1,8]]
    
    len <- as.integer(str_split(name_in[[1,13]], '=', simplify = T)[,2])
    
    chain <- case_when(
      grepl('TR[AD][VJ]',gene) | grepl('IG[KL][VJ]',gene) ~ 'A',
      grepl('TR[BG][VJ]',gene) | grepl('IGH[VJ]',gene) ~ 'B',
      .default = 'X'
    )
    
    species <- species_assigner(gene,speciesx)
    
    if (chain %in% c('A','B')){
      
      region <- gsub('-REGION','',name_in[[1,5]])
      seq <- as.character(in_AA)[[i]]
      if (region =='V'){
        CDR_info <- V_gene_parser(seq, species, len)
      } else {
        CDR_info <- J_gene_parser(seq, species, chain)
      } 
      
      
      out <- c(gene, 
               species, 
               chain, 
               region, 
               tolower(as.character(in_NT[[i]])),
               paste0('+',frame),
               seq,
               CDR_info[[1]],
               CDR_info[[2]])
      
      names(out) <- col_vals
      
      list_length <- length(paired_NT_AA) + 1
      paired_NT_AA[[list_length]] <- out
    }
    
  }
  
  
}

#bind back together
ref_db <- do.call(rbind, paired_NT_AA)

#runs through now. bit slow. 
#might prefilter to just VDJ to save some time.
#need to add D region
#need to duplicate the TRAVs with the '_gd' added to the species



