require(tidyverse)
require(Biostrings)

# define the species to process and some other info ====

white_list_species <- c(
  "Homo sapiens", 
  "Macaca mulatta_rheMacS_CGG.01",
  "Macaca mulatta_RUp15",
  "Macaca mulatta_AG07107",
  "Macaca mulatta_17573",
  "Macaca mulatta",
  "Bos taurus",
  "Bos taurus_Holstein",
  "Bos taurus_Hereford", 
  "Mus musculus_C57BL/10",
  "Mus musculus_C57BL/6",
  "Mus musculus_C57BL/6J",
  "Mus musculus_BALB/cJ",
  "Mus musculus_BALB/c",
  "Mus musculus_C3H",
  "Mus musculus_B10.A",
  "Mus musculus_A/J",
  "Mus musculus_C58",
  "Mus musculus_SK",
  "Mus musculus_DBA/2J",
  "Mus musculus_SJL"
  
  )

stat_df_cols <- c('accession','gene','species',
                  'functionality','label','nuc_start_stop',
                  'nuc_length','codon_start','5p_nt_add',
                  '3p_nt_sub','n_nt_errors','aa_length',
                  'aa_gap_length','partial','rev_complement')

col_vals <- c('id','organism','chain',
              'region','nucseq','frame',
              'aligned_protseq',
              'cdr_columns','cdrs','gene_family')


# download imgt refs ====

tmp <- tempdir()
AAseq_file <- paste(tmp, 'IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP', sep = "")

download.file("https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP", 
              AAseq_file, mode = "wb")

nucseq_file <-  paste(tmp, 'IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP', sep = "")
download.file("https://www.imgt.org/download/GENE-DB/IMGTGENEDB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP", 
              nucseq_file, mode = "wb")

# read in and filter ====

in_AA <- readAAStringSet(AAseq_file)
in_NT <- readDNAStringSet(nucseq_file)

# filter the entries on species, functionality, gene type

keep_list <- which(
  str_split(names(in_AA), pattern = '[|]',simplify = T)[,3] %in% white_list_species & #filter on species
    grepl('[VDJ]-REGION', str_split(names(in_AA), pattern = '[|]',simplify = T)[,5]) & #filter on gene regions
    grepl("partial", str_split(names(in_AA), pattern = '[|]',simplify = T)[,14]) == F & #filter out partial
    str_split(names(in_AA), pattern = '[|]',simplify = T)[,4] == 'F'  #filter on functional
    )

in_AA <-in_AA[keep_list]
in_NT <-in_NT[keep_list]

# stats ====

rough_stats_df <- as.data.frame(
  str_split(
    names(in_AA), 
    pattern = '[|]',
    simplify = T)
) %>%
  select(-V16) %>%
  mutate(V13 = as.integer(str_split(V13,'=',simplify = T)[,2]))

colnames(rough_stats_df) <- stat_df_cols

filtered_stats_df <- rough_stats_df %>%
  mutate(gene_family = str_extract(gene,'^[TI][RG][ABGDKLH][VJD]')) %>%
  mutate(species2 = case_when(
    grepl('Homo sapiens', species) ~ 'human',
    grepl('Mus musculus', species) ~ 'mouse',
    grepl('Macaca mulatta', species) ~ 'rhesus',
    grepl('Bos taurus',species) ~ 'bovine',
    .default = species))


# the max length per species/gene class varies. Here we're making a table of max
# lengths so we can add the gaps in

length_df <- filtered_stats_df %>%
  group_by(species2, gene_family)  %>%
  summarise(max_length = max(aa_gap_length)) %>% 
  pivot_wider( values_from = 'max_length', names_from = 'gene_family')

length_mat <- as.matrix(length_df[,2:ncol(length_df)])
rownames(length_mat) <- length_df[[1]]

# functions to parse and annotate entries ====
species_assigner <- function(gene, species){
  
  chain2 <- case_when(
    grepl('IG',gene) ~ '_ig',
    grepl('TR[DG]',gene) ~ '_gd',
    .default = '')
    

  species2 <- case_when(
    species =='Homo sapiens' ~'human',
    grepl('Mus musculus', species) ~ 'mouse',
    grepl('Macaca mulatta', species) ~ 'rhesus',
    grepl('Bos taurus',species) ~ 'bovine'
    )
  
  new_species <- paste0(species2, chain2)
  

  return(new_species)
}

V_gene_parser <- function(seq, gene_family, species){ # 
  
  species3 <- str_split(species,'_', simplify = T)[,1]
  end_length <- length_mat[[species3,gene_family]]
  
  cdr1 <- substr(seq, 27, 38) #IMGT gap settings
  cdr2 <-substr(seq, 56, 65) #IMGT gap settings
  cdr2.5 <- substr(seq, 81,86) #IMGT gap settings
  cdr3_start <- tail(str_locate_all(seq, 'C')[[1]][,1], n=1)
  cdr3 <- substr(seq, cdr3_start, nchar(seq))
  gaps_to_add <- (end_length-(cdr3_start-1))-nchar(cdr3)
  cdr3 <- paste0(cdr3, paste(rep('.',gaps_to_add),collapse = ''))
  
  CDR_columns <- paste(
    paste(27, 38, sep = '-'),
    paste(56, 65, sep = '-'),
    paste(81, 86, sep = '-'),
    paste(cdr3_start, end_length, sep = '-'),
    sep = ';'
  )
  
  if (grepl('ig',species)){
    CDRs <- paste(cdr1,cdr2,cdr2.5,cdr3,sep = ';')
  } else {
    CDRs <- paste(cdr1,cdr2,cdr3,sep = ';')
  }
  
  
  return(list(CDR_columns, CDRs))
} 

J_gene_parser <- function(seq, gene_family, species){

  speciesX <- str_split(species,'_', simplify = T)[,1]
  max_length <- length_mat[[speciesX, gene_family]]
  
  if (nchar(seq) < max_length){
    gaps_to_add <- max_length-nchar(seq)
    seq <- paste0(paste(rep('.',gaps_to_add),collapse = ''), seq)
  }
  
  ending <- CDR3ending(species, gene_family)
  seq_out <- subseq(seq, 1, ending)
  CDR_columns <- paste0('1-', ending)

  return(list(CDR_columns, seq_out))
}

D_gene_parser <- function(seq, gene_family, species){

  speciesX <- str_split(species,'_', simplify = T)[,1]
  max_length <- length_mat[[speciesX, gene_family]]
  
  if (nchar(seq) < max_length){
    gaps_to_add <- max_length-nchar(seq)
    seq <- paste0(seq, paste(rep('.',gaps_to_add),collapse = ''))
  }
  
  CDR_columns <- ''
  
  return(list(CDR_columns, seq))
}

# table defining the end bounds of the CDR3 in the J segments
# for each species/gene family. The position is fixed within a combo but differ across (like the lengths)
# Unlike the lengths, these are defined manually. 
# Adding a new species would require updating the case switch in the fxn below.

CDR3ending <- function(species, chain) {
  case_when(
    species =='bovine' & chain == 'IGHJ' ~ 7,
    species =='bovine' & chain %in% c('IGKJ','IGLJ') ~ 3,
    species =='bovine' & chain == 'TRAJ' ~ 12,
    species =='bovine' & chain %in% c('TRBJ','TRGJ') ~ 8,
    species =='bovine' & chain == 'TRDJ' ~ 9,
    species =='human' & chain == 'IGHJ' ~ 10,
    species =='human' & chain %in% c('IGKJ','IGLJ') ~ 3,
    species =='human' & chain == 'TRAJ' ~ 12,
    species =='human' & chain == 'TRBJ' ~ 8,
    species =='human' & chain == 'TRDJ' ~ 9,
    species =='human' & chain == 'TRGJ' ~ 11, 
    species =='rhesus' & chain == 'IGHJ' ~ 7,
    species =='rhesus' & chain %in% c('IGKJ','IGLJ') ~ 3,
    species =='rhesus' & chain == 'TRAJ' ~ 12,
    species =='rhesus' & chain == 'TRBJ' ~ 8,
    species =='rhesus' & chain == 'TRDJ' ~ 9,
    species =='rhesus' & chain == 'TRGJ' ~ 11, 
    species =='mouse' & chain == 'IGHJ' ~ 7,
    species =='mouse' & chain %in% c('IGKJ','IGLJ') ~ 3,
    species =='mouse' & chain == 'TRAJ' ~ 12,
    species =='mouse' & chain == 'TRBJ' ~ 7,
    species =='mouse' & chain == 'TRDJ' ~ 9,
    species =='mouse' & chain == 'TRGJ' ~ 8
  )[1]
}

# run through the list and parse entries ====
paired_NT_AA <- list()
for (i in seq_along(in_AA)){

  name_in <- str_split(names(in_AA)[i], pattern = '[|]',simplify = T)

  gene <- name_in[[1,2]]
  speciesx <- name_in[[1,3]]
  frame <- name_in[[1,8]]
  gene_family <- str_extract(gene,'^[TI][RG][ABGDKLH][VJD]')
  len <- as.integer(str_split(name_in[[1,13]], '=', simplify = T)[,2])
  
  chain <- case_when(
    grepl('TR[AD][VDJ]',gene) | grepl('IG[KL][VDJ]',gene) ~ 'A',
    grepl('TR[BG][VDJ]',gene) | grepl('IGH[VDJ]',gene) ~ 'B',
    .default = 'X'
  )
  
  species <- species_assigner(gene,speciesx)
  
  if (chain %in% c('A','B')){
    
    region <- gsub('-REGION','',name_in[[1,5]])
    seq <- as.character(in_AA)[[i]]
    if (region =='V'){
      CDR_info <- V_gene_parser(seq, gene_family, species)
    } else if (region =='J') {
      CDR_info <- J_gene_parser(seq, gene_family, species)
    } else {
      CDR_info <- D_gene_parser(seq, gene_family, species)
    }
    
    
    out <- c(gene, 
             species, 
             chain, 
             region, 
             tolower(as.character(in_NT[[i]])),
             paste0('+',frame),
             seq,
             CDR_info[[1]],
             CDR_info[[2]],
             gene_family)
    
    names(out) <- col_vals
    
    list_length <- length(paired_NT_AA) + 1
    paired_NT_AA[[list_length]] <- out
    
    
  }
  
}

# bind back together
ref_db <- as.data.frame(do.call(rbind, paired_NT_AA))

# duplicate the TRAVs with the '_gd' added to the species

TRAVtogd <- ref_db %>% 
  filter(gene_family == 'TRAV') %>%
  mutate(chain = 'B') %>%
  mutate(organism = paste0(organism,'_gd'))

# full ref
full_ref_db <- bind_rows(ref_db, TRAVtogd) %>%
  arrange(organism, chain, region)

# write out
write_tsv(full_ref_db, 'combo_xcr.tsv')

ref_stats <- full_ref_db %>%
  group_by(organism, chain, region) %>%
  tally(name = 'count')

write_tsv(ref_stats, 'combo_xcr_stats.tsv')
