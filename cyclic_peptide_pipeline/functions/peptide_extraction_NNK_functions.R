#Extract the peptide from a dna_seq sequence. regexPattern should be the sequence right before the peptides. 

#' Extract the peptide sequence from a read
#'
#' @param dna_seq The sequence from which peptide is extracted
#' @param regexPattern The DNA sequence immediately before the peptide sequence. Should be long enough to be unique
#' @param pepSize The size of the peptide gene, in nt
#'
#' @return
#' @export
#'
#' @examples
extract.peptide <- function(dna_seq, regexPattern = "TGGCTTCATTGCGAGCAAT", pepSize = 24){
  pepPosition <- str_locate(dna_seq %>% as.character(), regexPattern)[,2]
  sensiblePos <- pepPosition+pepSize < 145 #Added this for weird case where the pattern is present near the end
  peptideList <- subseq(x = dna_seq[sensiblePos], start = pepPosition[sensiblePos]+1, end = pepPosition[sensiblePos]+pepSize)
  rm(pepPosition, sensiblePos)
  return(peptideList)
}

#Extract all peptides from the fastq of a run. Assume a fastq.gz file that has been merged by FLASH

extract.peptides.fastq = function(file_name, 
                                  destination, 
                                  front_pattern = "TGGCTTCATTGCGAGCAAT",
                                  back_pattern = "TGTCTGTCTTACGACA",
                                  short_peptide_size = 12,
                                  long_peptide_size = 24
) {
  # Load a fastQ.gz file with the peptides. The files are too big to be used
  # entirely, so I'll have to use FastqStremer. .gz compressed files actually process faster than .fastq files.
  
  # I first need to get a connection established.
  # open the connection
  reads_per_batch <- 10^6
  stream <- FastqStreamer(file_name, n = reads_per_batch)
  on.exit(close(stream))
  
  #This is just to track the speed of processing of the file.
  i = 1
  ptm = proc.time()
  
  dna_seq <- yield(stream)
  
  # This DF will be used to write out the total number of peptides produced at the end.
  tracking_df <- data.frame(n_reads = 0,
                            n_short_peptide = 0,
                            n_long_peptide = 0)
  
  # This bit below checks before the while loop that output doesn't exist to avoid appending to an existing file
  
  if (file.exists(file.path(destination, 
                        gsub(basename(file_name), 
                             pattern = ".fastq.gz", 
                             replacement = paste("_peptide", 
                                                 long_peptide_size, 
                                                 ".fa", sep=""))))){
    stop("The output file for large peptides already exists.")
  }
  
  if (file.exists(file.path(destination, 
                            gsub(basename(file_name), 
                                 pattern = ".fastq.gz", 
                                 replacement = paste("_peptide", 
                                                     short_peptide_size, 
                                                     ".fa", sep=""))))){
    stop("The output file for short peptides already exists.")
  }
  
  while(length(dna_seq) > 0){

    length_dna_seq <- length(dna_seq)
    tracking_df <- tracking_df %>% mutate(n_reads = n_reads+length_dna_seq) # Tracking purpose
    
    # Filter the one that have a weird size and take the DNAStringSet
    dna_seq <- dna_seq[between(width(dna_seq),150,210)] %>% 
      sread()

    #I need to sort the fwd from reverse read.    
    # Get all the sequence that match the fwd pattern for intein
    
    fwd_dna_seq <- dna_seq[grepl(pattern = front_pattern, dna_seq %>%
                                as.character())]
    
    # Get all the sequence that match the complement of the front pattern for intein and take their reverse
    # complement
    
    rev_dna_seq = dna_seq[grepl(pattern = front_pattern %>% 
                                DNAString() %>%
                                reverseComplement(), 
                              dna_seq %>%
                                as.character())] %>%
                    reverseComplement()
    
    # Recreate dna_seq, from fwd_dna_seq and rev_dna_seq
    
    dna_seq <- c(fwd_dna_seq, rev_dna_seq)
    rm(fwd_dna_seq, rev_dna_seq)
    
    #Remove those that do not contain back_pattern, as I won't be able to extract the peptide

    dna_seq = dna_seq[grepl(pattern = back_pattern, dna_seq %>%
                                  as.character())]
    
    # I can get the size of the peptide by looking at the
    # intein sequence 
    
    front_pos = str_locate(dna_seq %>% as.character(), front_pattern)[,2]+1
    
    back_pos = str_locate(dna_seq %>% as.character(), back_pattern)[,1]-1
    
    make_sense <- front_pos < back_pos #Check the back is after the front
    
    peptide_sequence = subseq(x = dna_seq[make_sense], 
                        start =  front_pos[make_sense], 
                        end = back_pos[make_sense])
    
    short_peptide_sequence = peptide_sequence[width(peptide_sequence) == short_peptide_size]  
    
    long_peptide_sequence = peptide_sequence[width(peptide_sequence) == long_peptide_size]
    
    rm(peptide_sequence, dna_seq, front_pos, back_pos, make_sense)
    
    # I'll store all the sequence in a ShortRead object
    
    long_peptide_sequence = ShortRead(sread = long_peptide_sequence, 
                                      id = BStringSet(rep(paste(basename(file_name),
                                                                long_peptide_size,"nt peptide "),
                                                                   length(long_peptide_sequence))))
    
    
    short_peptide_sequence = ShortRead(sread = short_peptide_sequence, 
                                      id = BStringSet(rep(paste(basename(file_name),
                                                                short_peptide_size,"nt peptide "),
                                                          length(short_peptide_sequence))))
    
    tracking_df <- tracking_df %>% mutate(n_short_peptide = n_short_peptide + length(short_peptide_sequence),
                                          n_long_peptide = n_long_peptide + length(long_peptide_sequence),)
    
    
    writeFasta(short_peptide_sequence, 
               file.path(destination, 
                         gsub(basename(file_name), 
                              pattern = ".fastq.gz", 
                              replacement = paste0("_peptide", 
                                                  short_peptide_size, 
                                                  ".fa"))), 
               mode = "a")  # The mode append it to a file if existing
    
    writeFasta(long_peptide_sequence, 
               file.path(destination, 
                         gsub(basename(file_name), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  long_peptide_size, 
                                                  ".fa", sep=""))), 
               mode = "a")
    dna_seq <- yield(stream)
    
    # Give progress info
    
    paste("Out of ", length_dna_seq," reads in batch ", i, " of file ",basename(file_name)," ",length(short_peptide_sequence)+length(long_peptide_sequence)," encoded for peptides.") %>% print()
    print(proc.time() - ptm)
    i = i+1
    gc()
  }
  #Load the next batch of reads
  
  tracking_df %>% mutate(file = basename(file_name))
  
  return(tracking_df)
  
}

# Older function, that doesn't assume the reads have been merged but are a mix of forward and reverse reads.

extract.peptides.fastq.unmerged = function(file_name, 
                                  destination, 
                                  front_pattern = "TTCATTGCGAGCAAT",
                                  back_pattern = "TGTCTGTCTTACGACA",
                                  short_peptide_size = 12,
                                  long_peptide_size = 24
) {
  # Load a fastQ.gz file with the peptides. The files are too big to be used
  # entirely, so I'll have to use FastqStremer. .gz compressed files actually process faster than .fastq files.
  
  # I first need to get a connection established.
  # open the connection
  reads_per_batch <- 10^6
  stream <- FastqStreamer(file_name, n = reads_per_batch)
  on.exit(close(stream))
  
  #This is just to track the speed of processing of the file.
  i = 1
  ptm = proc.time()
  
  dna_seq <- yield(stream)
  
  # This DF will be used to write out the total number of peptides produced at the end.
  tracking_df <- data.frame(n_reads = 0,
                            n_short_peptide = 0,
                            n_long_peptide = 0)
  
  # This bit below checks before the while loop that output doesn't exist to avoid appending to an existing file
  
  if (file.exists(file.path(destination, 
                            gsub(basename(file_name), 
                                 pattern = ".fastq.gz", 
                                 replacement = paste("_peptide", 
                                                     long_peptide_size, 
                                                     ".fa", sep=""))))){
    stop("The output file for large peptides already exists.")
  }
  
  if (file.exists(file.path(destination, 
                            gsub(basename(file_name), 
                                 pattern = ".fastq.gz", 
                                 replacement = paste("_peptide", 
                                                     short_peptide_size, 
                                                     ".fa", sep=""))))){
    stop("The output file for short peptides already exists.")
  }
  
  while(length(dna_seq) > 0){
    # Filter the one that are way too small (<145 bp out of 150, about 0.5% of seq for
    # NNK7)
    length_dna_seq <- length(dna_seq)
    tracking_df <- tracking_df %>% mutate(n_reads = n_reads+length_dna_seq) # Tracking purpose
    
    dna_seq = dna_seq[width(dna_seq) > 145]
    
    # Convert to a DNAStringSet
    
    dna_seq = sread(dna_seq)
    
    # Get all the sequence that match the fwd pattern for intein
    
    fwd_dna_seq <- dna_seq[grepl(pattern = front_pattern, dna_seq %>%
                                as.character())]
    
    # I will subset the NNK7 and NNK3
    
    fwd_short_seq = fwd_dna_seq[grepl(pattern = "^TGCA(.)*", fwd_dna_seq %>%
                                    as.character())]
    
    fwd_long_seq = fwd_dna_seq[grepl(pattern = "^AAAA(.)*", fwd_dna_seq %>%
                                   as.character())]
    
    # And get the sequence of the peptides
    
    fwdShortPep = extract.peptide(dna_seq = fwd_short_seq, 
                                  regexPattern = front_pattern,
                                  pepSize = short_peptide_size)
    
    fwdLongPep = extract.peptide(dna_seq = fwd_long_seq, 
                                 regexPattern = front_pattern,
                                 pepSize = long_peptide_size)
    
    # Get all the sequence that match the rev pattern for intein and take their reverse
    # complement
    
    rev_dna_seq = dna_seq[grepl(pattern = back_pattern %>% 
                               DNAString() %>%
                               reverseComplement(), 
                             dna_seq %>%
                               as.character())] %>%
      reverseComplement()
    
    #Remove those that do not contain fronPattern or back_pattern
    rev_dna_seq = rev_dna_seq[grepl(pattern = front_pattern, rev_dna_seq %>%
                                  as.character())]
    rev_dna_seq = rev_dna_seq[grepl(pattern = back_pattern, rev_dna_seq %>%
                                  as.character())]
    
    # I do not have a barcode but I can get the size of the peptide by looking at the
    # intein sequence anyway
    
    front_pos = str_locate(rev_dna_seq %>% as.character(), front_pattern)[,2]+1
    
    back_pos = str_locate(rev_dna_seq %>% as.character(), back_pattern)[,1]-1
    
    make_sense <- front_pos < back_pos #Check the back is after the front
    
    revPeptide = subseq(x = rev_dna_seq[make_sense], 
                        start =  front_pos[make_sense], 
                        end = back_pos[make_sense])
    
    revShortPep = revPeptide[width(revPeptide) == short_peptide_size]  
    
    revLongPep = revPeptide[width(revPeptide) == long_peptide_size]  
    
    # I'll store all the sequence in a ShortRead object
    
    fwdLongPep = ShortRead(sread = fwdLongPep, id = BStringSet(rep(paste(long_peptide_size,"AA peptide fw"),
                                                                   length(fwdLongPep))))
    
    revLongPep = ShortRead(sread = revLongPep, id = BStringSet(rep(paste(long_peptide_size,"AA peptide rv"),
                                                                   length(revLongPep))))
    
    longPep <- append(fwdLongPep, revLongPep)
    
    
    fwdShortPep = ShortRead(sread = fwdShortPep, id = BStringSet(rep(paste(short_peptide_size,"AA peptide fw"),
                                                                     length(fwdShortPep))))
    
    revShortPep = ShortRead(sread = revShortPep, id = BStringSet(rep(paste(short_peptide_size,"AA peptide rv"),
                                                                     length(revShortPep))))
    
    shortPep <- append(fwdShortPep, revShortPep)
    
    tracking_df <- tracking_df %>% mutate(n_short_peptide = n_short_peptide + length(shortPep),
                                          n_long_peptide = n_long_peptide + length(longPep),)
    
    
    writeFasta(shortPep, 
               file.path(destination, 
                         gsub(basename(file_name), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  short_peptide_size, 
                                                  ".fa", sep=""))), 
               mode = "a")  # The mode append it to a file if existing
    
    writeFasta(longPep, 
               file.path(destination, 
                         gsub(basename(file_name), 
                              pattern = ".fastq.gz", 
                              replacement = paste("_peptide", 
                                                  long_peptide_size, 
                                                  ".fa", sep=""))), 
               mode = "a")
    dna_seq <- yield(stream)
    
    # Give progress info
    
    paste("Out of ", length_dna_seq," reads in batch ", i, " of file ",basename(file_name)," ",length(shortPep)+length(longPep)," encoded for peptides.") %>% print()
    print(proc.time() - ptm)
    i = i+1
    gc()
  }
  #Load the next batch of reads
  
  tracking_df %>% mutate(file = basename(file_name))
  
  return(tracking_df)
  
}