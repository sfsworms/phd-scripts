primer_seq <- function(pep_seq){
  pep_seq <- pep_seq %>%
    DNAString() %>%
    reverseComplement() %>%
    as.character()
  prim_seq <- paste0("GAATGCCATATTCAACCGTCAGGATTTCGGTGTCGTAAGACAGACA",
                       pep_seq,
                       "GCAATTGCTCGCAATGAA")
  return(prim_seq)
}
