create_F1 <- function(dat, dat_quanti, pedi){
  
  dos <- biall2dos(dat[, -1])
  nl <- ncol(dos)

  dos_quanti <- biall2dos(dat_quanti[, -1])
  nl_quanti <- ncol(dos_quanti)

  #Combine neutral and quantitative loci - quantinemo data
  combined_nl <- nl + nl_quanti
  combined_founders <- array(0, dim = c(nft, combined_nl, 2))
  combined_dat <- cbind(dat[, -1], dat_quanti[, -1])

  tmp1 <- as.matrix(combined_dat %/% 10 - 1)
  tmp2 <- as.matrix(combined_dat %% 10 - 1)
  combined_founders[, , 1] <- tmp1
  combined_founders[, , 2] <- tmp2

  #F1 generation
  nfounders_combined <- rbind(combined_founders[, , 1], combined_founders[, , 2])

  genos_F1_combined <- drop.along.ped(pedi, founders.genotypes = nfounders_combined, nloc = combined_nl, maplength = 1000)

  #Extract neutral and quantitative genotypes
  genos_F1_neutral <- genos_F1_combined[, 1:nl, ]
  genos_F1_quanti <- genos_F1_combined[, (nl + 1):combined_nl, ]

  dos_F1_neutral <- genos_F1_neutral[, , 1] + genos_F1_neutral[, , 2]
  bed_F1_neutral <- as.bed.matrix(dos_F1_neutral)

  dos_F1_quanti <- genos_F1_quanti[, , 1] + genos_F1_quanti[, , 2]
  bed_F1_quanti <- as.bed.matrix(dos_F1_quanti)
  
  data_rep <- list(dos_F1_neutral = dos_F1_neutral,
                  dos_F1_quanti = dos_F1_quanti)

}