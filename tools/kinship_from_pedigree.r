###########################################################################
#' @title Kinship Matrix from Pedigree
#'
#' @description Generates a kinship matrix from a pedigree data frame. 
#' The matrix represents the relatedness between individuals based on shared parents.
#'
#' @usage kinship_from_pedigree(pedigree)
#'
#' @param pedigree A data frame with three columns: 
#' \itemize{
#'   \item{\code{id}: Individual IDs (unique for each individual).}
#'   \item{\code{sire}: IDs of the sire (father); NA for founders.}
#'   \item{\code{dam}: IDs of the dam (mother); NA for founders.}
#' }
#'
#' @return A symmetric kinship matrix where rows and columns are labeled with individual IDs. 
#' Diagonal elements represent self-relatedness (1), and off-diagonal elements represent relatedness between individuals.
#' 
#' @details The function computes kinship coefficients under the assumption of non-inbred founders.
#' Half-sibling relationships (shared sire or dam) yield a coefficient of \eqn{\frac{1}{4}}. 
#' If individuals share both parents (full siblings), their kinship coefficient is \eqn{\frac{1}{4} + \frac{1}{4} = \frac{1}{2}}.
#'
#' @examples
#'
#'
#' @author Isabela do O \email{isabela.doo@@unil.ch}
#' 
#' @references 
#' 
kinship_from_pedigree <- function(pedigree) {
  
  ids <- pedigree$id
  num_individuals <- length(ids)
  kinship_matrix <- matrix(0, nrow = num_individuals, ncol = num_individuals)
  rownames(kinship_matrix) <- ids
  colnames(kinship_matrix) <- ids
  
  #Self-relatedness = 1 
  diag(kinship_matrix) <- 1
  
  
  for (i in seq_len(num_individuals)) {
    for (j in seq_len(i - 1)) { 
      sire_i <- pedigree$sire[i]
      dam_i <- pedigree$dam[i]
      sire_j <- pedigree$sire[j]
      dam_j <- pedigree$dam[j]
      
      if (!is.na(sire_i) && !is.na(sire_j) && sire_i == sire_j) {
        kinship_matrix[i, j] <- kinship_matrix[j, i] <- kinship_matrix[i, j] + (1/4) #half sib
      }
      if (!is.na(dam_i) && !is.na(dam_j) && dam_i == dam_j) {
        kinship_matrix[i, j] <- kinship_matrix[j, i] <- kinship_matrix[i, j] + (1/4) #half sib
      }
      #If you're "halfsib" on both sides, then  you're full sib!
    }
  }
  
  return(kinship_matrix)
}