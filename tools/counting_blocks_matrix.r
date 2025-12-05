#######################################################################
#' Counting blocks of non-zero rows in a matrix
#'
#'
#' This function counts the number of blocks of non-zero rows in a given matrix.
#' A block is defined as a contiguous sequence of non-zero rows.
#'
#' @param mat A matrix where rows are considered for counting blocks
#'
#' @return data frame of the block number, the number of rows in the block, and the number of columns in the block
#'
counting_blocks_matrix <- function(mat) {
  n <- nrow(mat)
  blocks <- list()
  visited <- rep(FALSE, n)  # Track visited rows
  block_count <- 0

  for (i in seq_len(n)) {
    if (!visited[i]) {
      #Find the start of a new block
      block_count <- block_count + 1
      block_rows <- which(mat[i, ] != 0)

      #Check for contiguous rows forming a block
      block_size <- 1
      while ((i + block_size) <= n && 
             all(mat[i + block_size, block_rows] != 0) &&
             all(mat[i + block_size, -block_rows] == 0)) {
        block_size <- block_size + 1
      }
      
      #Mark the rows as visited
      visited[i:(i + block_size - 1)] <- TRUE
      
      #Store block information
      blocks[[block_count]] <- list(block = block_count, rows = block_size, cols = length(block_rows))
    }
  }

  
  block_df <- do.call(rbind, lapply(blocks, as.data.frame))
  return(block_df)
}