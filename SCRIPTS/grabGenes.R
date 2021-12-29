# DESCRIPTION:
#   function to find genes positioned in the given region; uses karyoploteR
#   TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, and regioneR
#' @param positions - data.frame with 'start' and 'end' columns (any other
#'   columns will be ignored)
#' @param chromosome - which chromosome? string in format "chrA", where 'A'
#'   needs to be one of 1-22 or X or Y
#' 

grabGenes <- function(
  positions,
  chromosome,
  asGRangesList = TRUE
){
  requireNamespace('org.Hs.eg.db', quietly = TRUE)
  requireNamespace('AnnotationDbi', quietly = TRUE)
  requireNamespace('Biobase', quietly = TRUE)
  
  genes_data_list <- purrr::map(1:nrow(positions), function(row){
    cur_row <- positions[row,]
    positions_ranges <- regioneR::toGRanges(
      data.frame(
        chr = chromosome,
        start = cur_row$start,
        end = cur_row$end
      )
    )
    pdf(NULL)
    kp_tmp <- karyoploteR::plotKaryotype(
      chromosome = chromosome,
      zoom = positions_ranges
    )
    genes.data <- suppressMessages(karyoploteR::makeGenesDataFromTxDb(
      TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
      karyoplot = kp_tmp
    ))
    if(!is.null(genes.data)){
      genes.data <- karyoploteR::addGeneNames(genes.data)
      out <- karyoploteR::mergeTranscripts(genes.data)
    } else {
      out <- NULL
    }
    dev.off()
    return(out)
  })
  
  if(!asGRangesList){
    genes_df <- purrr::map(genes_data_list, function(d){
      if(length(d$genes) == 0){
        return(NULL)
      }
      return(
        tibble(
          pos = GenomicRanges::start(d$genes),
          end = GenomicRanges::end(d$genes),
          gene_name = GenomicRanges::mcols(d$genes)$name
        )
      )
    }) %>% bind_rows() %>%
      distinct()
      return(genes_df)
  }
  return(genes_data_list)
}
