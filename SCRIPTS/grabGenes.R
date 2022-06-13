# DESCRIPTION:
#   function to find genes positioned in the given region; uses karyoploteR
#   TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, and regioneR
#' @param positions - data.frame with 'chromosome_name", 'start', and 'end'
#'   columns (all other columns will be ignored)
#' @param asGRanges - whether to return the results as GRanges (check
#'    regioneR::GRanges for help) - boolean, defaults to FALSE
#' 

grabGenes <- function(
  positions,
  asGRangesList = TRUE
){
  requireNamespace('org.Hs.eg.db', quietly = TRUE)
  requireNamespace('AnnotationDbi', quietly = TRUE)
  requireNamespace('Biobase', quietly = TRUE)
  
  # check the format of 'chromosome_name'
  positions_cp <- positions %>%
    mutate(chromosome_name = ifelse(
      grepl("chr", chromosome_name),
      yes = chromosome_name,
      no = paste0("chr", chromosome_name)
    ))

  genes_data_list <- purrr::map(1:nrow(positions_cp), function(row){
    cur_row <- positions_cp[row,]
    positions_ranges <- regioneR::toGRanges(
      data.frame(
        chr = cur_row$chromosome_name,
        start = cur_row$start,
        end = cur_row$end
      )
    )
    pdf(NULL)
    kp_tmp <- karyoploteR::plotKaryotype(
      chromosomes = cur_row$chromosome_name,
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
        as_tibble(d$genes) %>%
          dplyr::select(chr = seqnames, start, end, strand, gene_name = name)
      )
    }) %>% bind_rows() %>%
      distinct()
      return(genes_df)
  }
  return(genes_data_list)
}
