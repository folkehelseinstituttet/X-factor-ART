# DESCRIPTION:
# function to grab regulatory regions from ensembl - using biomaRt
#' @param positions - positions of markers for which to search the
#'    regions in the mart, given as a data.frame with columns named as
#'    filters that are given below
#' @param mart - which database to use - check description in biomaRt::getBM
#' @param attribs - which attributes to fetch - check description in biomaRt::getBM
#' @param filters - which filters to use - check description in biomaRt::getBM
#' @param asGRanges - whether to return the results as GRanges (check
#'    regioneR::GRanges for help) - boolean, defaults to FALSE
#' 

grabRegulRegions <- function(
  positions,
  mart,
  attribs,
  filters,
  asGRanges = FALSE
){
  if(!requireNamespace('biomaRt')){
    stop("Couldn't find 'biomaRt' package - please install and try once more.",
         call. = FALSE)
  }
  requireNamespace('biomaRt')
  
  if(!all(names(positions) %in% filters)){
    stop("Columns of 'positions' data.frame should be the same as filters!",
         call. = FALSE)
  }
  
  regulatory_regs <- purrr::map(seq_len(nrow(positions)), function(row){
    cur_cpg <- as.list(positions[row, ])
    cur_out <- biomaRt::getBM(
      attributes = attribs,
      filters = filters,
      values = cur_cpg,
      mart = mart
    )
    if(nrow(cur_out) == 0){
      return(NULL)
    }
    return(cur_out)
  })
  regulatory_regs <- dplyr::bind_rows(regulatory_regs)
  
  if(asGRanges){
    if(nrow(regulatory_regs) == 0){
      return(NULL)
    }
    regulatory_regs <- regioneR::toGRanges(
      as.data.frame(
        regulatory_regs %>%
          mutate(seqnames = paste0("chr", chromosome_name)) %>%
          dplyr::select(
            seqnames,
            start = chromosome_start,
            end = chromosome_end,
            everything()
          ) %>%
          distinct()
      )
    )
  }
  return(regulatory_regs)
}
