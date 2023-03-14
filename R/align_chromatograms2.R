#' @title Aligning peaks based on retention times
#'
#' @description Custom function for aligning chromatograms.
#' It is a wrapper around GCalignR::align_chromatograms function, thus it is
#' recommended to check its documentation (?GCalignR::align_chromatograms).
#'
#' Its main purpose is to facilitate the usage of the
#' function with lapply, and within portable scripts.
#'
#' It sets a seed before applying using GCalignR::align_chromatograms.
#'
#' The name of the alignment criteria parameters are more similar to their
#' description in the original paper on the GCalignR package.
#'
#' rt_col_name is set to "RT" by default, as it is the name given by
#' import_gcms_data.
#'
#' If the given data set contain only one sample, it returns a list with a
#' similar structure to that of the typical GCalignR::align_chromatograms
#' output. However, this list is not a GCalignR object.
#'
#' @param data2align
#' Dataset containing peaks that need to be aligned and matched.
#'
#' @param blanks
#' Character vector of names of negative controls.
#'
#' @param linear_shift_criteria
#' Numeric value indicating the size range within the linear shift of peaks can
#' occur during the initial full alignment step.
#'
#'
#' @param partial_alignment_threshold
#' Numeric value indicating the threshold for the partial peak alignment step.
#'
#' @param row_merging_threshold
#' Numeric value indicating the threshold for the merging rows step.
#'
#' @import tidyr
#' @import dplyr
#' @import GCalignR
#'
#' @examples
#'
#' library(tidyr)
#'
#' grouping_info <- grouping_info |>
#'   unite(group_label
#'         , where(is.factor)
#'         , sep = "_"
#'         , remove = FALSE)
#'
#' # Nest the data frames in sublists by group
#' samples_data_list <- mg_list(sample.info = grouping_info
#'                              , group.label = "group_label"
#'                              , samples.data.list = samples_data_list)
#'
#' # Aligning a single list of data frames
#' Winter_IW <- samples_data_list$`Winter_In-hive workers_A. m. mellifera`
#'
#' aligned_samples_data_list <- Winter_IW |>
#'   align_chromatograms2(blanks = NULL
#'                        , linear_shift_criteria = 0.02
#'                        , partial_alignment_threshold = 0.05
#'                        , row_merging_threshold = 0.15)
#'
#' # Using lapply to align all the data lists nested within a list
#' aligned_samples_data_list <- samples_data_list |>
#'   lapply(align_chromatograms2
#'          , blanks = NULL
#'          , linear_shift_criteria = 0.02
#'          , partial_alignment_threshold = 0.05
#'          , row_merging_threshold = 0.15)
#'
#' # Parallelized alignment of several data group lists nested within a list
#' library(doParallel)
#'
#' ## Define the number of cores that will be used for the parallel processes
#' no_cores <- 2
#'
#' ## Set a processing cluster
#' cl <- makeCluster(no_cores)
#'
#' ## Register the processing cluster to be used for parallel operations
#' registerDoParallel(cl)
#'
#' aligned_samples_data_list <- parLapply(cl
#'                                        , samples_data_list
#'                                        , align_chromatograms2
#'                                        , blanks = NULL
#'                                        , linear_shift_criteria = 0.02
#'                                        , partial_alignment_threshold = 0.05
#'                                        , row_merging_threshold = 0.15)
#'
#' stopCluster(cl)
#' registerDoSEQ()
#'
#' @export
align_chromatograms2 <- function(data2align
                                 , blanks = NULL
                                 , linear_shift_criteria
                                 , partial_alignment_threshold
                                 , row_merging_threshold){
  if (length(data2align) == 1) {
    nombre <- names(data2align)

    RT <- data2align[[1]] |>
      mutate("mean_RT" = get("RT")) |>
      select(contains("mean_RT"), contains("RT")) |>
      as.data.frame()
    colnames(RT) <- c("mean_RT", nombre)

    Area <- data2align[[1]] |>
      mutate("mean_RT" = get("RT")) |>
      select(contains("mean_RT"), contains("Area")) |>
      as.data.frame()
    colnames(Area) <- c("mean_RT", nombre)

    df <- list("RT" = RT, "Area" = Area)
  }

  if (length(data2align) > 1) {
    withr::local_seed(12345)
    df <- align_chromatograms(data = data2align
                              , rt_col_name = "RT"
                              , max_linear_shift = linear_shift_criteria
                              , max_diff_peak2mean = partial_alignment_threshold
                              , min_diff_peak2peak = row_merging_threshold
                              , blanks = blanks
    )
  }

  df
}
