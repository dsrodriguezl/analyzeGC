# samples_data_list ----
#' GC/MS integration data of honey bees cuticular hydrocarbon extracts.
#'
#' Each data frame corresponds tot he integration results for the CHC extract of
#' a honey bee worker.
#'
#' @format ##
#' A list with 20 tibble data frames, each with 2 columns:
#' \describe{
#'   \item{RT}{Retention time measured for each peak within a sample}
#'   \item{Area}{Area result of the integration of the peaks within a sample}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"samples_data_list"

# standards_data_list ----
#' GC/MS integration data of n-alkane analytical standard solutions
#' (04070-1ML and 04071-5ML, Sigma-Aldrich).
#'
#' @format ##
#' A list with 2 tibble data frames, each with 2 columns:
#' \describe{
#'   \item{RT}{Retention time measured for each peak within a sample}
#'   \item{Area}{Area result of the integration of the peaks within a sample}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"standards_data_list"

# grouping_info ----
#' Group membership information of honeybee samples in the `samples_data_list`
#' data object.
#'
#' @format ##
#' A tibble data frame with 20 rows and 4 columns:
#' \describe{
#'   \item{Season}{Season during which the samples were collected (all in
#'   Winter).}
#'   \item{Task}{The type of worker (in- or out-hive) the bees were when
#'   collected.}
#'   \item{Subspecies}{Subspecies of the bees}
#'   \item{Individual}{Individual ID number}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"grouping_info"

# aligned_samples_data_list ----
#' GC/MS data of honey bees cuticular hydrocarbons aligned by group
#' (In- and out-hive workers).
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{Winter_In-hive workers_A. m. mellifera}{Alignment of the GC/MS
#'   integration dara of the in-hive workers.}
#'   \item{Winter_Out-hive workers_A. m. mellifera}{Alignment of the GC/MS
#'   integration dara of the out-hive workers.}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"aligned_samples_data_list"

# aligned_standards ----
#' GC/MS data of n-alkane analytical standard solutions aligned
#'
#'
#' @format ##
#' A GCalignR object.
#'
#' @source Data obtained by the author of the package, for educational purposes.
"aligned_standards"

# corrected_samples_list ----
#' Corrected alignment of the honey bees cuticular hydrocarbons GC/MS data.
#'
#' @format ##
#' A list with two sublists:
#' \describe{
#'   \item{Winter_In-hive workers_A. m. mellifera}{
#'   List with the RT and Area data frames for the corrected alignment of
#'   the in-hive workers.
#'   }
#'   \item{Winter_Out-hive workers_A. m. mellifera}{
#'   List with the RT and Area data frames for the corrected alignment of
#'   the out-hive workers.
#'   }
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"corrected_samples_list"

# corrected_samples_list2 ----
#' Corrected alignment of the honey bees cuticular hydrocarbons GC/MS data, with
#' corrected mean RT.
#'
#'
#' @format ##
#' A list with two sublists:
#' \describe{
#'   \item{Winter_In-hive workers_A. m. mellifera}{
#'   List with the RT and Area tibble data frames for the corrected alignment
#'   with recalculated mean RT of the in-hive workers.
#'   }
#'   \item{Winter_Out-hive workers_A. m. mellifera}{
#'   List with the RT and Area tibble data frames for the corrected alignment
#'   with recalculated mean RT of the out-hive workers.
#'   }
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"corrected_samples_list2"

# comps_id_std ----
#' Compound identification of the n-alkanes in the analytical standard
#' solutions (04070-1ML and 04071-5ML, Sigma-Aldrich).
#'
#' @format ##
#' A tibble data frame with 59 rows and 5 columns:
#' \describe{
#'   \item{Peak}{Peak labels}
#'   \item{Compound}{ID of the compound(s) contained by a peak}
#'   \item{mean_RT}{Mean RT of the peak among samples}
#'   \item{H2005}{Sample corresponding to an n-alkane 04071-5ML analytical
#'   standard solution.}
#'   \item{L2005}{Sample corresponding to an n-alkane 04070-1ML analytical
#'   standard solution.}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"comps_id_std"

# std_info ----
#' Standards information data.
#'
#'
#' @format ##
#' A tibble data frame with 30 rows and 9 columns:
#' \describe{
#'   \item{Compound}{ID of the compound(s) contained by a peak}
#'   \item{mean_RT}{Mean RT of the peak among samples}
#'   \item{area}{Measured area under the curve of a GC/MS peak, as obtained
#'   from the peak integration.}
#'   \item{Chain.length}{Integer indicating the chain length of the hydrocarbon
#'   molecule.}
#'   \item{Class}{Hydrocarbons class see `help("get_hc_info")` for details}
#'   \item{Mod.position}{Position of unsaturations or methyl-branches in the
#'   hydrocarbon chain}
#'   \item{median_area}{Median of the integrated area for all peaks in the data
#'   frame}
#'   \item{area_correction }{Proportional correction to apply to the absolute
#'   abundance (area) of a peak to match the median}
#'   \item{corrected_area}{Area under the peak after correction}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"std_info"

# comps_id_list ----
#' Identification of compounds in the cuticular hydrocarbons of the honey bees
#'
#' @format ##
#' A list with two tibble data frames:
#' \describe{
#'   \item{Winter_In-hive workers_A. m. mellifera}{
#'   Compound identification for the cuticular hydrocarbon extracts of the
#'   in-hive workers
#'   }
#'   \item{Winter_Out-hive workers_A. m. mellifera}{
#'   Compound identification for the cuticular hydrocarbon extracts of the
#'   out-hive workers
#'   }
#'   Both tibble data frames contain the labels of each peak (Peak), the ID of
#'   the compound(s) contained by a peak (Compound), mean RT of the peak
#'   (mean_RT), and RT of the peaks for each sample (ID number of each sample).
#'
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"comps_id_list"

# comps_info_list ----
#' Information of hydrocarbons in the cuticular hydrocarbon extracts of the
#' honey bees.
#'
#' @format ##
#' A list with two tibble data frames:
#' \describe{
#'   \item{Winter_In-hive workers_A. m. mellifera}{
#'   Information of compounds in the cuticular hydrocarbon extracts of the
#'   in-hive workers}
#'   \item{Winter_Out-hive workers_A. m. mellifera}{
#'   Information of compounds in the cuticular hydrocarbon extracts of the
#'   out-hive workers}
#'   Both tibble data frames they contain the labels of each peak (Peak),
#'   the ID of the hydrocarbons in the corresponding peak (Compound), the
#'   hydrocarbon chain length (Chain.length), the hydrocarbon class (Class),
#'   and the position of unsaturations or methyl-branches in the hydrocarbon
#'   chain (Mod.position).
#'   These tibble data frames were generated by applying `get_hc_info` on the
#'   tibble data frames contained in `comps_id_list`.
#'
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"comps_info_list"

# adjusted_samples_list ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"adjusted_samples_list"

# unfiltered_samples_list ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"unfiltered_samples_list"

# filtered_samples_list ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"filtered_samples_list"

# filtered_samples_list2 ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"filtered_samples_list2"

# samples_plus_ri_list ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"samples_plus_ri_list"

# group_tables_list ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"group_tables_list"

# master_table ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"master_table"

# duplicated_compounds_presence ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"duplicated_compounds_presence"

# master_table2 ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"master_table2"

# group_tables_list2 ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"group_tables_list2"

# group_tables_list2 ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"master_table_reassembled"

# master_table_transformed ----
#' GC/MS data of honeybees cuticular hydrocarbons aligned by group
#'
#'
#' @format ##
#' A list with two GCalignR objects:
#' \describe{
#'   \item{}{}
#'   \item{}{}
#' }
#' @source Data obtained by the author of the package, for educational purposes.
"master_table_final"
