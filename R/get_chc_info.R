#' @title get the CHC information data frame
#'
#' @description Function to generate the data frame holding the compound
#' information for the hydrocarbons present in each sample, within a table
#'
#' @param comps.id
#' data frame with the compound information for the hydrocarbons found in each
#' sample.
#'
#' The data frame should have a **"Compound"** and a **"mean_RT"** columns, plus
#'  a column for every sample.
#'
#' The Compound column should have the name of the hydrocarbons with the
#' following format:
#' C*numeric chain length* _ *class abbreviation* _
#' *numeric position of unsaturations / methyl-branches*.
#'
#' As an example the entry for the compound 15-Me C20 would be C20_Me_15-.
#'
#' The function automatically extracts the information of the compound from the
#' entry in the **"Compound"** column, and replaces it by the correct name.
#'
#' If a peak contains more than one compound with unsaturations /
#' methyl-branches, they must be separated by a ";" (e.g. C20_Me_13-;15-
#' for a peak containing both 13-Me C20 and 15-Me C20).
#'
#'
#' ## Class abbreviations accepted by the function:
#' n-alkanes: ane
#' alkene: ene
#' alkadiene: diene
#' monomethyl-branched alkanes: Me
#' Dimethyl-branched alkanes: Dime
#' Trimethyl-branched alkanes: Trime
#'
#' The rest of the data frame should be just as it was extracted from the RT
#' entry of the aligned standards.
#'
#' @import dplyr
#' @import tidyr
#'
#' @export
get_chc_info <- function(comps.id) {
  # Extract compounds information from the comps table
  group.comps.info <- comps.id |>
    select(contains("Peak"), contains("Compound")) |>
    cbind(data.frame(row.names = rownames(comps.id)
                     , t(
                       as.data.frame(
                         strsplit(comps.id$Compound
                                  , "_")
                         )
                       )
                     )
          ) |>
    as_tibble()

  # cat('\n')
  # print("The compounds' information was extracted")
  # cat('\n')

  # Define columns' names
  colnames(group.comps.info) <-
    c(colnames(comps.id |>
                 select(contains("Peak"):contains("Compound")))
      , "Chain.length", "Class", "Mod.position")
  # cat('\n')
  # print("The columns of the comps_info data frame were renamed")
  # cat('\n')

  # Correct entries format
  # cat('\n')
  # print("Correcting the format of the entries in the data frame")
  # cat('\n')

  ## Compound names
  ### Store the compound names in a vector, where they will be altered into
  ### their final format
  comps.names <- group.comps.info$Compound

  ### Set the first part of the identified CHC names (full name in the case
  ### of the alkanes)
  comps.names[!is.na(comps.names)] <-
    paste0(group.comps.info |>
             filter(!is.na(get("Compound"))) |>
             pull("Mod.position")
           , group.comps.info |>
             filter(!is.na(get("Compound"))) |>
             pull("Class")
           , group.comps.info |>
             filter(!is.na(get("Compound"))) |>
             pull("Chain.length")) |>
    str_replace_all(paste(c("ane", "ene", "diene")
                          , collapse = '|')
                    , "NA") |>
    str_remove_all("NA")

  ### Finish formatting the names of the identified CHC
  comps.names[!is.na(comps.names)] <-
    paste(comps.names[!is.na(comps.names)]
          , group.comps.info |>
            filter(!is.na(get("Compound"))) |>
            pull("Class") |>
            str_replace(paste(c("ane"
                                , "Me"
                                , "Dime"
                                , "Trime")
                              , collapse = '|')
                        , "NA")
          , sep = ":") |>
    str_remove(":NA")

  ### Change the compound names in the comps_info data frame to
  ### the correct format names
  group.comps.info$Compound <- comps.names

  # cat('\n')
  # print("The names of the compoubnds were formated")
  # cat('\n')

  ## Chain length
  group.comps.info$Chain.length <- group.comps.info$Chain.length |>
    str_remove("C") |>
    as.integer()

  # cat('\n')
  # print("Format of Chain length was set as integer")
  # cat('\n')

  ## Class
  group.comps.info['Class'][group.comps.info['Class'] == "ane"] <-
    "Alkane"
  group.comps.info['Class'][group.comps.info['Class'] == "ene"] <-
    "Alkene"
  group.comps.info['Class'][group.comps.info['Class'] == "diene"] <-
    "Alkadiene"
  group.comps.info['Class'][group.comps.info['Class'] == "Me"] <-
    "Methyl"
  group.comps.info['Class'][group.comps.info['Class'] == "Dime"] <-
     "Dimethyl"
  group.comps.info['Class'][group.comps.info['Class'] == "Trime"] <-
    "Trimethyl"
  #  group.comps.info['Class'][group.comps.info['Class'] == "Tetrame"] <-
  #    "Tetramethyl"

  # cat('\n')
  # print("The class of the compounds was formatted")
  # cat('\n')

  group.comps.info |>
    as_tibble()
}
