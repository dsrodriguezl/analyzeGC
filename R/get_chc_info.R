#' @title
#'
#' @description
#'
#' @param group.comps.id
#'
#' @export
get_chc_info <- function(group.comps.id) {
  # Extract compounds information from the comps table
  group.comps.info <- group.comps.id |>
    select(Peak:Compound) |>
    cbind(data.frame(row.names = rownames(group.comps.id)
                     , t(
                       as.data.frame(
                         strsplit(group.comps.id$Compound
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
  colnames(group.comps.info) <- c(colnames(group.comps.id |>
                                             select(Peak:Compound)
  )
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
             filter(!is.na(Compound)) |>
             pull("Mod.position")
           , group.comps.info |>
             filter(!is.na(Compound)) |>
             pull("Class")
           , group.comps.info |>
             filter(!is.na(Compound)) |>
             pull("Chain.length")) |>
    str_replace_all(paste(c("ane", "ene", "diene")
                          , collapse = '|')
                    , "NA") |>
    str_remove_all("NA")

  ### Finish formatting the names of the identified CHC
  comps.names[!is.na(comps.names)] <-
    paste(comps.names[!is.na(comps.names)]
          , group.comps.info |>
            filter(!is.na(Compound)) |>
            pull(Class) |>
            str_replace(paste(c("ane", "Me"
                                , "Dime", "Trime")
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
  #  group.comps.info['Class'][group.comps.info['Class'] == "Dime"] <-
  #    "Dimethyl"
  #  group.comps.info['Class'][group.comps.info['Class'] == "Trime"] <-
  #    "Trimethyl"
  #  group.comps.info['Class'][group.comps.info['Class'] == "Tetrame"] <-
  #    "Tetramethyl"

  # cat('\n')
  # print("The class of the compounds was formatted")
  # cat('\n')

  group.comps.info |>
    as_tibble()
}
