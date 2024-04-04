load_all()

group_tables_list <- samples_plus_ri_list |>
  lapply(shape_group_table2)

master_table <- build_master_table2(group_tables_list)


grouping_info <- grouping_info |>
  unite(group_label
        , where(is.factor)
        , sep = "_"
        , remove = FALSE)

group.tables.list <- retrieve_group_tables2(group.label = "group_label"
                      , master.table = master_table
                      , grouping.info = grouping_info)

duplicated_compounds_presence <-
  retrieve_group_tables2(group.label = "group_label"
                        , master.table = master_table
                        , grouping.info = grouping_info) |>
  assess_duplicated_compounds2()


