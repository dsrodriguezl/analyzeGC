load_all()

df <- data.frame(
  id    = 1:5,
  group = c("A", "B", "C", "D", "E"),
  score = c(10, 20, 30, 40, 50),
  stringsAsFactors = FALSE
  )

rownames(df) <- letters[1:5]

df

df |>
  (function(df = df, origin = c("a", "group"), target = c("b", "group")){
  df[target[1], target[2]] <- df[origin[1], origin[2]]
  df[origin[1], origin[2]] <- NA
  df
})()


aligned_samples_data_list$`Winter_In-hive workers_A. m. mellifera`$aligned

