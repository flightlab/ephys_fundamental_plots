############################### package loading ################################
## Specify the packages you'll use in the script
packages <- c("tidyverse",
              "readxl",
              "zoo",
              "gridExtra",
              "R.matlab",
              "cowplot",
              "easystats",
              "circular",
              "splines",
              "MESS", ## area under curve
              "zoo" ## rolling means
)
## Now for each package listed, first check to see if the package is already
## installed. If it is installed, it's simply loaded. If not, it's downloaded
## from CRAN and then installed and loaded.
package.check <- lapply(packages,
                        FUN = function(x) {
                          if (!require(x, character.only = TRUE)) {
                            install.packages(x, dependencies = TRUE)
                            library(x, character.only = TRUE)
                          }
                        }
)

#### %not_in% ####
`%not_in%` <- Negate(`%in%`)

############################### metadata import ################################
## list all files
csv_files <-
  list.files("./data", pattern = ".csv",
             full.names = TRUE)
mat_files <-
  list.files("./data", pattern = ".mat",
             full.names = TRUE)
csv_file_info <-
  tibble(
    csv_files = csv_files,
    basename =  basename(csv_files) %>% str_remove(".csv"),
    basedate =  basename(csv_files) %>% str_sub(start = 1, end = 12)
  )
mat_file_info <-
  tibble(
    mat_files = mat_files,
    basename =  basename(mat_files) %>% str_remove(".mat"),
    basedate =  basename(mat_files) %>% str_sub(start = 1, end = 12)
  )

## now find which files have both CSVs and MATs
csv_mat_filejoin <-
  inner_join(csv_file_info, mat_file_info, by = "basename")

base_names <- csv_mat_filejoin$basename


########################## data import and processing ##########################
## now work on importing
gc()
metadata_sets <- NULL
meta_splits <- NULL
data_splits <- NULL

starttime <- Sys.time()
for (i in 1:nrow(csv_mat_filejoin)) {

  print(i)

  csv_data_sets <- NULL
  mat_data_sets <- NULL
  joined_data_sets <- NULL

  mat_import <-
    R.matlab::readMat(csv_mat_filejoin[i,"mat_files"])

  csv_data_sets[[i]] <-
    read_csv(csv_mat_filejoin[i,"csv_files"],
             show_col_types = FALSE) %>%
    rename(
      Spatial_Frequency = `Spatial Frequency`,
      Temporal_Frequency = `Temporal Frequency`
    )
  csv_data_sets[[i]]$SF_cpd[csv_data_sets[[i]]$Spatial_Frequency == 0.000668] <-
    2^-6
  csv_data_sets[[i]]$SF_cpd[csv_data_sets[[i]]$Spatial_Frequency == 0.001336] <-
    2^-5
  csv_data_sets[[i]]$SF_cpd[csv_data_sets[[i]]$Spatial_Frequency == 0.00267]  <-
    2^-4
  csv_data_sets[[i]]$SF_cpd[csv_data_sets[[i]]$Spatial_Frequency == 0.0053]   <-
    2^-3
  csv_data_sets[[i]]$SF_cpd[csv_data_sets[[i]]$Spatial_Frequency == 0.0106]   <-
    2^-2
  csv_data_sets[[i]]$SF_cpd[csv_data_sets[[i]]$Spatial_Frequency == 0.0212]   <-
    2^-1

  initial <- tibble(
    Trial = "initialization",
    Spatial_Frequency = csv_data_sets[[i]]$Spatial_Frequency[1],
    SF_cpd = csv_data_sets[[i]]$SF_cpd[1],
    Temporal_Frequency = csv_data_sets[[i]]$Temporal_Frequency[1],
    Direction  = csv_data_sets[[i]]$Direction[1],
    Time = 0.000
  )

  first_moving_mat <-
    mat_import[[stringr::str_which(names(mat_import), "Ch3")[1]]][[5]][,1][1]
  first_moving_csv <-
    csv_data_sets[[i]] %>%
    filter(Trial == "moving") %>%
    select(Time) %>%
    slice(1) %>%
    as.numeric()
  first_blank <-
    csv_data_sets[[i]] %>%
    filter(Trial == "blank") %>%
    select(Time) %>%
    slice(1) %>%
    as.numeric()
  first_mvbl_diff <- first_moving_csv - first_blank

  first_csv_tmp <-
    bind_rows(initial, csv_data_sets[[i]]) %>%
    ## Add the first event time to "Time" and subtract first_mvbl_diff (~2 secs)
    mutate(Time = Time + first_moving_mat - first_mvbl_diff - first_blank) %>%
    ## Make character version of Time for joining later
    mutate(Time_char = as.character(round(Time,3)))

  ## duplicate the initialization for ease of setting T0
  inception <-
    initial %>%
    mutate(Time_char = as.character(round(Time,3)))
  inception$Trial[1] <- "inception"

  first_csv <-
    bind_rows(inception, first_csv_tmp)
  first_csv$Stim_end <- c(first_csv$Time[-1], max(first_csv$Time) + 3)

  ## get final time
  final_time <- first_csv$Stim_end[nrow(first_csv)]

  # ## enforce timing offsets (if needed)
  # first_csv <-
  #   first_csv %>%
  #   mutate(
  #     Time = Time + 0.09,
  #     Time_char = as.character(Time)
  #   )
  # first_csv$Time[1] <- 0
  # first_csv$Time_char[1] <- "0"

  # if (i == 28) {
  #   ### FOR FILE 28 in newest batch ONLY ###
  #   pho <- read_csv("./NEWDATA2/file28_photod.csv", col_names = FALSE)
  #   spi <- read_csv("./NEWDATA2/file28_spikestransposed.csv", col_names = FALSE)
  #   pho_vec <- c(0,as_vector(pho[,1]))
  #   spi_vec <- c(0,as_vector(spi[,1]))
  #
  #   mat_data_sets[[i]] <-
  #     data.frame(
  #       Time =
  #         seq(
  #           from = 0, by = 0.001,
  #           length.out = 0 + length(spi_vec)
  #         ), #t(mat_import$tt),
  #       Spikes = spi_vec,
  #       Photod = pho_vec
  #     ) %>%
  #     as_tibble() %>%
  #     mutate(Time_char = as.character(Time))
  # } else {
  mat_data_sets[[i]] <-
    data.frame(
      Time =
        seq(
          from = 0, by = 0.001,
          length.out = 1 + length(t(mat_import$spikes))
        ), #t(mat_import$tt),
      Spikes = c(0, t(mat_import$spikes)),
      Photod = c(0, mat_import$photodiode)
    ) %>%
    as_tibble() %>%
    mutate(Time_char = as.character(Time)) %>%
    filter(Time <= final_time)
  #}


  joined_one_full <-
    mat_data_sets[[i]] %>%
    full_join(first_csv, by = "Time_char") %>%
    rename(Time_mat = Time.x,
           Time_csv = Time.y) %>%
    mutate(Time = as.numeric(Time_char)) %>%
    ## carry metadata forward
    mutate(
      Trial = zoo::na.locf(Trial, type = "locf"),
      Spatial_Frequency = zoo::na.locf(Spatial_Frequency, type = "locf"),
      SF_cpd = zoo::na.locf(SF_cpd, type = "locf"),
      Temporal_Frequency = zoo::na.locf(Temporal_Frequency, type = "locf"),
      Direction = zoo::na.locf(Direction, type = "locf"),
      Time_csv = zoo::na.locf(Time_csv, type = "locf"),
      Stim_end = zoo::na.locf(Stim_end, type = "locf")
    ) %>%
    ## Calculate velocity
    mutate(
      Speed = Temporal_Frequency/SF_cpd,
      Log2_Speed = log2(Speed)
    ) #%>%
  ## enforce timing offsets (if needed)
  # mutate(
  #   Time_mat = Time_mat + 0.09,
  #   Time_char = as.character(Time_mat),
  #   Time = Time + 0.09
  # )

  metadata_one_full <-
    first_csv %>%
    mutate(
      Speed = Temporal_Frequency/SF_cpd,
      Log2_Speed = log2(Speed),
      Stim_end_diff = c(0, diff(Stim_end))
    ) #%>%
  # ## enforce timing offsets (if needed)
  # mutate(
  #   Time = Time + 0.09,
  #   Time_char = as.character(Time)
  # )

  ## quality control
  stimtime_diffs <- round(metadata_one_full$Stim_end_diff)[-c(1:2)]
  stimtime_reps <- length(stimtime_diffs)/3
  stimtime_expectation <- rep(c(1,1,3), stimtime_reps)
  if (!all(stimtime_diffs == stimtime_expectation)) {
    print("stimtime issue; investigate")
  }

  mark_for_removal <-
    which(round(metadata_one_full$Stim_end_diff) %not_in% c(1, 3))
  if (any(mark_for_removal == 1 | mark_for_removal == 2)) {
    mark_for_removal <- mark_for_removal[mark_for_removal > 2]
  }
  if (length(mark_for_removal) > 0) {
    metadata_sets[[i]] <-
      metadata_one_full %>%
      filter(Time < metadata_one_full[mark_for_removal[1],]$Time)
    joined_data_sets[[i]] <-
      joined_one_full %>%
      filter(Time_mat < metadata_one_full[mark_for_removal[1],]$Time)
  } else {
    metadata_sets[[i]] <-
      metadata_one_full
    joined_data_sets[[i]] <-
      joined_one_full
  }

  meta_splits[[i]] <-
    metadata_sets[[i]] %>%
    filter(!Trial == "inception") %>%
    filter(!Trial == "initialization") %>%
    group_by(Spatial_Frequency, Temporal_Frequency, Direction) %>%
    group_split()

  data_splits[[i]] <-
    joined_data_sets[[i]] %>%
    filter(!Trial == "inception") %>%
    filter(!Trial == "initialization") %>%
    group_by(Spatial_Frequency, Temporal_Frequency, Direction) %>%
    group_split()

  rm(
    first_csv, inception, initial, mat_import, first_csv_tmp, metadata_one_full,
    joined_one_full, joined_data_sets, csv_data_sets, mat_data_sets
  )
  message("File ", i, ": ", csv_mat_filejoin[i,"basename"], " imported")
  gc()
}

endtime <- Sys.time()
endtime - starttime

## tidy up
gc()

names(metadata_sets) <- csv_mat_filejoin$basename#base_names
names(meta_splits)   <- csv_mat_filejoin$basename#base_names
names(data_splits)   <- csv_mat_filejoin$basename#base_names

################################# preprocessing ################################

metadata_combos <- NULL
for (i in 1:length(metadata_sets)) {
  metadata_combos[[i]] <-
    metadata_sets[[i]] %>%
    distinct(SF_cpd, Temporal_Frequency, Speed, Direction) %>%
    arrange(Direction) %>%
    arrange(desc(SF_cpd)) %>%
    mutate(
      plot_order = 1:length(meta_splits[[i]]),
      name = paste(Direction, "Deg,",  Speed, "Deg/s")
    )
}
names(metadata_combos) <- csv_mat_filejoin$basename #base_names

#### __SET BIN SIZE HERE ####
bin_size = 10 ## 10 or 100

slice_size = NULL
slicemin = NULL
slicemax = NULL
if (bin_size == 10){
  slice_size <- 501
  slicemin <- 202
  slicemax <- 498
} else if (bin_size == 100){
  slice_size <- 51
  slicemin <- 21
  slicemax <- 49
} else {
  print("bin_size is non-standard")
}


#### __reorganize in stim order ####
all_replicate_data_reorganized <-
  vector(mode = "list", length = length(meta_splits))
name_sets <-
  vector(mode = "list", length = length(meta_splits))
gc()

starttime <- Sys.time()
for (i in 1:length(meta_splits)){ # i = number of files
  print(i)
  replicate_data_reorganized <- NULL
  for (j in 1:length(meta_splits[[i]])) { # j = {direction,speed}
    d <- data_splits[[i]][[j]]
    m <- meta_splits[[i]][[j]] %>%
      group_by(Trial) %>%
      mutate(Replicate = row_number())

    name_sets[[i]][[j]] <-
      paste(m$Direction[1], "Deg,",  m$Speed[1], "Deg/s")

    replicates_ordered <- NULL
    for (k in 1:max(m$Replicate)){
      tmp <-
        m %>%
        filter(Replicate == k)

      if (nrow(tmp) == 3 ) {
        doot <-
          d %>%
          filter(Time >= min(tmp$Time)) %>%
          filter(Time <= max(tmp$Stim_end))
        doot$bin <-
          rep(1:ceiling(nrow(doot)/bin_size), each = bin_size)[1:nrow(doot)]

        replicates_ordered[[k]] <-
          doot %>%

          ## IF YOU ARE BINNING, RUN THIS:
          group_by(bin) %>%
          summarise(
            Trial = first(Trial),
            Time_bin_mid = mean(Time_mat),
            Time_bin_begin = min(Time_mat),
            Time_bin_end = max(Time_mat),
            Spike_rate = sum(Spikes)/(max(Time_mat) - min(Time_mat)),
            Photod_mean = mean(Photod)
          ) %>%
          mutate(
            Time_stand = Time_bin_mid - min(Time_bin_mid),
            Blank_end = tmp$Stim_end[1] - min(Time_bin_mid),
            Static_end = tmp$Stim_end[2] - min(Time_bin_mid),
            SF_cpd = m$SF_cpd[1],
            Temporal_Frequency = m$Temporal_Frequency[1],
            Speed = m$Speed[1],
            Direction = m$Direction[1],
            Replicate = k
          ) %>%
          select(Speed, SF_cpd, Temporal_Frequency, Direction,
                 everything()) %>%
          filter(Time_stand >= 0) %>%
          filter(bin < slice_size + 1)

        ## IF YOU ARE NOT BINNING, RUN THIS:
        # mutate(
        #   Time_stand = Time_mat - min(Time_mat),
        #   Time_begin = min(Time_mat),
        #   Time_end = max(Time_mat),
        #   Blank_end = tmp$Stim_end[1] - min(Time_mat),
        #   Static_end = tmp$Stim_end[2] - min(Time_mat),
        #   Replicate = k
        # ) %>%
        #   select(Speed, SF_cpd, Temporal_Frequency, Direction,
        #          everything()) %>%
        #   filter(Time_stand >= 0)

      }
    }

    replicate_data_reorganized[[j]] <-
      replicates_ordered %>%
      bind_rows()

    all_replicate_data_reorganized[[i]][[j]] <-
      replicate_data_reorganized[[j]]

    rm(replicates_ordered)

    gc()
  }
}
endtime <- Sys.time()
endtime - starttime
gc()

for (i in 1:length(all_replicate_data_reorganized)) {
  for (j in 1:length(all_replicate_data_reorganized[[i]])) {
    names(all_replicate_data_reorganized[[i]])[[j]] <- name_sets[[i]][[j]]
  }
}
names(all_replicate_data_reorganized) <- csv_mat_filejoin$basename #base_names

## check to see if all averaged replicate sets are in the same order with
## respect to direction and speed

names_grid <- matrix(
  ncol = length(all_replicate_data_reorganized),
  nrow = length(all_replicate_data_reorganized)
)
for (i in 1:length(all_replicate_data_reorganized)) {
  for (j in 1:length(all_replicate_data_reorganized)) {
    names_grid[i,j] <-
      identical(names(all_replicate_data_reorganized[[i]]),
                names(all_replicate_data_reorganized[[j]])
      )
  }
}

################################### export csv #################################
## Export a csv for each session with all data organized

export_path <- "./reformatted_data_csv/"
#condition <-  "_unbinned"
condition <-  "_binsize10"
#condition <-  "_binsize100"

for (i in 1:length(all_replicate_data_reorganized)) {
  print(i)
  dat <-
    all_replicate_data_reorganized[[i]] %>%
    bind_rows()
  write_csv(
    dat,
    file =
      paste0(
        export_path,
        names(all_replicate_data_reorganized)[i],
        condition,
        ".csv"
      )
  )
  rm(dat)
}

gc()

################################## quality check ###############################
## To guard against spurious matlab pulses, ensure that each replicate consists
## of:
## 1) 1-sec blank
## 2) 1-sec stationary
## 3) 3-sec moving

converted_csv_mats <-
  csv_mat_filejoin %>%
  filter(basename %in% unbinned_basenames)

csv_assessment <- NULL
for (i in 1:nrow(converted_csv_mats)) {

  print (i)

  csv_read <- NULL
  # mat_data_sets <- NULL
  # joined_data_sets <- NULL

  ## import the csv
  csv_read <-
    read_csv(converted_csv_mats[i,"csv_files"],
             show_col_types = FALSE) %>%
    rename(
      Spatial_Frequency = `Spatial Frequency`,
      Temporal_Frequency = `Temporal Frequency`
    )
  csv_read$SF_cpd[csv_read$Spatial_Frequency == 0.000668] <-
    2^-6
  csv_read$SF_cpd[csv_read$Spatial_Frequency == 0.001336] <-
    2^-5
  csv_read$SF_cpd[csv_read$Spatial_Frequency == 0.00267]  <-
    2^-4
  csv_read$SF_cpd[csv_read$Spatial_Frequency == 0.0053]   <-
    2^-3
  csv_read$SF_cpd[csv_read$Spatial_Frequency == 0.0106]   <-
    2^-2
  csv_read$SF_cpd[csv_read$Spatial_Frequency == 0.0212]   <-
    2^-1

  ## get diffs in time
  time_diffs <-
    diff(csv_read$Time) %>%
    ## round to nearest second
    round()

  if (any(time_diffs %not_in% c(1, 3))) {
    print("fuckery")
    csv_assessment[i] <- "investigate"
  } else(
    csv_assessment[i] <- "all clear"
  )

  rm(csv_read, time_diffs)

}

################################## export matching #############################

unbinned_original_filelist <-
  list.files("./data_csv_exports/", pattern = "_unbinned.csv",
             full.names = TRUE)
unbinned_original_basenames <-
  unbinned_original_filelist %>%
  str_remove("./data_csv_exports/") %>%
  str_remove("_unbinned.csv")

unbinned_newer_filelist <-
  list.files("./reformatted_data_csv/", pattern = "_unbinned.csv",
             full.names = TRUE)
unbinned_newer_basenames <-
  unbinned_newer_filelist %>%
  str_remove("./reformatted_data_csv/") %>%
  str_remove("_unbinned.csv")

## first check
identical(unbinned_original_basenames, unbinned_newer_basenames)

file_checks <- NULL
## second check
for (i in 1:length(unbinned_original_filelist)) {

  print(i)

  original <-
    read_csv(unbinned_original_filelist[i],
             show_col_types = FALSE) %>%
    select(-any_of(c("bin")))

  newer <-
    read_csv(unbinned_newer_filelist[i],
             show_col_types = FALSE) %>%
    select(-any_of(c("bin")))

  file_checks[[i]] <- all.equal(original, newer)

  if(!identical(original,newer)) {
    paste0("file #", i, ", ", unbinned_original_basenames[i], " not identical")
  }

  rm(original, newer)

}
