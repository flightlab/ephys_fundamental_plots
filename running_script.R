#### 2.1 R packages & versioning ####

## Specify the packages you'll use in the script
packages <- c("tidyverse",
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


#### 2.2 %not_in% ####
`%not_in%` <- Negate(`%in%`)


#### 5 Data wrangling ####

#### 5.1.1 Identify files to import ####

## List all files of each file type
csv_files <-
  list.files("./data", pattern = ".csv",
             full.names = TRUE)
mat_files <-
  list.files("./data", pattern = ".mat",
             full.names = TRUE)

## Generate metadata tibbles for each file type
csv_file_info <-
  tibble(
    csv_files = csv_files,
    ## Extract the basename by removing the file extension
    basename =  basename(csv_files) %>% str_remove(".csv"),
    ## NOTE: PLEASE ADJUST THE FOLLOWING LINE TO BE ABLE TO EXCTRACT OUT THE
    ## DATE BASED ON YOUR NAMING CONVENTION
    basedate =  basename(csv_files) %>% str_sub(start = 1, end = 12)
  )
mat_file_info <-
  tibble(
    mat_files = mat_files,
    ## Extract the basename by removing the file extension
    basename =  basename(mat_files) %>% str_remove(".mat"),
    ## NOTE: AGAIN, PLEASE ADJUST THE FOLLOWING LINE TO BE ABLE TO EXCTRACT OUT
    ## THE DATE BASED ON YOUR NAMING CONVENTION
    basedate =  basename(mat_files) %>% str_sub(start = 1, end = 12)
  )

## Matchmake between .MAT data and .CSV log files
csv_mat_filejoin <-
  inner_join(csv_file_info, mat_file_info, by = "basename") %>%
  ## OPTIONAL STEP: remove any rows where either the .MAT or .CSV is missing
  drop_na()

## Store a vector of basenames in the environment. This will become useful later
base_names <- csv_mat_filejoin$basename


#### 5.1.2 Data import and preliminary labeling ####

## Set up empty vectors that will collect sets of replicates that we will be
## splitting up
metadata_sets <- NULL
meta_splits <- NULL
data_splits <- NULL
gc()

starttime <- Sys.time() ## Optional, will help you assess run time
for (i in 1:nrow(csv_mat_filejoin)) {

  ## Which file # are we working on?
  print(i)

  ## Set up temporary objects in which to eventually write data
  csv_data_sets <- NULL
  mat_data_sets <- NULL
  joined_data_sets <- NULL

  ## Import the matlab file. This may take some time
  mat_import <-
    R.matlab::readMat(csv_mat_filejoin[i,"mat_files"])

  ## Read in the corresponding csv log file
  csv_data_sets[[i]] <-
    read_csv(as.character(csv_mat_filejoin[i,"csv_files"]),
             show_col_types = FALSE) %>%
    ## Rename columns for convenience
    rename(
      Spatial_Frequency = `Spatial Frequency`,
      Temporal_Frequency = `Temporal Frequency`,
      Cycles_Per_Pixel = `Cycles Per Pixel`
    )

  ## The log file does not have time = 0, so set up a separate tibble to
  ## add this info in later. Some of the metadata will just be filler for now.
  initial <- tibble(
    Trial = "initialization",
    Spatial_Frequency =  csv_data_sets[[i]]$Spatial_Frequency[1],
    Cycles_Per_Pixel  =  csv_data_sets[[i]]$Cycles_Per_Pixel[1],
    Temporal_Frequency = csv_data_sets[[i]]$Temporal_Frequency[1],
    Direction  = csv_data_sets[[i]]$Direction[1],
    Time = 0.000
  )

  ## Find photodiode
  ## It is almost always in channel 2, but we should be sure to check before
  ## extracting automatically
  photod_default_channel <-
    mat_import[[stringr::str_which(names(mat_import), "Ch2")[1]]]
  if (!photod_default_channel[1][[1]][1] == "waveform") {
    warning("File ", i,": Photodiode channel identity uncertain")
  }

  ## Find spikes
  ## Similarly, spikes are almost always in channel 5, but we should check
  spikes_default_channel <-
    mat_import[[stringr::str_which(names(mat_import), "Ch5")[1]]]
  if("codes" %not_in% attributes(spikes_default_channel)$dimnames[[1]]) {
    warning("File ", i,": Sorted spikes channel identity uncertain")
  }
  ## If that worked, see if we can automatically determine the "times" and
  ## "codes" slot numbers
  times_location <-
    which(attributes(spikes_default_channel)$dimnames[[1]] == "times")
  codes_location <-
    which(attributes(spikes_default_channel)$dimnames[[1]] == "codes")

  ## Find matlab's stimulus change log
  stim_change_channel <-
    mat_import[[stringr::str_which(names(mat_import), "Ch3")[1]]]
  ## Each sweep should be 5 secs. We'll check that the median is 5
  ## If this results in an error, then the channel identity could be wrong, or
  ## there may have been an issue with sweep duration during the recording
  ## process
  if(!median(round(diff(stim_change_channel[[5]][,1]),0)) == 5) {
    warning("File ", i,": stim change channel identity uncertain")
  }

  ## Determine when the onset of motion occurred according to matlab
  first_moving_mat <-
    stim_change_channel[[5]][,1][1]
  ## Find the first "moving" phase in the log file
  first_moving_csv <-
    csv_data_sets[[i]] %>%
    filter(Trial == "moving") %>%
    select(Time) %>%
    slice(1) %>%
    as.numeric()
  ## Find the first "blank" phase in the log file
  first_blank <-
    csv_data_sets[[i]] %>%
    filter(Trial == "blank") %>%
    select(Time) %>%
    slice(1) %>%
    as.numeric()
  ## Compute the difference between these two
  first_mvbl_diff <- first_moving_csv - first_blank

  ## Check to see if the final row of the metadata is "moving" and truncate
  ## if not
  ## This can effectively be done by truncating after the final "moving" phase
  max_moving_sweep <-
    max(which(csv_data_sets[[i]]$Trial == "moving"))

  first_csv_tmp <-
    bind_rows(initial, csv_data_sets[[i]]) %>%
    ## Add 1 to max moving sweep since we tacked on "initial" in previous step
    ## Then slice to restrict any extraneous partial sweeps
    slice_head(n = (max_moving_sweep + 1)) %>%
    ## Add the first event time to "Time" and subtract first_mvbl_diff (~2 secs)
    ## What this does is shift the log csv's time stamping to match the matlab
    ## file's stim change channel's time stamping
    mutate(Time = Time + first_moving_mat - first_mvbl_diff - first_blank) %>%
    ## Make character version of Time for joining later
    ## This will be crucial for _join functions
    mutate(Time_char = as.character(round(Time,3)))

  ## Duplicate the initialization for ease of setting T0
  inception <-
    initial %>%
    mutate(Time_char = as.character(round(Time,3)))
  inception$Trial[1] <- "inception"

  ## Bind the initialization rows
  first_csv <-
    bind_rows(inception, first_csv_tmp)
  ## Compute stimulus end times
  first_csv$Stim_end <- c(first_csv$Time[-1], max(first_csv$Time) + 3)

  ## Get final time
  final_time <- first_csv$Stim_end[nrow(first_csv)]

  ## Extract photodiode data
  ## First generate a time sequence to match to the photodiode trace
  Time_vec <- seq(
    from = 0.0000,
    by = 1 / 25000,
    length.out = length(photod_default_channel[9][[1]][, 1])
  )
  ## The key thing is to get a character version of time from this
  Time_char_vec <- as.character(round(Time_vec, 3))

  ## Grab the photodiode data
  photod_full <-
    tibble(Photod =
             photod_default_channel[9][[1]][, 1])
  ## Add numeric time
  photod_full$Time <-
    seq(
      from = 0.0000,
      by = 1 / 25000,
      length.out = nrow(photod_full)
    )
  options(scipen = 999)
  photod_full <-
    photod_full %>%
    ## Add the character time
    add_column(Time_char = Time_char_vec) %>%
    ## Use the charcter time to define a group
    group_by(Time_char) %>%
    ## Then average the photodiode within
    summarise(Photod = mean(Photod)) %>%
    ungroup() %>%
    mutate(Time = round(as.numeric(Time_char), 3)) %>%
    arrange(Time) %>%
    filter(Time <= final_time)

  ## Extract all spike data
  all_spike_dat <-
    tibble(
      Time =
        spikes_default_channel[times_location][[1]][, 1],
      code =
        spikes_default_channel[codes_location][[1]][, 1]) %>%
    ## Characterize time, for purposes of joining later
    mutate(Time_char = as.character(round(Time, 3)))

  ## How many distinct neurons are there?
  n_cells <- sort(unique(all_spike_dat$code))

  if(length(n_cells) > 1) { ## if there's more than one distinct neuron
    all_spike_dat_tmp <-
      all_spike_dat %>%
      ## Group by identity of spiking neuron
      group_by(code) %>%
      ## Split into separate dfs, one per neuron
      group_split()

    ## Additional cells are labeled as "Spike_n"
    all_cells <- NULL
    for (j in n_cells) {
      #print(j)
      new_name = paste0("Spikes_", j)
      all_cells[[j+1]] <-
        all_spike_dat_tmp[[j+1]]

      ## Consolidate to 3 decimal places
      all_cells[[j+1]] <-
        all_cells[[j+1]] %>%
        group_by(Time_char) %>%
        summarise(code = mean(code)) %>%
        mutate(code = ceiling(code)) %>%
        ungroup() %>%
        mutate(Time = round(as.numeric(Time_char), 3)) %>%
        arrange(Time) %>%
        filter(Time <= final_time)

      names(all_cells[[j+1]])[match("code", names(all_cells[[j+1]]))] <-
        new_name
      ## Replace "j" with 1 to indicate presence/absence of spike rather than
      ## cell identity
      all_cells[[j+1]][new_name] <- 1

      ## If the identity is 1, replace "Spikes_1" with just "Spikes"
      if (new_name == "Spikes_1") {
        names(all_cells[[j+1]])[match(new_name, names(all_cells[[j+1]]))] <-
          "Spikes"
      }
    }

    ## Consolidate
    all_spike_dat <-
      all_cells %>%
      ## Tack on additional spike columns
      reduce(full_join, by = "Time_char") %>%
      arrange(Time_char) %>%
      ## Remove time.n columns but
      ## Do not remove Time_char
      select(-contains("Time.")) %>%
      ## Regenerate numeric time
      mutate(
        Time = as.numeric(Time_char)
      ) %>%
      select(Time, Time_char, Spikes, everything()) %>%
      filter(Time <= final_time)

  } else { ## If there's only 1 neuron
    all_spike_dat <-
      all_spike_dat %>%
      group_by(Time_char) %>%
      summarise(code = mean(code)) %>%
      mutate(code = ceiling(code)) %>%
      ungroup() %>%
      rename(Spikes = code) %>%
      mutate(Time = round(as.numeric(Time_char), 3)) %>%
      arrange(Time) %>%
      filter(Time <= final_time) %>%
      select(Time, Time_char, everything())
  }

  options(scipen = 999)
  mat_data_sets[[i]] <-
    ## Generate a time sequence from 0 to final_time
    tibble(
      Time = seq(from = 0, to = final_time, by = 0.001)
    ) %>%
    ## Character-ize it
    mutate(Time_char = as.character(round(Time, 5))) %>%
    ## Join in the photodiode data
    left_join(photod_full, by = "Time_char") %>%
    select(-Time.y) %>%
    rename(Time = Time.x) %>%
    arrange(Time) %>%
    ## Join in the spike data
    left_join(all_spike_dat, by = "Time_char") %>%
    select(-Time.y) %>%
    rename(Time = Time.x) %>%
    arrange(Time) %>%
    filter(Time <= final_time) %>%
    ## Replace NAs with 0 in the Spike columns only
    mutate(
      across(starts_with("Spikes"), ~replace_na(.x, 0))
    )

  ## Merge the matlab data with the metadata
  joined_one_full <-
    mat_data_sets[[i]] %>%
    ## Join by the character version of time NOT the numerical!!
    full_join(first_csv, by = "Time_char") %>%
    ## Rename columns for clarity of reference
    rename(Time_mat = Time.x,
           Time_csv = Time.y) %>%
    ## Convert character time to numeric time
    mutate(Time = round(as.numeric(Time_char), 3)) %>%
    ## Carry metadata forward
    mutate(
      Trial = zoo::na.locf(Trial, type = "locf"),
      Spatial_Frequency = zoo::na.locf(Spatial_Frequency, type = "locf"),
      Cycles_Per_Pixel  = zoo::na.locf(Cycles_Per_Pixel, type = "locf"),
      Temporal_Frequency = zoo::na.locf(Temporal_Frequency, type = "locf"),
      Direction = zoo::na.locf(Direction, type = "locf"),
      Time_csv = zoo::na.locf(Time_csv, type = "locf"),
      Stim_end = zoo::na.locf(Stim_end, type = "locf")
    ) %>%
    ## Calculate velocity
    mutate(
      Speed = round(Temporal_Frequency/Spatial_Frequency, 0),
      Log2_Speed = log2(Speed)
    )

  ## Add info to metadata
  metadata_one_full <-
    first_csv %>%
    mutate(
      Speed = round(Temporal_Frequency/Spatial_Frequency, 0),
      Log2_Speed = log2(Speed),
      Stim_end_diff = c(0, diff(Stim_end))
    )

  ## Some quality control checks
  ## What are the stim time differences?
  stimtime_diffs <- round(metadata_one_full$Stim_end_diff)[-c(1:2)]
  ## How many total reps were recorded?
  stimtime_reps <- length(stimtime_diffs)/3
  ## What do we expect the overall structure to look like?
  stimtime_expectation <- rep(c(1,1,3), stimtime_reps)
  ## Does reality match our expectations?
  if (!all(stimtime_diffs == stimtime_expectation)) {
    ## If you get this, investigate the file further and determine what went
    ## wrong
    print("stimtime issue; investigate")
  }

  ## Sometimes the final sweep gets carried for an indefinite amount of time
  ## before the investigator manually shuts things off. The following block
  ## truncates accordingly
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

  ## Organize the metadata for export in the R environment
  meta_splits[[i]] <-
    metadata_sets[[i]] %>%
    ## Get rid of the non-sweep info
    filter(!Trial == "inception") %>%
    filter(!Trial == "initialization") %>%
    ## Group by trial
    #group_by(Spatial_Frequency, Temporal_Frequency, Direction) %>%
    ## Split by trial group
    group_split(Spatial_Frequency, Temporal_Frequency, Direction)

  data_splits[[i]] <-
    joined_data_sets[[i]] %>%
    ## Get rid of the non-sweep info
    filter(!Trial == "inception") %>%
    filter(!Trial == "initialization") %>%
    ## Group by trial
    group_by(Spatial_Frequency, Temporal_Frequency, Direction) %>%
    ## Split by trial group
    group_split()

  ## Do some cleanup so large objects don't linger in memory
  rm(
    first_csv, inception, initial, mat_import, first_csv_tmp,
    photod_default_channel, spikes_default_channel, photod_full,
    all_spike_dat, all_spike_dat_tmp, all_cells,
    metadata_one_full, joined_one_full, joined_data_sets,
    csv_data_sets, mat_data_sets
  )
  message("File ", i, ": ", csv_mat_filejoin[i,"basename"], " imported")
  gc()
}

endtime <- Sys.time()
endtime - starttime ## Total elapsed time

## Tidy up how R has been using RAM by running garbage collection
gc()

## Name each data set according to the basename of the file
names(metadata_sets) <- csv_mat_filejoin$basename #base_names
names(meta_splits)   <- csv_mat_filejoin$basename #base_names
names(data_splits)   <- csv_mat_filejoin$basename #base_names

## Old method: Extract matlab spike and photodiode data
# #mat_data_sets[[i]] <-
# spikedat <-
#   data.frame(
#     Time =
#       seq(
#         from = 0, by = 0.001,
#         length.out = 1 + length(t(mat_import$spikes))
#       ), #t(mat_import$tt),
#     Spikes = c(0, t(mat_import$spikes)),
#     Photod = c(0, mat_import$photodiode)
#   ) %>%
#   as_tibble() %>%
#   mutate(Time_char = as.character(Time)) %>%
#   filter(Time <= final_time)


## Get organized lists of stimuli that were used
## This will ultimately be used for arranging data by stimuli in a sensible
## order
metadata_combos <- NULL
for (i in 1:length(metadata_sets)) {
  metadata_combos[[i]] <-
    metadata_sets[[i]] %>%
    ## Get unique stimulus parameters
    distinct(Spatial_Frequency, Temporal_Frequency, Speed, Direction) %>%
    arrange(Direction) %>%
    ## Sort by SF (smallest to largest)
    arrange(desc(Spatial_Frequency)) %>%
    mutate(
      ## Set a "plot_order" which provides a running order of stimuli
      plot_order = 1:length(meta_splits[[i]]),
      name = paste(Direction, "Deg,",  Speed, "Deg/s")
    )
}
names(metadata_combos) <- csv_mat_filejoin$basename #base_names


#### 5.2 Organizing replicates (required) and binning (optional) ####
#### __SET BIN SIZE HERE ####

## Set bin size here
## Units are in ms (e.g. 10 = 10ms)
bin_size = 100 ## 10 or 100 or 1 (1 = "unbinned")

slice_size = NULL
slicemin = NULL
slicemax = NULL
condition = NULL
if (bin_size == 10){
  slice_size <- 501
  slicemin <- 202
  slicemax <- 498
  condition <- "_binsize10"
} else if (bin_size == 100){
  slice_size <- 51
  slicemin <- 21
  slicemax <- 49
  condition <- "_binsize100"
} else if (bin_size == 1){
  slice_size <- NULL
  slicemin <- NULL
  slicemax <- NULL
  condition <- "_unbinned"
} else {
  stop("bin_size is non-standard")
}


all_replicate_data_reorganized <-
  vector(mode = "list", length = length(meta_splits))
name_sets <-
  vector(mode = "list", length = length(meta_splits))
gc()

starttime <- Sys.time()
for (i in 1:length(meta_splits)){

  ## i = file number
  print(i)

  ## We'll need to collect data at a per-stimulus level and on a per-replicate
  ## level within the per-stimulus level
  ## "j" will be used to designate a unique stimulus
  ## We'll first create an empty object in which to collect stimulus-specific
  ## data
  replicate_data_reorganized <- NULL
  ## For each of j unique stimuli...
  for (j in 1:length(meta_splits[[i]])) { # j = {direction,speed}
    ## Isolate the j-th data
    d <- data_splits[[i]][[j]]
    ## And the j-th log data
    m <- meta_splits[[i]][[j]] %>%
      group_by(Trial) %>%
      ## Label separate replicates
      mutate(Replicate = row_number())

    ## Extract a stimulus label to a name_set that will be used later
    name_sets[[i]][[j]] <-
      paste(m$Direction[1], "Deg,",  m$Speed[1], "Deg/s")

    ## Set up a temporary object to deal with per-replicate data
    replicates_ordered <- NULL
    ## "k" will be used to designate replicate number
    for (k in 1:max(m$Replicate)){
      tmp <-
        m %>%
        filter(Replicate == k)

      ## If you have a complete replicate (i.e., blank, stationary, moving)
      if (nrow(tmp) == 3 ) {
        ## Grab the specific data
        doot <-
          d %>%
          filter(Time >= min(tmp$Time)) %>%
          filter(Time <= max(tmp$Stim_end))
        ## Add bin information
        doot$bin <-
          rep(1:ceiling(nrow(doot)/bin_size), each = bin_size)[1:nrow(doot)]

        if (bin_size == 1) {
          ## IF YOU ARE NOT BINNING, RUN THIS:
          replicates_ordered[[k]] <-
            doot %>%
            mutate(
              ## Construct a standardized time within the sweep
              Time_stand = Time_mat - min(Time_mat),
              ## When does the sweep begin
              Time_begin = min(Time_mat),
              ## When does the sweep end
              Time_end = max(Time_mat),
              ## Delineate the end of the blank phase
              Blank_end = tmp$Stim_end[1] - min(Time_mat),
              ## Delineate the end of the stationary phase
              Static_end = tmp$Stim_end[2] - min(Time_mat),
              ## Label the replicate number
              Replicate = k
            ) %>%
            ## Bring stim info to first few columns
            select(Speed, Spatial_Frequency, Temporal_Frequency, Direction,
                   everything()) %>%
            ## Just in case there some hang over
            filter(Time_stand >= 0)
        } else { ## IF YOU ARE BINNING, RUN THIS:

          ## First grab time and meta info
          time_and_meta <-
            doot %>%
            ## WITHIN EACH BIN:
            group_by(bin) %>%
            summarise(
              ## Label the trial
              Trial = first(Trial),
              ## Midpoint of bin
              Time_bin_mid = mean(Time_mat),
              ## Bin beginning
              Time_bin_begin = min(Time_mat),
              ## Bin end
              Time_bin_end = max(Time_mat),
              ## Spike rate = sum of spikes divided by elapsed time
              Spike_rate = sum(Spikes) / (max(Time_mat) - min(Time_mat)),
              Photod_mean = mean(Photod)
            )

          ## Now deal with Spike_N columns
          hold_spike_n <-
            doot %>%
            select(starts_with("Spikes_")) %>%
            add_column(bin = doot$bin) %>%
            add_column(Time_mat = doot$Time_mat) %>%
            group_by(bin) %>%
            summarise(across(starts_with("Spikes_"),
                             ~ sum(.x) / (max(Time_mat) - min(Time_mat))))

          ## Put them together
          replicates_ordered[[k]] <-
            time_and_meta %>%
            left_join(hold_spike_n, by = "bin") %>%
            mutate(
              ## Add in metadata (following same definitions above)
              Time_stand = Time_bin_mid - min(Time_bin_mid),
              Blank_end = tmp$Stim_end[1] - min(Time_bin_mid),
              Static_end = tmp$Stim_end[2] - min(Time_bin_mid),
              Spatial_Frequency = m$Spatial_Frequency[1],
              Temporal_Frequency = m$Temporal_Frequency[1],
              Speed = m$Speed[1],
              Direction = m$Direction[1],
              Replicate = k
            ) %>%
            ## Bring stim info to first few columns
            select(Speed, Spatial_Frequency, Temporal_Frequency, Direction,
                   everything()) %>%
            ## Just in case there some hang over
            filter(Time_stand >= 0) %>%
            filter(bin < slice_size + 1)

          rm(time_and_meta, hold_spike_n)
        }
      }
      }

    ## Now insert it within the collector of per-stimulus data
    replicate_data_reorganized[[j]] <-
      replicates_ordered %>%
      bind_rows()

    ## Now insert it within the overall data collector
    all_replicate_data_reorganized[[i]][[j]] <-
      replicate_data_reorganized[[j]]

    ## Toss out temporary objects and clean up
    rm(replicates_ordered, d, m, tmp)
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


#### 5.3 Data export ####

## Export a csv for each session with all data organized

## Declare export destination
export_path <- "./data/"
## The "condition" will be appended to the file name.

## Export each tibble within all_replicate_data_reorganized
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


#### 6 Raster and mean spike rate plots ####

#### 6.1 Data sets ####

## File paths and basenames of _unbinned.csv files
unbinned_filelist <-
  list.files("./data/", pattern = "_unbinned.csv",
             full.names = TRUE)
unbinned_basenames <-
  unbinned_filelist %>%
  str_remove("./data/") %>%
  str_remove("_unbinned.csv")

## File paths and basenames of _binsize10.csv files
bin10_filelist <-
  list.files("./data/", pattern = "_binsize10.csv",
             full.names = TRUE)
bin10_basenames <-
  bin10_filelist %>%
  str_remove("./data/") %>%
  str_remove("_binsize10.csv")

## File paths and basenames of _binsize100.csv files
bin100_filelist <-
  list.files("./data/", pattern = "_binsize100.csv",
             full.names = TRUE)
bin100_basenames <-
  bin100_filelist %>%
  str_remove("./data/") %>%
  str_remove("_binsize100.csv")



#### 6.2 Raster plot ####

## Read in the data
## Since there is only 1 example file, we'll simplify things by just
## subsetting unbinned_filelist
unbinned_data <-
  read_csv(unbinned_filelist[1], show_col_types = FALSE) %>%
  as_tibble()

## determine the max number of replicates
max_reps <- max(unbinned_data$Replicate)

## Generate the code for the ggplot and save it as "rasterplot"
rasterplot <-
  unbinned_data %>%
  ## Remove any rows where spiking does not occur in the Spikes column
  filter(Spikes == 1) %>%
  ## Convert Trial and Speed into factors and specify their level ordering
  ## This will make it easier to get the subplots in the order we want them
  mutate(
    Trial = factor(Trial,
                   levels = c("blank", "stationary", "moving")),
    Speed = factor(Speed,
                   levels = c(1024, 644, 407, 256, 128, 64, 32, 16, 8, 4))) %>%
  ggplot(aes(x = Time_stand, y = Replicate)) +
  ## The next three blocks will undershade each subplot according to stimulus
  ## phase (i.e., blank, stationary, moving)
  annotate("rect",
           xmin = 0, xmax = first(unbinned_data$Blank_end),
           ymin = 0.5, ymax = max_reps + 0.5,
           alpha = 0.075, color = NA, fill = "red") +
  annotate("rect",
           xmin = first(unbinned_data$Blank_end),
           xmax = first(unbinned_data$Static_end),
           ymin = 0.5, ymax = max_reps + 0.5,
           alpha = 0.075, color = NA, fill = "darkgoldenrod1") +
  annotate("rect",
           xmin = first(unbinned_data$Static_end), xmax = 5,
           ymin = 0.5, ymax = max_reps + 0.5,
           alpha = 0.075, color = NA, fill = "forestgreen") +
  ## Up to 10 replicates were used, so we will force the y-axis to go to 10
  scale_y_continuous(
    limits = c(0.5, max_reps + 0.5),
    expand = c(0, 0),
    breaks = c(max_reps/2, max_reps)
  ) +
  ## There are multiple ways to plot a spike event. Since 100% of the rows in
  ## this filtered data set are spike events, we can simply plot a symbol at
  ## each time (Time_stand) that appears in the data. The `|` symbol is a good
  ## choice.
  geom_point(pch = '|', size = 1.5) +
  xlab("Time (sec)") +
  ggtitle(paste0(unbinned_basenames[1], " raster")) +
  ## Use facet_grid() to create a grid of subplots. Rows will correspond to
  ## Speeds, and columns correspond to Directions
  facet_grid(rows = vars(Speed), cols = vars(Direction)) +
  theme_classic() +
  theme(legend.position = 'none',
        panel.spacing = unit(0.1, "lines"))


#### 6.3 Mean spike rate plots ####

## Bin size = 100

## Read in the data
## Since there is only 1 example file, we'll simplify things by just
## subsetting bin100_filelist
bin100_data <-
  read_csv(bin100_filelist[1], show_col_types = FALSE) %>%
  as_tibble()

## Compute SE and other metrics and add this to our data set
dataslices_100 <-
  bin100_data %>%
  ## Split by direction and speed, because we will use those to define each
  ## subplot
  group_split(Direction, Speed) %>%
  ## Group by time bin
  purrr::map(group_by, bin) %>%
  ## Within each time bin, compute the following:
  purrr::map(transmute,
             ## first() can be used for metadata such as Speed or Direction
             Speed = first(Speed),
             Direction = first(Direction),
             ## I generally compute the mean within each bin for the following:
             Time_stand = mean(Time_stand),
             Blank_end = mean(Blank_end),
             Static_end = mean(Static_end),
             Mean_spike_rate = mean(Spike_rate),
             ## To get SE, divide s.d. by the square root of sample size
             Spike_rate_SE = sd(Spike_rate)/sqrt(n()),
             Mean_photod_rate = mean(Photod_mean),
             ## SE of photodiode
             Photod_SE = sd(Photod_mean)/sqrt(n())
  ) %>%
  purrr::map(ungroup) %>%
  bind_rows() %>%
  ## We'll manually set the levels of the "Speed" column to ensure they plot in
  ## a desired order (fastest = highest, slowest = lowest)
  mutate(across(
    Speed,
    factor,
    levels = c("1024", "644", "407", "256", "128", "64", "32", "16", "8", "4")
  ))

## Generate the mean spike rate plot using ggplot
bin100_msr_plot <-
  dataslices_100 %>%
  ## The same code block can be used to generate either the mean spike rate
  ## (shown below) or photodiode trace (commented out)
  ggplot(aes(x = Time_stand,
             y = Mean_spike_rate #Mean_photod_rate
  )) +
  ## We'll actually start by placing red, yellow, and green vertical lines to
  ## distinguish between blank, stationary, and moving phases
  ## This comes first so that it is the bottom-most layer and doesn't obstruct
  ## the data
  geom_vline(xintercept = 0, col = "red") +
  geom_vline(xintercept = first(dataslices_100$Blank_end),
             col = "darkgoldenrod1") +
  geom_vline(xintercept = first(dataslices_100$Static_end),
             col = "forestgreen") +
  ## We'll use `geom_ribbon()` to shade in the SE traces
  geom_ribbon(aes(
    ymin = Mean_spike_rate - Spike_rate_SE,
    ymax = Mean_spike_rate + Spike_rate_SE
    # ymin = Mean_photod_rate - Photod_SE,
    # ymax = Mean_photod_rate + Photod_SE
  ),
  fill = "grey80") +
  ## `geom_line()` will be used to draw the mean spike rate itself on top of
  ## the SE traces
  geom_line(linewidth = 0.05) +
  ## Add a title to help us know what cell this is
  ggtitle(bin10_basenames[1]) +
  xlab("Time (sec)") +
  ylab("Spike rate (spikes/sec)") +
  ## To sub-plot by Speed and Direction, I typically use `facet_grid()`. This
  ## method allows me to explicitly declare what the row- and column-wise
  ## grouping variables are
  facet_grid(rows = vars(Speed), cols = vars(Direction)) +
  theme_classic()


#### 6.2.1 Export to PDF ####

## Use the `pdf()` function to start the graphics device driver for producing
## PDFs
## Aspects such as page size and centering mode can be adjusted
pdf(
  file = "./path/to/directory/filename.pdf",
  width = 10.5,
  height = 8,
  title = "2023-02-16_001 raster",
  paper = "USr",
  bg = "white",
  pagecentre = TRUE,
  colormodel = "srgb"
)
## Now add the plot to the PDF simply by calling plot()
plot(rasterplot)
## To declare an end to this PDF writing session, use `dev.off()`
dev.off()


#### 7 Polar direction tuning plots ####

#### 7.1 Data import and baseline rate measurement ####

## Read in the data
## Since there is only 1 example file, we'll simplify things by just
## subsetting bin10_filelist
bin10_data <-
  read_csv(bin10_filelist[1], show_col_types = FALSE) %>%
  as_tibble() %>%
  ## Again, we will set Speed as an ordered factor to help control plotting
  ## later on
  mutate(across(
    Speed,
    factor,
    levels = c("1024", "644", "407", "256", "128", "64", "32", "16", "8", "4")
  ))

## Extract all rows corresponding to our desired baseline epoch
baseline_df <-
  bin10_data %>%
  filter(Time_stand >= Blank_end + 0.5) %>% ## 0.5 sec after blank end
  filter(Time_stand <= Static_end - 0.05)   ## 0.05 sec before static end

## Compute the mean baseline
global_base <- mean(baseline_df$Spike_rate)
## Compute the SE
global_base_se <-
  sd(baseline_df$Spike_rate) / sqrt(length(baseline_df$Spike_rate)#/45
  )

## Construct a summary data frame
baseline_summary_df <-
  baseline_df %>%
  group_by(Speed) %>%
  summarize(
    speed_specific_baseline = mean(Spike_rate),
    speed_specific_baseline_se = sd(Spike_rate)/sqrt(length(Spike_rate)),
    global_baseline = global_base,
    global_baseline_se = global_base_se
  )

## This tibble contains speed-specific baselines (and SE) along with the global
## mean baseline (and SE)
baseline_summary_df

#### __naive ####
## global baseline values, full 3-sec motion period
naive_df <-
  bin10_data %>%
  filter(Time_stand > Static_end) %>%
  group_by(Speed, Direction, Replicate) %>%
  summarize(
    mean_spike_rate = mean(Spike_rate)
  )
naive_360 <-
  naive_df %>%
  filter(Direction == 0) %>%
  transmute(
    Speed = Speed,
    Direction = 360,
    Replicate = Replicate,
    mean_spike_rate = mean_spike_rate)
naive_vecsum_df <-
  naive_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  ## NOTE: FOLLOWING THE JNP AND CB PAPERS, WE ARE SUBTRACTING BASELINE HERE AND
  ## THEN IF ANY AVERAGED FIRING RATES ARE NEGATIVE, THEY ARE SHIFTED SO THE
  ## LOWEST ONE IS ZERO. THE GENERATED PLOTS STILL SHOW NON-BASELINE-SUBTRACTED
  ## FIRING RATES (AND BASELINE AS A RED RING), BUT COMPUTATION OF VECTOR SUM
  ## AND SI HAVE BEEN BASELINE SUBTRACTED (AND THE ENTIRE CURVE IS SHIFTED
  ## UPWARDS IF ANY PART IS NEGATIVE)
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      x = cos(Direction * pi / 180) * mean_spike_rate,
      y = sin(Direction * pi / 180) * mean_spike_rate,
      Speed = first(Speed)
  ) %>%
  map(summarise,
      x = mean(x),
      y = mean(y),
      Speed = first(Speed)) %>%
  map(transmute,
      vector_sum = (atan2(y, x) * 180 / pi) %% 360,
      Speed = first(Speed)
  ) %>%
  bind_rows()
naive_si_df<-
  naive_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      a = (sin(Direction * pi / 180) * mean_spike_rate),
      b = (cos(Direction * pi / 180) * mean_spike_rate),
      c = mean(mean_spike_rate),
      Speed = first(Speed)
  ) %>%
  map(summarise,
      a = mean(a),
      b = mean(b),
      c = mean(c),
      Speed = first(Speed)) %>%
  map(transmute,
      si = sqrt(a ^ 2 + b ^ 2) / c,
      Speed = first(Speed)
  ) %>%
  bind_rows()

naive_data <-
  naive_df %>%
  bind_rows(naive_360) %>%
  left_join(naive_vecsum_df, by = "Speed") %>%
  left_join(naive_si_df, by = "Speed") %>%
  drop_na(mean_spike_rate)
naive_max_y <- max(naive_data$mean_spike_rate)

naive <-
  naive_data %>%
  left_join(baseline_summary_df, by = "Speed") %>%
  ggplot(aes(x = Direction, y = mean_spike_rate)) +
  geom_ribbon(aes(
    x = Direction,
    ymin = global_baseline - global_baseline_se,
    ymax = global_baseline + global_baseline_se
  ),
  fill = "red") +
  stat_smooth(method = "glm", formula = y ~ ns(x,8), size = 0.3) +
  geom_point() +
  geom_label(aes(label = round(si, 2) , x = 315, y = naive_max_y * 1.2),
             size = 3) +
  geom_vline(aes(xintercept = vector_sum), colour = "grey30",
             size = 0.75) +
  coord_polar(direction = 1, start = pi/2) +
  scale_x_continuous(
    breaks = c(0, 90, 180, 270),
    expand = c(0, 0),
    limits = c(0, 360)
  )+
  facet_grid(rows = vars(Speed)) +
  ggtitle("full 3-sec motion") +
  theme_minimal()

#### __cb_it ####
## global baseline values, 40-200 msec motion period
cbit_df <-
  bin10_data %>%
  filter(Time_stand >= Static_end + 0.04) %>%
  filter(Time_stand <= Static_end + 0.2) %>%
  group_by(Speed, Direction, Replicate) %>%
  summarize(
    mean_spike_rate = mean(Spike_rate)
  )
cbit_360 <-
  cbit_df %>%
  filter(Direction == 0) %>%
  transmute(
    Speed = Speed,
    Direction = 360,
    Replicate = Replicate,
    mean_spike_rate = mean_spike_rate)
cbit_vecsum_df <-
  cbit_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      x = cos(Direction * pi / 180) * mean_spike_rate,
      y = sin(Direction * pi / 180) * mean_spike_rate,
      Speed = first(Speed)
  ) %>%
  map(summarise,
      x = mean(x),
      y = mean(y),
      Speed = first(Speed)) %>%
  map(transmute,
      vector_sum = (atan2(y, x) * 180 / pi) %% 360,
      Speed = first(Speed)
  ) %>%
  bind_rows()
cbit_si_df <-
  cbit_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      a = (sin(Direction * pi / 180) * mean_spike_rate),
      b = (cos(Direction * pi / 180) * mean_spike_rate),
      c = mean(mean_spike_rate),
      Speed = first(Speed)
  ) %>%
  map(summarise,
      a = mean(a),
      b = mean(b),
      c = mean(c),
      Speed = first(Speed)) %>%
  map(transmute,
      si = sqrt(a ^ 2 + b ^ 2) / c,
      Speed = first(Speed)
  ) %>%
  bind_rows()

cbit_data <-
  cbit_df %>%
  bind_rows(cbit_360) %>%
  left_join(cbit_vecsum_df, by = "Speed") %>%
  left_join(cbit_si_df, by = "Speed") %>%
  drop_na(mean_spike_rate)
cbit_max_y <- max(cbit_data$mean_spike_rate)

cb_it <-
  cbit_data %>%
  left_join(baseline_summary_df, by = "Speed") %>%
  ggplot(aes(x = Direction, y = mean_spike_rate)) +
  geom_ribbon(aes(
    x = Direction,
    ymin = global_baseline - global_baseline_se,
    ymax = global_baseline + global_baseline_se
  ),
  fill = "red") +
  stat_smooth(method = "glm", formula = y ~ ns(x,8), size = 0.3) +
  geom_point() +
  geom_label(aes(label = round(si, 2) , x = 315, y = cbit_max_y * 1.2),
             size = 3) +
  geom_vline(aes(xintercept = vector_sum), colour = "grey30",
             size = 0.75) +
  coord_polar(direction = 1, start = pi/2) +
  scale_x_continuous(
    breaks = c(0, 90, 180, 270),
    expand = c(0, 0),
    limits = c(0, 360)
  )+
  facet_grid(rows = vars(Speed)) +
  ggtitle("40-200 msec motion") +
  theme_minimal()

#### __cb_ss ####
## global baseline values, second half of motion period
cbss_df <-
  bin10_data %>%
  filter(Time_stand > Static_end + 1.5) %>%
  filter(Time_stand < 4.95) %>%
  group_by(Speed, Direction, Replicate) %>%
  summarize(
    mean_spike_rate = mean(Spike_rate)
  )
cbss_360 <-
  cbss_df %>%
  filter(Direction == 0) %>%
  transmute(
    Speed = Speed,
    Direction = 360,
    Replicate = Replicate,
    mean_spike_rate = mean_spike_rate)
cbss_vecsum_df <-
  cbss_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      x = cos(Direction * pi / 180) * mean_spike_rate,
      y = sin(Direction * pi / 180) * mean_spike_rate,
      Speed = first(Speed)
  ) %>%
  map(summarise,
      x = mean(x),
      y = mean(y),
      Speed = first(Speed)) %>%
  map(transmute,
      vector_sum = (atan2(y, x) * 180 / pi) %% 360,
      Speed = first(Speed)
  ) %>%
  bind_rows()
cbss_si_df <-
  cbss_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      a = (sin(Direction * pi / 180) * mean_spike_rate),
      b = (cos(Direction * pi / 180) * mean_spike_rate),
      c = mean(mean_spike_rate),
      Speed = first(Speed)
  ) %>%
  map(summarise,
      a = mean(a),
      b = mean(b),
      c = mean(c),
      Speed = first(Speed)) %>%
  map(transmute,
      si = sqrt(a ^ 2 + b ^ 2) / c,
      Speed = first(Speed)
  ) %>%
  bind_rows()

cbss_data <-
  cbss_df %>%
  bind_rows(cbss_360) %>%
  left_join(cbss_vecsum_df, by = "Speed") %>%
  left_join(cbss_si_df, by = "Speed") %>%
  drop_na(mean_spike_rate)
cbss_max_y <- max(cbss_data$mean_spike_rate)


cb_ss <-
  cbss_data %>%
  left_join(baseline_summary_df, by = "Speed") %>%
  ggplot(aes(x = Direction, y = mean_spike_rate)) +
  geom_ribbon(aes(
    x = Direction,
    ymin = global_baseline - global_baseline_se,
    ymax = global_baseline + global_baseline_se
  ),
  fill = "red") +
  stat_smooth(method = "glm", formula = y ~ ns(x,8), size = 0.3) +
  geom_point() +
  geom_label(aes(label = round(si, 2) , x = 315, y = cbss_max_y * 1.2),
             size = 3) +
  geom_vline(aes(xintercept = vector_sum), colour = "grey30",
             size = 0.75) +
  coord_polar(direction = 1, start = pi/2) +
  scale_x_continuous(
    breaks = c(0, 90, 180, 270),
    expand = c(0, 0),
    limits = c(0, 360)
  )+
  facet_grid(rows = vars(Speed)) +
  ggtitle("steady state motion") +
  theme_minimal()


#### __arb ####
## global baseline values, 0-500 msec motion period
arb_df <-
  bin10_data %>%
  filter(Time_stand >= Static_end + 0.001) %>%
  filter(Time_stand <= Static_end + 0.500) %>%
  group_by(Speed, Direction, Replicate) %>%
  summarize(
    mean_spike_rate = mean(Spike_rate)
  )
arb_360 <-
  arb_df %>%
  filter(Direction == 0) %>%
  transmute(
    Speed = Speed,
    Direction = 360,
    Replicate = Replicate,
    mean_spike_rate = mean_spike_rate)
arb_vecsum_df <-
  arb_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      x = cos(Direction * pi / 180) * mean_spike_rate,
      y = sin(Direction * pi / 180) * mean_spike_rate,
      Speed = first(Speed)
  ) %>%
  map(summarise,
      x = mean(x),
      y = mean(y),
      Speed = first(Speed)) %>%
  map(transmute,
      vector_sum = (atan2(y, x) * 180 / pi) %% 360,
      Speed = first(Speed)
  ) %>%
  bind_rows()
arb_si_df <-
  arb_df %>%
  ungroup() %>%
  drop_na(mean_spike_rate) %>%
  group_split(Speed) %>%
  map(group_by, Direction) %>%
  map(summarize,
      mean_spike_rate = mean(mean_spike_rate) - global_base,
      Speed = first(Speed)) %>%
  map(mutate,
      mean_spike_rate =
        case_when(min(mean_spike_rate) < 0 ~ mean_spike_rate + abs(min(mean_spike_rate)),
                  TRUE ~ mean_spike_rate)) %>%
  map(transmute,
      a = (sin(Direction * pi / 180) * mean_spike_rate),
      b = (cos(Direction * pi / 180) * mean_spike_rate),
      c = mean(mean_spike_rate),
      Speed = first(Speed)
  ) %>%
  map(summarise,
      a = mean(a),
      b = mean(b),
      c = mean(c),
      Speed = first(Speed)) %>%
  map(transmute,
      si = sqrt(a ^ 2 + b ^ 2) / c,
      Speed = first(Speed)
  ) %>%
  bind_rows()

arb_data <-
  arb_df %>%
  bind_rows(arb_360) %>%
  left_join(arb_vecsum_df, by = "Speed") %>%
  left_join(arb_si_df, by = "Speed") %>%
  drop_na(mean_spike_rate)
arb_max_y <- max(arb_data$mean_spike_rate)

arb <-
  arb_data %>%
  left_join(baseline_summary_df, by = "Speed") %>%
  ggplot(aes(x = Direction, y = mean_spike_rate)) +
  geom_ribbon(aes(
    x = Direction,
    ymin = global_baseline - global_baseline_se,
    ymax = global_baseline + global_baseline_se
  ),
  fill = "red") +
  stat_smooth(method = "glm", formula = y ~ ns(x,8), size = 0.3) +
  geom_point() +
  geom_label(aes(label = round(si, 2) , x = 315, y = arb_max_y * 1.2),
             size = 3) +
  geom_vline(aes(xintercept = vector_sum), colour = "grey30",
             size = 0.75) +
  coord_polar(direction = 1, start = pi/2) +
  scale_x_continuous(
    breaks = c(0, 90, 180, 270),
    expand = c(0, 0),
    limits = c(0, 360)
  )+
  facet_grid(rows = vars(Speed)) +
  ggtitle("0-500 msec motion") +
  theme_minimal()


#### __polar cowplot ####
cow_polar <-
  plot_grid(naive, cb_it, cb_ss, arb,
            nrow = 1)

pdf(file =
      paste0("./",
             bin10_basenames[i],
             "_polar_set.pdf"),
    width = 22, height = 17,
    pagecentre = TRUE, colormodel = "srgb")
plot(cow_polar)
dev.off()
