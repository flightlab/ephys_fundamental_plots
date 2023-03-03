############################### package loading ################################
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
      Temporal_Frequency = `Temporal Frequency`
    )
  ## Set up a `SF_cpd` column that tranlates SFs to cycles per degree
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

  ## The log file does not have time = 0, so set up a separate tibble to
  ## add this info in later. Some of the metadata will just be filler for now.
  initial <- tibble(
    Trial = "initialization",
    Spatial_Frequency = csv_data_sets[[i]]$Spatial_Frequency[1],
    SF_cpd = csv_data_sets[[i]]$SF_cpd[1],
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
  if(!median(round(diff(stim_change_channel[[5]][,1]),1)) == 5) {
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
    ## Join in the spike data
    left_join(all_spike_dat, by = "Time_char") %>%
    select(-Time.y) %>%
    rename(Time = Time.x) %>%
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
    )

  ## Add info to metadata
  metadata_one_full <-
    first_csv %>%
    mutate(
      Speed = Temporal_Frequency/SF_cpd,
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
    all_spike_dat, all_spike_dat_tmp, first_cell, all_cells,
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

## Set bin size here
## Units are in ms (e.g. 10 = 10ms)
bin_size = 10 ## 10 or 100 or 1 (1 = "unbinned")

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
            select(Speed, SF_cpd, Temporal_Frequency, Direction,
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
              SF_cpd = m$SF_cpd[1],
              Temporal_Frequency = m$Temporal_Frequency[1],
              Speed = m$Speed[1],
              Direction = m$Direction[1],
              Replicate = k
            ) %>%
            ## Bring stim info to first few columns
            select(Speed, SF_cpd, Temporal_Frequency, Direction,
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

################################### export csv #################################
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

################################# raster plot ##################################
## read in data
unbinned_data <-
  read_csv("./data/04132022_009m_unbinned.csv") %>%
  as_tibble()

## Generate the ggplot
rasterplot <-
  unbinned_data %>%
  ## Remove any rows where spiking does not occur in the Spikes column
  filter(Spikes == 1) %>%
  ## Convert Trial and Speed into factors and specify their level ordering
  ## This will make it easier to get the subplots in the order we want them
  mutate(Trial = factor(Trial, levels = c("blank", "stationary", "moving")),
         Speed = factor(Speed, levels = c(1024, 256, 32, 4, 0.5, 0.062))) %>%
  ggplot(aes(x = Time_stand, y = Replicate)) +
  ## The next three blocks will undershade each subplot according to stimulus
  ## phase (i.e., blank, stationary, moving)
  annotate("rect",
           xmin = 0, xmax = first(unbinned_data$Blank_end),
           ymin = 0.5, ymax = 10.5,
           alpha = 0.1, color = NA, fill = "red") +
  annotate("rect",
           xmin = first(unbinned_data$Blank_end),
           xmax = first(unbinned_data$Static_end),
           ymin = 0.5, ymax = 10.5,
           alpha = 0.1, color = NA, fill = "darkgoldenrod1") +
  annotate("rect",
           xmin = first(unbinned_data$Static_end), xmax = 5,
           ymin = 0.5, ymax = 10.5,
           alpha = 0.1, color = NA, fill = "forestgreen") +
  ## Up to 10 replicates were used, so we will force the y-axis to go to 10
  scale_y_continuous(
    limits = c(0.5, 10.5),
    expand = c(0, 0),
    breaks = c(5, 10)
  ) +
  ## There are multiple ways to plot a spike event. Since 100% of the rows in
  ## this filtered data set are spike events, we can simply plot a symbol at
  ## each time (Time_stand) that appears in the data. The `|` symbol is a good
  ## choice.
  geom_point(pch = '|', size = 1.5) +
  xlab("Time (sec)") +
  ggtitle(paste0("04132022_009m raster")) +
  ## Use facet_grid() to create a grid of subplots. Rows will correspond to
  ## Speeds, and columns correspond to Directions
  facet_grid(rows = vars(Speed), cols = vars(Direction)) +
  theme_classic() +
  theme(legend.position = 'none',
        panel.spacing = unit(0.1, "lines"))

rasterplot

## Should you elect to export this as a PDF
pdf(
  file = "./path/to/directory/filename.pdf",
  width = 10.5,
  height = 8,
  title = "04132022_009m raster",
  paper = "USr",
  bg = "white",
  pagecentre = TRUE,
  colormodel = "srgb"
)
plot(rasterplot)
dev.off()
