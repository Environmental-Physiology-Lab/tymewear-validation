## bind branch list, write to excel, return tracked file path
write_xlsx_target <- function(..., file_path) {
    df <- dplyr::bind_rows(...)
    writexl::write_xlsx(df, path = file_path)
    return(file_path)
}


## find file path
return_path <- function(folder = ".", id, device) {
    files <- list.files(folder, full.names = TRUE)
    matched <- files[grepl(id, files) & grepl(device, files)]

    if (length(matched) == 0) {
        stop("No file found matching id '", id, "' and device '", device, "'")
    }
    if (length(matched) > 1) {
        warning("Multiple files matched. Returning first match.")
    }

    return(matched[1L])
}


create_combined_data <- function(
    id = "BP17",
    parvo_data = targets::tar_read(parvo_data_BP17),
    tyme_data = targets::tar_read(tyme_data_BP17)
) {
    ## manual time sync correction for tymewear drift
    if (id == "BP17") {
        start_time <- min(parvo_data$timestamp)

        tyme_data <- tyme_data |>
            mutate(
                elapsed = as.numeric(timestamp - start_time, units = "mins"),
                timestamp = timestamp +
                    case_when(
                        between(elapsed, 22, 31) ~ 101,
                        between(elapsed, 35.6, 41.5) ~ 80,
                        between(elapsed, 45, 51.5) ~ 71,
                        between(elapsed, 56.1, 62.6) ~ 49,
                        between(elapsed, 67.8, 74.3) ~ 45,
                        .default = 0
                    ),
                elapsed = NULL
            )
    } else if (id == "BP10") {
        tyme_data$timestamp <- tyme_data$timestamp + (70.78 * 60)
    } else if (!id %in% c(sprintf("TW%02d", c(1:6)), "BP12")) {
        ## correct starting tymewear timestamp
        ## don't need this step for the early participants
        ## but probably will need it for future participants
        diff <- parvo_data$timestamp[1L] - tyme_data$timestamp[1L]
        tyme_data$timestamp <- tyme_data$timestamp + diff
    }

    ## convert to data.table and set keys for rolling join
    parvo_dt <- as.data.table(parvo_data)
    tyme_dt <- as.data.table(tyme_data)
    setkey(parvo_dt, timestamp)
    setkey(tyme_dt, timestamp)
    ## rolling join, find nearest match by bxb sample
    combined_data <- tyme_dt[parvo_dt, roll = "nearest"] |>
        arrange(timestamp) |>
        mutate(time = as.numeric(timestamp - timestamp[1L])) |>
        relocate(time) |>
        select(-timestamp)

    ## manual outlier rejection for one participant Parvo data
    if (id == "TW06") {
        combined_data <- mnirs::replace_mnirs(
            combined_data,
            nirs_channels = c(RR, VE, Vt, VO2kg),
            time_channel = time,
            outlier_cutoff = 3,
            width = 120,
            method = "none"
        )
    }

    return(combined_data)
}

## return data frame with all stages of cross-correlation matched BxB,
## optionally smoothed to a single named level
create_long_data <- function(
    id = "TW01",
    parvo_events = targets::tar_read(parvo_events_BP17),
    combined_data = targets::tar_read(combined_data_BP17),
    ## single named smoothing level, e.g. c(width = 1), c(span = 5), c(bin = 15)
    ## NULL returns unsmoothed BxB
    smooth = NULL
) {
    span <- c(-360, 120)

    purrr::imap(parvo_events$time, \(t_end, .i) {
        df <- combined_data |>
            ## select target interval range WITH ONSET & REC KINETICS
            filter(between(time, t_end + span[1L], t_end + span[2L])) |>
            ## global outlier detection on errant breaths
            mnirs::replace_mnirs(
                nirs_channels = c(br, ve, vt, RR, VE, Vt, VO2kg),
                time_channel = time,
                outlier_cutoff = 3,
                span = 30,
                method = "median"
            )

        ## cross-correlation bxb alignment (signed lag: +lag / -lead)
        ccf_values <- ccf(df$RR, df$br, lag.max = 500, plot = FALSE)
        t_lag <- ccf_values$lag[which.max(ccf_values$acf)]

        df <- df |>
            mutate(across(c(br, ve, vt), \(.x) data.table::shift(.x, t_lag)))

        ## bin averages over fixed time windows, else moving average filter
        if (isTRUE(names(smooth) == "bin")) {
            df <- df |>
                ## breath duration weights before overwriting per-sample time
                mutate(
                    dt = c(diff(time)[1L], diff(time)),
                    time = ceiling(time / smooth) * smooth,
                ) |>
                summarise(
                    .by = time,
                    across(c(br, ve, vt, RR, VE, Vt, VO2kg), \(.x) {
                        weighted.mean(.x, dt, na.rm = TRUE)
                    })
                ) |>
                mutate(
                    ## convert blank values to NA
                    across(c(br, ve, vt, RR, VE, Vt, VO2kg), \(.x) {
                        replace(.x, !is.finite(.x), NA_real_)
                    }),
                )
        } else if (!is.null(smooth)) {
            df <- rlang::inject(filter_mnirs(
                data = df,
                method = "moving_average",
                partial = TRUE,
                na.rm = TRUE,
                !!!smooth ## inject window
            ))
        }

        df |>
            ## select 5-min interval range omitting first & last 15-sec
            filter(between(time, t_end + span[1L], t_end + span[2L])) |>
            ## add identifying info and refactor time to zero
            mutate(
                id = id,
                stage = .i,
                position = factor(
                    if_else(stage %in% c(2, 4), "RP", "UP"),
                    levels = c("UP", "RP")
                ),
                smooth = if (is.null(smooth)) {
                    NULL
                } else {
                    paste(names(smooth), "=", smooth)
                },
                time = time - (t_end - 300), ## 0 == interval start time
            ) |>
            relocate(id, stage, any_of("smooth"))
    }) |>
        ## create single data frame for all participant samples
        purrr::list_rbind()
}
