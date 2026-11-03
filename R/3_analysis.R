## per-yvar tymewear column & display labels for analysis output
analysis_vars <- list(
    ## units are plotmath fragments, rendered via axis_label()
    RR = list(
        tyme = "br",
        unit_parvo = "breaths%.%min^'-1'",
        unit_tyme = "breaths%.%min^'-1'",
        unit_small = "'br/min'",
        label = "Respiratory Frequency"
    ),
    ## tymewear reports VE in uncalibrated volume units, not litres
    VE = list(
        tyme = "ve",
        unit_parvo = "L%.%min^'-1'",
        unit_tyme = "vol%.%min^'-1'",
        unit_small = "'vol/min'",
        label = "Minute Ventilation"
    )
)


## reshape paired device columns to long format for mixed models
## RR & VE pivoted separately then bound, to keep one row per device-sample
create_device_long <- function(df) {
    bind_cols(
        pivot_longer(df, c(br, RR), names_to = "device", values_to = "RR") |>
            select(-c(ve, VE)),
        pivot_longer(df, c(ve, VE), names_to = "device", values_to = "VE") |>
            select(VE)
    ) |>
        mutate(
            device = factor(
                if_else(device == "br", "Tymewear", "Parvo"),
                levels = c("Parvo", "Tymewear")
            ),
            stage = factor(stage),
            position = as_factor(if_else(stage %in% c(2, 4), "RP", "UP")),
            ## device-position combination, for position-specific smooths
            dp = interaction(device, position, drop = TRUE),
        )
}


## drop participants whose chest strap required shortening
filter_excluded_ids <- function(df, exclude = FALSE) {
    if (exclude) {
        filter_out(df, id %in% exclude_ids)
    } else {
        df
    }
}


## repeated measures correlation, accounting for within-participant covariance
fit_rmcorr <- function(df, yvar = c("RR", "VE")) {
    yvar <- match.arg(yvar)
    rmcorr::rmcorr(
        participant = "id",
        measure1 = yvar,
        measure2 = analysis_vars[[yvar]]$tyme,
        dataset = df
    )
}


## bland-altman limits of agreement, nested within participant
## x/y flipped so delta = Parvo - Tymewear
fit_agreement <- function(df, yvar = c("RR", "VE"), log_tf = FALSE) {
    yvar <- match.arg(yvar)
    tyme <- analysis_vars[[yvar]]$tyme

    SimplyAgree::agreement_limit(
        data = drop_na(df, all_of(c(yvar, tyme))),
        x = tyme,
        y = yvar,
        id = "id",
        data_type = "nest",
        loa_calc = "blandaltman",
        ## VE shows proportional bias across the measurement range; the
        ## log-ratio fits stay flat so the sympercent bias is a single value
        prop_bias = yvar == "VE" && !log_tf,
        log_tf = log_tf,
        log_tf_display = "sympercent"
    )
}


## TOST equivalence of a contrast against +/-3% log-ratio bounds
## (van Hooren et al 2023), 90% CI matching TOST at alpha = 0.05. lme4 does not data-mask, so build formula by name
## contrast = "device" tests Tymewear vs Parvo; "position" tests RP vs UP with
## the device effect and its interaction held in the model
fit_equivalence <- function(
    df_long,
    yvar = c("RR", "VE"),
    bound = log(1.03),
    contrast = c("device", "position")
) {
    yvar <- match.arg(yvar)
    contrast <- match.arg(contrast)

    ## position contrast allows the time-course to differ per device-position,
    ## so the smooth varies by dp; device contrast holds one shared time term
    if (contrast == "device") {
        model <- lmerTest::lmer(
            as.formula(str_glue(
                "log({yvar}) ~ device + time + (1 | id / stage)"
            )),
            data = df_long
        )

        ## mean within-participant log-ratio for the contrast of interest
        eq <- marginaleffects::avg_comparisons(
            model,
            variables = "device",
            conf_level = 0.90,
            equivalence = c(-bound, bound)
        )
    } else {
        model <- mgcv::gam(
            as.formula(str_glue(
                "log({yvar}) ~ s(time, by = dp, k = 10) + device * position +
                s(id, bs = 're')"
            )),
            data = df_long
        )

        eq <- marginaleffects::avg_comparisons(
            model,
            variables = "position",
            by = "device",
            hypothesis = ~reference,
            conf_level = 0.90,
            equivalence = c(-bound, bound)
        )
    }

    return(list(model = model, eq = eq, bound = bound, contrast = contrast))
}


## naive correlation ignoring repeated measures, for contrast with rmcorr
fit_cor <- function(df, yvar = c("RR", "VE")) {
    yvar <- match.arg(yvar)
    cor.test(df[[yvar]], df[[analysis_vars[[yvar]]$tyme]])
}


## intraclass correlation treating device as the repeated item
## long data carries both RR & VE columns, so drop incomplete cells first
fit_icc <- function(df_long, yvar = c("RR", "VE")) {
    yvar <- match.arg(yvar)
    SimplyAgree::reli_aov(
        measure = yvar,
        item = "device",
        id = "id",
        data = as.data.frame(drop_na(df_long, all_of(yvar))),
        other_ci = TRUE
    )
}


## absolute & percent error, plus log-ratio CV partitioned per device
## superseded by rmcorr & bland-altman, retained for reference
summarise_error <- function(df, yvar = c("RR", "VE")) {
    yvar <- match.arg(yvar)
    tyme <- analysis_vars[[yvar]]$tyme

    ## log-ratio needs strictly positive pairs, and >= 2 samples per
    ## participant to centre within-participant
    cv <- df |>
        filter(
            is.finite(.data[[yvar]]),
            is.finite(.data[[tyme]]),
            .data[[yvar]] > 0,
            .data[[tyme]] > 0
        ) |>
        mutate(log_ratio = log(.data[[yvar]] / .data[[tyme]])) |>
        filter(.by = id, n() >= 2) |>
        mutate(.by = id, centred_ratio = log_ratio - mean(log_ratio)) |>
        summarise(
            n_participants = n_distinct(id),
            n_pairs = n(),
            ## residual sd of the log ratio, df corrected for participant means
            sd_log_difference = sqrt(
                sum(centred_ratio^2) / (n_pairs - n_participants)
            ),
            ## typical within-participant disagreement between devices
            cv_between_devices = 100 * sqrt(exp(sd_log_difference^2) - 1),
            ## per-device share, assuming equal independent error
            cv_each_device = 100 * sqrt(exp(sd_log_difference^2 / 2) - 1),
            ## positive -> Parvo higher than Tymewear
            mean_bias = 100 * (exp(mean(log_ratio)) - 1)
        )

    error <- df |>
        summarise(
            MAE = mean(abs(.data[[tyme]] - .data[[yvar]]), na.rm = TRUE),
            MAPE = mean(
                abs(.data[[tyme]] - .data[[yvar]]) / .data[[yvar]] * 100,
                na.rm = TRUE
            )
        )

    return(bind_cols(cv, error))
}
