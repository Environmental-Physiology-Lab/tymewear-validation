## na.rm overrides to silence warnings
geom_point <- \(...) ggplot2::geom_point(na.rm = TRUE, ...)
geom_line <- \(...) ggplot2::geom_line(na.rm = TRUE, ...)

## opaque points over gridlines: outlined ring, white knockout, translucent fill
## mapping = NULL inherits the plot-level fill aesthetic
geom_point_stack <- function(mapping = NULL, sizes = c(1.4, 1.2, 2)) {
    return(list(
        geom_point(fill = "white", size = sizes[1], shape = 21, stroke = 1),
        geom_point(colour = "white", size = sizes[2]),
        geom_point(
            mapping,
            size = sizes[3],
            shape = 21,
            stroke = 0,
            alpha = 0.6
        )
    ))
}

## axis title with plotmath units, e.g. "Parvo VE (L·min^-1)"
## text is split on spaces so plotmath treats each word as a separate token
axis_label <- function(text, unit) {
    parse(
        text = str_glue(
            "bold({str_replace_all(text, ' ', '~')}~'('*{unit}*')')"
        )
    )
}

## common plot settings
## a function, not a list: ggplot layer objects carry environments that hash
## unstably, invalidating unrelated targets whenever this file changes
plot_settings <- function() {
    list(
        theme_epl(),
        scale_x_continuous(
            name = "Time (mm:ss)",
            breaks = breaks_timespan(n = 4),
            labels = format_hmmss,
            expand = expansion(0.01)
        ),
        scale_colour_epl(
            name = NULL,
            aesthetics = c("colour", "fill")
        ),
        guides(
            colour = guide_legend(
                override.aes = list(size = 4, shape = 19, alpha = 1)
            )
        ),
        geom_point_stack()
    )
}

## per-yvar tyme column, axis limits, and label
yvar_labels <- list(
    RR = list(
        yvar_tyme = "br",
        ylims = c(0, 75),
        ylab = "Respiratory Frequency (brpm)"
    ),
    Vt = list(
        yvar_tyme = "vt",
        ylims = c(0, 6),
        ylab = "Tidal Volume (L)"
    ),
    VE = list(
        yvar_tyme = "ve",
        ylims = c(0, 260),
        ylab = expression(bold(Ventilation ~ '(' * vol %.% min^'-1' * ')'))
    )
)

## save plot to file & return file path for {targets} pipeline
save_plot_target <- function(plot, file_path, width = 180, height = 180) {
    ggsave(
        file_path,
        plot = plot,
        width = width,
        height = height,
        dpi = 300,
        units = "mm",
        device = "png",
        bg = "white",
        create.dir = TRUE
    )
    return(file_path)
}

plot_session <- function(
    id = "TW04",
    combined_data = tar_read(combined_data_TW05),
    parvo_events = tar_read(parvo_events_TW05),
    yvar = c("RR", "Vt", "VE")
) {
    yvar <- match.arg(yvar)
    labels <- yvar_labels[[yvar]]

    ## reshape paired Tymewear/Parvo columns to long format
    df <- combined_data |>
        rename(
            Tymewear = any_of(labels$yvar_tyme),
            Parvo = any_of(yvar)
        ) |>
        pivot_longer(
            cols = c(Tymewear, Parvo),
            names_to = "name",
            values_to = "value",
        )

    interval_times <- c(
        vapply(parvo_events$time, \(x) x - c(360, 300, 0, -60), numeric(4))
    )

    plot <- ggplot(df, aes(x = time, y = value, colour = name, fill = name)) +
        plot_settings() +
        scale_y_continuous(
            name = labels$ylab,
            breaks = scales::pretty_breaks(n = 6),
            expand = expansion(0.01)
        ) +
        coord_cartesian(xlim = c(NA, NA), ylim = labels$ylims) +
        geom_vline(xintercept = interval_times, linetype = "dotted")

    return(save_plot_target(
        plot,
        figures_path(yvar, str_glue("session_{id}_{yvar}.png")),
        height = 180 * 2 / 3
    ))
}

plot_stages <- function(
    id = "TW04",
    long_id_data = tar_read(long_id_data_TW04),
    yvar = c("RR", "Vt", "VE")
) {
    yvar <- match.arg(yvar)
    labels <- yvar_labels[[yvar]]

    ## reshape paired Tymewear/Parvo columns to long format
    df <- long_id_data |>
        rename(
            Tymewear = any_of(labels$yvar_tyme),
            Parvo = any_of(yvar)
        ) |>
        pivot_longer(
            cols = c(Tymewear, Parvo),
            names_to = "name",
            values_to = "value",
        )

    plot <- ggplot(df, aes(x = time, y = value, colour = name, fill = name)) +
        plot_settings() +
        labs(
            title = paste("Participant:", id),
            y = labels$ylab
        ) +
        facet_wrap(~stage, ncol = 2) +
        theme(
            plot.title = element_text(face = "bold"),
            panel.grid.major.y = element_line(),
            strip.background = element_blank(),
            strip.text.x = element_blank()
        ) +
        scale_y_continuous(
            breaks = scales::pretty_breaks(n = 6),
            expand = expansion(0.01)
        ) +
        coord_cartesian(xlim = c(NA, NA), ylim = labels$ylims) +
        geom_vline(xintercept = c(0, 300), linetype = "dotted")

    return(save_plot_target(
        plot,
        figures_path(yvar, str_glue("stage_{id}_{yvar}.png"))
    ))
}


## analysis plots ==============================================
## NB. geom_line is shadowed above with na.rm = TRUE, so line layers below are
## qualified with ggplot2:: where the unmodified geom is wanted.

## per-yvar axis limits, annotation placement, and rounding.
## VE_excl shares VE units & labels, differing only in the equivalence x-range
analysis_plot_params <- list(
    RR = list(
        yvar = "RR",
        ylab = "bold(italic('f')['R'])",
        rmcorr_lims = c(0, 70),
        loa_digits = 1,
        eq_xlim = c(-8, 8)
    ),
    VE = list(
        yvar = "VE",
        ylab = "bold(dot(V)['E'])",
        rmcorr_lims = c(0, 270),
        loa_digits = 0,
        eq_xlim = c(-13, 50)
    ),
    VE_excl = list(
        yvar = "VE",
        ylab = "bold(dot(V)['E'])",
        rmcorr_lims = c(0, 170),
        loa_digits = 0,
        eq_xlim = c(-13, 30)
    ),
    ## position contrast, equivalence plot only
    RR_pos = list(yvar = "RR", eq_xlim = c(-8, 8)),
    VE_pos = list(yvar = "VE", eq_xlim = c(-10, 8))
)

## equivalence decision region fills & labels, shared across yvars
## NA x anchors to the lower equivalence bound, filled in per yvar
eq_labels <- tribble(
    ~x , ~y  , ~hjust , ~label              ,
    NA , 2.8 , -0.1   , "Non-inferiority →" ,
     0 , 1.8 ,  0.5   , "Equivalence"       ,
     0 , 0.8 ,  1.1   , "← Inferiority"     ,
     0 , 0.8 , -0.1   , "Superiority →"     ,
)

eq_fills <- c(
    "Non-inferiority" = "#00468Bff",
    "Equivalence" = "#ff7f00",
    "Inferiority" = "#db5555",
    "Superiority" = "#42B540FF"
)

## decision regions bounded by the equivalence limit b
eq_regions <- function(b) {
    return(tribble(
        ~xmin , ~xmax , ~ymin , ~ymax , ~fill             ,
        -b    , Inf   ,     2 ,     3 , "Non-inferiority" ,
        -b    , b     ,     1 ,     2 , "Equivalence"     ,
        -Inf  ,     0 ,     0 ,     1 , "Inferiority"     ,
            0 , Inf   ,     0 ,     1 , "Superiority"
    ))
}

## beeswarm & interval markers, repeated across the equivalence plots
quasirandom_stack <- function(diff) {
    swarm <- \(...) {
        ggbeeswarm::geom_quasirandom(
            data = tibble(),
            aes(x = diff, y = 1.5),
            orientation = "y",
            cex = 2,
            width = 1 / 3,
            ...
        )
    }

    return(list(
        swarm(shape = 21, stroke = 1),
        swarm(colour = "#ffbf7f"),
        swarm(fill = "black", shape = 21, stroke = 0, alpha = 0.01)
    ))
}

## per-participant panels for a single interval, both devices overlaid
## exclude_ids flags excluded participants with a red panel border
plot_descriptive <- function(
    df = tar_read(results_data),
    yvar = c("RR", "VE"),
    exclude_ids = NULL
) {
    yvar <- match.arg(yvar)
    tyme <- analysis_vars[[yvar]]$tyme

    ylab <- switch(
        yvar,
        RR = expression(bold(
            italic('f')['R'] ~ '(' * breaths %.% min^'-1' * ')'
        )),
        VE = expression(bold(dot(V)['E'] ~ '(' * L %.% min^'-1' * ')'))
    )
    ybreaks <- switch(yvar, RR = seq(0, 60, 20), VE = seq(0, 200, 50))

    ## flag excluded-participant panels with red border + faint tint
    highlight <- if (!is.null(exclude_ids)) {
        geom_rect(
            data = ~ distinct(filter(.x, id %in% exclude_ids), id),
            aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            colour = "red",
            fill = "red",
            alpha = 0.05,
            linewidth = 1
        )
    }

    p <- filter(df, stage == 4) |>
        ggplot(aes(time)) +
        facet_wrap(~id) +
        labs(x = "Time (mm:ss)", y = ylab) +
        mnirs::theme_mnirs(border = "full") +
        theme(
            strip.text = element_blank(),
            panel.grid.major.y = element_line(
                linewidth = 0.333,
                colour = "grey80"
            ),
            ## nudge outermost x-axis labels inboard to avoid panel overhang
            axis.text.x = element_text(
                colour = "grey40",
                hjust = c(0.3, 0.5, 0.5, 0.6)
            )
        ) +
        scale_x_continuous(
            breaks = seq(0, 360, 120),
            labels = c("0", "2:00", "4:00", "6:00")
        ) +
        scale_y_continuous(
            limits = c(0, NA),
            breaks = ybreaks,
            expand = expansion(c(0, 0.03))
        ) +
        scale_colour_mnirs(name = NULL) +
        highlight +
        geom_line(
            aes(y = .data[[tyme]], colour = "Tymewear"),
            linewidth = 0.5
        ) +
        geom_line(aes(y = .data[[yvar]], colour = "Parvo"), linewidth = 0.5)

    return(p)
}


## repeated measures correlation, per-participant fitted slopes
plot_rmcorr <- function(
    fit = tar_read(rmcorr_fit_RR),
    yvar = c("RR", "VE", "VE_excl")
) {
    yvar <- match.arg(yvar)
    params <- analysis_plot_params[[yvar]]
    vars <- analysis_vars[[params$yvar]]

    data <- tibble(fit$model$model) |>
        mutate(fitted = fit$model$fitted.values)

    label <- str_glue(
        "rmcorr = {r} [{ci}], {p}",
        r = mnirs:::signif_trailing(fit$r, 2),
        ci = paste(mnirs:::signif_trailing(fit$CI, 2), collapse = ', '),
        p = mnirs:::signif_pvalue(fit$p)
    )

    p <- ggplot(data, aes(x = Measure1, y = Measure2, fill = Participant)) +
        mnirs::theme_mnirs() +
        theme(
            legend.position = "none",
            plot.title = element_text(face = "bold"),
            plot.title.position = "plot",
            # plot.margin = margin_auto(10, 5, 10, 10),
        ) +
        coord_cartesian(
            xlim = params$rmcorr_lims,
            ylim = params$rmcorr_lims
        ) +
        scale_x_continuous(
            name = axis_label(str_glue("Parvo {params$ylab}"), vars$unit_parvo),
            expand = expansion(0)
        ) +
        scale_y_continuous(
            name = axis_label(
                str_glue("Tymewear {params$ylab}"),
                vars$unit_tyme
            ),
            expand = expansion(0)
        ) +
        geom_point_stack(aes(fill = Participant), c(2.4, 2.3, 3.0)) +
        geom_line(aes(y = fitted), colour = "grey20", linewidth = 1.8) +
        geom_line(aes(y = fitted, colour = Participant), linewidth = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
        annotate(
            "text",
            x = -Inf,
            y = Inf,
            label = label,
            size = 5,
            hjust = -0.015,
            vjust = 1
        )

    return(p)
}


## bland-altman difference against mean, with bias and limits of agreement
plot_agreement <- function(
    fit = tar_read(agree_fit_VE_excl),
    yvar = c("RR", "VE", "VE_excl")
) {
    yvar <- match.arg(yvar)
    params <- analysis_plot_params[[yvar]]
    vars <- analysis_vars[[params$yvar]]

    data <- tibble(plot(fit)$data)

    ## prop_bias fits return LoA at the min/mean/max of the device means, so
    ## bias & limits slope across the panel; flat fits carry a single "overall"
    sloped <- is.numeric(fit$loa$avg)

    loa <- fit$loa |>
        slice(c(1, 3)) |>
        select(avg, lower_loa, bias, upper_loa) |>
        pivot_longer(-avg, names_to = "line", values_to = "y") |>
        mutate(
            line = factor(line, c("lower_loa", "bias", "upper_loa")),
            vjust = if_else(line == "lower_loa", 1.2, -0.2),
            ## geom_labelline() drops the hyphen-minus when tracing text along a
            ## path, so labels are parsed as plotmath, which renders the sign as
            ## a minus glyph. The value stays a quoted literal to keep the
            ## trailing zeros from signif_trailing(), with the sign split out in
            ## front of it so plotmath treats it as an operator, not text.
            label_expr = str_replace(
                str_glue(
                    "'{mnirs:::signif_trailing(y, params$loa_digits)}'",
                    "~{vars$unit_small}"
                ),
                "^'-",
                "-'"
            )
        )

    ## spread flat limits across the x-range so both cases draw as labelled
    ## line segments
    if (!sloped) {
        loa <- cross_join(select(loa, -avg), tibble(avg = range(data$mean)))
    }

    ## geom_labelline draws one label per group: sloped limits are labelled at
    ## both ends (lines first, then a text-only overlay), flat at the right end
    ends <- tibble(
        hjust = if (sloped) c(0.02, 0.98) else 0.98,
        text_only = if (sloped) c(FALSE, TRUE) else FALSE
    )

    ## one layer per line and end, because parse = TRUE makes geomtextpath apply
    ## the first group's vjust to every group sharing a layer
    loa_layers <- loa |>
        cross_join(ends) |>
        mutate(
            .by = c(line, hjust),
            ## label the end of the line nearest the requested hjust
            label = label_expr[which.max(if_else(hjust > 0.5, avg, -avg))]
        ) |>
        group_split(line, hjust) |>
        lapply(\(d) geomtextpath::geom_labelline(
            data = d,
            aes(x = avg, y = y, group = line, label = label),
            parse = TRUE,
            linetype = "dashed",
            linewidth = 0.7,
            size = 4.8,
            hjust = d$hjust[1],
            vjust = d$vjust[1],
            text_smoothing = 0,
            text_only = d$text_only[1],
            ## semi-opaque box separates labels from the point cloud
            boxcolour = NA,
            boxlinewidth = 0,
            fill = alpha("white", 0.8),
            label.padding = unit(0.1, "lines")
        ))

    ggplot(data, aes(x = mean, y = delta)) +
        labs(subtitle = "") +
        mnirs::theme_mnirs() +
        theme(
            legend.position = "none",
            plot.title = element_text(face = "bold"),
            plot.title.position = "plot",
            axis.title.y = element_text(hjust = 1),
        ) +
        coord_cartesian(
            # xlim = c(0, max(data$mean) + min(data$mean)),
            ylim = max(abs(data$delta)) * c(-1, 1)
        ) +
        scale_x_continuous(
            name = axis_label(
                str_glue("Means of {params$ylab}"),
                vars$unit_tyme
            ),
            breaks = scales::breaks_pretty(n = 5)
        ) +
        scale_y_continuous(
            name = axis_label(
                str_glue("Differences between {params$ylab}"),
                vars$unit_tyme
            ),
            breaks = scales::breaks_pretty(n = 7)
        ) +
        geom_point_stack(aes(fill = id), c(2.4, 2.3, 3.0)) +
        geom_hline(yintercept = 0, linetype = "dotted") +
        loa_layers
}


## TOST decision regions with the contrast and paired differences
plot_equivalence <- function(
    fit = tar_read(eq_fit_RR),
    df = tar_read(analysis_data_RR),
    yvar = c("RR", "VE", "VE_excl", "RR_pos", "VE_pos")
) {
    yvar <- match.arg(yvar)
    params <- analysis_plot_params[[yvar]]
    vars <- analysis_vars[[params$yvar]]
    eq <- fit$eq
    by_position <- identical(fit$contrast, "position")

    ## equivalence bound as symmetric percent
    b <- 100 * fit$bound

    ## individual paired differences on the symmetric-percent scale
    ## device pairs row-wise; position pairs participant-device log-means
    ## across the racing & upright stages
    diff <- if (by_position) {
        df |>
            summarise(
                .by = c(id, position),
                across(all_of(c(vars$tyme, params$yvar)), \(x) {
                    mean(log(x), na.rm = TRUE)
                })
            ) |>
            pivot_longer(-c(id, position), names_to = "device") |>
            pivot_wider(names_from = position, values_from = value) |>
            mutate(diff = 100 * (RP - UP)) |>
            pull(diff)
    } else {
        100 * log(df[[vars$tyme]] / df[[params$yvar]])
    }

    regions <- eq_regions(b)
    ## non-inferiority arrow anchors on the lower bound
    labels <- mutate(eq_labels, x = coalesce(x, -b))

    ## RR spans a narrow axis, so its baseline stops at the bounds and breaks
    ## are set explicitly; VE runs full width with default breaks
    segment <- if (yvar == "RR") {
        geom_segment(
            aes(x = -b + 0.1, xend = b - 0.1, y = y),
            data = tibble(y = c(1, 2)),
            colour = "grey10"
        )
    } else {
        geom_segment(
            aes(x = -Inf, xend = Inf, y = y),
            data = tibble(y = c(1, 2)),
            colour = "grey10"
        )
    }

    title <- if (by_position) {
        expression(paste(
            bold("racing"),
            " vs reference ",
            bold("upright"),
            " position"
        ))
    } else {
        expression(paste(
            bold("Tymewear"),
            " vs criterion ",
            bold("ParvoMedics")
        ))
    }

    p <- ggplot() +
        mnirs::theme_mnirs() +
        labs(
            title = title,
            subtitle = bquote(paste(
                bold(.(str_glue("±{round(b)}% equivalence bounds"))),
                " on the log-ratio scale"
            )),
            x = str_glue("{vars$label} difference (symmetric %)"),
            y = NULL
        ) +
        theme(
            plot.subtitle = element_text(size = rel(0.9)),
            legend.position = "none",
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()
        ) +
        coord_cartesian(xlim = params$eq_xlim, ylim = c(0, 3)) +
        scale_x_continuous(expand = expansion(0.03)) +
        scale_y_continuous(expand = expansion(0)) +
        scale_fill_manual(values = eq_fills) +
        geom_rect(
            data = regions,
            aes(
                xmin = xmin,
                xmax = xmax,
                ymin = ymin,
                ymax = ymax,
                fill = fill
            ),
            alpha = 0.5
        ) +
        geom_vline(
            xintercept = c(-b, 0, b),
            linetype = c("dashed", "solid", "dashed")
        ) +
        segment +
        # quasirandom_stack(diff) +
        geom_label(
            data = labels,
            aes(x = x, y = y, label = label),
            size = 5,
            linewidth = 0.3,
            hjust = labels$hjust
        ) +
        geom_errorbar(
            aes(
                xmin = 100 * eq$conf.low,
                xmax = 100 * eq$conf.high,
                y = 1.5
            ),
            linewidth = 1,
            width = 0.1
        ) +
        geom_point(aes(x = 100 * eq$estimate, y = 1.5), size = 1)

    return(p)
}


## population-level time-course from the position gam, per device & position
## model is fit on log(yvar), so predictions back-transform to observed units
plot_position <- function(
    fit = tar_read(eq_pos_fit_RR_pos),
    df_long = tar_read(results_long),
    yvar = c("RR_pos", "VE_pos")
) {
    yvar <- match.arg(yvar)
    params <- analysis_plot_params[[yvar]]
    vars <- analysis_vars[[params$yvar]]

    ## random effects held at their modes, giving population-level curves
    pred <- marginaleffects::predictions(
        fit$model,
        newdata = marginaleffects::datagrid(
            model = fit$model,
            time = mnirs:::seq_range(df_long$time, length.out = 100),
            position = unique,
            device = unique
        ),
        transform = exp
    ) |>
        tibble()

    p <- ggplot(pred, aes(x = time, y = estimate)) +
        aes(colour = device, fill = device, linetype = position) +
        labs(subtitle = "", linetype = NULL) +
        theme(
            ## legends stacked inside plot area, device above position
            legend.position = "inside",
            legend.position.inside = c(0.99, 1.1),
            legend.justification = c(1, 1),
            legend.direction = "horizontal",
            legend.box = "vertical",
            legend.box.just = "right",
            legend.box.spacing = unit(0, "pt"),
            legend.spacing.y = unit(0, "pt"),
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_blank(),
            legend.key = element_blank(),
            panel.grid.major.y = element_line()
        ) +
        scale_x_continuous(
            name = "Time (mm:ss)",
            breaks = breaks_timespan(n = 4),
            labels = format_hmmss,
            expand = expansion(0.02)
        ) +
        scale_y_continuous(
            name = axis_label(vars$label, vars$unit_tyme),
            limits = c(0, NA),
            breaks = scales::pretty_breaks(n = 6),
            expand = expansion(0)
        ) +
        scale_colour_mnirs(name = "Device", aesthetics = c("colour", "fill")) +
        scale_linetype_manual(
            name = "Position",
            values = c(UP = "solid", RP = "22")
        ) +
        guides(
            colour = guide_legend(order = 1, nrow = 1),
            linetype = guide_legend(
                order = 2,
                nrow = 1,
                override.aes = list(colour = "grey30")
            )
        ) +
        geom_point(
            data = df_long,
            aes(y = .data[[params$yvar]]),
            size = 1.2,
            alpha = 0.1,
            shape = 16,
            show.legend = FALSE
        ) +
        geom_ribbon(
            aes(ymin = conf.low, ymax = conf.high),
            colour = NA,
            alpha = 0.2,
            show.legend = FALSE
        ) +
        geom_line(linewidth = 0.9)

    return(p)
}
