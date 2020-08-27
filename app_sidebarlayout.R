for (p in c("shiny", "dplyr", "tidyr", "gtools", "tibble", "stringr", "ggplot2", "DT"))
    library(p, character.only = TRUE)

dirichlet_diff <- function(k, d, prop = 1, n = 1000) { 
    # k: "shrinkage" factor
    # d: data frame with row names giving the level names, and each columns the cohort-wise counts
    # prop: proportion of levels to be different (default 100% = 1)
    # n: number of samples (larger n => more precise estimates and longer run-time)
    
    diri_tibble <- function(a)
        # Create a tibble with n draws from a distribution with the alpha parameter = counts per level of categorical variable
        as_tibble(rdirichlet(n, a), .name_repair = ~ str_replace_all(row.names(d), " ", "_"))
    
    # Posterior dirichlet distributions (dirichlets prior + multinomial likelihoods = dirichlets posterior)
    post_params <- d + k
    post_samples <- bind_rows(lapply(names(post_params), 
                                     function(.) mutate(diri_tibble(post_params[[.]]), cohort = .)))
    
    # Make all pairwise comparisons
    d <- combn(names(post_params), 2, simplify = FALSE) %>%
        # For each pair, find difference between each draw for each level, give appropriate label
        lapply(function(.) mutate(select(filter(post_samples, cohort == .[[1]]), -cohort) - 
                                      select(filter(post_samples, cohort == .[[2]]), -cohort),
                                  comparison = paste(., collapse = "_vs_"))) %>%
        bind_rows() %>%
        gather(level, p_diff, -comparison) # Make long format
    
    q <- group_by(d, comparison, level) %>%
        do(tibble(quantile = paste0("p", c(2.5, 50, 97.5)), 
                  value = quantile(.$p_diff, c(0.025, 0.5, 0.975)))) %>%
        spread(quantile, value)
    
    q_significant <- group_by(q, comparison) %>%
        filter(mean(p2.5 > 0 | p97.5 < 0) >= prop) %>%
        filter(p2.5 > 0 | p97.5 < 0)
    
    list(quantiles = mutate(q, k = k),
         different_cohorts = mutate(q_significant, k = k),
         comparisons = mutate(d, k = k),
         post_samples = mutate(post_samples, k))
}

default_theme <- theme_minimal() + 
    theme(strip.background = element_rect(fill = grey(0.95), size = 0),
          panel.spacing.x = grid::unit(2, "lines"),
          panel.spacing.y = grid::unit(1, "lines"))
theme_set(default_theme)

# Define UI for application that draws a histogram
ui <- fluidPage(
    titlePanel("Dirichlet-distributed multinomial variables"),
    hr(),
    sidebarLayout(
        sidebarPanel(
            radioButtons("obs", "Variables", choices = list("4 levels, pairwise similar" = "levels_similar",
                                                               "4 levels, all dissimilar" = "levels_dissimilar",
                                                               "Age, pairwise similar" = "age_similar",
                                                               "Age, all dissimilar" = "age_dissimilar",
                                                               "Custom:" = "custom")),
            textAreaInput("custom_obs", "Observed counts", value = "A = 30 65 10 30\nB = 10 20 30 50\nC = 40 60 10 25\nD = 20 55 50 95", 
                          width = "100%", height = "13ch", label = NULL),
            hr(),
            textInput("k_input", "k value(s)", value = "10 50", width = "100%"),
            sliderInput("prop", "Fraction threshold", value = 100, min = 0, max = 100, step = 5, width = "100%", post = "%", ticks = FALSE),
            sliderInput("n_samples", "Number of samples", value = 1000, min = 500, max = 5000, step = 500, width = "100%", ticks = FALSE),
            hr(),
            h4("The data"),
            tableOutput("parsed_dat"),
            width = 3
        ),
        mainPanel(
            tabsetPanel(
                tabPanel("Posterior distributions", uiOutput("post_dist")),
                tabPanel("Distributions of differences", uiOutput("post_diff")),
                tabPanel("Comparisons as numbers", 
                         br(),
                         p("The table shows the median (95% CrI) for differences between the distributions of each level."), 
                         div(dataTableOutput("diff_as_table"), style = "font-size:75%"))
            )
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    continuous_data <- function(...) {
        lapply(list(...), 
               function(x) table(cut(rnorm(100, x[2], x[3]), breaks = seq(0, 200, by = 5), labels = FALSE)) %>%
                   tibble(cohort = LETTERS[x[1]], grp = names(.), n = .)) %>%
            bind_rows() %>%
            spread(cohort, n, fill = 0) %>%
            mutate(grp = seq(0, 100, by = 5)[as.integer(grp)],
                   grp = sprintf("%i-%i", grp, grp + 4)) %>%
            column_to_rownames("grp")
    }
    cat_data <- function(...) {
        o <- str_split(..., "\n")[[1]] %>% 
            lapply(function(.) str_split(., pattern = " = ")[[1]]) %>% 
            lapply(function(.) tibble(!!sym(.[1]) := as.integer(str_split(.[2], " ")[[1]]))) %>%
            data.frame()
        row.names(o) <- paste0("level_", 1:nrow(o))
        return(o)
    }
    
    dat <- reactive({
        switch(input$obs,
               "levels_similar" = cat_data("A = 130 65 10 30\nB = 10 120 30 50\nC = 40 60 110 25\nD = 4 55 50 95"),
               "levels_dissimilar" = cat_data("A = 30 65 10 30\nB = 20 25 15 5\nC = 40 60 10 25\nD = 100 20 75 30"),
               "age_similar" = continuous_data(c(1, 70, 10), c(2, 50, 7), c(3, 65, 10), c(4, 55, 9)),
               "age_dissimilar" = continuous_data(c(1, 80, 5), c(2, 65, 5), c(3, 50, 5), c(4, 35, 5)),
               cat_data(input$custom_obs))
    })
    
    k <- reactive({
        as.numeric(str_split(input$k_input, " ")[[1]])
    })
    
    dir_diff <- reactive({
        res <- list()
        for (k in k())
          res[[paste(k)]] <- dirichlet_diff(k, dat(), input$prop / 100, input$n_samples)
        return(res)
    })
    
    output$parsed_dat <- renderTable({ addmargins(as.matrix(dat())) }, rownames = TRUE, digits = 0)
    
    # Posterior distributions plot
    post_samples <- reactive({
        bind_rows(lapply(names(dir_diff()), function(.) dir_diff()[[.]]$post_samples)) %>%
            gather(level, value, -k, -cohort)
    })
    
    output$post_dist_plot <- renderPlot({
        ggplot(post_samples(), aes(x = value, colour = level, group = level)) +
            stat_density(geom = "line", position = "identity") +
            scale_x_continuous(labels = scales::percent) +
            facet_grid(cohort ~ k) +
            labs(x = "Proportion", y = "Density",
                 caption = "Columns: different values of k; rows: cohorts") 
    })
    
    output$post_dist <- renderUI({
        plotOutput("post_dist_plot", height = n_distinct(post_samples()$cohort) * 200, width = "100%")
    })
    
    comparisons_dat <- reactive({
        bind_rows(lapply(names(dir_diff()), function(.) dir_diff()[[.]]$comparisons))
    })
    
    # Diff. between posterior distributions plot
    output$post_diff_plot <- renderPlot({
        ggplot(comparisons_dat(), aes(x = p_diff, colour = factor(level))) +
            geom_vline(aes(xintercept = 0), linetype = 2, size = 0.2) +
            stat_density(geom = "line", position = "identity") +
            facet_grid(str_replace_all(comparison, "_", " ") ~ k) +
            labs(x = "Difference", y = "Density", 
                 caption = "Columns: different values of k; rows: pairwise comparisons between cohorts") +
            guides(colour = guide_legend(title = "Level"))
    })
    
    output$post_diff <- renderUI({
        plotOutput("post_diff_plot", height = n_distinct(comparisons_dat()$comparison) * 200, width = "100%")
    })
    
    # Diff. as table
    output$diff_as_table <- renderDataTable({
        d <- dir_diff()
        all <- lapply(names(d), function(.) d[[.]]$quantiles) %>%
            bind_rows() %>%
            select(comparison, k, level)
        different <- lapply(names(d), function(.) d[[.]]$different_cohorts) %>%
            bind_rows() %>%
            transmute(k, level, median_95int = sprintf("%0.2f (%0.2f, %0.2f)", p50, p2.5, p97.5))
        left_join(all, different, by = c("comparison", "k", "level")) %>% 
            spread(comparison, median_95int) %>%
            arrange(k, level)  %>%
            rename_all(~ str_replace_all(., "_", " "))
    }, options = list(paging = FALSE, searching = FALSE))
}

# Run the application 
shinyApp(ui = ui, server = server)

