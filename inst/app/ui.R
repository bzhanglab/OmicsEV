
library(shinyjs)
library("shiny")
library(shinydashboard)
library(shinythemes)
library(shinysky)
library(bsplus)

library(parallel)
ncpus <- detectCores()

actionButton <-
    function(inputId,
             label,
             btn.style = "" ,
             css.class = "") {
        if (btn.style %in% c("primary",
                             "info",
                             "success",
                             "warning",
                             "danger",
                             "inverse",
                             "link")) {
            btn.css.class <- paste("btn", btn.style, sep = "-")
        } else
            btn.css.class = ""

        tags$button(
            id = inputId,
            type = "button",
            class = paste("btn action-button", btn.css.class, css.class, collapse =
                              " "),
            label
        )
    }

shinyUI(dashboardPage(
    skin = "green",

    dashboardHeader(
        title = "OmicsEV",
        titleWidth = 300,
        dropdownMenuOutput("taskMenu")
    ),

    dashboardSidebar(
        width = 300,
        selectInput(
            "dataType",
            width = 300,
            label = "Data type:",
            choices = c("protein" = "protein",
                        "gene" = "gene"),
            selected = "protein",
            multiple = FALSE
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "The data type for input datasets: protein or gene", placement = "left")
            ),

        textAreaInput(
            inputId = "data_dir",
            label = "Input datasets",
            resize = "vertical"
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "A folder contains multiple omics datasets in tsv format. Row: genes/proteins, column: samples", placement = "left")
            ),

        textAreaInput(
            inputId = "sample_list",
            label = "Sample list file",
            resize = "vertical"
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "A file in tsv format contains sample class, batch and order information", placement = "left")
            ),

        textAreaInput(
            inputId = "x2",
            label = "x2",
            resize = "vertical"
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "A tsv format file which contains protein or gene expression data. When the value of parameter is not NULL, the data in this file will be used to perform correlation analysis between this dataset and the dataset in data_dir.", placement = "left")
            ),

        textAreaInput(
            inputId = "class_color",
            label = "Class color:",
            resize = "vertical"
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "A tsv format file which contains color for each class", placement = "left")
            ),


        numericInput(
            "missing_value_cutoff",
            "Missing value cutoff:",
            0.5,
            min = 0,
            max = 1,
            step = 0.1
        ),

        numericInput(
            "cpu",
            "CPUs:",
            ncpus,
            min = 0,
            max = ncpus,
            step = 1
        ),

        textAreaInput(
            inputId = "use_class",
            label = "Class for correlation:",
            resize = "vertical"
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "The class of samples which will be used to perform correlation analysis", placement = "left")
            ),

        textAreaInput(
            inputId = "ml_class",
            label = "Class for prediction:",
            resize = "vertical"
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "The class of samples which will be used for phenotype prediction", placement = "left")
            ),

        textAreaInput(
            inputId = "out_dir",
            label = "Output folder",
            resize = "vertical"
        ) %>%
            shinyInput_label_embed(
                shiny_iconlink() %>%
                    bs_embed_tooltip(title = "Output folder.", placement = "left")
            ),

        splitLayout(
            cellArgs = list(style = "padding: 6px"),
            actionButton("start", "Start", "info"),
            actionButton("stop", "Stop", "info")
        ),
        br()

    ),


    dashboardBody(
        busyIndicator(text = "It's running ...",
                      img = "shinysky/busyIndicator/ajaxloaderq.gif", wait =
                          1000),
        fluidRow(box(
            width = NULL,
            height = NULL,
            solidHeader = TRUE,
            uiOutput("bodyUI")
        ))
    )
))
