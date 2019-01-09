

library(shinyjs)
library(shiny)
library(OmicsEV)


options(shiny.maxRequestSize = 1000 * 1024 ^ 2)


shinyServer(function(input, output) {
    vdata <- reactiveValues(data = list(result = NULL))

    observeEvent(input$start, {
        showNotification(
            paste(date(), " -> Your task is running ...", sep = ""),
            type = "message",
            duration = NULL,
            closeButton = TRUE
        )

        ## input parameters
        vdata$data$dataType <- input$dataType
        vdata$data$data_dir <- input$data_dir
        vdata$data$sample_list <- input$sample_list
        vdata$data$x2 <- input$x2
        vdata$data$missing_value_cutoff <- input$missing_value_cutoff
        vdata$data$cpu <- input$cpu

        if(!is.null(input$use_class)){
            vdata$data$use_class <- unlist(strsplit(x = input$use_class,split = ","))
        }else{
            vdata$data$use_class <- NULL
        }

        vdata$data$ml_class <- input$ml_class
        vdata$data$out_dir <- input$out_dir
        vdata$data$class_color <- input$class_color

        tryCatch({
            html_report <- run_omics_evaluation(
                data_dir = vdata$data$data_dir,
                sample_list = vdata$data$sample_list,
                x2 = vdata$data$x2,
                class_color = vdata$data$class_color,
                use_class = vdata$data$use_class,
                cpu = vdata$data$cpu,
                missing_value_cutoff = vdata$data$missing_value_cutoff,
                ml_class = vdata$data$ml_class,
                out_dir = vdata$data$out_dir
            )

            if (file.exists(html_report)) {
                vdata$data$result <- html_report
            } else{
                vdata$data$result <- NULL

            }
            showNotification(
                paste(
                    date(),
                    " -> Congratulations. Your task is finished! Now you can go to te result page to see the result.",
                    sep = ""
                ),
                type = "message",
                duration = NULL,
                closeButton = TRUE
            )

        },
        error = function(e) {
            vdata$data$result <- NULL
            showNotification(
                paste(
                    date(),
                    " -> Sorry! Your task is finished with error!",
                    sep = ""
                ),
                type = "error",
                duration = NULL,
                closeButton = TRUE
            )
            showModal(
                modalDialog(
                    title = "Error! Please check the input data and parameters! You can open an issue at github for help (https://github.com/bzhanglab/OmicsEV/issues).",
                    sprintf("caught Error: %s", e),
                    easyClose = TRUE,
                    footer = NULL,
                    fade = TRUE
                )
            )

        })
    })

    observeEvent(input$stop, {
        vdata$data = list(result = NULL)
    })


    output$bodyUI <- renderUI({
        if (is.null(vdata$data$result)) {
            br()
            h5("Welcome!")

        } else{
            p(
                "The evaluation analysis is done. Here is the HTML-based report:",
                vdata$data$result
            )

        }

    })

})
