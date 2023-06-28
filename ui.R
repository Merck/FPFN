
library(shiny)
library(shinythemes)
library(data.table)
library(ggplot2)
library(scales)
library(cowplot)
library(egg)

# Reference: Daniel I. S. Rosenbloom, Julie Dudášová, Casey Davis, Radha A. Railkar, Nitin Mehrotra, Jeffrey R. Sachs (2023).
# "Replicate testing of clinical endpoints can prevent no-go decisions for beneficial vaccines"
# This program is released under the GNU GPLv3 license. Copyright © 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

shinyUI(navbarPage(fluid = TRUE, theme = shinytheme("simplex"),
                   title = "Replicate Assay Testing Calculator",
                   fluidRow(column(12, 'Reference: Daniel I. S. Rosenbloom, Julie Dudášová, Casey Davis, Radha A. Railkar, Nitin Mehrotra, Jeffrey R. Sachs (2023).
                                     "Replicate testing of clinical endpoints can prevent no-go decisions for beneficial vaccines"')),
                   fluidRow(br()),
                   sidebarLayout(
                     sidebarPanel(
                       h3("Settings"),
                       numericInput("userEff", "True vaccine efficacy", 0.8, min = 0, max = 1, step = 0.05),
                       wellPanel(
                         fluidRow(
                           column(6, numericInput("userP", "True incidence in control arm, per year", 0.02, min = 0, max = 1, step = 0.01)),
                           column(6, numericInput("userDurationBtwnSamples", "Duration btwn samples, years (set to 1 if considering only a single sample)", 1, min = 0, max = 1e5, step = 0.25))
                         )
                       ),
                       wellPanel(
                         fluidRow(
                           column(6, numericInput("userFP", "Main / confirmatory assay FP rate", 0.03, min = 0, max = 1, step = 0.01)),
                           column(6, numericInput("userFN", "Main / confirmatory assay FN rate", 0.2, min = 0, max = 1, step = 0.01))
                         )
                       ),
                       wellPanel(
                         fluidRow(
                           column(6, numericInput("userFPi", "Initial assay FP rate", 0.12, min = 0, max = 1, step = 0.01)),
                           column(6, numericInput("userFNi", "Initial assay FN rate", 0.05, min = 0, max = 1, step = 0.01))
                         )
                       ),
                       sliderInput("userMaxN", "Maximum number of replicates to consider", 5, min = 3, max = 11, step = 1),
                       checkboxInput("userShowEvens", "Consider an even number of replicates?", FALSE),
                       width = 4
                     ), # sidebarPanel
                     
                     mainPanel(
                       fluidRow(title = 'Output',
                                h3('Output'),
                                conditionalPanel(
                                  condition = "input.userSHOWUNCERTAINTY == false",
                                  wellPanel(plotOutput("plot1")),
                                ),
                                conditionalPanel(
                                  condition = "input.userSHOWUNCERTAINTY == true",
                                  wellPanel(plotOutput("plot2")),
                                ),
                                wellPanel(
                                  conditionalPanel(
                                    condition = "input.userSHOWUNCERTAINTY == false",
                                    downloadButton("downloadDataNoUnc", "Download Output")
                                  ),
                                  conditionalPanel(
                                    condition = "input.userSHOWUNCERTAINTY == true",
                                    downloadButton("downloadDataWithUnc", "Download Output")
                                  )
                                ),
                                checkboxInput("userSHOWUNCERTAINTY", "Consider uncertainty in parameters?", FALSE),
                                conditionalPanel(condition = "input.userSHOWUNCERTAINTY == true",
                                                 title = 'Uncertainty settings',
                                                 wellPanel(
                                                   h3('Uncertainty settings'),
                                                   wellPanel(fluidRow(
                                                     column(2, h4('Main / confirmatory assay uncertainty')),
                                                     column(5, numericInput("userBetaN_FP", "Number of controls from which main assay FP rate was estimated (# neg. control observations)", 200, min = 0, step = 50)),
                                                     column(5, numericInput("userBetaN_FN", "Number of controls from which main assay FN rate was estimated (# pos. control observations)", 200, min = 0, step = 50))
                                                   )),
                                                   wellPanel(fluidRow(
                                                     column(2, h4('Initial assay uncertainty')),
                                                     column(5, numericInput("userBetaN_FPi", "Number of controls from which init. assay FP rate was estimated (# neg. control observations)", 200, min = 0, step = 50)),
                                                     column(5, numericInput("userBetaN_FNi", "Number of controls from which init. assay FN rate was estimated (# pos. control observations)", 200, min = 0, step = 50))
                                                   )),
                                                   wellPanel(fluidRow(
                                                     column(2, h4('Infection incidence uncertainty')),
                                                     column(5, numericInput("userBetaN_P", "Effective N from which incidence was estimated (e.g., # person-years)", 1000, min = 0, step = 50))
                                                   )),
                                                   wellPanel(fluidRow(
                                                     column(2, h4('Other settings')),
                                                     column(3, numericInput("userAlpha", "Alpha for CIs", 0.05, min = 0, max = 1, step = 0.01)),
                                                     column(3, numericInput("userNumSamples", "Number of samples to draw from uncertainty distributions", 2000, min = 500, max = 1e5, step = 500)),
                                                     column(3, numericInput("userSeed", "Random seed", 1, min = 1, max = 1e8, step = 1))
                                                   ))
                                                 )
                                ) # conditionalPanel
                       ), # fluidRow
                     ) # mainPanel
                   ), # sidebarLayout
                   fluidRow(column(12, 'This program is released under the GNU GPLv3 license. Copyright © 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.')),
                   fluidRow(br())
))
