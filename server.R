
library(shiny)
library(shinyBS)
library(ggplot2)
library(scales)
library(cowplot)
library(egg)

# Reference: Daniel I. S. Rosenbloom, Julie Dudášová, Casey Davis, Radha A. Railkar, Nitin Mehrotra, Jeffrey R. Sachs (2023).
# "Replicate testing of clinical endpoints can prevent no-go decisions for beneficial vaccines"
# This program is released under the GNU GPLv3 license. Copyright © 2023 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

##########################################################################################
##### Functions for calculating effective FP rate, FN rate, and dilution in efficacy #####
##########################################################################################

# Dilution function
# p = true incidence (per sample)
# fp = false positive rate
# fn = false negative rate
dilutionFunc = function(p, fp, fn){
  oip = p*(1-fn) + (1-p)*fp
  return(fp / oip)
}
dilutionFuncApprox = function(p, fp, fn){
  fp / (p + fp)
}

###############################################################################
### Defining false positive and negative rates for various assay strategies ###
###############################################################################

######## Replicate assay strategy:
#### False positives:
fpReplicate = function(n, minN, fp){
  pbinom(minN-1, n, fp, lower.tail = FALSE) # P(positives > minN-1)
}

fpMajority = function(n, fp){
  fpReplicate(n, floor(n/2+1), fp)
}

fpReplicateApprox = function(n, minN, fp){
  choose(n, minN)*fp^minN
}

fpMajorityApprox = function(n, fp){
  fpReplicateApprox(n, floor(n/2+1), fp)
}

#### False negatives:
fnReplicate = function(n, minN, fn){
  pbinom(n-minN, n, fn, lower.tail = FALSE) # P(positives > n-minN)
}

fnMajority = function(n, fn){
  fnReplicate(n, floor(n/2+1), fn)
}

fnReplicateApprox = function(n, minN, fn){
  choose(n, n-minN+1)*fn^(n-minN+1)
}

fnMajorityApprox = function(n, fn){
  fnReplicateApprox(n, floor(n/2+1), fn)
}


######## Conditional assay strategy:
#### False positives:
fpConditional = function(n, minN, fp1, fp2){
  fp1 * fpReplicate(n-1, minN-1, fp2)
}

fpMajorityConditional = function(n, fp1, fp2){
  fpConditional(n, floor(n/2+1), fp1, fp2)
}

fpConditionalApprox = function(n, minN, fp1, fp2){
  fp1 * fpReplicateApprox(n-1, minN-1, fp2)
}

fpMajorityConditionalApprox = function(n, fp1, fp2){
  fpConditionalApprox(n, floor(n/2+1), fp1, fp2)
}

#### False negatives:
fnConditional = function(n, minN, fn1, fn2){
  fn1 + (1-fn1)*fnReplicate(n-1, minN-1, fn2)
}

fnMajorityConditional = function(n, fn1, fn2){
  fnConditional(n, floor(n/2+1), fn1, fn2)
}

fnConditionalApprox = function(n, minN, fn1, fn2){
  fn1 + (1-fn1)*fnReplicateApprox(n-1, minN-1, fn2)
}

fnMajorityConditionalApprox = function(n, fn1, fn2){
  fnConditionalApprox(n, floor(n/2+1), fn1, fn2)
}


########################
###   Shiny Server   ###
########################

shinyServer(function(input, output, session){
  
  # Set up seed and N
  reactive({set.seed(input$userSeed)})
  allN = reactive({seq(1, input$userMaxN, ifelse(input$userShowEvens, 1, 2))})
  
  # Calculate per-sample incidence (userP is incidence per year; userDurationBtwnSamples is in years)
  perSampleIncidence = reactive({1 - (1 - input$userP)^input$userDurationBtwnSamples})
  
  ##### First show calculations without uncertainty #####
  
  dtCalc = reactive({
    dtCalc_Setup = CJ(
      p = perSampleIncidence(),
      fp = input$userFP,
      fpi = input$userFPi,
      fn = input$userFN,
      fni = input$userFNi,
      n = allN(),
      isConditional = c(0,1,2)
    )
    dtCalc_Setup[, netFP :=
                   ifelse(
                     isConditional == 0,
                     fpMajority(n, fp),
                     ifelse(
                       isConditional == 1,
                       fpMajorityConditional(n, fp, fp),
                       ifelse(
                         isConditional == 2,
                         fpMajorityConditional(n, fpi, fp),
                         NA)))]
    
    dtCalc_Setup[, netFN :=
                   ifelse(
                     isConditional == 0,
                     fnMajority(n, fn),
                     ifelse(
                       isConditional == 1,
                       fnMajorityConditional(n, fn, fn),
                       ifelse(
                         isConditional == 2,
                         fnMajorityConditional(n, fni, fn),
                         NA)))]
    
    dtCalc_Setup[, dil := dilutionFunc(p, netFP, netFN)]
    dtCalc_Setup[, eff := input$userEff]
    dtCalc_Setup[, netEff := eff * (1-dil)]
    
    dtCalc_Setup[, isConditional := factor(isConditional,
                                           labels = c(
                                             'Majority rule\n',
                                             'Confirmatory\nmajority rule\n',
                                             'Confirmatory\nmajority rule,\ndifferent initial\nassay\n'
                                           ))]
    dtCalc_Setup
  })
  
  
  
  ##### Next show calculations with uncertainty #####
  
  dtCalcUncertainty = reactive({
    if(input$userSHOWUNCERTAINTY){
      set.seed(input$userSeed)
      dtCalcUncertainty_Setup = CJ(
        sampleI = 1:input$userNumSamples,
        n = allN(),
        isConditional = c(0,1,2)
      )
      
      sampleDT = data.table(
        sampleI = 1:input$userNumSamples,
        pAnnual = rbeta(input$userNumSamples, 1 + input$userP*input$userBetaN_P, 1 + (1-input$userP)*input$userBetaN_P),
        fp = rbeta(input$userNumSamples, 1 + input$userFP*input$userBetaN_FP, 1 + (1-input$userFP)*input$userBetaN_FP),
        fn = rbeta(input$userNumSamples, 1 + input$userFN*input$userBetaN_FN, 1 + (1-input$userFN)*input$userBetaN_FN),
        fpi = rbeta(input$userNumSamples, 1 + input$userFPi*input$userBetaN_FPi, 1 + (1-input$userFPi)*input$userBetaN_FPi),
        fni = rbeta(input$userNumSamples, 1 + input$userFNi*input$userBetaN_FNi, 1 + (1-input$userFNi)*input$userBetaN_FNi)
      )
      sampleDT[, p := 1 - (1 - pAnnual)^input$userDurationBtwnSamples] # Convert each draw of pAnnual to per-sample time basis
      
      dtCalcUncertainty_Setup = merge(dtCalcUncertainty_Setup, sampleDT, by='sampleI')
      
      dtCalcUncertainty_Setup[, netFP :=
                                ifelse(
                                  isConditional == 0,
                                  fpMajority(n, fp),
                                  ifelse(
                                    isConditional == 1,
                                    fpMajorityConditional(n, fp, fp),
                                    ifelse(
                                      isConditional == 2,
                                      fpMajorityConditional(n, fpi, fp),
                                      NA)))]
      
      dtCalcUncertainty_Setup[, netFN :=
                                ifelse(
                                  isConditional == 0,
                                  fnMajority(n, fn),
                                  ifelse(
                                    isConditional == 1,
                                    fnMajorityConditional(n, fn, fn),
                                    ifelse(
                                      isConditional == 2,
                                      fnMajorityConditional(n, fni, fn),
                                      NA)))]
      
      dtCalcUncertainty_Setup[, dil := dilutionFunc(p, netFP, netFN)]
      dtCalcUncertainty_Setup[, eff := input$userEff]
      dtCalcUncertainty_Setup[, netEff := eff * (1-dil)]
      
      dtCalcUncertainty_Setup[, isConditional := factor(isConditional,
                                                        labels = c(
                                                          'Majority rule\n',
                                                          'Confirmatory\nmajority rule\n',
                                                          'Confirmatory\nmajority rule,\ndifferent initial\nassay\n'
                                                        ))]
    }else{
      # if input$userSHOWUNCERTAINTY == FALSE
      dtCalcUncertainty_Setup = NULL
    }
    dtCalcUncertainty_Setup
  }) # dtCalcUncertainty = reactive({
  
  # Render plot without uncertainty:
  output$plot1 = renderPlot({
    FPGraph = ggplot(data = dtCalc(),
                     aes(x = n, group = isConditional, color = isConditional, shape = isConditional)
    ) +
      geom_line(aes(y=netFP), size=2, alpha=0.7) +
      geom_point(aes(y=netFP), size=4) +
      theme_bw() +
      xlab(NULL) +
      ylab('Overall false\npositive rate') +
      theme(text = element_text(size=15), axis.title = element_text(size=18), axis.title.y = element_text(angle=0, vjust=0.5),
            legend.position = "none", axis.text.x = element_blank()) +
      labs(color = 'Strategy', shape = 'Strategy') +
      scale_color_brewer(palette = 'Set1') +
      scale_x_continuous(breaks = 1:20, minor_breaks = NULL) +
      scale_y_log10()
    
    FNGraphWithLegend = ggplot(data = dtCalc(),
                               aes(x = n, group = isConditional, color = isConditional, shape = isConditional)
    ) +
      geom_line(aes(y=netFN), size=2, alpha=0.7) +
      geom_point(aes(y=netFN), size=4) +
      theme_bw() +
      xlab(NULL) +
      ylab('Overall false\nnegative rate') +
      theme(text = element_text(size=15), axis.title = element_text(size=18), axis.title.y = element_text(angle=0, vjust=0.5),
            axis.text.x = element_blank()) +
      labs(color = 'Strategy', shape = 'Strategy') +
      scale_color_brewer(palette = 'Set1') +
      scale_x_continuous(breaks = 1:20, minor_breaks = NULL) +
      scale_y_continuous()
    
    FNGraph = FNGraphWithLegend + theme(legend.position = 'none')
    
    EfficacyGraph = ggplot(data = dtCalc(),
                           aes(x = n, group = isConditional, color = isConditional, shape = isConditional)
    ) +
      geom_line(aes(y=netEff), size=2, alpha=0.7) +
      geom_point(aes(y=netEff), size=4) +
      theme_bw() +
      xlab('# Replicate assays') +
      ylab('Expected\nobserved\nefficacy') +
      theme(text = element_text(size=15), axis.title = element_text(size=18), axis.title.y = element_text(angle=0, vjust=0.5),
            legend.position = "none") +
      labs(color = 'Strategy', shape = 'Strategy') +
      scale_color_brewer(palette = 'Set1') +
      scale_x_continuous(breaks = 1:20, minor_breaks = NULL) +
      scale_y_continuous(labels = scales::percent_format(1), limits = c(NA, input$userEff))
    
    
    grobs = ggplotGrob(FNGraphWithLegend)$grobs
    mylegend = grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
    
    # build grid without legends
    pgrid = plot_grid(FPGraph, FNGraph, EfficacyGraph, ncol = 1, align = 'v')
    # add legend
    plot_grid(pgrid, mylegend, ncol = 2, rel_widths = c(0.8, .2))
  })
  
  # Render plot with uncertainty:
  output$plot2 = renderPlot({
    if(input$userSHOWUNCERTAINTY){
      dtCalcUncertaintySummary = dtCalcUncertainty()[,
                                                     list(
                                                       nSamples = .N,
                                                       
                                                       netFP_median = median(netFP),
                                                       netFP_hi = quantile(netFP, 1 - input$userAlpha/2),
                                                       netFP_lo = quantile(netFP, input$userAlpha/2),
                                                       
                                                       netFN_median = median(netFN),
                                                       netFN_hi = quantile(netFN, 1 - input$userAlpha/2),
                                                       netFN_lo = quantile(netFN, input$userAlpha/2),
                                                       
                                                       netEff_median = median(netEff),
                                                       netEff_hi = quantile(netEff, 1 - input$userAlpha/2),
                                                       netEff_lo = quantile(netEff, input$userAlpha/2)
                                                     ),
                                                     by = list(isConditional, n)]
      
      # Remove ending newline from strip text:
      dtCalcUncertaintySummary[, isConditionalStripText := factor(isConditional,
                                                                  labels = c(
                                                                    'Majority rule',
                                                                    'Confirmatory\nmajority rule',
                                                                    'Confirmatory\nmajority rule,\ndifferent initial\nassay'
                                                                  ))]
      
      UncertaintyFPGraph = ggplot(data = dtCalcUncertaintySummary,
                                  aes(x = n, group = isConditional, color = isConditional, fill = isConditional, shape = isConditional)
      ) +
        geom_ribbon(aes(ymin=netFP_lo, ymax=netFP_hi), alpha=0.4, color=NA) +
        geom_line(aes(y=netFP_median), size=2, alpha=0.7) +
        geom_point(aes(y=netFP_median), size=4) +
        theme_bw() +
        xlab(NULL) +
        ylab('Overall false\npositive rate') +
        theme(text = element_text(size=15), axis.title = element_text(size=16), axis.title.y = element_text(angle=0, vjust=0.5),
              legend.position = "none", axis.text.x = element_blank()) +
        labs(color = 'Strategy', shape = 'Strategy', fill = 'Strategy') +
        scale_color_brewer(palette = 'Set1') +
        scale_fill_brewer(palette = 'Set1') +
        scale_x_continuous(breaks = 1:20, minor_breaks = NULL) +
        scale_y_log10() +
        facet_grid(cols = vars(isConditionalStripText))
      
      UncertaintyFNGraphWithLegend = ggplot(data = dtCalcUncertaintySummary,
                                            aes(x = n, group = isConditional, color = isConditional, fill = isConditional, shape = isConditional)
      ) +
        geom_ribbon(aes(ymin=netFN_lo, ymax=netFN_hi), alpha=0.4, color=NA) +
        geom_line(aes(y=netFN_median), size=2, alpha=0.7) +
        geom_point(aes(y=netFN_median), size=4) +
        theme_bw() +
        xlab(NULL) +
        ylab('Overall false\nnegative rate') +
        theme(text = element_text(size=15), axis.title = element_text(size=16), axis.title.y = element_text(angle=0, vjust=0.5),
              axis.text.x = element_blank(), strip.background = element_blank(), strip.text = element_blank()) +
        labs(color = 'Strategy', shape = 'Strategy', fill = 'Strategy') +
        scale_color_brewer(palette = 'Set1') +
        scale_fill_brewer(palette = 'Set1') +
        scale_x_continuous(breaks = 1:20, minor_breaks = NULL) +
        scale_y_continuous() +
        facet_grid(cols = vars(isConditionalStripText))
      
      UncertaintyFNGraph = UncertaintyFNGraphWithLegend + theme(legend.position = 'none')
      
      UncertaintyEfficacyGraph = ggplot(data = dtCalcUncertaintySummary,
                                        aes(x = n, group = isConditional, color = isConditional, fill = isConditional, shape = isConditional)
      ) +
        geom_ribbon(aes(ymin=netEff_lo, ymax=netEff_hi), alpha=0.4, color=NA) +
        geom_line(aes(y=netEff_median), size=2, alpha=0.7) +
        geom_point(aes(y=netEff_median), size=4) +
        theme_bw() +
        xlab('# Replicate assays') +
        ylab('Expected\nobserved\nefficacy') +
        theme(text = element_text(size=15), axis.title = element_text(size=16), axis.title.y = element_text(angle=0, vjust=0.5),
              legend.position = "none", strip.background = element_blank(), strip.text = element_blank()) +
        labs(color = 'Strategy', shape = 'Strategy', fill = 'Strategy') +
        scale_color_brewer(palette = 'Set1') +
        scale_fill_brewer(palette = 'Set1') +
        scale_x_continuous(breaks = 1:20, minor_breaks = NULL) +
        scale_y_continuous(labels = scales::percent_format(1), limits = c(NA, input$userEff)) +
        facet_grid(cols = vars(isConditionalStripText))
      
      
      grobsUncertainty = ggplotGrob(UncertaintyFNGraphWithLegend)$grobs
      mylegendUncertainty = grobsUncertainty[[which(sapply(grobsUncertainty, function(x) x$name) == "guide-box")]]
      
      # build grid without legends
      pgridUncertainty = plot_grid(UncertaintyFPGraph, UncertaintyFNGraph, UncertaintyEfficacyGraph,
                                   ncol = 1, align = 'v', rel_heights = c(1.35, 0.75, 1.1))
      # add legend
      plot_grid(pgridUncertainty, mylegendUncertainty, ncol = 2, rel_widths = c(0.85, .15))
    }else{
      # if input$userSHOWUNCERTAINTY == FALSE
      NULL
    }
  })
  
  # Generate data for download
  dtCalcForDownload = reactive({
    if(input$userSHOWUNCERTAINTY){
      thisDT = dtCalcUncertainty()
      # Remove '\n' and stray spaces from sampling strategy, fix up order:
      thisDT[, isConditional := gsub('\n', ' ', isConditional)]
      thisDT[isConditional == 'Majority rule ', isConditional := 'Majority rule']
      thisDT[isConditional == 'Confirmatory majority rule ', isConditional := 'Confirmatory majority rule']
      thisDT[isConditional == 'Confirmatory majority rule, different initial assay ', isConditional := 'Confirmatory majority rule, different initial assay']
      thisDT[, isConditional := factor(isConditional,
                                       levels = c(
                                         'Majority rule',
                                         'Confirmatory majority rule',
                                         'Confirmatory majority rule, different initial assay'
                                       ))]
      # Rename columns:
      names(thisDT)[[1]] = 'Simulation #'
      names(thisDT)[[2]] = '# Replicate Assays'
      names(thisDT)[[3]] = 'Sampling Strategy'
      names(thisDT)[[4]] = 'True incidence, per year'
      names(thisDT)[[5]] = 'FP rate (main or confirmatory assay)'
      names(thisDT)[[6]] = 'FN rate (main or confirmatory assay)'
      names(thisDT)[[7]] = 'FP rate (initial assay)'
      names(thisDT)[[8]] = 'FN rate (initial assay)'
      names(thisDT)[[9]] = 'TIP (true infected positives, per sampling occasion)'
      names(thisDT)[[10]] = 'Effective FP rate'
      names(thisDT)[[11]] = 'Effective FN rate'
      names(thisDT)[[12]] = 'Efficacy dilution'
      names(thisDT)[[13]] = 'True efficacy'
      names(thisDT)[[14]] = 'Observed (diluted) efficacy'
      # Create additional column:
      thisDT[, `Time between sampling occasions` := input$userDurationBtwnSamples]
      # Reorder:
      setcolorder(thisDT,
                  c(
                    "Sampling Strategy",
                    "# Replicate Assays",
                    "Simulation #",
                    "True incidence, per year",
                    "Time between sampling occasions",
                    "TIP (true infected positives, per sampling occasion)",
                    "FP rate (main or confirmatory assay)",
                    "FN rate (main or confirmatory assay)",
                    "FP rate (initial assay)",
                    "FN rate (initial assay)",
                    "Effective FP rate",
                    "Effective FN rate",
                    "Efficacy dilution",
                    "True efficacy",
                    "Observed (diluted) efficacy"
                  ))
      setorderv(thisDT,
                c(
                  "Sampling Strategy",
                  "# Replicate Assays",
                  "Simulation #"
                ))
    }else{
      # input$userSHOWUNCERTAINTY == FALSE
      thisDT = dtCalc()
      # Remove '\n' and stray spaces from sampling strategy, fix up order:
      thisDT[, isConditional := gsub('\n', ' ', isConditional)]
      thisDT[isConditional == 'Majority rule ', isConditional := 'Majority rule']
      thisDT[isConditional == 'Confirmatory majority rule ', isConditional := 'Confirmatory majority rule']
      thisDT[isConditional == 'Confirmatory majority rule, different initial assay ', isConditional := 'Confirmatory majority rule, different initial assay']
      thisDT[, isConditional := factor(isConditional,
                                       levels = c(
                                         'Majority rule',
                                         'Confirmatory majority rule',
                                         'Confirmatory majority rule, different initial assay'
                                       ))]
      # Rename columns:
      names(thisDT)[[1]] = 'TIP (true infected positives, per sampling occasion)'
      names(thisDT)[[2]] = 'FP rate (main or confirmatory assay)'
      names(thisDT)[[3]] = 'FP rate (initial assay)'
      names(thisDT)[[4]] = 'FN rate (main or confirmatory assay)'
      names(thisDT)[[5]] = 'FN rate (initial assay)'
      names(thisDT)[[6]] = '# Replicate Assays'
      names(thisDT)[[7]] = 'Sampling Strategy'
      names(thisDT)[[8]] = 'Effective FP rate'
      names(thisDT)[[9]] = 'Effective FN rate'
      names(thisDT)[[10]] = 'Efficacy dilution'
      names(thisDT)[[11]] = 'True efficacy'
      names(thisDT)[[12]] = 'Observed (diluted) efficacy'
      # Create additional column:
      thisDT[, `Time between sampling occasions` := input$userDurationBtwnSamples]
      thisDT[, `True incidence, per year` := input$userP]
      # Reorder:
      setcolorder(thisDT,
                  c(
                    "Sampling Strategy",
                    "# Replicate Assays",
                    "True incidence, per year",
                    "Time between sampling occasions",
                    "TIP (true infected positives, per sampling occasion)",
                    "FP rate (main or confirmatory assay)",
                    "FN rate (main or confirmatory assay)",
                    "FP rate (initial assay)",
                    "FN rate (initial assay)",
                    "Effective FP rate",
                    "Effective FN rate",
                    "Efficacy dilution",
                    "True efficacy",
                    "Observed (diluted) efficacy"
                  ))
      setorderv(thisDT,
                c(
                  "Sampling Strategy",
                  "# Replicate Assays"
                ))
    }
  })
  
  output$downloadDataNoUnc = downloadHandler(
    filename = 'VEDilutionNoUncertainty.txt',
    content = function(file){write.table(dtCalcForDownload(), file, row.names=F, sep='\t', quote=F)}
  )
  output$downloadDataWithUnc = downloadHandler(
    filename = 'VEDilutionWithUncertainty.txt',
    content = function(file){write.table(dtCalcForDownload(), file, row.names=F, sep='\t', quote=F)}
  )
  
}) # shinyServer


