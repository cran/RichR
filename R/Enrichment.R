#' @title Enrichment
#' @description Given a list of genes associated to diseases and a background list, it calculates the diseases enrichment. It calculates both pvalues from proportion test and Fisher's exact test. Adjusted BH pvalues are returned.
#' @param Background The background list of genes. In generela is the list of genes expressed.
#' @param Genes2Dis A data.frame with the Gene names and the Diseases. The package has two inbuild lists. g2d_clean and g2d_complete. Both lists contains data from 4 publications (GS2D, DisGeNET, Berto2016 and PsyGeNET), however in the clean similar names were treated as the same disease.
#' @param Genes A list of genes to test for enrichment
#'
#' @return a list contating the enrichment of diseases
#' @export
#'
#' @examples
#' data('g2d')
#'
#' g2d_clean = g2d$clean
#'
#' # The user can choose a particular disorder, or use the whole disease set to compare to
#'
#' g2d_ASD = subset(g2d_clean, g2d_clean$Disease %in% c('AUTISM'))
#' Enrichment(Background = g2d_clean$Gene.symbol,
#' Genes2Dis = g2d_ASD,
#' Genes = g2d_ASD$Gene.symbol[1:100])
#'
#'
#'
#' @author Deisy Morselli Gysi <deisy at bioinf.uni-leipzig.de>
#' @importFrom plyr join
#' @importFrom metap sumlog
#' @importFrom reshape2 dcast
#' @importFrom  magrittr "%>%" "%<>%"
#' @importFrom stats fisher.test na.exclude p.adjust prop.test

Enrichment =
  function(Background,
           Genes2Dis = g2d$clean,
           Genes){

    `%>%` <- magrittr::`%>%`
    `%<>%` <- magrittr::`%<>%`

    Genes2Dis = data.frame(Gene.symbol = Genes2Dis$Gene.symbol, Disease = Genes2Dis$Disease) %>% unique()
    Genes2Dis = subset(Genes2Dis, Genes2Dis$Gene.symbol %in% Background) %>% unique()
    Full = suppressMessages(unique(reshape2::dcast(formula = Gene.symbol~Disease, data = Genes2Dis )))
    Full = suppressMessages( plyr::join(data.frame(Gene.symbol = Background, BG = 'BG'), Full, type = 'full'))
    Full %<>% unique()
    CAT = ifelse(Full$Gene.symbol %in% Genes, 'Yes', 'No')
    CAT =   factor(CAT , c('No', 'Yes'))
    Full[is.na(Full)]<-'No'

    BACK = stats::aggregate(formula = Gene.symbol~Disease, data = Genes2Dis, FUN = length)
    GA = data.frame(Gene.symbol = Genes)
    Genes_a = suppressMessages( plyr::join(GA, Genes2Dis, type = 'left'))
    if(nrow(na.exclude(Genes_a)) > 0 ){
      OBS = stats::aggregate(formula = Gene.symbol~Disease, data = Genes_a, FUN = length)
    }
    else{
      OBS = BACK
      OBS$Gene.symbol = 0
    }
    names(OBS)[2] = 'OBS'

    EXP_OBS_CAT = suppressMessages( plyr::join(BACK, OBS))
    EXP_OBS_CAT[is.na(EXP_OBS_CAT)]<-0
    #########################################
    ######### Fisher's exact test & Proportion
    #########################################
    Out = matrix(NA, ncol = 3, nrow = (nrow(EXP_OBS_CAT))) %>% as.data.frame()
    names(Out) = c('Disease', 'Fisher', 'Prop')
    for( i in 1: nrow(EXP_OBS_CAT)){

      X = table(CAT, Full[,i+2]) %>% fisher.test()
      X = X$p.value
      Y = table(CAT, Full[,i+2]) %>% prop.test()
      Y = Y$p.value

      Out$Disease[i] = as.character(EXP_OBS_CAT$Disease)[i]
      Out$Fisher[i]  = X
      Out$Prop[i]    = Y
    }

    Out$FisherAdj = p.adjust(Out$Fisher, method = 'BH')
    Out$PropAdj = p.adjust(Out$Prop, , method = 'BH')

    Out$weight =   apply(Out[, 2:4], 1, FUN = function(x){ y = metap::sumlog (x); return(y$p)})

    Out = suppressMessages( plyr::join(EXP_OBS_CAT, Out))


    return(Out)
  }
