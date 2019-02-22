#' @title g2d
#' @docType data
#' @description g2d is a list with two data.frames based on the following  4 publications
#' (GS2D, Fontaine (2016) <doi:10.18547/gcb.2016.vol2.iss1.e33>,
#' DisGeNET, Pinero (2016) <doi:10.1093/nar/gkw943>,
#' Berto2016, Berto (2016) <doi:10.3389/fgene.2016.00031> and
#' PsyGeNET, Sacristan (2015) <doi:10.1093/bioinformatics/btv301>).
#' Those lists were combined and manually curated to have matching disease names.
#' The first list, clean, contains the curated data, the list complete contains complete data. In the former disease names might not match.
"g2d"
utils::globalVariables(names = 'g2d')
