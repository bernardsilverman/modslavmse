#' UK data
#'
#' 
#'  Data from the UK 2013 strategic assessment
#' 
#' This is a table of six lists used in the research \href{https://tinyurl.com/ydegfjaw}{published by the Home Office} as part of the strategy leading
#'  to the Modern Slavery Act 2015. The data are considered in six lists, labelled as follows:  LA--Local authorities; NG--Non-government organisations;
#' PF--Police forces; GO--Government organisations; GP--General public; NCA--National Crime Agency.  Each of the first six columns in the data frame
#'  corresponds to one of these lists.
#'   Each of the rows of the data frame corresponds to a possible combination of lists, with value 1 in the relevant
#'  column if the list is in that particular combination. 
#'  The last column of the data frame states the number of cases observed in that particular combination of lists. 
#'  Combinations of lists for which zero cases are observed are omitted. 
#'
#' @references \url{https://www.gov.uk/government/publications/modern-slavery-an-application-of-multiple-systems-estimation}
#' 
"UKdat"
#' UK data five list version
#'
#' UK data consolidated into five lists
#'
#' This reduces the UK data \code{\link{UKdat}} into five lists, constructed by combining
#'  the PF and NCA lists into a single PFNCA list
#'
#' 
"UKdat_5"
#' UK data four list version
#'
#' UK data consolidated into four lists
#'
#' This reduces the UK data \code{\link{UKdat}} into four lists, constructed by combining
#'  the PF and NCA lists into a single PFNCA list and by omitting the GP list
#'
#'
"UKdat_4"
#' Netherlands data
#' 
#'  Identified victims in The Netherlands, 2010-15
#' 
#' This is a table corresponding to six lists of victims in The Netherlands. The data are considered in six lists, labelled as follows: P = National Police; K = Border Police;
#' I = Inspectorate SZW (Ministry of Social Aairs and Employment); 
#' R = regional coordinators; O = residential treatment centers and shelters; Z = others
#' (for example, ambulatory care centers, organizations providing legal services,
#' Immigration and Naturalization Service). Each of the first six columns in the data frame
#'  corresponds to one of these lists.
#'   Each of the rows of the data frame corresponds to a possible combination of lists, with value 1 in the relevant
#'  column if the list is in that particular combination. 
#'  The last column of the data frame states the number of cases observed in that particular combination of lists. 
#'  Combinations of lists for which zero cases are observed are omitted. 
#'
#' @references J. J. van Dijk, M. Cruy, P. G. M. van der Heijden, and S. L. J. Kragten-Heerdink (2017). Monitoring Target 16.2 of the United Nations Sustainable Development Goals; a multiple systems estimation of the numbers of presumed human tracking victims in the Netherlands in 2010-2015 by year, age, gender, form of exploitation and nationality. Research Brief, available at \url{https://tinyurl.com/y9mpkach}.
#' 
"Ned"
#' Netherlands data five list version
#'
#' Netherlands data consolidated into five lists
#'
#' This reduces the Netherlands data \code{\link{Ned}} into five lists, constructed by combining
#'  the two smallest lists I and O into a single list.
#'
"Ned_5"
#' Kosovo data
#'
#' Data on 4400 observed killings in the Kosovo war between 20 March and 22 June 1999
#'
#' These data (processed from \code{\link[LCMCR]{kosovo_aggregate}}) give the numbers of cases on each possible combination of four lists. 
#' The lists are labelled as follows: EXH = exhumations; ABA = American Bar Association Central and East European Law Initiative; OSCE = Organization for Security and Cooperation in Europe; HRW = Human Rights Watch.
#' All 15 combinations have a nonzero count.
#' @references Ball, P., W. Betts, F. Scheuren, J. Dudukovich, and J. Asher (2002). Killings and Refugee Flow in Kosovo March-June 1999. American Association for the Advancement of Science. A Report to the International Criminal Tribunal for the Former Yugoslavia.
"Kosovo"
#' New Orleans data
#'
#' Victims related to modern slavery and human trafficking in New Orleans
#'
#' These data are collected into 8 lists.  For reasons of confidentiality the labels of the lists are anonymised.  
#' 
#' @references K. Bales, L. Murphy and B. W. Silverman (2018). How many trafficked and enslaved people are there in New Orleans? Unpublished report.
"NewOrl"
#' New Orleans data five list version
#'
#' New Orleans data consolidated into five lists
#'
#' This reduces the New Orleans data \code{\link{NewOrl}} into five lists, constructed by combining
#'  the four smallest lists B, E, F and G into a single list.
#'
#' 
"NewOrl_5"