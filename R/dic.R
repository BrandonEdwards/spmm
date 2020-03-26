#' Deviance Information Criterion
#'
#' \code{dic} returns the deviance information criterion for the model
#'   and allows you to specify which DIC to return (based on the 8
#'   DIC equations for mixture models proposed by Celeux et al. 2006).
#'
#' @export
#'

dic <- function(mod = NULL,
                form = NULL)
{
  switch(form,
         "1" = {
           print("DIC 1")
         },
         "2" = {
           print("DIC 2")
         },
         "3" = {
           print("DIC 3")
         },
         "4" = {
           print("DIC 4")
         },
         "5" = {
           print("DIC 5")
         },
         "6" = {
           print("DIC 6")
         },
         "7" = {
           print("DIC 7")
         },
         "8" = {
           print("DIC 8")
         }
  )
}
