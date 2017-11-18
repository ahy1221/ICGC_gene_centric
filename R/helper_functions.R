#' Get all alteration binary matrics
#' 
#' @param gt gene centric binary table
#' @param alterations A vector of alterations
#' @param alterations.names Name vector for alterations

get_all_alterations <- function(gt, alterations, alteration.names ) {
  if(length(alterations) != length(alteration.names)) {
    ERROR = sprintf("The types of alterations [ %i ] is not equal to the alterations names [ %i ]", 
                     length(alterations), length(alteration.names) )
    stop(ERROR)
  }
  alteration.list = list()
  for( idx in seq.int(length(alterations))) {
    alter =  alterations[idx]
    obj.name = alteration.names[idx]
    INFO = sprintf("[INFO] Getting %s to %s ...\n", alter, obj.name )
    cat(INFO)
    alteration.list[[obj.name]] = get_alteration(gt, alter)
    
  }
  
  return(alteration.list)
  
}
