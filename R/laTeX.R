

#' write_latex_table
#' 
#' dumps a data frame as a LaTeX table
#'
#' @param df (data.frame)
#' @param fname (character)
#'
#' @return
#' @export
#'
#' @examples
write_latex_table <- function(df, fname, 
                              path="/Users/brianconroy/Documents/research/dataInt/output/"){
  
  df[,] <- lapply(df[, ], as.character)
  tab <- "\\begin{table}[h]
  \\begin{center}
  \\begin{tabular}{l*{5}{c}r}\n"
  header <- paste(names(df), collapse=" & ")
  tab <- paste(tab, header, "\\\\ \n \\hline \n")
  for (i in 1:nrow(df)){
    row_ <- paste(df[i,], collapse=" & ")
    row_ <- paste(row_, '\\\\\n')
    tab <- paste(tab, row_, sep=" ")
  }
  end <- "\\hline
  \\end{tabular}
  \\caption[nrgk]
  {
  nrgk
  }
  \\label{tab:nrgk}
  \\end{center}
  \\end{table}"
  tab <- paste(tab, end, sep=" ")
  path <- paste(path, fname, sep="")
  write(tab, path)

}
