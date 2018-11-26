
# df <- data.frame(matrix(c(1,1,1,1), nrow=2))
# names(df) <- c("ckya", "blyat")
# fname = "blyat.txt"


write_latex_table <- function(df, fname){
  
  tab <- "\\begin{table}[h]
  \\begin{center}
  \\begin{tabular}{l*{5}{c}r}\n"
  header <- paste(names(df), collapse=" & ")
  tab <- paste(tab, header, "\\\\ \n \\hline \n")
  for (row in df){
    row_ <- paste(row, collapse=" & ")
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
  path <- paste("/Users/brianconroy/Documents/research/dataInt/output/", fname, sep="")
  write(tab, path)

}
