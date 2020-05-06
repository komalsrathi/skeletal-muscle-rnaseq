# Author: Komal S. Rathi
# Function: Annotate and filter limma output
# Date: 04/02/2020

annotate.limma <- function(x, foldchange) {
  
  # select specific columns
  x <- x %>%
    rownames_to_column("gene_symbol") %>%
    select(gene_symbol, matches("_logFC"), P.Value, adj.P.Val, AveExpr)

  # filter by logFC
  varname1 <- paste0(gsub('_logFC','',colnames(x)[2]),'_DEGAnnot')
  varname2 <- paste0(gsub('_logFC','',colnames(x)[3]),'_DEGAnnot')
  x <- x %>%
    #filter(abs(!!as.symbol(colnames(x)[2])) > foldchange) %>%
    mutate(!!varname1 := ifelse(!!as.symbol(colnames(x)[2]) > foldchange, "Up", "Down"),
           !!varname2 := ifelse(!!as.symbol(colnames(x)[3]) > foldchange, "Up", "Down"))
  
  # return
  return(x)
}