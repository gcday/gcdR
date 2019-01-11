#' Renames human genes to mouse 
#'
#'
#' @param genes list of human genes
#' 
#' @return list containing mouse homologs of human genes
#'
#' @examples
#' convertHumanGeneList(genes)
#'
#' @export
convertHumanGeneList <- function(genes){
  require("biomaRt")
  human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  # human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
 
  mouse = useMart("ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = genes, 
                   mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, 
                   uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}

#' Renames human genes to mouse 
#'
#'
#' @param pathways.list list of human genes
#' 
#' @return list containing mouse homologs of human genes
#'
#' @examples
#' convertHumanGeneList(genes)
#'
#' @export
convertHumanPathwayList <- function(pathways.list){
	require("biomaRt")
  human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  new.pathways = list()
	for (i in 1:length(pathways.list)) {
    set.name <- names(pathways.list)[i]
    message(set.name)
    pathways <- pathways.list[[set.name]]
    new.pathways[[set.name]] = list()
    for (pathway in names(pathways)) {
    	genesV2 = getLDS(attributes = c("hgnc_symbol"), 
                   filters = "hgnc_symbol", 
                   values = pathways[[pathway]], 
                   mart = human, 
                   attributesL = c("mgi_symbol"), 
                   martL = mouse, 
                   uniqueRows=T)
    	humanx <- unique(genesV2[, 2])
    	new.pathways[[set.name]][[pathway]] = humanx
    }
  }
  return(new.pathways)
}

#' Renames mouse genes to human 
#'
#'
#' @param genes list of mouse genes
#' 
#' @return list containing human corresponding genes
#'
#' @examples
#' convertMouseGeneList(genes)
#'
#' @export
convertMouseGeneList <- function(de.table){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  new.names <- getLDS(attributes = c("mgi_symbol"),
                      filters = "mgi_symbol", 
                      values = de.table$gene, mart = mouse,
                      attributesL = c("hgnc_symbol"),
                      martL = human, 
                      uniqueRows=T)
  renamed.table <- dplyr::left_join(de.table, new.names, by = c("gene" = "MGI.symbol"))
  renamed.table <- dplyr::filter(renamed.table, !is.na(HGNC.symbol))
  renamed.table$mouse_gene <- renamed.table$gene
  renamed.table$gene <- renamed.table$HGNC.symbol
  return(renamed.table)
}