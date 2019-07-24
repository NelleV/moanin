library(topGO)


#' Find enriched GO terms
#'
#' @param assignments boolean named vector
#' @param gene_id_to_go Gene ID to GO object from topGO
#' @param ontology string, optional, default: BP
#'	specficies which ontology to use. Can be 'BP', 'CC', or 'NF'
#' @param weighted, boolean, optional, default: FALSE
#'	Whether to use the weighted algorithm or not.
#' @param node_size integer, optional, default: 10
#'	Consider only GO terms with node_size number of genes.
#' @param ontology: BP, CC, NF
#' @export
find_enriched_go_terms = function(assignments, gene_id_to_go,
				  ontology="BP", 
				  weighted=FALSE,
				  node_size=10){
    gene_names = names(assignments)
    assignments = as.numeric(assignments)
    names(assignments) = gene_names
    if(!(ontology %in% c("BP", "CC", "NF"))){
	error_message = paste(
	    "moanin::find_enriched_go_terms: Ontology should be 'BP', 'CC',",
	    "or 'NF'. Ontology provided is",
	    ontology)
	stop(error_message)
    }


    getTopDiffGenes = function(data, cutOff=NULL){
      return(data < 0.5)
    }

    GOdata = new("topGOdata", ontology=ontology, allGenes=assignments,
		 nodeSize=node_size,
		 geneSel=getTopDiffGenes,
		 annot=annFUN.gene2GO, gene2GO=gene_id_to_go)

    if(weighted){
        resultFisher = topGO::runTest(GOdata, algorithm="classic", statistic="fisher")
    }else{
        resultFisher = topGO::runTest(GOdata, algorithm="weight", statistic="fisher")
    }
    n_nodes = length(resultFisher@score)

    allRes = topGO::GenTable(
	GOdata, 
	resultFisher=resultFisher,
	orderBy="resultFisher", ranksOf="resultFisher",
	topNodes=n_nodes)
        
    # P-value correct
    allRes[, "resultFisher_padj"] = stats::p.adjust(allRes$resultFisher, method="BH")
    wh = which(allRes[, "resultFisher"] <= 0.05)

    allRes = allRes[wh,]
    wh = which(apply(allRes[, c("Significant", "Expected")],
    	             1, function(x){x["Significant"] > x["Expected"]}))
    allRes = allRes[wh,]

    return(allRes)
}


#' Create the Gene to GO Term mapping
#'
#' @param genes dataframe, with `gene_col` and a column corresponding to the
#'  `go_id`
#' @param gene_col the column of the genes data frame that contains the
#'  correct gene reference. By default, is refseq_mrna
#' @export
create_go_term_mapping = function(genes, gene_col="refseq_mrna"){
    gene_id_go_mapping = NULL
    gene_names = unique(genes[, "refseq_mrna"])

    i = 1
    for(gene in gene_names){
	go_terms = genes[genes[, gene_col] == gene, "go_id"]
	if(length(go_terms) != 0){
	    gene_id_go_mapping$gene = go_terms
	    names(gene_id_go_mapping)[i] = gene
	    i = i + 1
	}
    }
    return(gene_id_go_mapping)
}
