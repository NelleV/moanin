
#' Find enriched GO terms
#'
#' @param assignments boolean named vector determining the gene subset to be
#'   tested for enrichment of GO terms. The names of the vector should be the
#'   gene names. Elements with TRUE will consist of the gene cluster.
#' @param gene_id_to_go List giving the Gene ID to GO object required for topGO
#'   (see \code{\link[pkg:topGO]{topGOdata-class}}). 
#'   \code{create_go_term_mapping} can construct such a list from a
#'   data-frame.
#'@param ontology string, optional, default: BP. specficies which ontology to
#'  use (passed to \code{ontology} argument in creating a new \code{topGOdata}
#'  object). Can be 'BP', 'CC', or 'NF'. See
#'  \code{\link[pkg:topGO]{topGOdata-class}}.
#'@param weighted boolean, optional, default: FALSE. Whether to use the
#'  weighted algorithm or not in \code{\link[pkg:topGO]{runTest}}.
#'@param node_size integer, optional, default: 10. Consider only GO terms with
#'  node_size number of genes, passed to \code{nodeSize} argument of
#'  \code{\link[pkg:topGO]{topGOdata-class}}
#'@details \code{find_enriched_go_terms} is a wrapper for running a GO
#'  enrichment analysis via the package \code{topGO}. This function creates a
#'  \code{\link[pkg:topGO]{topGOdata-class}} object, runs the function
#'  \code{\link[pkg:topGO]{runTest}} to test for enrichment using the
#'  \code{statistic="fisher"} option, and then runs
#'  \code{\link[pkg:topGO]{GenTable}}. This function then does some post-processing
#'  of the results, returning only GO terms that satisfy:
#'  \enumerate{
#'  \item{BH adjusted p-values less than 0.05 using
#'  \code{\link[stats]{p.adjust}}}
#'  \item{GO terms are \emph{enriched}, i.e. the number of genes from the GO 
#'  term found in the subset is greater than expected}
#'  }
#' @return Returns results in the format of \code{\link[pkg:topGO]{GenTable}}.
#' @aliases create_go_term_mapping
#' @examples 
#' data(exampleData)
#' head(testGenesGO) #gives the mapping of genes to GO
#' geneId2Go <- create_go_term_mapping(testGenesGO)
#' #create fake assignment of genes to group based on TRUE/FALSE values
#' inGroup=rep(FALSE,nrow(testData))
#' inGroup[1:10]=TRUE
#' names(inGroup)<-names(geneId2Go)
#' find_enriched_go_terms(inGroup,geneId2Go)
#' @seealso \code{\link{create_go_term_mapping}}, \code{\link[KEGGprofile]{find_enriched_pathway}}, \code{\link[pkg:topGO]{GenTable}}, \code{\link[pkg:topGO]{runTest}}, \code{\link[pkg:topGO]{topGOdata-class}}, \code{\link[stats]{p.adjust}} 
#' @export
#' @importFrom topGO annFUN.gene2GO runTest GenTable
find_enriched_go_terms <- function(assignments, gene_id_to_go,
                                  ontology="BP", 
                                  weighted=FALSE,
                                  node_size=10){
    if(is.null(names(assignments))) stop("assignments must be a named vector,",
    "where names match names of gene_id_to_go ")
    gene_names <- names(assignments)
    assignments <- as.numeric(assignments)
    names(assignments) <- gene_names
    if(!(ontology %in% c("BP", "CC", "NF"))){
        error_message <- paste(
            "moanin::find_enriched_go_terms: Ontology should be 'BP', 'CC',",
            "or 'NF'. Ontology provided is",
            ontology)
        stop(error_message)
    }
    
    # function to specify which genes are interesting based on the gene scores.
    # It should be present iff the allGenes object is of type numeric.
    getTopDiffGenes <- function(data, cutOff=NULL){
        return(data < 0.5)
    }
    
    GOdata <- new("topGOdata", ontology=ontology, allGenes=assignments,
                 nodeSize=node_size,
                 geneSel=getTopDiffGenes,
                 annot=topGO::annFUN.gene2GO, gene2GO=gene_id_to_go)
    
    if(weighted){
        resultFisher <- topGO::runTest(GOdata, algorithm="classic", 
                                      statistic="fisher")
    }else{
        resultFisher <- topGO::runTest(GOdata, algorithm="weight", 
                                      statistic="fisher")
    }
    n_nodes <- length(resultFisher@score)
    
    allRes <- topGO::GenTable(
        GOdata, 
        resultFisher=resultFisher,
        orderBy="resultFisher", ranksOf="resultFisher",
        topNodes=n_nodes)
    
    # P-value correct
    allRes[, "resultFisher_padj"] <- stats::p.adjust(allRes$resultFisher, 
                                                    method="BH")
    wh <- which(allRes[, "resultFisher"] <= 0.05)
    
    allRes <- allRes[wh,]
    wh <- which(apply(allRes[, c("Significant", "Expected")],
                     1, function(x){x["Significant"] > x["Expected"]}))
    allRes <- allRes[wh,]
    
    return(allRes)
}


#' Create the Gene to GO Term mapping
#'
#' @param genes dataframe, with two required columns. The first gives the gene
#'   names, with column name by the argument \code{gene_col}. The other column
#'   must be named "go_id" and give the genes GO id. Genes will have multiple GO
#'   id that they map to, and each go mapping of a gene is a separate row. Thus
#'   genes will be in multiple rows of the input.
#' @param gene_col the name of the column of the \code{genes} data frame that
#'   contains the correct gene reference. By default, is "refseq_mrna".
#' @return \code{create_go_term_mapping} returns a list giving the gene to GO id
#'   in the format required by \code{\link[pkg:topGO]{topGOdata-class}}.
#' @rdname find_enriched_go_terms
#' @export
create_go_term_mapping <- function(genes, gene_col="refseq_mrna"){
    gene_id_go_mapping <- NULL
    gene_names <- unique(genes[, "refseq_mrna"])
    
    i <- 1
    for(gene in gene_names){
        go_terms <- genes[genes[, gene_col] == gene, "go_id"]
        if(length(go_terms) != 0){
            gene_id_go_mapping$gene <- go_terms
            names(gene_id_go_mapping)[i] <- gene
            i <- i + 1
        }
    }
    return(gene_id_go_mapping)
}
