suppressPackageStartupMessages({
  library(rtracklayer)
  library(BiocGenerics)
  library(IRanges)
  library(GenomicRanges)
  library(S4Vectors)
  library(Gviz)
  library(ensembldb)
  library(SingleCellExperiment)
})


#' Custom coverage plot
#'
#' Plot the coverage in a region of the genome, for a set of bigwig files.
#'
#' @param se A \code{SummarizedExperiment} object.
#' @param rows Selected rows (i.e., features).
#' @param columns Selected columns (i.e., samples).
#' @param bigwig_files Paths to bigwig files.
#' @param bigwig_names Names (e.g. sample IDs) for bigwig files.
#' @param bigwig_conditions Groups for bigwig files, used for coloring the 
#' coverage plots.
#' @param id_field The rowData column to use as feature id. If omitted, 
#' row.names will be used.
#' @param showgene Gene ID to show.
#' @param granges optional path to .rds file containing GRanges object with gene
#'  annotations.
#' @param db The name of the EnsDb package to use as gene annotation. Ignored if
#' `granges` is given.
#' @param chr Chromosome to show. Ignored if \code{showgene} is not "".
#' @param start,end Start and end position of the region to show. Ignored if
#'   \code{showgene} is not "".
#' @param othersDense Logical; whether to collapse transcripts for genes others
#' than the one of interest (default TRUE)
#'
#' @return A gene coverage plot
#'
#' @author Charlotte Soneson, adapted by Pierre-Luc Germain
CUSTOM_GVIZ <- function(se, rows, columns, bigwig_files="", bigwig_names="",
                        bigwig_conditions="", showgene="", id_field="", 
                        granges="", db="", chr="", start="", end="",
                        othersDense=TRUE ) {
  options(ucscChromosomeNames = FALSE)

  ## ---------------------------------------------------------------------- ##
  ## Pre-flight checks
  ## ---------------------------------------------------------------------- ##
  ## Must have at least one of bigwig_files or gene db
  if (bigwig_files == "" && db == "" && granges=="") {
    return(NULL)
  }
  
  ## If no names are given, assign names to bigwig files
  if (bigwig_files != "" && bigwig_names == "") {
    bigwig_names <- paste(paste0("S", seq_along(strsplit(bigwig_files, ",")[[1]])),
                          collapse = ",")
  }
  
  ## If granges file does not exist, don't show annotation
  if (!file.exists(granges)) {
    granges <- ""
  }
  
  ## If granges not given and a db is defined, try to load it; if unable ignore
  if (granges == "" && db != "") {
    db <- tryCatch({
      library(db, character.only = TRUE)
      db
    }, error=function(e){
      warning(e)
      ""
    })
  }
  
  ## If no db given, the viewing region must be set
  if (granges == "" && db == "" && (chr == "" || start == "" || end == "")) {
    return(NULL)
  }
  
  ## Convert start and end positions to numeric values
  if (start != "") start <- as.numeric(start)
  if (end != "") as.numeric(end)

  ## If rows given, overwrite any provided showgene
  if (length(rows)==1) {
    if(id_field==""){
      showgene <- rows
    }else if(id_field %in% colnames(rowData(se))){
      showgene <- rowData(se)[rows,id_field]
      if(is.na(showgene)) return(NULL)
    }
  }
  
  ## Strip version number from the gene of interest if it exists
  showgene <- gsub("\\.[0-9]+$", "", showgene)

  if (showgene == "" && (chr == "" || is.na(start) || is.na(end))) {
    return(NULL)
  }

  sel.gr <- NULL
  if (showgene == ""){
    sel.gr <- GRanges(seqnames=chr, IRanges(start=start, end=end), strand="*")
  }
  
  ## ---------------------------------------------------------------------- ##
  ## Prepare the annotation
  ## ---------------------------------------------------------------------- ##
  if (granges != "") {
    gt <- .geneTracksFromGR(showgene, granges, sel.gr, othersDense=othersDense)
  } else {
    gt <- .geneTracksFromDb(showgene, sel.gr, db, othersDense = othersDense)
  }
  tracks <- c(Gviz::GenomeAxisTrack(), gt[[1]], gt[[2]])

  ## ---------------------------------------------------------------------- ##
  ## Set title and viewing range
  ## ---------------------------------------------------------------------- ##
  if (showgene != ""){
    # exit if the feature of interest was not found
    if(is.null(gt[[1]]) || length(gt[[1]])==0) return(NULL)
    chr <- seqnames(gt[[1]])[1]
    start <- min(start(gt[[1]]))
    end <- max(end(gt[[1]]))
    start <- round(start - 0.15*(end - start))
    end <- round(end + 0.05*(end - start))
  }
  plot_title <- paste0(showgene, ifelse(showgene=="",""," - "), 
                       chr, ":", start, "-", end)

  ## ---------------------------------------------------------------------- ##
  ## Prepare bigWig files
  ## ---------------------------------------------------------------------- ##
  if (bigwig_files != "") 
    tracks <- c( .bwTracks(bigwig_files, bigwig_names, bigwig_conditions),
                 tracks, )
  
  ## Plot tracks
  Gviz::plotTracks(tracks, chromosome = chr, from = start,
                   to = end, main = plot_title,
                   transcriptAnnotation = "transcript",
                   min.width = 0, min.distance = 0, collapse = FALSE)
}


#' .geneTracksFromDb
#' 
#' Create Gviz gene tracks from an EnsDb and either a feature name or a GRanges
#'
#' @param feature The id/name of the feature
#' @param sel.gr A selection GRanges, ignored if `features`!=""
#' @param db An `EnsDb`
#' @param fillBy Exon fill color options, currently only accepts 
#' "protein_coding" (default, colors transcripts differently depending on 
#' whether they are protein-coding or not) or "" (default color for all 
#' transcripts)
#' @param othersDense Logical; whether to collapse transcripts for genes others
#' than the one of interest (default TRUE)
#'
#' @return A list of length 2 with (eventually) GeneRegionTracks
.geneTracksFromDb <- function( feature="", sel.gr=NULL, db, 
                               fillBy=c("protein_coding",""), 
                               othersDense=!is.null(feature) ){
  if(is.null(sel.gr) && (is.null(feature) || feature=="")) 
    return(list(NULL,NULL))
  if(is.character(db)){
    if(!exists(db)) return(list(NULL,NULL))
    db <- get(db)
  }
  colFn <- switch( match.arg(fillBy),
                   protein_coding=function(e) 
                     c("gray80","darkblue")[(e$tx_biotype=="protein_coding") + 1],
                   function(e) NULL )
  
  if(feature != ""){
    # get feature of interest
    id_fields <- c("symbol","gene_id","gene_name","tx_id")
    filt <- as.formula( paste0("~ ", paste0(id_fields, "=='", 
                                            feature,"'", collapse=" | ")) )
    e <- exons(db, filter=filt, columns=c(id_fields, "tx_biotype"))
    e$transcript <- e$tx_id
    tr1 <- GeneRegionTrack(e, name=feature, col="gray80", fill=colFn(e), 
                           transcriptAnnotation="transcript", col.title="black")
    
    # get other features in the region
    filt <- GRangesFilter( GRanges(seqnames(e)[1], 
                                   IRanges(min(start(e)), max(end(e)))) )
  }else{
    tr1 <- NULL
    filt <- GRangesFilter( sel.gr )
  }
  
  e2 <- exons(db, filter=filt, columns=c(id_fields, "tx_biotype"))
  e2 <- e2[!(e2$tx_id %in% e$tx_id)]
  e2$transcript <- e2$tx_id
  
  tr2 <- GeneRegionTrack(e, name="others", col="gray80", fill=colFn(e2), 
                         transcriptAnnotation="transcript", col.title="black",
                         stacking=ifelse(othersDense,"dense","pack"))
  list(tr1, tr2)
}

#' .geneTracksFromGR
#' 
#' Create Gviz gene tracks from an annotation GRanges object and either a 
#' feature name or a selection GRanges
#'
#' @param showgene The id/name of the feature
#' @param granges Path to .rds file containing GRanges object with gene annotations
#' @param sel.gr A selection GRanges, ignored if `features`!=""
#' @param othersDense Logical; whether to collapse transcripts for genes others
#' than the one of interest (default TRUE)
#'
#' @return A list of length 2 with (eventually) GeneRegionTracks
.geneTracksFromGR <- function( showgene="", granges, sel.gr=NULL, 
                               othersDense=TRUE ){
  if(is.null(sel.gr) && (is.null(showgene) || showgene==""))
    return(list(NULL,NULL))
  ## Read the GRanges object
  if (caching$granges == granges && !is.null(caching$gr0)) {
    gr0 <- caching$gr0
  } else {
    caching$gr0 <- readRDS(granges)
    caching$granges <- granges
    gr0 <- caching$gr0
  }
  
  ## If a gene has been defined (either via rows or via showgene), set the
  ## viewing range accordingly
  if (showgene != "") {
    gr <- BiocGenerics::subset(gr0, tolower(gene) == tolower(showgene) |
                                 tolower(gene_name) == tolower(showgene))
    ## Select only one gene if there are many with the same name
    gr <- BiocGenerics::subset(gr, gene == gene[1])
    ## Other features in the region
    sel.gr <- GRanges( seqnames=GenomeInfoDb::seqnames(gr)[1],
                       ranges=IRanges::IRanges(start=start(gr), end=end(gr)),
                       strand = "*" )
    gro <- gr0[IRanges::overlapsAny(gr0, sel.gr)]
    gro <- gro[!(S4Vectors::`%in%`(gro, gr))]
    grtr <- GeneRegionTrack(gr, showId = TRUE, col = NULL, fill = "gray80",
                            name = showgene, col.title = "black", 
                            transcriptAnnotation="transcript")
  } else {
    gr <- NULL
    gro <- gr0[IRanges::overlapsAny(gr0, sel.gr)]
  }
  
  grtr2 <- GeneRegionTrack(gro, showId = TRUE, col = "black", fill = "white",
                           name = "others", col.title = "black", 
                           transcriptAnnotation="transcript", 
                           stacking=ifelse(othersDense,"dense","pack"))  
  list(grtr, grtr2)
}

.bwTracks <- function(bigwig_files, bigwig_names, bigwig_conditions){
  ## Reformat bigWig file paths and names (provided to the function as
  ## character strings)
  bigwig_files <- strsplit(bigwig_files, ",")[[1]]
  bigwig_names <- strsplit(bigwig_names, ",")[[1]]
  if (bigwig_conditions != "") {
    bigwig_conditions <- strsplit(bigwig_conditions, ",")[[1]]
    names(bigwig_conditions) <- bigwig_names
  }
  names(bigwig_files) <- bigwig_names
  
  ## ---------------------------------------------------------------------- ##
  ## Define colors if bigwig_conditions is provided
  ## ---------------------------------------------------------------------- ##
  ## Define colors for coverage tracks
  color_list <- rep(c("#DC050C", "#7BAFDE", "#B17BA6", "#F1932D", "#F7EE55",
                      "#90C987", "#777777", "#E8601C", "#1965B0", "#882E72",
                      "#F6C141", "#4EB265", "#CAEDAB"),
                    ceiling(length(unique(bigwig_conditions))/13))
  
  if (length(bigwig_conditions) > 1 || bigwig_conditions != "") {
    usecol <- color_list[match(bigwig_conditions,
                               unique(bigwig_conditions))]
  } else {
    usecol <- rep("gray", length(bigwig_files))
  }
  names(usecol) <- bigwig_names
  
  ## ------------------------------------------------------------------ ##
  ## Show only selected sample(s)
  ## ------------------------------------------------------------------ ##
  ## If columns is specified, subset bigwig files
  if (!is.null(columns)) {
    bigwig_files <- bigwig_files[columns]
    bigwig_conditions <- bigwig_conditions[columns]
    usecol <- usecol[columns]
  }
  
  ## ------------------------------------------------------------------ ##
  ## Prepare final plot
  ## ------------------------------------------------------------------ ##
  ## Set up coverage tracks
  lapply(seq_along(bigwig_files), function(i) {
    assign(paste0("covtr", i),
           Gviz::DataTrack(range = bigwig_files[i],
                           type = "histogram",
                           name = names(bigwig_files)[i],
                           col.title = "black",
                           fill = usecol[i],
                           col = usecol[i],
                           col.histogram = usecol[i],
                           fill.histogram = usecol[i]))
  })
}

# Set up a cache for the GRanges object
caching <- new.env()
