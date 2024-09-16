#' @export
load_fragments <- function(ds,
                           sampleId) {
  path <- ds@fragments[[sampleId]]
  fragments <- readRDS(path)
  return(fragments)
}

#' @importFrom GenomicRanges strand
#' @export
get_gene_region_strand <- function(regions,
                                   gene_name) {
  as.vector(strand(regions[regions$gene_name == gene_name]))
}

#' @export
get_gene_region_gr <- function(regionSet,
                               gene_name) {
  return(regionSet[regionSet$gene_name == gene_name])
}

#' @importFrom GenomicRanges findOverlaps
#' @export
get_gene_region_fragments <- function(fragments, regionSet, gene_name) {
  gene_gr <- get_gene_region_gr(regionSet, gene_name)
  overlaps <- findOverlaps(fragments, gene_gr, ignore.strand = TRUE)
  query_idx <- overlaps@from
  return(fragments[query_idx])
}

#' @importFrom GenomicRanges coverage start end seqnames
#' @importFrom tibble as_tibble
#' @export
get_gene_region_coverage <- function(fragments, regionSet, gene_name) {
  gene_gr <- get_gene_region_gr(regionSet, gene_name)
  gene_chr <- seqnames(gene_gr)
  start <- start(gene_gr)
  end <- end(gene_gr)
  gene_region_pos <- start(gene_gr):end(gene_gr)
  gene_region_fragments <- get_gene_region_fragments(fragments, regionSet, gene_name)
  gene_region_cov <- coverage(gene_region_fragments)[[gene_chr]][gene_region_pos]
  coverage_tib <- as_tibble(cbind(gene_region_pos, as.vector(gene_region_cov)))
  colnames(coverage_tib) <- c("position", "coverage")
  return(coverage_tib)
}

#' @importFrom GenomicRanges start end
#' @export
get_gene_tss <- function(gene_name, tssgr) {
  # if (genome != "ChrAccRAnnotationHg38") {
  #   message("Unsupported genome specification: ", genome)
  #   stop()
  # }
  # if (!require(genome, quietly = TRUE)) {
  #   install.packages("https://muellerf.s3.amazonaws.com/data/ChrAccR/data/annotation/ChrAccRAnnotationHg38_0.0.1.tar.gz")
  # }
  # library(ChrAccRAnnotationHg38)
  # tssgr <- get("getGeneAnnotation", asNamespace(genome))(anno = "gencode_coding", type = "tssGr")
  gene_tssgr <- tssgr[tssgr$gene_name == gene_name]
  if (length(gene_tssgr) != 1 || start(gene_tssgr) != end(gene_tssgr)) {
    message("Something wrong with TSS.")
    stop()
  }
  gene_tss <- start(gene_tssgr)
  return(gene_tss)
}

#' @importFrom tibble as_tibble
#' @import rlang
#' @importFrom dplyr mutate select group_by summarize
#' @importFrom tidyr unnest
#' @export
pool_gene_region_coverage.DsATAC = function(ds,
                                     gene_name,
                                     groupby,
                                     regions) {
  metadata = as_tibble(ds@sampleAnnot, rownames = "sampleid")
  metadata = metadata %>%
    mutate(region_coverage = lapply(sampleid, function(s, r, g, ds) {
                      message("Loading fragments for ", s)
                      frags = load_fragments(ds, s)
                      cov = get_gene_region_coverage(frags, r, g)
                      return(cov)
                      }, regions, gene_name, ds))
  metadata = metadata %>%
    select(sampleid, {{ groupby }}, region_coverage) %>%
    unnest(cols  = c(region_coverage))
  metadata = metadata %>%
    group_by({{ groupby }}, position) %>%
    summarize(coverage = mean(coverage))
  return(metadata)
}

#' @importFrom tibble as_tibble
#' @import rlang
#' @importFrom dplyr mutate select group_by summarize
#' @importFrom tidyr unnest
#' @export
pool_gene_region_coverage.list = function(ds,
                                          fragments_list,
                                          gene_name,
                                          groupby,
                                          regions) {
  metadata = as_tibble(ds@sampleAnnot, rownames = "sampleid")
  metadata = metadata %>%
    mutate(region_coverage = lapply(sampleid, function(s, r, g, ds) {
                      frags = fragments_list[[s]]
                      cov = get_gene_region_coverage(frags, r, g)
                      return(cov)
                      }, regions, gene_name, ds))
  metadata = metadata %>%
    select(sampleid, {{ groupby }}, region_coverage) %>%
    unnest(cols  = c(region_coverage))
  metadata = metadata %>%
    group_by({{ groupby }}, position) %>%
    summarize(coverage = mean(coverage))
  return(metadata)
}
