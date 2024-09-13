#' @importFrom Rsamtools BamFile
#' @export
check_seqnames = function(bam_paths, granges) {
  for (bam_path in bam_paths) {
    bam = BamFile(bam_path)
    bam_seqnames = seqnames(seqinfo(bam))
    granges_seqnames = levels(seqnames(granges))
    missing_in_bam = granges_seqnames[which(!granges_seqnames %in% bam_seqnames)]
    missing_in_granges = bam_seqnames[which(!bam_seqnames %in% granges_seqnames)]
    n_missing_in_bam = length(missing_in_bam)
    n_missing_in_granges = length(missing_in_granges)
    if (n_missing_in_bam > 0) {
      message("Warning: There are ", n_missing_in_bam, " reference names in GRanges that are missing in the BAM.")
    }
    if (n_missing_in_granges > 0) {
      message("Warning: There are ", n_missing_in_granges, " reference names in BAM that are missing in the Granges.")
    }
  }
}
