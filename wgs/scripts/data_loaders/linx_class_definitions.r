## Class definitions for storing LINX output data as R objects


#####################
# CLASS DEFINITIONS #
#####################


LINX_FUSION <- setClass("LINX_FUSION",  slots = list(
  fivePrimeBreakendId="numeric",
  threePrimeBreakendId="numeric",
  name="character",
  reported="logical",
  reportedType="character",
  phased="character",
  likelihood="character",
  chainLength="numeric",
  chainLinks="numeric",
  chainTerminated="logical",
  domainsKept="character",
  domainsLost="character",
  skippedExonsUp="numeric",
  skippedExonsDown="numeric",
  fusedExonUp="numeric",
  fusedExonDown="numeric",
  geneStart="character",
  geneContextStart="character",
  transcriptStart="character",
  geneEnd="character",
  geneContextEnd="character",
  transcriptEnd="character",
  junctionCopyNumber="numeric"))

LINX_BREAKEND <- setClass("LINX_BREAKEND",  slots = list(
  id="numeric",
  svId="numeric",
  isStart="logical",
  gene="character",
  transcriptId="character",
  canonical="logical",
  geneOrientation="character",
  disruptive="logical",
  reportedDisruption="logical",
  undisruptedCopyNumber="numeric",
  regionType="character",
  codingContext="character",
  biotype="character",
  exonicBasePhase="numeric",
  nextSpliceExonRank="numeric",
  nextSpliceExonPhase="numeric",
  nextSpliceDistance="numeric",
  totalExonCount="numeric",
  type="character",
  chromosome="character",
  orientation="numeric",
  strand="numeric",
  chrBand="character",
  exonUp="numeric",
  exonDown="numeric",
  junctionCopyNumber="numeric",
  fusions="list"))


LINX_CLUSTER <- setClass("LINX_CLUSTER",  slots = list(
  clusterId="numeric",
  category="character",
  synthetic="logical",
  resolvedType="character",
  clusterCount="numeric",
  clusterDesc="character"))

LINX_SV <- setClass("LINX_SV", slots = list(
  vcfId="character",
  svId="numeric",
  clusterId="numeric",
  clusterReason="character",
  fragileSiteStart="logical",
  fragileSiteEnd="logical",
  isFoldback="logical",
  lineTypeStart="character",
  lineTypeEnd="character",
  junctionCopyNumberMin="numeric",
  junctionCopyNumberMax="numeric",
  geneStart="character",
  geneEnd="character",
  replicationTimingStart="numeric",
  replicationTimingEnd="numeric",
  localTopologyIdStart="numeric",
  localTopologyIdEnd="numeric",
  localTopologyStart="character",
  localTopologyEnd="character",
  localTICountStart="numeric",
  localTICountEnd="numeric",
  breakends="list",
  clusters="list",
  vcf_info="data.frame"))