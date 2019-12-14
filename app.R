suppressPackageStartupMessages({
  library(iSEE)
  library(SingleCellExperiment)
})

source("custom.R")

sce <- readRDS("shiny_sce.rds")$sce_gene
rowData(sce) <- tidyr::unnest(as.data.frame(rowData(sce)))

cdp <- customDataPlotDefaults(sce, 1)
cdp$Function <- "CUSTOM_GVIZ"
cdp$Arguments <- c("db EnsDb.Hsapiens.v86\nid_field gene_id")
cdp$RowSource <- "Row data plot 1"

rdp <- rowDataPlotDefaults(sce, 5)
rdp$YAxis <- "edgeR.cellineN61311.cellineN052611.mlog10PValue"
rdp$XAxis <- "Row data"
rdp$XAxisRowData <- "edgeR.cellineN61311.cellineN052611.logFC"

initialPanels <- DataFrame(
  Name=c("Row data plot 1", "Custom data plot 1", "Row statistics table 1"),
  Width=c(4L, 8L, 4L),
  Height=c(500L, 600L, 400L)
)

app <- iSEE::iSEE(
  se=sce,
  rowDataArgs = rdp,
  rowStatArgs = rowStatTableDefaults(sce, 5),
  initialPanels=initialPanels,
  customDataArgs = cdp,
  customDataFun = list(CUSTOM_GVIZ = CUSTOM_GVIZ)
)


