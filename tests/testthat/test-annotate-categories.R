# The ORF categories of annotate_ORFs(). tx_category() replaced a cascade of if
# where the last true one won; it must give the same label for all inputs.

test_that("tx_category() gives each ORF_category_Tx label", {
  # Annotated CDS from 101 to 203 in the transcript, with the stop codon:
  # ann_sta = 101, ann_sto = 200. The ORF (sta, sto) is without the stop codon.
  cases <- rbind(
    c(101, 200, "ORF_annotated"),
    c( 50, 200, "N_extension"),
    c(131, 200, "N_truncation"),
    c( 20,  79, "uORF"),          # ends before the CDS start
    c( 20, 120, "overl_uORF"),    # ends in the CDS
    c( 20, 101, "overl_uORF"),    # ends on the CDS start
    c( 50, 260, "NC_extension"),
    c(131, 260, "overl_dORF"),    # starts in the CDS
    c(200, 260, "overl_dORF"),    # starts on ann_sto
    c(231, 290, "dORF"),          # starts after ann_sto
    c(131, 170, "nested_ORF"),
    c(101, 170, "C_truncation"),
    c(101, 260, "C_extension"))
  for (k in seq_len(nrow(cases))) {
    expect_identical(tx_category(as.numeric(cases[k, 1]), as.numeric(cases[k, 2]), 101, 200),
                     cases[k, 3], info = paste(cases[k, 1:2], collapse = " "))
  }
})

# The cascade that tx_category() replaced (R/annotate_orfs.R L786-797 @ 6abfe9d)
old_tx_category <- function(sta, sto, ann_sta, ann_sto) {
  category <- NA
  if(sto==ann_sto){
    if(sta==ann_sta){category<-"ORF_annotated"}
    if(sta<ann_sta){category<-"N_extension"}
    if(sta>ann_sta){category<-"N_truncation"}
  }
  if(sto!=ann_sto){
    if(sta<ann_sta & sto<ann_sto){category<-"overl_uORF"}
    if(sta<ann_sta & sto<ann_sta){category<-"uORF"}
    if(sta<ann_sta & sto>ann_sto){category<-"NC_extension"}
    if(sta>ann_sta & sto>ann_sto){category<-"overl_dORF"}
    if(sta>ann_sto & sto>ann_sto){category<-"dORF"}
    if(sta>ann_sta & sto<ann_sto){category<-"nested_ORF"}
    if(sta==ann_sta & sto<ann_sto){category<-"C_truncation"}
    if(sta==ann_sta & sto>ann_sto){category<-"C_extension"}
  }
  category
}

test_that("tx_category() equals the old cascade for every order of the 4 values", {
  # Values 1 to 5 give every order of sta, sto, ann_sta and ann_sto, with ties,
  # also those that a real ORF or CDS cannot have (sta>sto, ann_sta>ann_sto)
  grid <- expand.grid(sta = 1:5, sto = 1:5, ann_sta = 1:5, ann_sto = 1:5)
  new <- mapply(tx_category, grid$sta, grid$sto, grid$ann_sta, grid$ann_sto)
  old <- mapply(old_tx_category, grid$sta, grid$sto, grid$ann_sta, grid$ann_sto)
  expect_identical(new, old)
  expect_false(anyNA(new))
  # sta from compatible_ORF_id_tr is numeric, not integer
  expect_identical(tx_category(20, 79L, 101L, 200L), "uORF")
})
