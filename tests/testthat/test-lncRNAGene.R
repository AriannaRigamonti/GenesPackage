# Test lncRNAGene class creation and methods
test_that("LncRNAGene class works", {
  lncrna_gene <- GenesPackage::createLncRNAGene(2L, "SYMBOL_LNC", "LncRNA Name", "LncRNA Description", "chr1", 1, 1000, "+", list(), "lncrna1", "RNA_SEQUENCE")
  expect_s4_class(lncrna_gene, "lncRNAGene")
  expect_equal(GenesPackage::getLncRNAID(lncrna_gene), "lncrna1")
  expect_equal(GenesPackage::getRNASequence(lncrna_gene), "RNA_SEQUENCE")

  # Test slot types
  expect_type(GenesPackage::getID(lncrna_gene), "integer")
  expect_type(GenesPackage::getSymbol(lncrna_gene), "character")
  expect_type(GenesPackage::getName(lncrna_gene), "character")
  expect_type(GenesPackage::getDescription(lncrna_gene), "character")
  expect_s4_class(GenesPackage::getStructure(lncrna_gene), "GRanges")
  expect_type(GenesPackage::getProduct(lncrna_gene), "list")
  expect_type(GenesPackage::getLncRNAID(lncrna_gene), "character")
  expect_type(GenesPackage::getRNASequence(lncrna_gene), "character")
})

# Test LncRNAGene constructor handles errors
test_that("LncRNAGene constructor handles errors", {
  expect_error(GenesPackage::createLncRNAGene(2L, "SYMBOL_LNC", "LncRNA Name", "LncRNA Description", "chr1", 1, 1000, "+", list(), "lncrna1", 123), "invalid object for slot \"RNASequence\" in class \"lncRNAGene\": got class \"numeric\", should be or extend class \"character\"")
})
