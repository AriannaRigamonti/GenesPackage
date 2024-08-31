# Test siRNAGene class creation and methods
test_that("SiRNAGene class works", {
  sirna_gene <- GenesPackage::createSiRNAGene(4L, "SYMBOL_SI", "siRNA Name", "siRNA Description", "chr1", 1, 1000, "+", list(), "sirna1", "SIRNA_SEQ")
  expect_s4_class(sirna_gene, "siRNAGene")
  expect_equal(GenesPackage::getSiRNAID(sirna_gene), "sirna1")
  expect_equal(GenesPackage::getSiRNASequence(sirna_gene), "SIRNA_SEQ")

  # Test slot types
  expect_type(GenesPackage::getID(sirna_gene), "integer")
  expect_type(GenesPackage::getSymbol(sirna_gene), "character")
  expect_type(GenesPackage::getName(sirna_gene), "character")
  expect_type(GenesPackage::getDescription(sirna_gene), "character")
  expect_s4_class(GenesPackage::getStructure(sirna_gene), "GRanges")
  expect_type(GenesPackage::getProduct(sirna_gene), "list")
  expect_type(GenesPackage::getSiRNAID(sirna_gene), "character")
  expect_type(GenesPackage::getSiRNASequence(sirna_gene), "character")
})
