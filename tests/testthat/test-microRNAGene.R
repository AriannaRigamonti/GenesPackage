# Test microRNAGene class creation and methods
test_that("MicroRNAGene class works", {
  mirna_gene <- GenesPackage::createMicroRNAGene(3L, "SYMBOL_MIR", "MicroRNA Name", "MicroRNA Description", "chr1", 1, 1000, "+", list(), "mirna1", "SEED_SEQ")
  expect_s4_class(mirna_gene, "microRNAGene")
  expect_equal(GenesPackage::getMicroRNAID(mirna_gene), "mirna1")
  expect_equal(GenesPackage::getMicroRNASequence(mirna_gene), "SEED_SEQ")

  # Test slot types
  expect_type(GenesPackage::getID(mirna_gene), "integer")
  expect_type(GenesPackage::getSymbol(mirna_gene), "character")
  expect_type(GenesPackage::getName(mirna_gene), "character")
  expect_type(GenesPackage::getDescription(mirna_gene), "character")
  expect_s4_class(GenesPackage::getStructure(mirna_gene), "GRanges")
  expect_type(GenesPackage::getProduct(mirna_gene), "list")
  expect_type(GenesPackage::getMicroRNAID(mirna_gene), "character")
  expect_type(GenesPackage::getMicroRNASequence(mirna_gene), "character")
})

# Test MicroRNAGene constructor handles errors
test_that("MicroRNAGene constructor handles errors", {
  expect_error(GenesPackage::createMicroRNAGene(3L, "SYMBOL_MIR", "MicroRNA Name", "MicroRNA Description", "chr1", 1, 1000, "+", list(), "mirna1", 123), "invalid object for slot \"microRNASequence\" in class \"microRNAGene\": got class \"numeric\", should be or extend class \"character\"")
})
