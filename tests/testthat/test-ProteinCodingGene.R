# Test ProteinCodingGene class creation and methods
test_that("ProteinCodingGene class works", {
  protein_gene <- GenesPackage::createProteinCodingGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list(), "protein1", "SEQUENCE")
  expect_s4_class(protein_gene, "ProteinCodingGene")
  expect_equal(GenesPackage::getProteinID(protein_gene), "protein1")
  expect_equal(GenesPackage::getProteinSequence(protein_gene), "SEQUENCE")

  # Test slot types
  expect_type(GenesPackage::getID(protein_gene), "integer")
  expect_type(GenesPackage::getSymbol(protein_gene), "character")
  expect_type(GenesPackage::getName(protein_gene), "character")
  expect_type(GenesPackage::getDescription(protein_gene), "character")
  expect_s4_class(GenesPackage::getStructure(protein_gene), "GRanges")
  expect_type(GenesPackage::getProduct(protein_gene), "list")
  expect_type(GenesPackage::getProteinID(protein_gene), "character")
  expect_type(GenesPackage::getProteinSequence(protein_gene), "character")
})

# Test ProteinCodingGene constructor handles errors
test_that("ProteinCodingGene constructor handles errors", {
  expect_error(GenesPackage::createProteinCodingGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list(), "protein1", 123), "invalid object for slot \"proteinSequence\" in class \"ProteinCodingGene\": got class \"numeric\", should be or extend class \"character\"")
})
