# Test lengthProduct handles edge cases
test_that("lengthProduct handles edge cases correctly", {
  protein_gene <- GenesPackage::createProteinCodingGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list(), "protein1", "SEQUENCE")
  expect_equal(GenesPackage::computeGeneLength(protein_gene), 1000)

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
