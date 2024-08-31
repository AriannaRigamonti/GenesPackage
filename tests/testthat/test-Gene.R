# Test Gene class creation and methods
test_that("Gene class works", {
  gene <- GenesPackage::createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())
  expect_s4_class(gene, "Gene")
  expect_equal(GenomicRanges::seqnames(GenesPackage::getStructure(gene)), GenomicRanges::seqnames(GenomicRanges::GRanges("chr1", IRanges::IRanges(1, 1000), "+")))
  expect_equal(GenesPackage::getID(gene), 1)
  expect_equal(GenesPackage::getSymbol(gene), "SYMBOL")
  expect_equal(GenesPackage::getName(gene), "Gene Name")
  expect_equal(GenesPackage::getDescription(gene), "Description")
})

# Test accessors and setters
test_that("Accessors and setters work correctly", {
  # Create a gene object
  gene <- GenesPackage::createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())

  # Test getID and setID
  GenesPackage::setID(gene) <- 2L
  expect_equal(GenesPackage::getID(gene), 2L)
  expect_type(GenesPackage::getID(gene), "integer")

  # Test getSymbol and setSymbol
  GenesPackage::setSymbol(gene) <- "NEW_SYMBOL"
  expect_equal(GenesPackage::getSymbol(gene), "NEW_SYMBOL")
  expect_type(GenesPackage::getSymbol(gene), "character")

  # Test getName and setName
  GenesPackage::setName(gene) <- "New Gene Name"
  expect_equal(GenesPackage::getName(gene), "New Gene Name")
  expect_type(GenesPackage::getName(gene), "character")

  # Test getDescription and setDescription
  GenesPackage::setDescription(gene) <- "New Description"
  expect_equal(GenesPackage::getDescription(gene), "New Description")
  expect_type(GenesPackage::getDescription(gene), "character")

  # Test getStructure and setStructure
  new_structure <- GenomicRanges::GRanges(seqnames = "chr2", ranges = IRanges::IRanges(start = 100, end = 2000), strand = "-")
  GenesPackage::setStructure(gene) <- new_structure
  expect_equal(GenomicRanges::seqnames(GenesPackage::getStructure(gene)), GenomicRanges::seqnames(new_structure))

  # Test getProduct and setProduct
  new_product <- list(func = "Enzyme", pathway = "Metabolic")
  GenesPackage::setProduct(gene) <- new_product
  expect_equal(GenesPackage::getProduct(gene), new_product)
})

# Test gene-specific methods
test_that("Gene-specific methods work correctly", {
  # Create a gene object
  gene <- GenesPackage::createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", list())

  # Test computeGeneLength
  expect_equal(GenesPackage::computeGeneLength(gene), 1000)
})

# Test constructors handle errors correctly
test_that("Constructors handle errors correctly", {
  expect_error(GenesPackage::createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", "start", 1000, "+", list()), "'start' and 'end' must be numeric vectors")
  expect_error(GenesPackage::createGene(1L, "SYMBOL", "Gene Name", "Description", "chr1", 1, 1000, "+", "not a list"), "invalid class")
})
