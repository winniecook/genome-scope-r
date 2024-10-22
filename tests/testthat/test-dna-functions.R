library(testthat)
library(Biostrings)

# Test data
test_sequence <- "ATGCTAGCTAGCTAG"
invalid_sequence <- "ATGCTX"

test_that("DNA validation works", {
  expect_true(is_valid_dna(test_sequence))
  expect_false(is_valid_dna(invalid_sequence))
  expect_true(is_valid_dna(""))
})

test_that("Nucleotide counting works", {
  counts <- count_nucleotides(test_sequence)
  expect_equal(counts["A"], 3)
  expect_equal(counts["T"], 4)
  expect_equal(counts["G"], 4)
  expect_equal(counts["C"], 4)
})

test_that("GC content calculation is correct", {
  gc <- calculate_gc_content(test_sequence)
  expect_equal(gc, 53.33, tolerance = 0.01)
})

test_that("Complementary sequence is correct", {
  comp <- get_complementary_sequence(test_sequence)
  expect_equal(comp, "TACGATCGATCGATC")
})

test_that("ORF finding works", {
  seq_with_orf <- "ATGCTAGCTAGCTAG"  # Contains start codon
  orfs <- find_orfs(seq_with_orf)
  expect_type(orfs, "list")
  expect_true(length(orfs) > 0)
})

test_that("DNA properties calculation works", {
  props <- calculate_dna_properties(test_sequence)
  expect_type(props, "list")
  expect_true(all(c("length", "molecular_weight", "melting_temp", "complexity") %in% names(props)))
})

test_that("Repeating pattern detection works", {
  seq_with_repeats <- "ATGATGATG"
  patterns <- find_repeating_patterns(seq_with_repeats)
  expect_true(nrow(patterns) > 0)
  expect_true("ATG" %in% patterns$pattern)
})

test_that("Structure prediction works", {
  struct <- predict_dna_structure(test_sequence)
  expect_type(struct, "list")
  expect_true(all(c("stacking_energy", "bendability", "stability") %in% names(struct)))
})

test_that("Visualization functions work", {
  plots <- create_advanced_plots(test_sequence)
  expect_type(plots, "list")
  expect_true(all(c("distribution", "gc_skew", "cumulative") %in% names(plots)))
})

test_that("Circular plot creation works", {
  plot <- create_circular_plot(test_sequence)
  expect_s3_class(plot, "ggplot")
})