context("tkmeans package testing suite")

test_that("Basic functionality not compromised", {
  set.seed(123)
  iris_mat <- as.matrix(iris[,1:4])
  expect_equal_to_reference(tkmeans_lowmem(iris_mat, 2 , 0.1, 1, 10, 0.001), "test_one.rds")
  expect_equal_to_reference(tkmeans_lowmem(iris_mat, 3 , 0.2, 1, 10, 0.001), "test_two.rds")
  expect_equal_to_reference(tkmeans_lowmem(iris_mat, 1 , 0.1, 1, 10, 0.001), "test_three.rds")
})
