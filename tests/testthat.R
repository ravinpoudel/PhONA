library(testthat)
library(PhONA)

test_check("PhONA")

test_that("Test sqrt_newton: positive numeric",{
  expected <- 30
  actual <- getMean(10,20, 30, 40 , 50)
  expect_equal(expected, actual)
})
