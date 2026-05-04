test_that(".data_bearing_classes is the canonical 4-class set", {
  expect_setequal(
    Nestimate:::.data_bearing_classes,
    c("netobject", "netobject_group", "mcml", "htna")
  )
})

test_that(".is_data_bearing recognises every member", {
  ne <- structure(list(), class = c("netobject", "cograph_network"))
  ng <- structure(list(), class = "netobject_group")
  mc <- structure(list(), class = "mcml")
  ht <- structure(list(), class = c("htna", "netobject", "cograph_network"))

  expect_true(Nestimate:::.is_data_bearing(ne))
  expect_true(Nestimate:::.is_data_bearing(ng))
  expect_true(Nestimate:::.is_data_bearing(mc))
  expect_true(Nestimate:::.is_data_bearing(ht))
})

test_that(".is_data_bearing rejects non-members", {
  expect_false(Nestimate:::.is_data_bearing(1:5))
  expect_false(Nestimate:::.is_data_bearing(data.frame(x = 1)))
  expect_false(Nestimate:::.is_data_bearing(matrix(0, 2, 2)))
  expect_false(Nestimate:::.is_data_bearing(list(a = 1)))
  expect_false(Nestimate:::.is_data_bearing(structure(list(), class = "tna")))
})

test_that(".data_bearing_class returns the canonical token", {
  ne <- structure(list(), class = c("netobject", "cograph_network"))
  ng <- structure(list(), class = "netobject_group")
  mc <- structure(list(), class = "mcml")
  ht <- structure(list(), class = c("htna", "netobject", "cograph_network"))

  expect_identical(Nestimate:::.data_bearing_class(ne), "netobject")
  expect_identical(Nestimate:::.data_bearing_class(ng), "netobject_group")
  expect_identical(Nestimate:::.data_bearing_class(mc), "mcml")
  # htna inherits netobject; htna must win.
  expect_identical(Nestimate:::.data_bearing_class(ht), "htna")
})

test_that(".data_bearing_class returns NA for non-members", {
  expect_true(is.na(Nestimate:::.data_bearing_class(1:5)))
  expect_true(is.na(Nestimate:::.data_bearing_class(structure(list(), class = "tna"))))
})
