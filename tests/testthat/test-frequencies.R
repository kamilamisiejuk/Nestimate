# ---- frequencies() tests ----

test_that("frequencies works with long format data", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 2, 2, 2),
    Time = c(1, 2, 3, 1, 2, 3),
    Action = c("A", "B", "A", "B", "A", "B"),
    stringsAsFactors = FALSE
  )

  freq <- frequencies(long_data, action = "Action", id = "Actor")

  expect_true(is.matrix(freq))
  expect_equal(nrow(freq), 2)
  expect_equal(ncol(freq), 2)
  expect_equal(rownames(freq), c("A", "B"))
  expect_equal(colnames(freq), c("A", "B"))
  expect_equal(freq["A", "B"], 2L)
  expect_equal(freq["B", "A"], 2L)
  expect_equal(freq["A", "A"], 0L)
  expect_equal(freq["B", "B"], 0L)
})

test_that("frequencies works with wide format data", {
  wide_data <- data.frame(
    T1 = c("A", "B"),
    T2 = c("B", "A"),
    T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  freq <- frequencies(wide_data, format = "wide")

  expect_true(is.matrix(freq))
  expect_equal(freq["A", "B"], 2L)
  expect_equal(freq["B", "A"], 2L)
})

test_that("frequencies auto-detects format", {
  long_data <- data.frame(
    Actor = c(1, 1, 2, 2),
    Time = c(1, 2, 1, 2),
    Action = c("X", "Y", "Y", "X"),
    stringsAsFactors = FALSE
  )
  freq <- frequencies(long_data, action = "Action", id = "Actor")
  expect_equal(freq["X", "Y"], 1L)
  expect_equal(freq["Y", "X"], 1L)

  wide_data <- data.frame(T1 = c("X", "Y"), T2 = c("Y", "X"),
                           stringsAsFactors = FALSE)
  freq2 <- frequencies(wide_data)
  expect_equal(freq2["X", "Y"], 1L)
  expect_equal(freq2["Y", "X"], 1L)
})

test_that("frequencies handles multiple ID columns", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 1),
    Group = c("G1", "G1", "G2", "G2"),
    Time = c(1, 2, 1, 2),
    Action = c("A", "B", "B", "A"),
    stringsAsFactors = FALSE
  )

  freq <- frequencies(long_data, action = "Action", id = c("Actor", "Group"))
  expect_equal(freq["A", "B"], 1L)
  expect_equal(freq["B", "A"], 1L)
})

test_that("frequencies handles NAs", {
  wide_data <- data.frame(
    T1 = c("A", "B"), T2 = c("B", NA), T3 = c(NA, NA),
    stringsAsFactors = FALSE
  )
  freq <- frequencies(wide_data, format = "wide")
  expect_equal(freq["A", "B"], 1L)
  expect_equal(sum(freq), 1L)

  long_data <- data.frame(
    Actor = c(1, 1, 1), Time = c(1, 2, 3),
    Action = c("A", NA, "B"), stringsAsFactors = FALSE
  )
  freq2 <- frequencies(long_data, action = "Action", id = "Actor")
  expect_equal(sum(freq2), 0L)
})

test_that("frequencies with cols parameter selects specific columns", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"), T3 = c("C", "C"),
    extra = c("ignore", "this"), stringsAsFactors = FALSE
  )
  freq <- frequencies(wide_data, cols = c("T1", "T2", "T3"), format = "wide")
  expect_true("C" %in% rownames(freq))
  expect_equal(freq["A", "B"], 1L)
  expect_equal(freq["B", "C"], 1L)
})

test_that("frequencies returns integer matrix", {
  wide_data <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"),
                           stringsAsFactors = FALSE)
  freq <- frequencies(wide_data, format = "wide")
  expect_true(is.integer(freq))
})

test_that("frequencies long format without id treats all as one sequence", {
  long_data <- data.frame(
    Time = 1:4, Action = c("A", "B", "A", "B"), stringsAsFactors = FALSE
  )
  freq <- frequencies(long_data, action = "Action")
  expect_equal(freq["A", "B"], 2L)
  expect_equal(freq["B", "A"], 1L)
})

test_that("frequencies errors on missing columns", {
  long_data <- data.frame(x = 1:3, y = c("A", "B", "C"))
  expect_error(frequencies(long_data, action = "Action", format = "long"),
               "not found")
  expect_error(frequencies(long_data, action = "y", id = "missing", format = "long"),
               "not found")
})

test_that("convert_sequence_format frequency format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"), T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "frequency")

  expect_true(is.data.frame(result))
  expect_true("id" %in% names(result))
  expect_true("A" %in% names(result))
  expect_true("B" %in% names(result))
  # Row 1 (id=1): A, B, A => A=2, B=1
  expect_equal(result[result$id == 1, "A"], 2L)
  expect_equal(result[result$id == 1, "B"], 1L)
  # Row 2 (id=2): B, A, B => B=2, A=1
  expect_equal(result[result$id == 2, "A"], 1L)
  expect_equal(result[result$id == 2, "B"], 2L)
})

test_that("convert_sequence_format onehot format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "B"), T3 = c("A", "B"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "onehot")

  expect_true(is.data.frame(result))
  # Row 1 (id=1): A and B present
  expect_equal(result[result$id == 1, "A"], 1L)
  expect_equal(result[result$id == 1, "B"], 1L)
  # Row 2 (id=2): only B present
  expect_equal(result[result$id == 2, "A"], 0L)
  expect_equal(result[result$id == 2, "B"], 1L)
})

test_that("convert_sequence_format edgelist format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "edgelist")

  expect_true(is.data.frame(result))
  expect_true(all(c("id", "from", "to") %in% names(result)))
  expect_equal(nrow(result), 2)
  # Row 1: A->B, Row 2: B->A
  expect_true(any(result$from == "A" & result$to == "B"))
  expect_true(any(result$from == "B" & result$to == "A"))
})

test_that("convert_sequence_format follows format from wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "follows")

  expect_true(is.data.frame(result))
  expect_true(all(c("id", "act", "follows") %in% names(result)))
  expect_equal(nrow(result), 2)
  # B follows A, A follows B
  expect_true(any(result$act == "B" & result$follows == "A"))
  expect_true(any(result$act == "A" & result$follows == "B"))
})

test_that("convert_sequence_format frequency format from long", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 2, 2),
    Time = c(1, 2, 3, 1, 2),
    Action = c("A", "B", "A", "B", "B"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(long_data, action = "Action",
                                    id_col = "Actor", time = "Time",
                                    format = "frequency")

  expect_true(is.data.frame(result))
  expect_true("Actor" %in% names(result))
  # Actor 1: A=2, B=1
  expect_equal(result[result$Actor == 1, "A"], 2L)
  expect_equal(result[result$Actor == 1, "B"], 1L)
  # Actor 2: A=0, B=2
  expect_equal(result[result$Actor == 2, "A"], 0L)
  expect_equal(result[result$Actor == 2, "B"], 2L)
})

test_that("convert_sequence_format edgelist format from long", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 2, 2),
    Time = c(1, 2, 3, 1, 2),
    Action = c("A", "B", "A", "B", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(long_data, action = "Action",
                                    id_col = "Actor", time = "Time",
                                    format = "edgelist")

  expect_true(all(c("Actor", "from", "to") %in% names(result)))
  # Actor 1: A->B, B->A  Actor 2: B->A
  expect_equal(nrow(result), 3)
})

test_that("convert_sequence_format handles NAs in wide", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", NA), T3 = c(NA, NA),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, id_col = "id",
                                    format = "frequency")

  # id=1: A=1, B=1   id=2: B=1
  expect_equal(result[result$id == 1, "A"], 1L)
  expect_equal(result[result$id == 1, "B"], 1L)
  expect_equal(result[result$id == 2, "B"], 1L)
})

test_that("convert_sequence_format with seq_cols selects specific columns", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    extra = c("X", "Y"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(wide_data, seq_cols = c("T1", "T2"),
                                    id_col = "id", format = "frequency")

  # "extra" column values should NOT appear as state columns
  expect_false("X" %in% names(result))
  expect_false("Y" %in% names(result))
})

test_that("convert_sequence_format does not guess id_col; explicit id_col honored (A11-F01 / Codex)", {
  wide_data <- data.frame(
    Student = c(1, 2),
    T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  # id_col = NULL: NO value-based id inference (the old "first column = id"
  # guess silently dropped the first time point). Every column is a state,
  # so Student's values 1/2 are states and nothing is dropped.
  el_noid <- convert_sequence_format(wide_data, format = "edgelist")
  expect_equal(nrow(el_noid), 4L)
  expect_true(any(c("1", "2") %in% c(el_noid$from, el_noid$to)))

  # A genuine leading id now requires passing id_col explicitly.
  el_id <- convert_sequence_format(wide_data, id_col = "Student",
                                   format = "edgelist")
  expect_equal(nrow(el_id), 2L)
  expect_false(any(c("1", "2") %in% c(el_id$from, el_id$to)))
  expect_true("Student" %in% names(el_id))
})

test_that("convert_sequence_format with multiple id_col from long", {
  long_data <- data.frame(
    Actor = c(1, 1, 1, 1),
    Group = c("G1", "G1", "G2", "G2"),
    Time = c(1, 2, 1, 2),
    Action = c("A", "B", "C", "A"),
    stringsAsFactors = FALSE
  )

  result <- convert_sequence_format(long_data, action = "Action",
                                    id_col = c("Actor", "Group"),
                                    time = "Time", format = "edgelist")

  expect_true(all(c("Actor", "Group", "from", "to") %in% names(result)))
  expect_equal(nrow(result), 2)
})

test_that("convert_sequence_format errors on missing columns", {
  data <- data.frame(x = 1:3, y = c("A", "B", "C"))

  expect_error(
    convert_sequence_format(data, seq_cols = c("missing1", "missing2")),
    "not found"
  )
})


# ---- .build_pairs_long() coverage ----

test_that(".build_pairs_long errors on missing action column", {
  df <- data.frame(actor = 1:3, time = 1:3, stringsAsFactors = FALSE)
  expect_error(
    Nestimate:::.build_pairs_long(df, action = "Action", id = NULL, time = "time"),
    "not found"
  )
})

test_that(".build_pairs_long errors on missing ID column", {
  df <- data.frame(actor = 1:3, Action = c("A", "B", "C"), stringsAsFactors = FALSE)
  expect_error(
    Nestimate:::.build_pairs_long(df, action = "Action", id = c("actor", "missing_id"), time = "time"),
    "not found"
  )
})

test_that(".build_pairs_long with multiple ID columns uses interaction grouping", {
  df <- data.frame(
    actor  = c(1, 1, 2, 2),
    group  = c("g1", "g1", "g1", "g1"),
    Action = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )
  result <- Nestimate:::.build_pairs_long(df, action = "Action",
                                          id = c("actor", "group"), time = "time")
  expect_s3_class(result, "data.frame")
  expect_true(all(c("from", "to") %in% names(result)))
  expect_equal(nrow(result), 2L)
})

test_that(".build_pairs_long with NULL id treats all rows as single sequence", {
  df <- data.frame(Action = c("A", "B", "A"), stringsAsFactors = FALSE)
  result <- Nestimate:::.build_pairs_long(df, action = "Action", id = NULL, time = "time")
  expect_equal(nrow(result), 2L)
  expect_equal(result$from, c("A", "B"))
  expect_equal(result$to, c("B", "A"))
})

test_that(".build_pairs_long with single-row group returns empty data.frame", {
  df <- data.frame(
    actor  = c(1, 2),
    Action = c("A", "B"),
    stringsAsFactors = FALSE
  )
  result <- Nestimate:::.build_pairs_long(df, action = "Action", id = "actor", time = "time")
  expect_equal(nrow(result), 0L)
  expect_true(all(c("from", "to") %in% names(result)))
})

test_that(".build_pairs_long sorts by time column when present", {
  df <- data.frame(
    actor  = c(1, 1, 1),
    time   = c(3, 1, 2),
    Action = c("C", "A", "B"),
    stringsAsFactors = FALSE
  )
  result <- Nestimate:::.build_pairs_long(df, action = "Action", id = "actor", time = "time")
  expect_equal(result$from, c("A", "B"))
  expect_equal(result$to, c("B", "C"))
})

# ---- .build_pairs_wide() coverage ----

test_that(".build_pairs_wide errors on missing columns", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  expect_error(
    Nestimate:::.build_pairs_wide(df, id = NULL, cols = c("T1", "T3")),
    "not found"
  )
})

test_that(".build_pairs_wide errors when fewer than 2 state columns", {
  df <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  expect_error(
    Nestimate:::.build_pairs_wide(df, id = NULL, cols = c("T1")),
    "At least 2"
  )
})

test_that(".build_pairs_wide with NULL cols uses all columns", {
  df <- data.frame(T1 = c("A", "B"), T2 = c("B", "A"), stringsAsFactors = FALSE)
  result <- Nestimate:::.build_pairs_wide(df, id = NULL, cols = NULL)
  expect_equal(nrow(result), 2L)
})

test_that(".build_pairs_wide with id excludes id from state cols", {
  df <- data.frame(id = c(1, 2), T1 = c("A", "B"), T2 = c("B", "A"),
                   stringsAsFactors = FALSE)
  result <- Nestimate:::.build_pairs_wide(df, id = "id", cols = NULL)
  # Only T1, T2 used, not id
  expect_true(all(result$from %in% c("A", "B")))
  expect_equal(nrow(result), 2L)
})

test_that(".build_pairs_wide row with all NAs returns empty data.frame", {
  df <- data.frame(T1 = c(NA, "A"), T2 = c(NA, "B"), stringsAsFactors = FALSE)
  result <- Nestimate:::.build_pairs_wide(df, id = NULL, cols = c("T1", "T2"))
  # First row has all NAs, contributes 0; second row A->B
  expect_equal(nrow(result), 1L)
  expect_equal(result$from, "A")
  expect_equal(result$to, "B")
})

# ---- .fmt_edgelist and .fmt_follows single-row edge cases ----

test_that("convert_sequence_format edgelist handles single-action sequences", {
  # Actor with only one action -> no transitions -> contributes 0 rows
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"),
    T2 = c("B", NA),
    stringsAsFactors = FALSE
  )
  result <- convert_sequence_format(wide_data, id_col = "id", format = "edgelist")
  # Actor 2 has only B (T2 is NA) -> 0 edges
  expect_equal(nrow(result), 1L)
  expect_equal(result$from, "A")
  expect_equal(result$to, "B")
})

test_that("convert_sequence_format follows handles single-action sequences", {
  wide_data <- data.frame(
    id = c(1, 2),
    T1 = c("A", "B"),
    T2 = c("B", NA),
    stringsAsFactors = FALSE
  )
  result <- convert_sequence_format(wide_data, id_col = "id", format = "follows")
  expect_equal(nrow(result), 1L)
  expect_equal(result$act, "B")
  expect_equal(result$follows, "A")
})
