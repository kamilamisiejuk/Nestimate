# Regression tests for FIX6 (A11 data-prep / conversion / wtna findings).
#
# A11-F01 (CRITICAL) convert_sequence_format(id_col=NULL) on wide V1,V2,V3
#                     data must NOT consume the first state column as id.
# A11-F02 (MEDIUM)    long_to_wide(fill_na=) must actually change output.
# A11-F03 (MEDIUM)    wtna(window_size=) must reject invalid values.
# A11-F04 (LOW)       roxygen window_size default doc says 3 (no behavior
#                     change) — covered by the unchanged-default assertion.

# ---------------------------------------------------------------------------
# A11-F01 : convert_sequence_format wide no-id must not drop first time point
# ---------------------------------------------------------------------------

test_that("A11-F01: convert_sequence_format @example produces correct counts", {
  # The function's own documented @example (V1,V2,V3 wide, no id column).
  seqs <- data.frame(
    V1 = c("A", "B", "A"), V2 = c("B", "A", "C"), V3 = c("A", "C", "B"),
    stringsAsFactors = FALSE
  )

  fq <- convert_sequence_format(seqs, format = "frequency")
  # Row 1 = A,B,A -> A appears twice, B once, C zero. The first 'A' (V1)
  # was previously silently dropped (counted A=1).
  expect_equal(fq$A[1], 2L)
  expect_equal(fq$B[1], 1L)
  expect_equal(fq$C[1], 0L)
  expect_equal(fq$A[2], 1L)
  expect_equal(fq$C[2], 1L)
  # No state column may be emitted as a pseudo-id.
  expect_false("V1" %in% names(fq))

  el <- convert_sequence_format(seqs, format = "edgelist")
  # 3 sequences x 2 transitions each = 6 edges (was 3: every sequence lost
  # its first transition).
  expect_equal(nrow(el), 6L)
})

test_that("A11-F01: frequencies() and convert_sequence_format() agree on no-id wide data", {
  seqs <- data.frame(
    V1 = c("A", "B", "A"), V2 = c("B", "A", "C"), V3 = c("A", "C", "B"),
    stringsAsFactors = FALSE
  )

  fr <- unclass(frequencies(seqs, format = "wide"))
  el <- convert_sequence_format(seqs, format = "edgelist")

  states <- sort(unique(c(el$from, el$to)))
  em <- table(factor(el$from, levels = states),
              factor(el$to, levels = states))
  em <- matrix(as.integer(em), nrow = length(states),
               dimnames = list(states, states))

  fr_m <- matrix(as.integer(fr), nrow = nrow(fr),
                 dimnames = dimnames(fr))
  # Two exported functions must produce the identical transition matrix on
  # the same canonical wide input.
  expect_equal(em, fr_m)
})

test_that("A11-F01: agreement holds on simulate_sequences() canonical shape", {
  seqs <- simulate_sequences(n_actors = 6, n_states = 4,
                             seq_length = 8, seed = 42)
  expect_true(all(grepl("^T[0-9]+$", names(seqs))))  # no id column

  fr <- unclass(frequencies(seqs, format = "wide"))
  el <- convert_sequence_format(seqs, format = "edgelist")
  # Total edges = sum of all transition counts in frequencies().
  expect_equal(nrow(el), sum(fr))

  states <- sort(unique(c(el$from, el$to)))
  em <- table(factor(el$from, levels = states),
              factor(el$to, levels = states))
  em <- matrix(as.integer(em), nrow = length(states),
               dimnames = list(states, states))
  fr_m <- matrix(as.integer(fr), nrow = nrow(fr), dimnames = dimnames(fr))
  expect_equal(em, fr_m)
})

test_that("A11-F01: a genuine leading id is used only when passed explicitly (Codex)", {
  # Codex P2: value-based id inference is unsound (it drops first events
  # when a start-only state does not recur). The contract is now: id_col is
  # NEVER guessed from data. A genuine leading id must be passed explicitly.
  wide_data <- data.frame(
    Student = c(101, 202), T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )

  # id_col = NULL: no guessing -> Student is a state, nothing dropped.
  el_noid <- convert_sequence_format(wide_data, format = "edgelist")
  expect_equal(nrow(el_noid), 4L)                       # 2 rows x 2 transitions
  expect_true(any(c("101", "202") %in% c(el_noid$from, el_noid$to)))

  # Explicit id_col: Student is the identifier, never counted as a state.
  result <- convert_sequence_format(wide_data, id_col = "Student",
                                    format = "frequency")
  expect_true("Student" %in% names(result))
  expect_equal(sort(result$Student), c(101, 202))
  expect_false("101" %in% names(result))
  expect_false("202" %in% names(result))
  expect_equal(result$A[result$Student == 101], 1L)
  expect_equal(result$B[result$Student == 101], 1L)
})

test_that("A11-F01: explicit id_col on wide data is unchanged", {
  wide_data <- data.frame(
    id = c(1, 2), T1 = c("A", "B"), T2 = c("B", "A"),
    stringsAsFactors = FALSE
  )
  r1 <- convert_sequence_format(wide_data, id_col = "id", format = "frequency")
  expect_true("id" %in% names(r1))
  expect_equal(r1$A[r1$id == 1], 1L)

  # seq_cols pinned, id_col = NULL: no column is guessed as the id (Codex).
  # A row-index id is synthesized; only the pinned seq_cols are states, so
  # neither 'grp' nor 'extra'/'Z' is treated as an id or a state.
  wide2 <- data.frame(
    grp = c("g1", "g2"), T1 = c("A", "B"), T2 = c("B", "A"),
    extra = c("Z", "Z"), stringsAsFactors = FALSE
  )
  r2 <- convert_sequence_format(wide2, seq_cols = c("T1", "T2"),
                                format = "frequency")
  expect_false("grp" %in% names(r2))   # not auto-detected as id
  expect_false("Z" %in% names(r2))     # extra not counted as a state
})

test_that("A11-F01: long-format path unaffected by the wide-id fix", {
  long <- data.frame(
    Actor = rep(1:2, each = 3), Time = rep(1:3, 2),
    Action = c("A", "B", "C", "B", "A", "C"),
    stringsAsFactors = FALSE
  )
  el <- convert_sequence_format(long, action = "Action", id_col = "Actor",
                                time = "Time", format = "edgelist")
  expect_true(all(c("Actor", "from", "to") %in% names(el)))
  expect_equal(nrow(el), 4L)  # 2 sequences x 2 transitions
})

# ---------------------------------------------------------------------------
# A11-F02 : long_to_wide(fill_na=) must control NA padding of ragged rows
# ---------------------------------------------------------------------------

test_that("A11-F02: long_to_wide fill_na TRUE vs FALSE differ on ragged data", {
  long <- data.frame(
    Actor = c("a", "a", "a", "b", "b"),
    Time = c(1, 2, 3, 1, 2),
    Action = c("X", "Y", "Z", "P", "Q"),
    stringsAsFactors = FALSE
  )

  w_true <- long_to_wide(long, id_col = "Actor", fill_na = TRUE)
  w_false <- long_to_wide(long, id_col = "Actor", fill_na = FALSE)

  # The argument must now have an observable effect.
  expect_false(identical(w_true, w_false))

  # TRUE pads the short sequence's missing time point with NA.
  expect_true("V3" %in% names(w_true))
  expect_true(any(is.na(w_true$V3)))

  # FALSE does not fill missing time points with NA -> ragged tail dropped,
  # no NA-padded cells emitted.
  v_false <- grep("^V[0-9]+$", names(w_false), value = TRUE)
  expect_false("V3" %in% names(w_false))
  expect_false(any(is.na(as.matrix(w_false[, v_false, drop = FALSE]))))
})

test_that("A11-F02: fill_na is a no-op when all sequences are equal length", {
  long <- data.frame(
    Actor = c("a", "a", "b", "b"),
    Time = c(1, 2, 1, 2),
    Action = c("X", "Y", "P", "Q"),
    stringsAsFactors = FALSE
  )
  expect_identical(
    long_to_wide(long, id_col = "Actor", fill_na = TRUE),
    long_to_wide(long, id_col = "Actor", fill_na = FALSE)
  )
})

test_that("A11-F02: fill_na on group_regulation_long ragged sequences differs", {
  d <- group_regulation_long
  # Cap to a few actors for speed; sequences are naturally ragged.
  ids <- head(unique(d$Actor), 8)
  sub <- d[d$Actor %in% ids, c("Actor", "Time", "Action")]

  w_true <- long_to_wide(sub, id_col = "Actor", time_col = "Time",
                         action_col = "Action", fill_na = TRUE)
  w_false <- long_to_wide(sub, id_col = "Actor", time_col = "Time",
                          action_col = "Action", fill_na = FALSE)

  expect_false(identical(w_true, w_false))
  expect_gt(ncol(w_true), ncol(w_false))  # padding adds trailing V-columns
  v_false <- grep("^V[0-9]+$", names(w_false), value = TRUE)
  expect_false(any(is.na(as.matrix(w_false[, v_false, drop = FALSE]))))
})

test_that("A11-F02: fill_na is validated as logical scalar", {
  long <- data.frame(
    Actor = c("a", "a"), Time = c(1, 2), Action = c("X", "Y"),
    stringsAsFactors = FALSE
  )
  expect_error(long_to_wide(long, id_col = "Actor", fill_na = "yes"))
  expect_error(long_to_wide(long, id_col = "Actor", fill_na = c(TRUE, FALSE)))
})

# ---------------------------------------------------------------------------
# A11-F03 : wtna(window_size=) must reject invalid values cleanly
# ---------------------------------------------------------------------------

test_that("A11-F03: wtna rejects invalid window_size", {
  oh <- data.frame(
    A = c(1, 0, 0, 1, 0, 0, 1, 0, 0),
    B = c(0, 1, 0, 0, 1, 0, 0, 1, 0),
    C = c(0, 0, 1, 0, 0, 1, 0, 0, 1)
  )
  # Non-integer previously produced a distinct garbage matrix.
  expect_error(wtna(oh, window_size = 2.7))
  # Zero / negative previously fell silently to the non-windowed branch.
  expect_error(wtna(oh, window_size = 0))
  expect_error(wtna(oh, window_size = -5))
  expect_error(wtna(oh, window_size = NA))
  expect_error(wtna(oh, window_size = c(2, 3)))
  expect_error(wtna(oh, window_size = "x"))
})

test_that("A11-F03: valid window_size values still work and store integer", {
  oh <- data.frame(
    A = c(1, 0, 0, 1, 0, 0, 1, 0, 0),
    B = c(0, 1, 0, 0, 1, 0, 0, 1, 0),
    C = c(0, 0, 1, 0, 0, 1, 0, 0, 1)
  )
  expect_s3_class(wtna(oh, window_size = 1), "netobject")
  expect_s3_class(wtna(oh, window_size = 2L), "netobject")
  expect_s3_class(wtna(oh, window_size = 4), "netobject")

  net2 <- wtna(oh, window_size = 2)
  expect_type(net2$params$window_size, "integer")
  expect_equal(net2$params$window_size, 2L)

  # Numeric whole number (2.0) is accepted (it equals its integer cast).
  expect_s3_class(wtna(oh, window_size = 2.0), "netobject")
})

test_that("A11-F03: window_size guard also covers method='both' / actor path", {
  oh <- data.frame(
    A = c(1, 0, 0, 1, 0, 0, 1, 0, 0),
    B = c(0, 1, 0, 0, 1, 0, 0, 1, 0),
    C = c(0, 0, 1, 0, 0, 1, 0, 0, 1),
    actor = rep(1:3, each = 3)
  )
  expect_error(wtna(oh, method = "both", window_size = 2.7, actor = "actor"))
  expect_error(wtna(oh, method = "both", window_size = 0, actor = "actor"))
  expect_s3_class(
    wtna(oh, method = "both", window_size = 2, actor = "actor"),
    "wtna_mixed"
  )
})

# ---------------------------------------------------------------------------
# A11-F04 : default window_size value is unchanged (doc-only correction)
# ---------------------------------------------------------------------------

test_that("A11-F04: window_size default value remains 3 (behavior unchanged)", {
  expect_equal(as.integer(formals(wtna)$window_size), 3L)
  expect_equal(as.integer(formals(prepare_onehot)$window_size), 3L)
})
