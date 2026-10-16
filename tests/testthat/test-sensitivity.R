# cm_sensitivity() (v0.6 item 3): one input at a time, entered and unknown
# inputs together. Reference values: inst/validation/reference_v06_sensitivity.py
# (an independent implementation, built on reference_v041.py).

check_sens <- function(s, ref) {
  got <- s$inputs[match(ref$input, s$inputs$input), ]
  expect_equal(got$input, ref$input)
  expect_equal(got$low, ref$low, tolerance = 1e-8)
  expect_equal(got$high, ref$high, tolerance = 1e-8)
  expect_equal(got$total_low, ref$total_low, tolerance = 1e-8)
  expect_equal(got$total_high, ref$total_high, tolerance = 1e-8)
  expect_equal(got$rank_shift, as.integer(ref$rank_shift))
}

test_that("five diseases, yield: every entered and unknown input", {
  ref <- data.frame(
    input = c("associations:1", "associations:10", "associations:11", "associations:12", "associations:2",
      "associations:3", "associations:4", "associations:5", "associations:6", "associations:7",
      "associations:8", "associations:9", "impacts:1", "impacts:2", "impacts:3", "impacts:4",
      "impacts:5", "interaction:LAM;MAS", "interaction:LAM;MET", "interaction:LAM;RP",
      "interaction:LAM;SCK", "interaction:MAS;MET", "interaction:MAS;RP", "interaction:MAS;SCK",
      "interaction:MET;RP", "interaction:MET;SCK", "interaction:SCK;RP", "prob:LAM", "prob:MAS",
      "prob:MET", "prob:RP", "prob:SCK"),
    low = c(1.523431907, 0.8, 1.2, 1.089672122, 1.813333093, 1.7, 4.444909324, 1.303640693, 0.3032958536,
      0.06, 0.02951931524, 1.14375, 3.9, 2.459358591, 4.57238383, 2, 3.36, -1.142092503,
      -2.051800706, -1.334268368, -0.690384066, -1.142092503, -1.142092503, -0.690384066,
      -1.334268368, -0.690384066, -0.690384066, 0.2258400792, 0.2532976995, 0.07951931524, 0.048,
      0.2889901591),
    high = c(2.126777039, 1.2, 1.8, 2.064841299, 2.669118, 2.3, 8.099152845, 1.731920151, 0.3987094972, 0.1,
      0.07048068476, 1.715625, 5.7, 4.247952171, 6.604533661, 3, 5.04, 1.142092503, 2.051800706,
      1.334268368, 0.690384066, 1.142092503, 1.142092503, 0.690384066, 1.334268368, 0.690384066,
      0.690384066, 0.2752502837, 0.3526640095, 0.1204806848, 0.072, 0.4110098409),
    total_low = c(2.7374769391, 2.7258954580, 2.7249009231, 2.7133281669, 2.7341293969, 2.7611108593, 2.7308579309,
      2.7361749104, 2.7858991261, 2.7243143085, 2.7333048458, 2.7264219828, 2.5488705008,
      2.5564420965, 2.6472054929, 2.5737764431, 2.6808765169, 2.7987000239, 2.7721942440,
      2.7347798141, 2.7769032363, 2.7481588982, 2.7344149085, 2.7789167852, 2.7313689489,
      2.7399297519, 2.7294177810, 2.6297217671, 2.5836661981, 2.6452958809, 2.6909570913,
      2.6434775010),
    total_high = c(2.6898771848, 2.7028247406, 2.7039030347, 2.7133281669, 2.6930467856, 2.6742524172, 2.6964119103,
      2.6903748432, 2.6471469843, 2.7026139740, 2.6942363569, 2.7028991419, 2.8777858330,
      2.8902412828, 2.7779655809, 2.8528798907, 2.7457798169, 2.6279563099, 2.6544620898,
      2.6918765197, 2.6497530975, 2.6784974356, 2.6922414253, 2.6477395486, 2.6952873848,
      2.6867265819, 2.6972385528, 2.8032643657, 2.8648595191, 2.7847065780, 2.7370049952,
      2.7830712399),
    rank_shift = c(0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 2, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1),
    stringsAsFactors = FALSE)
  res <- deconflate(read_five("impacts_yield.csv"), n_draws = 0)
  s <- cm_sensitivity(res)
  expect_equal(nrow(s$inputs), 32)
  expect_setequal(s$inputs$input, ref$input)
  check_sens(s, ref)
  expect_equal(sum(s$inputs$kind == "unknown"), 10)
  expect_true(all(diff(s$inputs$swing) <= 0))
  expect_equal(s$inputs$range[s$inputs$input == "impacts:2"], "1 SD")
  expect_equal(s$inputs$range[s$inputs$input == "impacts:4"], "20%")
  expect_equal(s$inputs$range[s$inputs$input == "interaction:LAM;MAS"], "50% of smaller impact")
  # A positive interaction lowers the aggregate (crude rows already contain it).
  i <- s$inputs[s$inputs$input == "interaction:LAM;MAS", ]
  expect_lt(i$total_high, i$total_low)
  b <- s$by_disease[s$by_disease$input == "associations:6", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(4.0623500989, 2.3953504886, 4.4997733092, 1.5441366842, 2.6509204875), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(4.1409642151, 2.1957415971, 4.5390161793, 1.2224299147, 2.6842938072), tolerance = 1e-7)
  b <- s$by_disease[s$by_disease$input == "impacts:2", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(4.2053548637, 1.4094772786, 4.5625321808, 1.4863870406, 2.7217311559), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(3.9888587823, 3.2705521797, 4.4726262944, 1.2616666072, 2.6085518658), tolerance = 1e-7)
  b <- s$by_disease[s$by_disease$input == "interaction:LAM;MAS", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(4.4659045995, 2.6308025223, 4.52459816, 1.4039253073, 2.667072882), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(3.741298224, 1.9375674914, 4.5159544244, 1.3576109568, 2.6700005901), tolerance = 1e-7)
  expect_equal(s$central$first, "LAM")
  expect_null(s$variance)
  # Printing: the 10 largest swings by default, all with n = Inf.
  out <- capture.output(print(s))
  expect_true(any(grepl("32 inputs varied one at a time (22 entered, 10 unknown); the 10 largest swings",
                        out, fixed = TRUE)))
  expect_true(any(grepl("22 more: print(x, n = Inf)", out, fixed = TRUE)))
  expect_false(any(grepl("more: print", capture.output(print(s, n = Inf)), fixed = TRUE)))
  # Only some inputs, or only entered or unknown ones.
  expect_equal(nrow(cm_sensitivity(res, inputs = "unknown")$inputs), 10)
  expect_equal(cm_sensitivity(res, keys = "impacts:2")$inputs$total_low, ref$total_low[ref$input == "impacts:2"],
               tolerance = 1e-8)
  # A larger relative change for fixed inputs.
  s3 <- cm_sensitivity(res, keys = "impacts:4", change = 0.3)
  expect_equal(c(s3$inputs$low, s3$inputs$high), c(1.75, 3.25))
  expect_error(cm_sensitivity(res, change = 1.5), "`change`")
})

test_that("an unknown pair and a three-way input on request", {
  ref <- data.frame(
    input = c("associations:5", "interaction:MAS;RP", "or_ratio:MAS:RP|MET", "pair:MAS:RP", "prob:RP"),
    low = c(1.303640693, -1.165722591, 0.5, 0.5, 0.048),
    high = c(1.731920151, 1.165722591, 2, 3, 0.072),
    total_low = c(2.7497650550, 2.7474092703, 2.7380197878, 2.7661960726, 2.7021235722),
    total_high = c(2.7051051669, 2.7075506979, 2.7161790700, 2.6788683505, 2.7541980233),
    rank_shift = c(0, 0, 0, 0, 0),
    stringsAsFactors = FALSE)
  f <- function(x) extdata("five_diseases", x)
  a <- utils::read.csv(f("associations.csv"), stringsAsFactors = FALSE)
  a <- a[!(a$disease == "MAS" & a$associated_with == "RP"), ]
  m <- suppressMessages(cm_read_inputs(diseases = f("diseases.csv"), associations = a,
                                       impacts = f("impacts_yield.csv")))
  res <- deconflate(m, n_draws = 0)
  s <- cm_sensitivity(res, three_way = list(c("MAS", "RP", "MET")))
  check_sens(s, ref)
  expect_equal(sum(startsWith(s$inputs$input, "pair:")), 1)
  expect_equal(s$inputs$range[s$inputs$input == "pair:MAS:RP"], "OR 0.5-3")
  # Without the request there is no three-way input; a range for one adds it.
  expect_false(any(startsWith(cm_sensitivity(res, inputs = "unknown")$inputs$input, "or_ratio:")))
  s2 <- cm_sensitivity(res, inputs = "unknown", keys = "or_ratio:MAS:RP|MET",
                       ranges = data.frame(input = "or_ratio:MAS:RP|MET", low = 0.5, high = 2))
  expect_equal(s2$inputs$total_low, ref$total_low[ref$input == "or_ratio:MAS:RP|MET"], tolerance = 1e-8)
})

test_that("event impacts: the attributable risk", {
  ref <- data.frame(
    input = c("associations:1", "impacts:2", "impacts:6", "interaction:MAS;SCK", "prob:SCK", "risk"),
    low = c(1.523431907, 1.28, 2.651125738, 0.67, 0.2889901591, 0.2),
    high = c(2.126777039, 1.92, 3.862510123, 1.5, 0.4110098409, 0.3),
    total_low = c(0.0566747952, 0.0346313819, 0.0699610480, 0.0612290519, 0.0521847390, 0.0460563230),
    total_high = c(0.0512502159, 0.0735413268, 0.0230942168, 0.0480079240, 0.0559435004, 0.0611045322),
    rank_shift = c(1, 2, 3, 1, 0, 0),
    stringsAsFactors = FALSE)
  res <- deconflate(read_five("culling.csv"), event_model = TRUE, overall_risk = 0.25, n_draws = 0)
  s <- cm_sensitivity(res, keys = ref$input)
  check_sens(s, ref)
  expect_equal(s$metric, "attributable risk")
  b <- s$by_disease[s$by_disease$input == "associations:1", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(1.1435840429, 1.2850012134, 1.07939418, 1.1753285705, 1.5079654093), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(0.9686719157, 1.1753065228, 1.0815728713, 1.1754736328, 1.5109634587), tolerance = 1e-7)
  b <- s$by_disease[s$by_disease$input == "impacts:2", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(0.7836219658, 0.7908043991, 1.093726854, 1.1939858574, 1.5211664044), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(1.4188766301, 1.8094723829, 1.0692259503, 1.1607983388, 1.5051629259), tolerance = 1e-7)
  b <- s$by_disease[s$by_disease$input == "impacts:6", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(1.5887339043, 1.5923944549, 1.081391554, 1.1933605254, 1.5015804707), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(0.3796246009, 0.7347592958, 1.0924912003, 1.1747437965, 1.5450727819), tolerance = 1e-7)
  b <- s$by_disease[s$by_disease$input == "interaction:MAS;SCK", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(1.0313850922, 1.4153437133, 1.0829187652, 1.349194821, 1.509611151), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(1.1096262842, 1.0872341619, 1.0782286718, 1.0187415687, 1.5108515498), tolerance = 1e-7)
  b <- s$by_disease[s$by_disease$input == "prob:SCK", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(1.0536329028, 1.2267289676, 1.0809874904, 1.1774616015, 1.5096580714), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(1.0668020417, 1.2401401896, 1.0807316305, 1.1717170112, 1.5090049964), tolerance = 1e-7)
  b <- s$by_disease[s$by_disease$input == "risk", ]
  expect_equal(b$disease, c("LAM", "MAS", "MET", "SCK", "RP"))
  expect_equal(b$adjusted_low, c(1.0622490761, 1.2338662246, 1.0496800172, 1.1804229425, 1.6114684064), tolerance = 1e-7)
  expect_equal(b$adjusted_high, c(1.0586820329, 1.2339152589, 1.1079231168, 1.1698545053, 1.444891139), tolerance = 1e-7)
  # The first-ranked disease changes when the MAS hazard ratio is low.
  expect_equal(s$inputs$first[s$inputs$input == "impacts:2"], "SCK (low)")
  expect_true(any(grepl("First place changes", capture.output(print(s)))))
})

test_that("variance shares and rank probabilities from the draws", {
  res <- deconflate(read_five("impacts_yield.csv"), n_draws = 60, seed = 3)
  # Rank probabilities: each disease's places over the accepted draws.
  r <- res$ranks
  v <- res$draws$values
  C <- v[, paste0("contribution:", r$disease)]
  first <- r$disease[apply(C, 1, which.max)]
  expect_equal(r$rank_1, as.numeric(table(factor(first, levels = r$disease))) / nrow(C))
  expect_equal(unname(rowSums(r[, grep("^rank_", names(r))])), rep(1, nrow(r)))
  expect_equal(sum(r$mean_rank), sum(seq_len(nrow(r))))
  expect_true(any(grepl("Rank by contribution over the draws", capture.output(print(summary(res))))))
  # Variance shares: squared standardised coefficients of a linear fit of
  # the total on the normal scores of the drawn inputs.
  s <- cm_sensitivity(res, inputs = "entered")
  U <- res$draws$uniforms
  fixed <- vapply(colnames(U), function(k) identical(res$draws$specs[[k]]$type, "fixed"), logical(1))
  Z <- stats::qnorm(U[, !fixed, drop = FALSE])
  y <- v[, "total"]
  fit <- stats::lm(y ~ Z)
  src2 <- (stats::coef(fit)[-1] * apply(Z, 2, stats::sd) / stats::sd(y))^2
  names(src2) <- colnames(Z)
  expect_equal(s$r_squared, summary(fit)$r.squared, tolerance = 1e-10)
  expect_equal(s$inputs$share[s$inputs$input == "impacts:2"], unname(src2["impacts:2"]), tolerance = 1e-10)
  expect_equal(sum(s$variance$share), sum(src2), tolerance = 1e-10)
  # Fixed inputs have no share.
  expect_true(is.na(s$inputs$share[s$inputs$input == "impacts:4"]))
  expect_true(any(grepl("share of the variance", capture.output(print(s)))))
  # Without draws: no shares and no ranks.
  r0 <- deconflate(read_five("impacts_yield.csv"), n_draws = 0)
  expect_null(r0$ranks)
  expect_null(cm_sensitivity(r0, keys = "impacts:1")$variance)
})

test_that("cm_threshold() takes a result and the same input keys", {
  m <- example_supplement()
  res <- deconflate(m, n_draws = 0)
  a <- cm_threshold(res, "pair:d2:d3", c(0.1, 100), conclusion = "rank", n_grid = 21)
  b <- cm_threshold(m, "pair:d2:d3", c(0.1, 100), conclusion = "rank", n_grid = 21)
  expect_equal(a$thresholds, b$thresholds)
  expect_error(cm_threshold(res, "pair:d2:d3", c(0.1, 100), overall_risk = 0.2), "settings come from the result")
  # Keys of earlier versions say what to write.
  expect_error(cm_threshold(m, "assoc:d2:d3", c(0.1, 100)), "write 'pair:<d1>:<d2>'", fixed = TRUE)
  expect_error(cm_threshold(m, "inter:d1:d2", c(0, 1)), "write 'interaction:<d1>;<d2>'", fixed = TRUE)
  expect_error(cm_threshold(m, "three:d1:d2:d3", c(0.5, 2)), "or_ratio")
  expect_error(cm_threshold(m, "impact:d1", c(0, 1)), "write 'impacts:<row>'", fixed = TRUE)
  expect_error(cm_threshold(m, "foo:d1", c(0, 1)), "Unknown input")
  # Event impacts: the settings (here the baseline risk) come from the result.
  ev <- deconflate(read_five("culling.csv"), event_model = TRUE, baseline_risk = 0.18, n_draws = 0)
  th <- cm_threshold(ev, "risk", c(0.05, 0.6), conclusion = "total", target = 0.05, n_grid = 11)
  th2 <- cm_threshold(read_five("culling.csv"), "risk", c(0.05, 0.6), conclusion = "total", target = 0.05,
                      event_model = TRUE, baseline_risk = 0.18, n_grid = 11)
  expect_equal(th$thresholds$threshold, th2$thresholds$threshold)
})

test_that("the result keeps what cm_sensitivity() needs", {
  res <- deconflate(example_supplement(), n_draws = 0)
  expect_true(all(c("plan", "risk", "dots") %in% names(res$settings)))
  old <- res
  old$settings <- NULL
  expect_error(cm_sensitivity(old), "does not carry the settings")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  expect_invisible(plot(cm_sensitivity(res)))
})
