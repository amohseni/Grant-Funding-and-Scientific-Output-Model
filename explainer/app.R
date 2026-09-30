# Funding the Gap: an interactive guide to
# "A model of optimal science funding: targeting the capability-resource gap"
# (Mohseni, DeDeo, Zollman). Presentation layer only. The analytic parts (production
# function, gap rule, the two heuristics, the certain-grant-versus-lottery comparison)
# are computed live in base R; every simulation result is read from the paper's own
# figure data in data/ (produced by the model repository's analysis scripts).

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(ggplot2)
})

# ------------------------------------------------------------------ model pieces
A_CONST <- 1
lam <- function(K, R) A_CONST * K * R / (K + R)            # expected research output
rpareto <- function(n, shape, xmin = 1) xmin * (1 - runif(n))^(-1 / shape)
mean_pareto <- function(shape, xmin = 1) xmin * shape / (shape - 1)

draw_population <- function(n, alpha, seed) {
  set.seed(seed)
  data.frame(K = rpareto(n, alpha), R0 = rpareto(n, alpha))
}

gap_rule <- function(K, R0, B) {
  if (B <= 0) return(rep(0, length(K)))
  spend <- function(c) sum(pmax(c * K - R0, 0)) - B
  hi <- (B + sum(R0)) / min(K) + 1
  c <- uniroot(spend, c(0, hi), tol = 1e-10)$root
  list(g = pmax(c * K - R0, 0), c = c)
}

# ------------------------------------------------------------------ data
rd <- function(f) read.table(file.path("data", f), header = TRUE)
D_review  <- rd("fig4-value.dat")
D_rounds  <- rd("fig5-rounds.dat")
D_seed    <- rd("fig10-floor-cost.dat")
D_lot_h   <- rd("fig11-heavy.dat")
D_lot_e   <- rd("fig11-even.dat")
D_grid    <- rd("fig16-grid.dat")
read_contours <- function(f) {
  d <- read.table(file.path("data", f), header = TRUE, blank.lines.skip = FALSE, fill = TRUE)
  d$seg <- cumsum(is.na(d$x))
  d <- d[!is.na(d$x), ]
  d$b <- 10^d$x
  d
}
C_targ <- read_contours("fig16-targ.dat")
C_rev  <- read_contours("fig16-rev.dat")
gini_of_alpha <- function(a) 1 / (2 * a - 1)
FIELDS <- c("heavy-tailed" = 1.3, "intermediate" = 2, "evenly spread" = 3.5)

# bilinear interpolation on the regime grid (log b x gini)
grid_value <- function(b, gini, col) {
  lb <- log10(b); lbs <- sort(unique(D_grid$logb)); gs <- sort(unique(D_grid$gini))
  lb <- min(max(lb, min(lbs)), max(lbs)); gini <- min(max(gini, min(gs)), max(gs))
  i <- max(which(lbs <= lb)); i2 <- min(i + 1, length(lbs)); j <- max(which(gs <= gini)); j2 <- min(j + 1, length(gs))
  v <- function(ii, jj) D_grid[[col]][abs(D_grid$logb - lbs[ii]) < 1e-6 & abs(D_grid$gini - gs[jj]) < 1e-6][1]
  tx <- if (i2 == i) 0 else (lb - lbs[i]) / (lbs[i2] - lbs[i])
  ty <- if (j2 == j) 0 else (gini - gs[j]) / (gs[j2] - gs[j])
  (1 - tx) * (1 - ty) * v(i, j) + tx * (1 - ty) * v(i2, j) + (1 - tx) * ty * v(i, j2) + tx * ty * v(i2, j2)
}

# ------------------------------------------------------------------ plotting
th <- function(base = 13) {
  theme_minimal(base_size = base, base_family = "") +
    theme(panel.grid = element_blank(),
          axis.line = element_line(colour = "black", linewidth = 0.35),
          axis.ticks = element_line(colour = "black", linewidth = 0.35),
          axis.title = element_text(size = base - 1),
          plot.title = element_text(size = base, face = "plain", hjust = 0, lineheight = 1.05),
          plot.subtitle = element_text(size = base - 2, colour = "grey30"),
          strip.text = element_text(size = base - 2, hjust = 0),
          strip.background = element_blank(),
          legend.position = "none",
          plot.margin = margin(8, 24, 8, 8))
}
GREY <- "grey45"

# ------------------------------------------------------------------ text helpers
p <- function(...) tags$p(...)
lead <- function(...) tags$p(class = "lead-text", ...)
notice <- function(...) tags$div(class = "notice", tags$span(class = "notice-label", "What to notice"), tags$p(...))
why <- function(...) tags$div(class = "why", tags$span(class = "notice-label", "Why"), tags$p(...))
side <- function(...) tags$div(class = "side", ...)
# Figures render at a fixed pixel size chosen for their content (a fixed aspect ratio and a width that keeps
# plot text in proportion to the page text); CSS scales the image down, never up, on narrow screens.
fig <- function(id, w, h) tags$div(class = "fig", style = sprintf("max-width:%dpx", w),
                                   tags$div(style = sprintf("width:100%%; aspect-ratio:%d/%d", w, h),
                                            plotOutput(id, width = "100%", height = "100%")))

css <- "
html { font-size: 18px; }
body { font-family: Georgia, 'Times New Roman', serif; color: #111; margin: 10pt; background: #fff; }
/* header: title on its own line, navigation on the line below, one hairline under both */
.navbar { border-bottom: 1px solid #111; padding: 1.4rem 0 0 0; margin-bottom: 1.6rem; background: #fff !important; }
.navbar > .container-fluid { flex-direction: column; align-items: flex-start; padding: 0 0.75rem; }
.navbar-brand { white-space: normal; padding: 0; margin: 0 !important; }
.brand { display: flex; flex-direction: column; }
.brand-title { font-size: 2.05rem; line-height: 1.2; font-weight: 400; letter-spacing: -0.005em; color: #111; max-width: 34ch; }
.brand-sub { font-size: 1rem; color: #666; margin-top: 0.45rem; font-style: italic; }
.navbar-nav { flex-direction: row !important; flex-wrap: wrap; gap: 0.35rem 1.6rem; margin: 0 0 0.9rem 0; }
.navbar .nav-link { padding: 0.15rem 0 !important; font-size: 0.95rem; color: #666 !important; border-bottom: 1px solid transparent; }
.navbar .nav-link:hover { color: #111 !important; }
.navbar .nav-link.active { color: #111 !important; border-bottom: 1px solid #111; }
.navbar-toggler, .navbar-toggle { display: none !important; }
.bslib-page-navbar > .container-fluid, .tab-content > .container-fluid, .container-fluid { border-top: none !important; }
.navbar-collapse { display: flex !important; padding-top: 25pt; }
/* figures: fixed proportions, capped width, scale down only */
.fig { width: 100%; margin: 0.2rem 0 1.2rem 0; }
.fig .shiny-plot-output img { display: block; }
@media (min-width: 992px) { .fig { margin-top: 0.4rem; } }
/* narrow screens: the figure comes first, then the text and controls */
@media (max-width: 991px) { .bslib-grid > .bslib-grid-item:nth-child(2) { order: -1; } .side { margin-top: 0.6rem; } }
/* text */
.side, .main-text, .readout, .notice, .why { max-width: 60ch; }
.lead-text { font-size: 1.12rem; line-height: 1.55; }
.side p, .main-text p { line-height: 1.55; }
.notice, .why { border-left: 1px solid #111; padding: 0.2rem 0 0.2rem 0.9rem; margin: 1.1rem 0; }
.notice-label { display: block; font-size: 0.72rem; letter-spacing: 0.1em; text-transform: uppercase; color: #666; margin-bottom: 0.25rem; }
.readout { font-size: 1rem; line-height: 1.5; border-top: 1px solid #111; padding: 0.7rem 0 0.2rem 0; margin-top: 1rem; }
.readout .num { font-weight: 600; }
.control-label, label { font-size: 0.92rem; color: #333; }
.small-note { font-size: 0.85rem; color: #666; line-height: 1.45; }
h2 { font-weight: 400; font-size: 1.5rem; margin: 0.2rem 0 0.7rem 0; }
h3 { font-weight: 400; font-size: 1.15rem; margin-top: 1.3rem; }
.takeaway { margin: 0.7rem 0; padding-left: 0.9rem; border-left: 1px solid #999; }
/* controls, kept quiet */
.irs--shiny .irs-bar, .irs--shiny .irs-single, .irs--shiny .irs-from, .irs--shiny .irs-to { background: #111; border-color: #111; }
.irs--shiny .irs-handle { border-color: #111; }
.irs--shiny .irs-min, .irs--shiny .irs-max, .irs--shiny .irs-grid-text { color: #888; }
.btn-outline-dark { border-radius: 0; }
.form-control { border-radius: 0; border-color: #bbb; }
.form-check-input:checked { background-color: #111; border-color: #111; }
"

# ================================================================== UI
ui <- page_navbar(
  title = tags$div(class = "brand", tags$span(class = "brand-title", "A Model of Optimal Science Funding: Targeting the Capability-Resource Gap"), tags$span(class = "brand-sub", "An interactive guide to the paper by Mohseni, DeDeo, and Zollman")),
  theme = bs_theme(version = 5, bg = "#ffffff", fg = "#111111", primary = "#111111",
                   base_font = font_collection("Georgia", "Times New Roman", "serif")),
  header = tags$head(tags$style(HTML(css))),
  fillable = FALSE,

  # ---------------------------------------------------------------- 1 start
  nav_panel("Start",
    layout_columns(col_widths = c(7, 5),
      div(class = "main-text",
        h2("What this guide covers"),
        lead("A research funder divides a fixed budget among researchers whose abilities and needs it cannot observe.
              The paper asks how it should do so, and what that implies for peer review, seed grants, and lotteries."),
        p("The paper's main result is that the allocation that produces the most research output gives each researcher
           the gap between their capability and their resources. What a funder should do about peer review, seed grants,
           and lotteries then depends on two features of its field: how large the funder's budget is relative to the
           resources researchers already have, and how unequal capability is among them."),
        p("This guide presents that result in seven steps. Each step has one figure you can adjust and a short explanation
           of what the figure shows. The figures on the model and the gap rule are computed as you move the controls.
           The figures on review, seed grants, lotteries, and the map of funders use the simulation results reported in
           the paper."),
        tags$ol(class = "main-text",
          tags$li("The model"),
          tags$li("The gap rule"),
          tags$li("What track records show"),
          tags$li("What peer review is worth"),
          tags$li("Seed grants and lotteries"),
          tags$li("Where your program stands"),
          tags$li("Summary and limits")),
        p(class = "small-note", "Mohseni, DeDeo, and Zollman, A model of optimal science funding: targeting the capability-resource gap (manuscript, 2026).
                                 The code and data are in the project repository, linked on the last page.")),
      div())),

  # ---------------------------------------------------------------- 2 model
  nav_panel("1 The model",
    layout_columns(col_widths = c(4, 8),
      side(
        h2("The model"),
        p("Each researcher has a capability ", em("K"), ", what they know and can do, and resources ", em("R"),
          ", what they have to work with. A researcher's expected research output is proportional to the harmonic mean of
           the two. Output rises when either input rises, but it is limited by whichever input is scarcer. A grant adds to
           a researcher's resources in the round it is given. Doing research raises the researcher's capability in later
           rounds."),
        sliderInput("m_K", "Capability K", min = 1, max = 30, value = 10, step = 0.5),
        sliderInput("m_R", "Resources R", min = 0.5, max = 30, value = 3, step = 0.5),
        notice("Hold capability fixed and raise resources. Output rises quickly at first and then levels off as it
                approaches the ceiling set by capability. Now hold resources at a small value and raise capability.
                Output barely changes, because resources are the scarce input."),
        uiOutput("m_readout")),
      div(class = "figcol", fig("m_plot", 640, 420)))),

  # ---------------------------------------------------------------- 3 gap rule
  nav_panel("2 The gap rule",
    layout_columns(col_widths = c(4, 8),
      side(
        h2("The gap rule"),
        p("Suppose the funder knows every researcher's capability and resources. The allocation that produces the most
           research output this round gives each researcher ", em("g = max(cK − R, 0)"), ". Each researcher has a target
           level of resources equal to ", em("c"), " times their capability. A researcher below the target receives the
           difference; a researcher at or above the target receives nothing. The budget determines ", em("c"), ".
           Grants are larger for researchers with more capability and smaller for researchers with more resources."),
        sliderInput("g_b", "Budget, as a fraction of the field's total baseline resources", min = 0.02, max = 2, value = 0.1, step = 0.02),
        sliderInput("g_alpha", "Capability inequality (Pareto tail; smaller is more unequal)", min = 1.2, max = 5, value = 2, step = 0.1),
        div(class = "d-flex gap-2 align-items-end",
            numericInput("g_seed", "Population seed", value = 7, min = 1, max = 9999, width = "9rem"),
            actionButton("g_redraw", "New population", class = "btn btn-outline-dark btn-sm mb-3")),
        radioButtons("g_show", "Show grants from", choices = c("the gap rule", "funding the track record", "funding the under-resourced"), selected = "the gap rule"),
        notice("The line is the funding frontier, R = cK. Every researcher to its left is funded and moved onto the line.
                Every researcher to its right receives nothing. As the budget rises, the line rotates toward the horizontal
                and more researchers are funded. The two other rules are common alternatives. Funding the track record
                gives grants to the most productive researchers, many of whom already have ample resources. Funding the
                under-resourced gives grants to the researchers with the fewest resources, many of whom have low capability.
                Either rule can produce less research output than dividing the budget equally."),
        uiOutput("g_readout")),
      div(class = "figcol", fig("g_plot", 640, 480)))),

  # ---------------------------------------------------------------- 4 records
  nav_panel("3 What track records show",
    layout_columns(col_widths = c(4, 8),
      side(
        h2("What track records show"),
        p("A real funder observes research output, not capability or resources. Expected output is a single number that
           depends on both inputs, so two very different researchers can have the same record: one with high capability
           and few resources, another with modest capability and ample resources. The gap rule would treat them very
           differently. The first may have the largest gap in the field; the second may have no gap at all."),
        sliderInput("r_Ks", "Capabilities to compare", min = 2, max = 40, value = c(3, 30), step = 1),
        sliderInput("r_y", "An observed output level", min = 0.5, max = 8, value = 2.5, step = 0.1),
        notice("Upper figure: when resources are small relative to capability, researchers of very different capability
                produce nearly the same output, so a small grant tells the funder little about capability. Lower figure:
                every point on the curve produces the same expected output. A record alone does not say where on the curve a researcher
                is."),
        why("On a small grant, output is close to proportional to resources and nearly independent of capability. Only
             grants comparable in size to capability itself make records informative about capability. For this reason a
             funder that relies on records alone stays far from the gap rule for many rounds, as the next page shows.")),
      div(class = "figcol", fig("r_plot1", 640, 300), fig("r_plot2", 640, 300)))),

  # ---------------------------------------------------------------- 5 review
  nav_panel("4 What peer review is worth",
    layout_columns(col_widths = c(4, 8),
      side(
        h2("What peer review is worth"),
        p("Peer review of a proposal gives the funder information about capability that is less mixed up with resources
           than the track record is. The paper models a review score as the researcher's capability plus random error,
           observed once. The value of review is measured as the share of the shortfall that review removes: the funder
           without review falls short of the complete-information optimum by some amount, and review recovers part of
           that amount."),
        sliderInput("v_tau", "Review noise (standard deviation of the score's error)", min = 0.05, max = 20, value = 1, step = 0.05),
        uiOutput("v_readout"),
        notice("Review is worth most where capability is heavy-tailed, and it keeps most of its value as the noise grows.
                Where capability is evenly spread, the same review score is worth much less, and its value falls quickly
                as the noise grows."),
        why("Where capability is heavy-tailed, most of the output a funder can add comes from getting resources to the
             few researchers of unusually high capability. Those researchers stand far apart from everyone else, so even a
             noisy score identifies them. Where capability is evenly spread, no researcher matters much more than another,
             so no score adds much."),
        p(class = "small-note", "Simulation results from the paper (Fig. 2B and 2C). The three fields have the same mean
                                 capability. The lower panel shows the intermediate field.")),
      div(class = "figcol", fig("v_plot", 640, 330), fig("v_rounds", 640, 260)))),

  # ---------------------------------------------------------------- 6 seed grants and lotteries
  nav_panel("5 Seed grants and lotteries",
    layout_columns(col_widths = c(4, 8),
      side(
        h2("Seed grants and lotteries"),
        p("Uniform seed grants and lotteries are two ways of spreading funds that give up some of the value of targeting.
           How much they give up depends on the field. Where capability is evenly spread and the budget is large, they
           give up little. Where capability is heavy-tailed and the budget is small, they give up a great deal."),
        h3("Why a lottery costs output even among equals"),
        p("Take a pool of researchers with the same capability and resources, and a fixed sum to divide among them.
           A lottery gives full grants to a few of them. Equal division gives each of them a smaller grant with the same
           expected value. Because output rises with resources at a diminishing rate, equal division produces more
           expected output."),
        sliderInput("l_m", "Researchers in the pool", min = 2, max = 20, value = 10, step = 1),
        sliderInput("l_k", "Lottery winners", min = 1, max = 10, value = 3, step = 1),
        uiOutput("l_readout"),
        notice("Upper panel: the output lost to seed grants grows faster than the share of the budget given out as seed
                grants. Lower panel: in a heavy-tailed field with a small budget, selecting by review produces the most
                output and a lottery over everyone the least. In an evenly spread field with a large budget the order
                reverses, and uniform funding produces the most. Where a partial lottery loses output, most of the loss
                comes from ignoring the review ranking, and only a small part from chance."),
        p(class = "small-note", "Simulation results from the paper (Fig. 3).")),
      div(class = "figcol", fig("s_plot", 720, 320), fig("lot_plot", 720, 400)))),

  # ---------------------------------------------------------------- 7 regime map
  nav_panel("6 Where your program stands",
    layout_columns(col_widths = c(4, 8),
      side(
        h2("Where your program stands"),
        p("Enter a program's annual awards to a field and the field's annual research spending from all other sources.
           Then choose how unequal capability is in the field. The map shows how much targeting is worth at that point
           and how much review adds, and the text below gives the paper's reading of that region of the map."),
        div(class = "d-flex gap-2 flex-wrap mb-2",
            actionButton("pre_nih", "Biomedical agency", class = "btn btn-outline-dark btn-sm"),
            actionButton("pre_nsf", "Science agency", class = "btn btn-outline-dark btn-sm"),
            actionButton("pre_erc", "Excellence council", class = "btn btn-outline-dark btn-sm"),
            actionButton("pre_fdn", "Foundation", class = "btn btn-outline-dark btn-sm")),
        numericInput("w_awards", "Annual awards to the field (millions)", value = 100, min = 0.1),
        numericInput("w_other", "The field's other annual research spending (millions)", value = 10000, min = 0.1),
        radioButtons("w_field", "Capability in the field is", choices = names(FIELDS), selected = "heavy-tailed"),
        uiOutput("w_readout"),
        p(class = "small-note", "A round of the model corresponds to a year of awards, and the map covers two rounds, so the
                                 map's budget scale is twice the annual ratio. The vertical position is the funder's own estimate.
                                 The paper suggests three observable indicators: how unequal research output is in the field, how
                                 dispersed review scores are, and how strongly output responds to substantial pilot grants.")),
      div(class = "figcol", fig("w_plot", 720, 420)))),

  # ---------------------------------------------------------------- 8 summary
  nav_panel("7 Summary",
    layout_columns(col_widths = c(7, 5),
      div(class = "main-text",
        h2("Summary"),
        div(class = "takeaway", p("The allocation that produces the most research output gives each researcher the gap between their capability
                                    and their resources. Funding the track record and funding the under-resourced each account for only one of the two,
                                    and either can produce less output than equal division.")),
        div(class = "takeaway", p("Track records do not reveal the gap, because output depends on capability and resources together. A funder that
                                    relies on records alone stays far from the optimum for many rounds.")),
        div(class = "takeaway", p("Peer review gives the funder information about capability that is less mixed up with resources than the track
                                    record. Review is worth most where capability is heavy-tailed and the budget is small. There, even a fairly noisy
                                    review score recovers most of that value.")),
        div(class = "takeaway", p("Seed grants and lotteries lose little where review is worth little: where capability is evenly spread and the
                                    budget is large. Where a partial lottery does lose output, most of the loss comes from ignoring the review
                                    ranking, and only a small part from chance.")),
        div(class = "takeaway", p("Two features of a field determine which of these applies: the funder's budget relative to the resources
                                    researchers already have, and how unequal capability is. A small foundation in a field with unequal capability
                                    gains a great deal from targeting. A large agency in a field with evenly spread capability loses little by
                                    spreading funds.")),
        h3("What the model holds fixed"),
        p("The model has one funder. Researchers do not change their behavior in response to the funding rule. Research output
           is a single quantity, so the model places no value on exploring new directions or on diversity. Review scores are
           noisy but not biased. Grants are used up in the round they are given; whatever lasts from them is counted as
           capability. The funder's objective is total expected research output, so fairness and risk are outside the model.
           The numbers in this guide come from the paper's simulations at the stated settings. The paper's claims are the
           qualitative comparisons, which hold across the ranges of parameters the paper examines."),
        h3("Sources"),
        p("Mohseni, A., DeDeo, S., and Zollman, K. J. S. A model of optimal science funding: targeting the capability-resource gap. Manuscript, 2026."),
        p("Model, simulation code, and figure data: ", tags$a(href = "https://github.com/amohseni/Grant-Funding-and-Scientific-Output-Model", "github.com/amohseni/Grant-Funding-and-Scientific-Output-Model"), "."),
        p(class = "small-note", "Supported by the John Templeton Foundation.")),
      div()))
)

# ================================================================== SERVER
server <- function(input, output, session) {

  # ---- 1 model
  output$m_plot <- renderPlot({
    Rs <- seq(0.05, 30, length.out = 300)
    d <- data.frame(R = Rs, y = lam(input$m_K, Rs))
    ggplot(d, aes(R, y)) +
      geom_hline(yintercept = input$m_K, colour = GREY, linetype = "dotted") +
      annotate("text", x = 30, y = input$m_K, label = "ceiling = A·K", hjust = 1, vjust = -0.4, colour = GREY, size = 3.8) +
      geom_line(linewidth = 0.7) +
      geom_point(data = data.frame(R = input$m_R, y = lam(input$m_K, input$m_R)), size = 3) +
      scale_y_continuous(limits = c(0, 31), expand = c(0, 0)) + scale_x_continuous(expand = c(0.01, 0)) +
      labs(x = "resources R", y = "expected research output", title = sprintf("Expected research output against resources, at capability K = %.1f", input$m_K)) +
      th()
  }, res = 110)
  output$m_readout <- renderUI({
    y <- lam(input$m_K, input$m_R)
    div(class = "readout", HTML(sprintf("Expected output is <span class='num'>%.2f</span>; the ceiling set by capability is %.1f. The scarce input is <span class='num'>%s</span>. One more unit of the scarce input adds %.3f to expected output; one more unit of the other input adds %.3f.",
      y, input$m_K, if (input$m_R < input$m_K) "resources" else "capability",
      max(input$m_K^2, input$m_R^2) / (input$m_K + input$m_R)^2, min(input$m_K^2, input$m_R^2) / (input$m_K + input$m_R)^2)))
  })

  # ---- 2 gap rule
  pop <- reactive({ input$g_redraw; draw_population(60, input$g_alpha, input$g_seed + input$g_redraw) })
  alloc <- reactive({
    d <- pop(); B <- input$g_b * 60 * mean_pareto(input$g_alpha)
    gr <- gap_rule(d$K, d$R0, B); k <- max(1, sum(gr$g > 0))
    base <- lam(d$K, d$R0)
    tr <- rep(0, 60); tr[order(-base)[1:k]] <- B / k
    ur <- rep(0, 60); ur[order(d$R0)[1:k]] <- B / k
    un <- rep(B / 60, 60)
    out <- function(g) sum(lam(d$K, d$R0 + g)) - sum(base)
    list(d = d, B = B, c = gr$c, k = k, g = list("the gap rule" = gr$g, "funding the track record" = tr, "funding the under-resourced" = ur),
         gain = c(gap = out(gr$g), track = out(tr), under = out(ur), uniform = out(un)))
  })
  output$g_plot <- renderPlot({
    a <- alloc(); d <- a$d; g <- a$g[[input$g_show]]
    d$funded <- g > 0; d$R1 <- d$R0 + g
    xmax <- max(d$R1, d$R0) * 1.05; ymax <- max(d$K) * 1.05
    ggplot(d) +
      geom_abline(intercept = 0, slope = 1 / a$c, linewidth = 0.5) +
      geom_segment(data = d[d$funded, ], aes(x = R0, y = K, xend = R1, yend = K), arrow = arrow(length = unit(4, "pt"), type = "closed"), linewidth = 0.4, colour = GREY) +
      geom_point(aes(R0, K, shape = funded), size = 2.2, fill = "white") +
      scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1)) +
      annotate("text", x = min(xmax, ymax * a$c) * 0.55, y = min(xmax, ymax * a$c) * 0.55 / a$c, label = "frontier R = cK", hjust = 1.08, vjust = -0.3, size = 3.8, colour = GREY) +
      coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
      labs(x = "resources R", y = "capability K", title = sprintf("Grants under %s\nFilled points are funded; each arrow is a grant", input$g_show)) +
      th()
  }, res = 110)
  output$g_readout <- renderUI({
    a <- alloc(); gn <- a$gain
    f <- function(x) sprintf("%.2f", x)
    div(class = "readout", HTML(sprintf(
      "Research output added over no funding, for this population and budget:<br>the gap rule <span class='num'>%s</span> &nbsp;·&nbsp; equal division <span class='num'>%s</span> &nbsp;·&nbsp; funding the track record <span class='num'>%s</span> &nbsp;·&nbsp; funding the under-resourced <span class='num'>%s</span>.<br><span class='small-note'>The gap rule funds %d of the 60 researchers. Each of the other two rules divides the same budget equally among the %d researchers it selects.</span>",
      f(gn["gap"]), f(gn["uniform"]), f(gn["track"]), f(gn["under"]), a$k, a$k)))
  })

  # ---- 3 records
  output$r_plot1 <- renderPlot({
    Ks <- c(input$r_Ks[1], round(exp(mean(log(input$r_Ks))), 0), input$r_Ks[2])
    Rs <- seq(0.05, 12, length.out = 300)
    d1 <- do.call(rbind, lapply(Ks, function(k) data.frame(K = k, R = Rs, y = lam(k, Rs))))
    d1$lab <- ifelse(d1$R == max(Rs), sprintf("K = %g", d1$K), NA)
    ggplot(d1, aes(R, y, group = K)) + geom_line(linewidth = 0.6) +
      geom_text(aes(label = lab), na.rm = TRUE, hjust = -0.15, size = 3.8) +
      coord_cartesian(xlim = c(0, 14.5), ylim = c(0, max(d1$y) * 1.08), expand = FALSE) +
      labs(x = "resources R", y = "expected research output", title = "Expected output against resources, for three capabilities") + th()
  }, res = 110)
  output$r_plot2 <- renderPlot({
    y <- input$r_y
    Rc <- seq(y * 1.02, 40, length.out = 300)
    d2 <- data.frame(R = Rc, K = y * Rc / (Rc - y)); d2 <- d2[d2$K <= 40, ]
    ex <- data.frame(R = c(y * 1.25, min(36, y * 6)))
    ex$K <- y * ex$R / (ex$R - y); ex$lab <- c("high capability, few resources", "modest capability, ample resources")
    ggplot(d2, aes(R, K)) + geom_line(linewidth = 0.6) +
      geom_point(data = ex, size = 2.6) +
      geom_text(data = ex, aes(label = lab), hjust = 0, nudge_x = 1.2, vjust = c(0.3, -0.5), size = 3.8) +
      coord_cartesian(xlim = c(0, 40), ylim = c(0, 40), expand = FALSE) +
      labs(x = "resources R", y = "capability K", title = sprintf("Researchers whose expected output is %.1f", y)) + th()
  }, res = 110)

  # ---- 4 review
  output$v_plot <- renderPlot({
    d <- rbind(data.frame(tau = D_review$tau, s = D_review$s13, e = D_review$e13, field = "heavy-tailed"),
               data.frame(tau = D_review$tau, s = D_review$s20, e = D_review$e20, field = "intermediate"),
               data.frame(tau = D_review$tau, s = D_review$s35, e = D_review$e35, field = "evenly spread"))
    lab <- d[d$tau == min(d$tau), ]; lab$vj <- c(-0.35, 1.25, 0.5)
    ggplot(d, aes(tau, s, group = field)) +
      geom_vline(xintercept = input$v_tau, colour = GREY, linetype = "dotted") +
      geom_linerange(aes(ymin = s - e, ymax = s + e), colour = GREY, linewidth = 0.4) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
      geom_text(data = lab, aes(label = field, vjust = vj), hjust = 1.08, size = 3.6) +
      scale_x_log10(breaks = c(0.05, 0.2, 1, 5, 20), expand = expansion(mult = c(0.35, 0.05))) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      labs(x = "review noise (log scale)", y = "share of shortfall recovered", title = "Share of the shortfall that review recovers, by field") + th()
  }, res = 110)
  output$v_rounds <- renderPlot({
    d <- rbind(data.frame(round = D_rounds$round, s = D_rounds$sf4, who = "no review"), data.frame(round = D_rounds$round, s = D_rounds$sf5, who = "with review"))
    lab <- d[d$round == 20, ]
    ggplot(d, aes(round, s, group = who)) + geom_line(linewidth = 0.6) +
      geom_text(data = lab, aes(label = who), hjust = -0.1, size = 3.6) +
      scale_x_continuous(breaks = c(1, 5, 10, 15, 20), expand = expansion(mult = c(0.02, 0.25))) +
      scale_y_continuous(limits = c(0, 0.5), expand = c(0, 0)) +
      labs(x = "round", y = "shortfall (share of gain)", title = "Shortfall below the complete-information optimum,\nby round, in the intermediate field") + th()
  }, res = 110)
  output$v_readout <- renderUI({
    interp <- function(col) approx(log(D_review$tau), D_review[[col]], xout = log(input$v_tau), rule = 2)$y
    div(class = "readout", HTML(sprintf("At this level of noise, review recovers <span class='num'>%d%%</span> of the shortfall in the heavy-tailed field, <span class='num'>%d%%</span> in the intermediate field, and <span class='num'>%d%%</span> in the evenly spread field.",
      round(100 * interp("s13")), round(100 * interp("s20")), round(100 * interp("s35")))))
  })

  # ---- 5 seed grants and lotteries
  observeEvent(input$l_m, updateSliderInput(session, "l_k", max = input$l_m, value = min(input$l_k, input$l_m)))
  output$l_readout <- renderUI({
    K <- 10; R0 <- 2; B <- 10; m <- input$l_m; k <- min(input$l_k, m)
    certain <- m * (lam(K, R0 + B / m) - lam(K, R0))
    lottery <- k * (lam(K, R0 + B / k) - lam(K, R0))
    div(class = "readout", HTML(sprintf("A pool of %d researchers with equal capability and resources, and a fixed sum to divide. Equal division adds <span class='num'>%.2f</span> units of expected output. A lottery that gives full grants to %d of them adds <span class='num'>%.2f</span>. Equal division produces more whenever grants can be divided and output has diminishing returns to resources.", m, certain, k, lottery)))
  })
  output$s_plot <- renderPlot({
    d <- rbind(data.frame(x = D_seed$xseed, y = D_seed$hb01, f = "heavy-tailed, b = 0.2"),
               data.frame(x = D_seed$xseed, y = D_seed$hb05, f = "heavy-tailed, b = 1"),
               data.frame(x = D_seed$xseed, y = D_seed$hb1, f = "heavy-tailed, b = 2"),
               data.frame(x = D_seed$xseed, y = D_seed$db05, f = "intermediate, b = 1"))
    lab <- d[d$x == max(d$x), ]
    ggplot(d, aes(x, y, group = f)) + geom_line(linewidth = 0.6) + geom_point(size = 1.4) +
      geom_text(data = lab, aes(label = f), hjust = -0.05, size = 3.4) +
      scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75), limits = c(0, 1.25), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, 20), expand = c(0, 0)) +
      labs(x = "fraction of the budget given out as seed grants", y = "output lost (%)", title = "Research output lost to uniform seed grants,\nas a percent of the gain from funding") + th()
  }, res = 110)
  output$lot_plot <- renderPlot({
    mk <- function(D, f, u, l) rbind(data.frame(tau = D$tau, y = D$ss, s = "top fifth by review", field = f),
                                     data.frame(tau = D$tau, y = D$ws, s = "top two fifths by review", field = f),
                                     data.frame(tau = D$tau, y = D$sl, s = "partial lottery", field = f),
                                     data.frame(tau = D$tau, y = u, s = "uniform funding", field = f),
                                     data.frame(tau = D$tau, y = l, s = "lottery over all", field = f))
    d <- rbind(mk(D_lot_h, "heavy-tailed, tight budget", 0.445, 0.378), mk(D_lot_e, "evenly spread, ample budget", 0.831, 0.414))
    d$field <- factor(d$field, levels = c("heavy-tailed, tight budget", "evenly spread, ample budget"))
    d$s <- factor(d$s, levels = c("top fifth by review", "top two fifths by review", "partial lottery", "uniform funding", "lottery over all"))
    ggplot(d, aes(tau, y, group = s, linetype = s, shape = s, colour = s)) +
      geom_line(linewidth = 0.6) + geom_point(size = 1.9, fill = "white", na.rm = TRUE) +
      scale_shape_manual(values = c(16, 1, 17, NA, NA), name = NULL) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "dashed", "dotted"), name = NULL) +
      scale_colour_manual(values = c("black", "black", "black", GREY, GREY), name = NULL) +
      facet_wrap(~field) +
      scale_x_log10(breaks = c(0.3, 1, 3, 10), labels = c("0.3", "1", "3", "10"), expand = expansion(mult = c(0.06, 0.06))) +
      scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
      guides(shape = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2), colour = guide_legend(nrow = 2)) +
      labs(x = "review noise (log scale)", y = "share of the optimum's gain", title = "Gain of five funding schemes,\nas a share of the optimal allocation's gain") +
      th() + theme(legend.position = "bottom", legend.key.width = unit(26, "pt"), legend.text = element_text(size = 11), legend.margin = margin(4, 0, 0, 0), panel.spacing.x = unit(18, "pt"))
  }, res = 110)

  # ---- 6 regime map
  observeEvent(input$pre_nih, { updateNumericInput(session, "w_awards", value = 33100); updateNumericInput(session, "w_other", value = 22000) })
  observeEvent(input$pre_nsf, { updateNumericInput(session, "w_awards", value = 6700);  updateNumericInput(session, "w_other", value = 47000) })
  observeEvent(input$pre_erc, { updateNumericInput(session, "w_awards", value = 2300);  updateNumericInput(session, "w_other", value = 86000) })
  observeEvent(input$pre_fdn, { updateNumericInput(session, "w_awards", value = 100);   updateNumericInput(session, "w_other", value = 10000) })
  where <- reactive({
    req(input$w_awards > 0, input$w_other > 0)
    r <- input$w_awards / input$w_other; b <- min(max(2 * r, 0.01), 3)
    gini <- gini_of_alpha(FIELDS[[input$w_field]])
    list(r = r, b = b, gini = gini, targ = grid_value(b, gini, "targ"), rev = grid_value(b, gini, "rev"), clipped = 2 * r < 0.01 || 2 * r > 3)
  })
  output$w_plot <- renderPlot({
    w <- where()
    yb <- gini_of_alpha(FIELDS); ylabs <- names(FIELDS)
    pt <- data.frame(b = w$b, gini = w$gini)
    mkpanel <- function(C, title, ends) {
      lab <- do.call(rbind, lapply(split(C, C$seg), function(s) s[nrow(s), ]))
      ggplot(C, aes(b, y, group = seg)) + geom_path(linewidth = 0.5) +
        geom_text(data = lab[lab$b > 2.9, ], aes(label = level), hjust = -0.2, size = 3.3) +
        geom_point(data = pt, aes(b, gini), inherit.aes = FALSE, size = 3.6, shape = 21, fill = "black", colour = "white", stroke = 1.2) +
        scale_x_log10(breaks = c(0.01, 0.1, 1), limits = c(0.01, 3), expand = expansion(mult = c(0, 0.2))) +
        scale_y_continuous(breaks = yb, labels = ylabs, limits = c(0.1, 0.72), expand = c(0, 0)) +
        labs(x = "budget relative to field resources", y = if (ends) "capability inequality" else NULL, title = title) +
        th() + theme(axis.text.y = element_text(size = if (ends) 11 else 0), axis.title.x = element_text(size = 11))
    }
    p1 <- mkpanel(C_targ, "A  Value of targeting", TRUE)
    p2 <- mkpanel(C_rev, "B  Value of review", FALSE)
    grid::grid.newpage(); grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2, widths = grid::unit(c(1.12, 1), "null"))))
    print(p1, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)); print(p2, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  }, res = 110)
  output$w_readout <- renderUI({
    w <- where()
    regime <- if (w$targ >= 1) {
      "In this region targeting matters most. The optimal allocation adds at least as much research output as uniform funding adds in total. Review that identifies the few researchers of unusually high capability is worth most here, and even a fairly noisy review score recovers most of that value. Seed grants and lotteries give up the most here."
    } else if (w$targ >= 0.25) {
      "In this region targeting is worth a moderate amount. Review adds a real but limited share of research output, and seed grants and lotteries give up a corresponding share of the value of targeting."
    } else {
      "In this region targeting is worth little. The optimal allocation and equal division produce nearly the same research output. Review adds little however accurate it is, seed grants and lotteries give up little, and the best use of the budget in the model is small grants to many researchers."
    }
    div(class = "readout", HTML(sprintf(
      "Annual awards are <span class='num'>%.3g</span> of the field's other spending, so the program sits at budget scale <span class='num'>b = %.2g</span>%s, in a field with %s capability.<br>Value of targeting: <span class='num'>%.2f</span> times uniform funding's gain. Value of review at the default noise: <span class='num'>%.0f%%</span> of the research output that review-informed funding adds.<br><br>%s",
      w$r, w$b, if (w$clipped) " (at the edge of the map)" else "", input$w_field, w$targ, 100 * w$rev, regime)))
  })
}

shinyApp(ui, server)
