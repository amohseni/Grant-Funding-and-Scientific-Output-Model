# ============================================================
# sweep_T.R  —  T-round sweep manifest + runner (§5 manifest, §7 schema)
# ------------------------------------------------------------
# Drives run_simulation_T over the final parameter manifest and writes
# per-sweep RDS (summary + raw + raw-long schedule) plus plots.
#
# Depends on: model.R (the model engine — run_simulation_T etc.) and sweep.R
# (reuses summarize_sweep + plotting helpers + `%||%`). Source both first.
#
# Conventions / decisions (see PROGRESS.md):
#   * Budget: mean-E[R] normalization, default b=0.5. b-axes clamped to [0.1,1].
#   * Reproducibility: seed = trial index, SHARED across cells (common random
#     numbers). Deliberate deviation from spec §7's hash(cell_id) rule: CRN
#     across cells reduces the variance of cross-cell CONTRASTS (the horizon
#     trend PG(T) is exactly such a contrast). run_simulation_T sets its own
#     RNG from the seed arg, so runs are reproducible and worker-order-stable.
#   * Parallelism: parallel::mclapply (fork) with an explicit core count;
#     falls back to serial lapply if mc.cores == 1.
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(parallel)
})

# ---- Base params for sweeps (spec §2 defaults, our budget decision) ----
SWEEP_BASE_PARAMS_T <- list(
  T_rounds = 2, n = 50,
  k_min = 1.0, k_shape = 2.0, r_min = 1.0, r_shape = 2.0, rho_kr = 0,
  gamma = 1.0, epsilon = 0.1, b = 0.5, n_steps = 50,
  tau_r = 1.0, tau_k = 1.0, use_resource_signal = TRUE,
  n_pre_rounds = 0, x_seed = 0.25, M = 400, strategies = 1:9   # M bumped 200->400
)                                                              # (ESS headroom after SIR removal)

SWEEP_CORES <- max(1L, detectCores() - 1L)

# ============================================================
# MANIFEST (§5). Param name map: ε→epsilon, τ_K→tau_k, τ_R→tau_r,
# α_K→k_shape, ρ_c→rho_kr, T→T_rounds, n_pre→n_pre_rounds.
# b-axes clamped to [0.1, 1] per user directive (trims noted inline).
# ============================================================
SWEEP_CONFIGS_T <- list(

  # ---------- Tier 1: the horizon question (~55% of compute) ----------
  horizon_growth = list(
    name = "horizon_growth", tier = 1,
    description = "Horizon T vs knowledge-growth epsilon. Headline: does forward beat myopic as T grows, and how does growth modulate it?",
    grid_fn = function() expand.grid(T_rounds = 1:5, epsilon = c(0.05, 0.15, 0.3, 0.55, 0.85)),
    varied_params = c("T_rounds", "epsilon"),
    primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                        color_var = "epsilon", title = "Forward vs myopic (S8-S5) across horizon",
                        y_label = "S8 - S5 (expected output)"),
    secondary_plot = list(type = "heatmap", x_var = "T_rounds", y_var = "epsilon",
                          fill_var = "fwd_vs_myo_PG_mean", text_var = "fwd_vs_myo_PG_mean",
                          title = "S8 - S5 across horizon and growth", fill_label = "S8-S5")
  ),
  horizon_noise = list(
    name = "horizon_noise", tier = 1,
    description = "Horizon T vs grant-signal noise tau_k.",
    grid_fn = function() expand.grid(T_rounds = 1:5, tau_k = c(0.1, 0.3, 1, 3, 10)),
    varied_params = c("T_rounds", "tau_k"),
    primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                        color_var = "tau_k", title = "Forward vs myopic (S8-S5) across horizon",
                        y_label = "S8 - S5 (expected output)"),
    secondary_plot = list(type = "heatmap", x_var = "T_rounds", y_var = "tau_k",
                          fill_var = "signal_fwd_mean", text_var = "signal_fwd_mean",
                          title = "Forward signal value (S8-S7) across horizon and noise", fill_label = "S8-S7")
  ),
  horizon_scale = list(
    name = "horizon_scale", tier = 1,
    description = "Horizon T vs budget scale b (clamped to [0.1,1]; spec 0.05 dropped).",
    grid_fn = function() expand.grid(T_rounds = 1:5, b = c(0.1, 0.3, 0.5, 1.0)),
    varied_params = c("T_rounds", "b"),
    primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                        color_var = "b", title = "Forward vs myopic (S8-S5) across horizon",
                        y_label = "S8 - S5 (expected output)"),
    secondary_plot = list(type = "heatmap", x_var = "T_rounds", y_var = "b",
                          fill_var = "fwd_vs_myo_PG_mean", text_var = "fwd_vs_myo_PG_mean",
                          title = "S8 - S5 across horizon and budget", fill_label = "S8-S5")
  ),

  # ---------- Tier 2: T=2 cross-sections ----------
  signal_precision = list(
    name = "signal_precision", tier = 2,
    description = "Fine sweep of grant-signal noise tau_k at T=2.",
    grid_fn = function() data.frame(tau_k = c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1, 1.5, 2, 3, 5, 10, 20)),
    varied_params = c("tau_k"),
    primary_plot = list(type = "line", x_var = "tau_k", y_var = "signal_fwd_mean",
                        color_var = "tau_k", title = "Forward signal value (S8-S7) vs tau_k",
                        y_label = "S8 - S7 (expected output)"),
    secondary_plot = NULL
  ),
  signal_value = list(
    name = "signal_value", tier = 2,
    description = "Knowledge tail alpha_K vs grant-signal noise tau_k at T=2.",
    grid_fn = function() expand.grid(k_shape = c(1.3, 1.8, 2.5, 3.5), tau_k = c(0.3, 1, 3, 10)),
    varied_params = c("k_shape", "tau_k"),
    primary_plot = list(type = "heatmap", x_var = "tau_k", y_var = "k_shape",
                        fill_var = "signal_fwd_mean", text_var = "signal_fwd_mean",
                        title = "Forward signal value (S8-S7) across noise and inequality", fill_label = "S8-S7"),
    secondary_plot = list(type = "line", x_var = "tau_k", y_var = "signal_fwd_mean",
                          color_var = "k_shape", title = "Forward signal value vs tau_k by tail",
                          y_label = "S8 - S7")
  ),
  funder_scale = list(
    name = "funder_scale", tier = 2,
    description = "Budget scale b at T=2 (clamped to [0.1,1]; spec 1.5/2.5/5 and 0.05 dropped).",
    grid_fn = function() data.frame(b = c(0.1, 0.2, 0.3, 0.5, 0.75, 1.0)),
    varied_params = c("b"),
    primary_plot = list(type = "line", x_var = "b", y_var = "fwd_vs_myo_PG_mean",
                        color_var = "b", title = "Forward vs myopic (S8-S5) vs budget", y_label = "S8 - S5"),
    secondary_plot = NULL
  ),
  seed_value = list(
    name = "seed_value", tier = 2,
    description = "Budget b vs seed fraction x_seed at T=2 (b clamped to [0.1,1]; spec 2 dropped).",
    grid_fn = function() expand.grid(b = c(0.1, 0.3, 0.5, 1.0), x_seed = c(0.1, 0.25, 0.5, 0.75)),
    varied_params = c("b", "x_seed"),
    primary_plot = list(type = "heatmap", x_var = "b", y_var = "x_seed",
                        fill_var = "seed_fwd_mean", text_var = "seed_fwd_mean",
                        title = "Seed value under forward planner (S9-S7)", fill_label = "S9-S7"),
    secondary_plot = NULL
  ),
  alpha_regime = list(
    name = "alpha_regime", tier = 2,
    description = "Grant-signal noise tau_k vs growth epsilon: where does forward spend (b_idx / alpha)?",
    grid_fn = function() expand.grid(tau_k = c(0.3, 1, 3, 10), epsilon = c(0.05, 0.2, 0.5, 0.8)),
    varied_params = c("tau_k", "epsilon"),
    primary_plot = list(type = "heatmap", x_var = "tau_k", y_var = "epsilon",
                        fill_var = "b_idx_S8_mean", text_var = "b_idx_S8_mean",
                        title = "Forward schedule center-of-mass b_idx (S8)", fill_label = "b_idx"),
    secondary_plot = NULL
  ),

  # ---------- Tier 3 ----------
  pre_rounds = list(
    name = "pre_rounds", tier = 3,
    description = "Pre-round baseline observations n_pre at T=2.",
    grid_fn = function() data.frame(n_pre_rounds = c(0, 1, 5, 10, 20)),
    varied_params = c("n_pre_rounds"),
    primary_plot = list(type = "line", x_var = "n_pre_rounds", y_var = "signal_fwd_mean",
                        color_var = "n_pre_rounds", title = "Signal value vs pre-round observations",
                        y_label = "S8 - S7"),
    secondary_plot = NULL
  ),
  regime_map = list(
    name = "regime_map", tier = 3,
    description = "3D regime map: budget b x noise tau_k x tail alpha_K at T=2 (b clamped to [0.1,1]; spec 3 dropped).",
    grid_fn = function() expand.grid(b = c(0.1, 0.3, 1.0), tau_k = c(0.3, 1, 3), k_shape = c(1.5, 2, 3)),
    varied_params = c("b", "tau_k", "k_shape"),
    primary_plot = list(type = "heatmap", x_var = "tau_k", y_var = "b", facet_var = "k_shape",
                        fill_var = "signal_fwd_mean", text_var = "signal_fwd_mean",
                        title = "Forward signal value (S8-S7) regime map", fill_label = "S8-S7"),
    secondary_plot = NULL
  ),

  # ---------- Tier 4 ----------
  correlation = list(
    name = "correlation", tier = 4,
    description = "K-R Gaussian-copula correlation rho_c vs tail alpha_K at T=2, tau_k=0.3 (sharpest signal).",
    grid_fn = function() cbind(expand.grid(rho_kr = c(-0.5, 0, 0.5, 0.8), k_shape = c(1.3, 2, 3.5)), tau_k = 0.3),
    varied_params = c("rho_kr", "k_shape"),
    primary_plot = list(type = "heatmap", x_var = "rho_kr", y_var = "k_shape",
                        fill_var = "signal_fwd_mean", text_var = "signal_fwd_mean",
                        title = "Signal value vs K-R correlation and tail (tau_k=0.3)", fill_label = "S8-S7"),
    secondary_plot = list(type = "line", x_var = "rho_kr", y_var = "signal_fwd_mean",
                          color_var = "k_shape", title = "Signal value vs correlation by tail",
                          y_label = "S8 - S7")
  ),
  pop_size = list(
    name = "pop_size", tier = 4,
    description = "Population size n at T=2 (finite-n check).",
    grid_fn = function() data.frame(n = c(20, 50, 100, 200)),
    varied_params = c("n"),
    primary_plot = list(type = "line", x_var = "n", y_var = "fwd_vs_myo_PG_mean",
                        color_var = "n", title = "Forward vs myopic (S8-S5) vs population size",
                        y_label = "S8 - S5"),
    secondary_plot = NULL
  ),
  resource_noise = list(
    name = "resource_noise", tier = 4,
    description = "Resource-signal noise tau_r at T=2.",
    grid_fn = function() data.frame(tau_r = c(0.3, 1, 3)),
    varied_params = c("tau_r"),
    primary_plot = list(type = "line", x_var = "tau_r", y_var = "signal_fwd_mean",
                        color_var = "tau_r", title = "Signal value vs resource-signal noise",
                        y_label = "S8 - S7"),
    secondary_plot = NULL
  ),

  # ---------- Supplementary sweeps (Roadmap Stage 1.1; see NEXT_SIMULATIONS.md) ----------
  tail_map = list(
    name = "tail_map", tier = 2,
    description = "Knowledge tail k_shape vs resource tail r_shape at T=2, sharp signal (tau_k=0.3). NOTE: budget B∝E[R] co-varies with r_shape (see NEXT_SIMULATIONS.md A1).",
    grid_fn = function() cbind(expand.grid(k_shape = c(1.3, 2, 3.5), r_shape = c(1.3, 2, 3.5)), tau_k = 0.3),
    varied_params = c("k_shape", "r_shape"),
    primary_plot = list(type = "heatmap", x_var = "r_shape", y_var = "k_shape",
                        fill_var = "signal_fwd_mean", text_var = "signal_fwd_mean",
                        title = "Signal value vs knowledge and resource tails", fill_label = "S8-S7"),
    secondary_plot = NULL
  ),
  info_value = list(
    name = "info_value", tier = 1,
    description = "Information channel isolated: eps~0 (no compounding) across horizon T and signal noise tau_k.",
    grid_fn = function() cbind(expand.grid(T_rounds = c(2, 3, 4, 5), tau_k = c(0.3, 1, 3)), epsilon = 1e-4),
    varied_params = c("T_rounds", "tau_k"),
    primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                        color_var = "tau_k", title = "Information value (eps->0) across horizon",
                        y_label = "S8 - S5 (pure info value)"),
    secondary_plot = NULL
  ),
  horizon_long = list(
    name = "horizon_long", tier = 1,
    description = "Long horizon T=5..10 at eps in {0.3,0.85}: does the forward advantage saturate? Uses n_steps=400 (finer greedy granularity): the default n_steps=50 is too coarse beyond T=5 and inflates the advantage ~4x (see tests/horizon_long_convergence.R).",
    grid_fn = function() cbind(expand.grid(T_rounds = c(5, 6, 7, 8, 10), epsilon = c(0.3, 0.85)), n_steps = 400),
    varied_params = c("T_rounds", "epsilon"),
    primary_plot = list(type = "line", x_var = "T_rounds", y_var = "fwd_vs_myo_PG_mean",
                        color_var = "epsilon", title = "Forward advantage at long horizons",
                        y_label = "S8 - S5"),
    secondary_plot = NULL
  ),

  # ---------- Tier 1: the resource-poverty question (front-load reversal) ----------
  # Does a resource-poor community make the forward planner FRONT-load (invest early
  # to bootstrap growth + learning), reversing the paper's back-loading? Requires the
  # DECOUPLED budget (budget_ref="K"): under the default coupling B=b·n·E[R], driving
  # r_min->0 also zeroes the budget (a trivial null). With budget_ref="K" the funder
  # keeps a fixed purse (B=b·n·E[K], r_min-independent) while baseline resources vary.
  # Axes: baseline resources r_min (poor -> rich) x compounding rate epsilon (info-
  # dominated -> knowledge-dominated). Focal metric b_idx_S8 (schedule center-of-mass;
  # <0.5 front-load, >0.5 back-load) and its 0.5 contour = the reversal boundary.
  # T=5 fixed (back-loading is clearer at longer T); base tau_k=1, k_shape=2.
  resource_regime = list(
    name = "resource_regime", tier = 1,
    description = "Baseline resources r_min x compounding epsilon at fixed decoupled purse (budget_ref=K), T=5, smooth. Maps the front-load<->back-load reversal (b_idx_S8=0.5 contour).",
    grid_fn = function() cbind(
      expand.grid(r_min   = c(0.001, 0.03, 0.3, 1, 3),
                  epsilon = c(1e-4, 0.01, 0.03, 0.1, 0.3, 0.85)),
      T_rounds = 5, budget_ref = "K", n_steps = 50),
    varied_params = c("r_min", "epsilon"),
    primary_plot = list(type = "heatmap", x_var = "epsilon", y_var = "r_min",
                        fill_var = "b_idx_S8_mean", text_var = "b_idx_S8_mean",
                        title = "Forward schedule b_idx (S8): <0.5 front-load, >0.5 back-load",
                        fill_label = "b_idx"),
    secondary_plot = list(type = "line", x_var = "epsilon", y_var = "b_idx_S8_mean",
                          color_var = "r_min", title = "Schedule center-of-mass vs compounding, by baseline resources",
                          y_label = "b_idx (S8)")
  ),

  # ---------- Tier 1: the pure-exploration corner (paid information isolated) ----------
  # The resource_regime sweep held tau_k=1 — a readable FREE peer-review signal, which
  # violates the "you learn nothing without funding" premise. Here the free signal is
  # switched off (tau_k=100 ≈ uninformative vs K~Pareto(1.3)), the community is poor
  # (r_min=0.001 -> initial pubs ~0, uninformative), and the talent tail is heavy
  # (k_shape=1.3: identifying winners is most of the problem). All information must then
  # be PAID FOR via funded output (force D pure). Axes: tau_k (free-info kill switch) x
  # epsilon (growth force off->on). Watch b_idx_S8 (<0.5 = money front-loaded),
  # alpha_S8 (round-1 share), and the full alpha_t schedule in rawlong (distinguishes
  # true front-loading from "seed-and-harvest": positive early seed, mass late).
  exploration_corner = list(
    name = "exploration_corner", tier = 1,
    description = "Paid-information corner: free signal off (tau_k up to 100), poor community (r_min=0.001), heavy tail (k_shape=1.3), T=6, decoupled purse. Does the planner front-load money to buy learning?",
    grid_fn = function() cbind(
      expand.grid(tau_k   = c(1, 3, 10, 100),
                  epsilon = c(1e-4, 0.03, 0.1, 0.3)),
      r_min = 0.001, k_shape = 1.3, T_rounds = 6, budget_ref = "K"),
    varied_params = c("tau_k", "epsilon"),
    primary_plot = list(type = "heatmap", x_var = "epsilon", y_var = "tau_k",
                        fill_var = "b_idx_S8_mean", text_var = "b_idx_S8_mean",
                        title = "Schedule b_idx (S8) in the paid-info corner (r_min=0.001)",
                        fill_label = "b_idx"),
    secondary_plot = list(type = "line", x_var = "tau_k", y_var = "fwd_vs_myo_PG_mean",
                          color_var = "epsilon", title = "Planning value as the free signal degrades",
                          y_label = "S8 - S5")
  ),

  # Companion: fix eps~0 (pure information) and vary poverty x signal. Does the
  # exploration tilt require poverty, or only signal absence?
  exploration_poverty = list(
    name = "exploration_poverty", tier = 1,
    description = "Pure-info schedule tilt across baseline resources r_min x free-signal noise tau_k, at eps~0, T=6, heavy tail, decoupled purse.",
    grid_fn = function() cbind(
      expand.grid(r_min = c(0.001, 0.03, 0.3, 1),
                  tau_k = c(1, 100)),
      epsilon = 1e-4, k_shape = 1.3, T_rounds = 6, budget_ref = "K"),
    varied_params = c("r_min", "tau_k"),
    primary_plot = list(type = "line", x_var = "r_min", y_var = "b_idx_S8_mean",
                        color_var = "tau_k", title = "Pure-info schedule tilt vs baseline resources",
                        y_label = "b_idx (S8)"),
    secondary_plot = list(type = "line", x_var = "r_min", y_var = "fwd_vs_myo_PG_mean",
                          color_var = "tau_k", title = "Planning value vs baseline resources (eps~0)",
                          y_label = "S8 - S5")
  ),

  # Depth: paid information requires grants deep enough that output is TALENT-limited.
  # With lambda = K·g/(K+g), a small grant (g << K) pins output to the grant size
  # (lambda ~ g), so funded output reveals almost nothing about K — probe showed the
  # planner ties uniform funding in the poverty/no-signal corner at b=0.5 (per-round
  # per-capita grant ~0.33 << K~Pareto(1.3, min 1)). Talent-revealing output needs
  # g >~ K, i.e. per-capita depth 2bE[K]/T >~ E[K] -> b >~ T/2. Axis: budget scale b
  # from thin to deep x free-signal on/off, in the poverty corner at eps~0 (pure info).
  # Watch: out_S8 - out_S2 (does paid info beat uniform at depth?), alpha_S8, b_idx_S8.
  exploration_depth = list(
    name = "exploration_depth", tier = 1,
    description = "Funding depth b x free-signal tau_k in the poverty corner (r_min=0.001, eps~0, heavy tail, T=6): does deep funding make paid information work, and which way does it tilt the schedule?",
    grid_fn = function() cbind(
      expand.grid(b     = c(0.5, 1.5, 3, 6),
                  tau_k = c(1, 100)),
      epsilon = 1e-4, r_min = 0.001, k_shape = 1.3, T_rounds = 6, budget_ref = "K"),
    varied_params = c("b", "tau_k"),
    primary_plot = list(type = "line", x_var = "b", y_var = "b_idx_S8_mean",
                        color_var = "tau_k", title = "Schedule tilt vs funding depth (poverty, eps~0)",
                        y_label = "b_idx (S8)"),
    secondary_plot = list(type = "line", x_var = "b", y_var = "fwd_vs_myo_PG_mean",
                          color_var = "tau_k", title = "Planning value vs funding depth",
                          y_label = "S8 - S5")
  ),

  # ---------- Diagnostics (SWEEP_HANDOFF_2026-08-05, Package A) ----------
  # D3: the seed floor imposed on the SIGNAL strategies (S10 = S5+seed, S11 = S8+seed),
  # where targeting is worth most. Readouts: seed_sig_myo = S10-S5, seed_sig_fwd = S11-S8,
  # vs the family targeting value S5-S2 (P-A1: cost ~ kappa * targeting * floored share).
  D3_seed_signal = list(
    name = "D3_seed_signal", tier = 1,
    description = "Seed floor on the signal strategies (S10, S11) across b x x_seed, base and heavy+sharp regimes, T=2 smooth.",
    grid_fn = function() rbind(
      cbind(expand.grid(b = c(0.1, 0.3, 0.5, 1.0), x_seed = c(0.1, 0.25, 0.5, 0.75)),
            k_shape = 2.0, tau_k = 1.0),
      cbind(expand.grid(b = c(0.1, 0.3, 0.5, 1.0), x_seed = c(0.1, 0.25, 0.5, 0.75)),
            k_shape = 1.3, tau_k = 0.3)),
    varied_params = c("b", "x_seed", "k_shape"),
    primary_plot = list(type = "line", x_var = "x_seed", y_var = "seed_sig_fwd_mean",
                        color_var = "b", title = "Cost of the seed floor on the forward+signal family",
                        y_label = "S11 - S8"),
    secondary_plot = list(type = "line", x_var = "x_seed", y_var = "seed_sig_myo_mean",
                          color_var = "b", title = "Cost of the seed floor on the myopic+signal family",
                          y_label = "S10 - S5")
  ),

  # Tier B concentration gate (SWEEP_HANDOFF Package B3): one static round (T=1),
  # best static funder (S5), fine grid over the CES exponent. P-B2: Gini of the
  # optimal allocation monotone as ces_gamma falls toward Leontief, at every
  # (k_shape, b). Leontief endpoint runs separately (n_steps=800, degeneracy caveat).
  sigma_tierB = list(
    name = "sigma_tierB", tier = 1,
    description = "Concentration of the optimal static allocation (S5, T=1) across CES exponent x talent tail x budget; greedy n_steps=400.",
    grid_fn = function() cbind(
      expand.grid(ces_gamma = c(0, -0.25, -0.5, -1, -2, -3, -6, -12),
                  k_shape   = c(1.3, 2, 3.5),
                  b         = c(0.1, 0.5, 1.0)),
      T_rounds = 1),
    varied_params = c("ces_gamma", "k_shape", "b"),
    primary_plot = list(type = "line", x_var = "ces_gamma", y_var = "gini_g1_S5_mean",
                        color_var = "k_shape", title = "Concentration of optimal funding vs CES exponent",
                        y_label = "Gini of round-1 grants (S5)"),
    secondary_plot = list(type = "line", x_var = "ces_gamma", y_var = "coverage_S5_mean",
                          color_var = "k_shape", title = "Coverage vs CES exponent",
                          y_label = "share of researchers funded (S5)")
  ),

  # Back-loading attribution (SWEEP_HANDOFF Package C2): decouple free growth
  # (eps_free: K compounds from baseline resources regardless of grants) from paid
  # growth (eps_paid: the grant-caused increment). P-C1: pure free back-loads.
  # P-C2: pure paid FRONT-loads. Diagonal = coupled model (identity).
  bload_decouple = list(
    name = "bload_decouple", tier = 1,
    description = "Free vs paid growth decoupled (eps_free x eps_paid), T=5 smooth harmonic: which channel drives back-loading?",
    grid_fn = function() cbind(
      expand.grid(eps_free = c(0, 0.1, 0.3, 0.85),
                  eps_paid = c(0, 0.1, 0.3, 0.85)),
      T_rounds = 5, epsilon = 0.3),   # epsilon value irrelevant once both rates given;
                                      # fixed so the globals always engage off-diagonal
    varied_params = c("eps_free", "eps_paid"),
    primary_plot = list(type = "heatmap", x_var = "eps_paid", y_var = "eps_free",
                        fill_var = "b_idx_S8_mean", text_var = "b_idx_S8_mean",
                        title = "Schedule b_idx (S8): free vs paid growth",
                        fill_label = "b_idx"),
    secondary_plot = list(type = "line", x_var = "eps_paid", y_var = "fwd_vs_myo_PG_mean",
                          color_var = "eps_free", title = "Planning value across the decoupled grid",
                          y_label = "S8 - S5")
  ),

  # D4: persistent every-round floor in the no-signal myopic family. x_seed = 1 must
  # reproduce S2 exactly (P-A2 identity check on the implementation).
  D4_seed_persistent = list(
    name = "D4_seed_persistent", tier = 1,
    description = "Persistent every-round seed floor (seed_persistent=TRUE), S6 vs S4/S2, x_seed up to 1.0 (identity check), T=2 smooth.",
    grid_fn = function() rbind(
      cbind(expand.grid(x_seed = c(0.25, 0.5, 0.75, 0.9, 1.0)), k_shape = 2.0, tau_k = 1.0),
      cbind(expand.grid(x_seed = c(0.25, 0.5, 0.75, 0.9, 1.0)), k_shape = 1.3, tau_k = 0.3)),
    varied_params = c("x_seed", "k_shape"),
    primary_plot = list(type = "line", x_var = "x_seed", y_var = "seed_myo_mean",
                        color_var = "k_shape", title = "Cost of a persistent uniform floor (S6 - S4)",
                        y_label = "S6 - S4"),
    secondary_plot = NULL
  )
)

# ============================================================
# METRIC EXTRACTION (§7 schema)
# ============================================================
# Wide summary-row fields for one run: out/alpha per strategy, forward b_idx
# and gini, rho_s, plus the derived contrasts.
extract_metrics_T <- function(res) {
  out <- list()
  for (s in 1:11) {
    r <- res$strategies[[s]]
    out[[paste0("out_S", s)]]   <- if (is.null(r)) NA_real_ else r$total_expected
    out[[paste0("alpha_S", s)]] <- if (is.null(r) || is.null(r$alpha)) NA_real_ else r$alpha
  }
  for (s in c(7, 8, 9, 11)) {
    r <- res$strategies[[s]]
    out[[paste0("b_idx_S", s)]] <- if (is.null(r)) NA_real_ else r$b_idx
  }
  for (s in c(5, 7, 8, 9, 10, 11)) {
    r <- res$strategies[[s]]
    out[[paste0("gini_g1_S", s)]] <- if (is.null(r)) NA_real_ else r$gini_g1
  }
  # concentration metrics for the static Tier-B gate (round-1 grant vector, S5)
  r5 <- res$strategies[[5]]
  if (!is.null(r5)) {
    g1 <- r5$g_rounds[[1]]
    out$coverage_S5    <- mean(g1 > 1e-8)
    out$top10_share_S5 <- if (sum(g1) > 0) {
      k <- max(1L, floor(length(g1) / 10))
      sum(sort(g1, decreasing = TRUE)[seq_len(k)]) / sum(g1)
    } else NA_real_
  } else {
    out$coverage_S5 <- NA_real_; out$top10_share_S5 <- NA_real_
  }
  # D2: budget share of the S8 allocation displaced by the S11 floor (round 1)
  r8 <- res$strategies[[8]]; r11 <- res$strategies[[11]]
  out$disp_S11 <- if (!is.null(r8) && !is.null(r11)) {
    sum(abs(r11$g_rounds[[1]] - r8$g_rounds[[1]])) / res$params$B_total
  } else NA_real_
  out$rho_s <- res$rho_s
  # derived contrasts (exact names — downstream depends on them)
  o <- function(s) out[[paste0("out_S", s)]]
  out$fwd_vs_myo_P     <- o(7) - o(4)
  out$fwd_vs_myo_PG    <- o(8) - o(5)
  out$fwd_vs_myo_PS    <- o(9) - o(6)
  out$signal_myo       <- o(5) - o(4)
  out$signal_fwd       <- o(8) - o(7)
  out$seed_myo         <- o(6) - o(4)
  out$seed_fwd         <- o(9) - o(7)
  out$optimal_vs_naive <- o(8) - o(3)
  out$seed_sig_myo     <- o(10) - o(5)   # cost of the floor on the myopic+signal family
  out$seed_sig_fwd     <- o(11) - o(8)   # cost of the floor on the forward+signal family
  out
}

# Long-format schedule rows for one run: one row per (strategy, round).
extract_schedule_long <- function(res, cell_id, trial) {
  rows <- list()
  for (s in 1:11) {
    r <- res$strategies[[s]]
    if (is.null(r)) next
    Tr <- length(r$g_rounds)
    for (t in seq_len(Tr)) {
      g_t <- r$g_rounds[[t]]
      rows[[length(rows) + 1L]] <- data.frame(
        cell_id = cell_id, trial = trial, strategy = s, round = t,
        alpha_t = r$alpha_t[t], mean_g_t = mean(g_t), gini_g_t = gini(g_t)
      )
    }
  }
  do.call(rbind, rows)
}

# ============================================================
# RUNNER
# ============================================================
run_sweep_T <- function(sweep_name, seeds = 1:200, base_params = NULL,
                        cores = SWEEP_CORES) {
  config <- SWEEP_CONFIGS_T[[sweep_name]]
  if (is.null(config)) stop("Unknown sweep: ", sweep_name)
  grid <- config$grid_fn()
  bp   <- modifyList(SWEEP_BASE_PARAMS_T, base_params %||% list())

  # flatten (cell, seed) task list
  tasks <- list()
  for (ci in seq_len(nrow(grid))) {
    row_vals <- as.list(grid[ci, , drop = FALSE])
    for (sd in seeds) tasks[[length(tasks) + 1L]] <- list(cell_id = ci, seed = sd, row_vals = row_vals)
  }

  run_one <- function(task) {
    run_params <- modifyList(bp, task$row_vals)
    run_params$seed <- task$seed
    res <- do.call(run_simulation_T, run_params)
    wide <- as.data.frame(c(task$row_vals, list(cell_id = task$cell_id, seed = task$seed),
                            extract_metrics_T(res)), stringsAsFactors = FALSE)
    list(wide = wide, long = extract_schedule_long(res, task$cell_id, task$seed))
  }

  out <- if (cores > 1L) mclapply(tasks, run_one, mc.cores = cores)
         else lapply(tasks, run_one)

  # surface any worker errors instead of silently dropping
  errs <- vapply(out, function(x) inherits(x, "try-error") || is.null(x$wide), logical(1))
  if (any(errs)) stop(sprintf("%d/%d tasks failed in sweep %s", sum(errs), length(out), sweep_name))

  list(raw      = do.call(rbind, lapply(out, `[[`, "wide")),
       raw_long = do.call(rbind, lapply(out, `[[`, "long")))
}

# ============================================================
# DRIVER (checkpointed): skip sweeps whose summary already exists.
# ============================================================
sweep_one_T <- function(sweep_name, seeds = 1:200, base_params = NULL,
                        out_dir = "sweep_results/T_run", save = TRUE, cores = SWEEP_CORES,
                        resume = TRUE) {
  config <- SWEEP_CONFIGS_T[[sweep_name]]
  summ_path <- file.path(out_dir, paste0(sweep_name, "_summary.rds"))
  if (resume && file.exists(summ_path)) {
    cat(sprintf("[skip] %s (summary exists)\n", sweep_name)); return(invisible(NULL))
  }
  cat(sprintf("\n=== Sweep: %s (tier %d) ===\n%s\n", config$name, config$tier, config$description))
  t0 <- Sys.time()
  r  <- run_sweep_T(sweep_name, seeds = seeds, base_params = base_params, cores = cores)
  summ <- summarize_sweep(r$raw, config)
  summ <- summ[, !grepl("^cell_id_", names(summ)), drop = FALSE]   # drop noise agg
  dt <- as.numeric(Sys.time() - t0, units = "secs")
  cat(sprintf("  %d cells x %d seeds done in %.0fs (%.3fs/run)\n",
              nrow(config$grid_fn()), length(seeds), dt, dt / nrow(r$raw)))

  if (save) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(r$raw,      file.path(out_dir, paste0(sweep_name, "_raw.rds")))
    saveRDS(r$raw_long, file.path(out_dir, paste0(sweep_name, "_rawlong.rds")))
    saveRDS(summ,       summ_path)
    p1 <- tryCatch(plot_sweep_from_config(summ, config, "primary"),  error = function(e) NULL)
    p2 <- tryCatch(plot_sweep_from_config(summ, config, "secondary"), error = function(e) NULL)
    if (!is.null(p1)) ggsave(file.path(out_dir, paste0(sweep_name, "_primary.png")), p1, width = 8, height = 5.5, dpi = 120)
    if (!is.null(p2)) ggsave(file.path(out_dir, paste0(sweep_name, "_secondary.png")), p2, width = 8, height = 5.5, dpi = 120)
    cat(sprintf("  saved -> %s/\n", out_dir))
  }
  invisible(list(raw = r$raw, raw_long = r$raw_long, summary = summ, config = config, secs = dt))
}

main_sweep_T <- function(seeds = 1:200, base_params = NULL,
                         out_dir = "sweep_results/T_run", cores = SWEEP_CORES) {
  for (name in names(SWEEP_CONFIGS_T)) {
    sweep_one_T(name, seeds = seeds, base_params = base_params, out_dir = out_dir, cores = cores)
  }
  cat("\nAll T-round sweeps complete ->", out_dir, "\n")
}
