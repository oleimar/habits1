# functions to read data from hdf5 file created by simulation program

# return data table for "one-column" phenotype data
h5_dt <- function(hf_name) {
  require(hdf5r)
  require(data.table)
  f_h5 <- H5File$new(hf_name, mode = "r") # nolint: object_name_linter.
  lambda <- f_h5[["lambda"]][]
  thr1 <- f_h5[["thr1"]][]
  thr2 <- f_h5[["thr2"]][]
  t_thr <- f_h5[["t_thr"]][]
  Rew <- f_h5[["Rew"]][]
  R_ewm <- f_h5[["R_ewm"]][]
  R_tot <- f_h5[["R_tot"]][]
  choice <- as.integer(f_h5[["choice"]][] + 1)
  n_Rew <- f_h5[["n_Rew"]][]
  n_off <- f_h5[["n_off"]][]
  i_num <- as.integer(f_h5[["i_num"]][] + 1)
  g_num <- as.integer(f_h5[["g_num"]][] + 1)
  female <- f_h5[["female"]][]
  alive <- f_h5[["alive"]][]
  f_h5$close_all()
  data.table(
    lambda = lambda,
    thr1 = thr1, thr2 = thr2,
    t_thr = t_thr,
    Rew = Rew, R_ewm = R_ewm, R_tot = R_tot,
    choice = choice,
    n_Rew = n_Rew, n_off = n_off,
    i_num = i_num, g_num = g_num, 
    female = female, alive = alive
  )
}

# return matrix where each row is a parameter y_bar
h5_y_bar <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  y_bar <- t(f_h5[["y_bar"]][, ])
  f_h5$close_all()
  y_bar
}

# return matrix where each row is a parameter rho_hat
h5_rho_hat <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  rho_hat <- t(f_h5[["rho_hat"]][, ])
  f_h5$close_all()
  rho_hat
}

# return matrix where each row is a parameter tau_hat
h5_tau_hat <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  tau_hat <- t(f_h5[["tau_hat"]][, ])
  f_h5$close_all()
  tau_hat
}

# return matrix where each row is a parameter att
h5_att <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  att <- t(f_h5[["att"]][, ])
  f_h5$close_all()
  att
}

# return matrix where each row is a parameter hs
h5_hs <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  hs <- t(f_h5[["hs"]][, ])
  f_h5$close_all()
  hs
}

# return matrix where each row is a parameter n_chs
h5_n_chs <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  n_chs <- t(f_h5[["n_chs"]][, ])
  f_h5$close_all()
  n_chs
}

# return matrix where each row is a maternal gamete value
h5_mat_gam <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  mat_gam <- t(f_h5[["MatGam"]][, ])
  f_h5$close_all()
  mat_gam
}

# return matrix where each row is a paternal gamete value
h5_pat_gam <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  pat_gam <- t(f_h5[["PatGam"]][, ])
  f_h5$close_all()
  pat_gam
}

# return data table for learning history data
h5_hdt <- function(hf_name) {
  require(hdf5r)
  require(data.table)
  f_h5 <- H5File$new(hf_name, mode = "r")
  i_num <- as.integer(f_h5[["i_num"]][] + 1)
  g_num <- as.integer(f_h5[["g_num"]][] + 1)
  lph <- as.integer(f_h5[["lph"]][] + 1)
  tstep <- as.integer(f_h5[["tstep"]][] + 1)
  choice <- as.integer(f_h5[["choice"]][] + 1)
  n_off <- f_h5[["n_off"]][]
  Rew <- f_h5[["Rew"]][]
  R_ewm <- f_h5[["R_ewm"]][]
  hs <- f_h5[["hs"]][]
  w_tr <- f_h5[["w_tr"]][]
  regret <- f_h5[["regret"]][]
  f_h5$close_all()
  data.table(
    i_num = i_num, g_num = g_num, lph = lph, tstep = tstep,
    choice = choice, n_off = n_off, Rew = Rew, R_ewm = R_ewm, 
    hs = hs, w_tr = w_tr, regret = regret
  )
}
