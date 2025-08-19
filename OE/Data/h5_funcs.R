# functions to read data from hdf5 file created by simulation program

# return data table for "one-column" phenotype data
h5_dt <- function(hf_name) {
  require(hdf5r)
  require(data.table)
  f_h5 <- H5File$new(hf_name, mode = "r") # nolint: object_name_linter.
  w0 <- f_h5[["w0"]][]
  alph_w <- f_h5[["alph_w"]][]
  lambda <- f_h5[["lambda"]][]
  thr1 <- f_h5[["thr1"]][]
  thr2 <- f_h5[["thr2"]][]
  t_thr <- f_h5[["t_thr"]][]
  beta <- f_h5[["beta"]][]
  Rew <- f_h5[["Rew"]][]
  R_ewm <- f_h5[["R_ewm"]][]
  R_tot <- f_h5[["R_tot"]][]
  delt <- f_h5[["delt"]][]
  choice <- as.integer(f_h5[["choice"]][] + 1)
  n_Rew <- f_h5[["n_Rew"]][]
  n_off <- f_h5[["n_off"]][]
  i_num <- as.integer(f_h5[["i_num"]][] + 1)
  g_num <- as.integer(f_h5[["g_num"]][] + 1)
  female <- f_h5[["female"]][]
  alive <- f_h5[["alive"]][]
  f_h5$close_all()
  data.table(
    w0 = w0, alph_w= alph_w, 
    lambda = lambda,
    thr1 = thr1, thr2 = thr2,
    t_thr = t_thr, beta = beta,
    Rew = Rew, R_ewm = R_ewm, R_tot = R_tot,
    delt = delt, choice = choice,
    n_Rew = n_Rew, n_off = n_off,
    i_num = i_num, g_num = g_num, 
    female = female, alive = alive
  )
}

# return matrix where each row is a parameter w
h5_w <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  w <- t(f_h5[["w"]][, ])
  f_h5$close_all()
  w
}

# return matrix where each row is a parameter att
h5_att <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  att <- t(f_h5[["att"]][, ])
  f_h5$close_all()
  att
}

# # return matrix where each row is a parameter v_ewm
# h5_v_ewm <- function(hf_name) {
#   require(hdf5r)
#   f_h5 <- H5File$new(hf_name, mode = "r")
#   v_ewm <- t(f_h5[["v_ewm"]][, ])
#   f_h5$close_all()
#   v_ewm
# }

# # return matrix where each row is a parameter s_ewm
# h5_s_ewm <- function(hf_name) {
#   require(hdf5r)
#   f_h5 <- H5File$new(hf_name, mode = "r")
#   s_ewm <- t(f_h5[["s_ewm"]][, ])
#   f_h5$close_all()
#   s_ewm
# }

# return matrix where each row is a parameter n_feat
h5_n_feat <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  n_feat <- t(f_h5[["n_feat"]][, ])
  f_h5$close_all()
  n_feat
}

# return matrix where each row is a parameter l_feat
h5_l_feat <- function(hf_name) {
  require(hdf5r)
  f_h5 <- H5File$new(hf_name, mode = "r")
  l_feat <- t(f_h5[["l_feat"]][, ])
  f_h5$close_all()
  l_feat
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
  delt <- f_h5[["delt"]][]
  R_ewm <- f_h5[["R_ewm"]][]
  w <- f_h5[["w"]][]
  hs <- f_h5[["hs"]][]
  w_tr <- f_h5[["w_tr"]][]
  regret <- f_h5[["regret"]][]
  f_h5$close_all()
  data.table(
    i_num = i_num, g_num = g_num, lph = lph, tstep = tstep,
    choice = choice, n_off = n_off, Rew = Rew, delt = delt,
    R_ewm = R_ewm, w = w,
    hs = hs, w_tr = w_tr, regret = regret
  )
}
