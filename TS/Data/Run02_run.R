# Run EvoDom program and add some results to file,

source("h5_funcs.R")

reps <- 100

for (rep in 1:reps) {

  system("../build/EvoLearn Run02_inp.toml")

  d1 <- h5_dt("Run02_pop.h5")
  av_alive <- with(d1, mean(alive))
  av_lambda <- with(d1, mean(lambda))
  sd_lambda <- with(d1, sd(lambda))
  av_thr1 <- with(d1, mean(thr1))
  sd_thr1 <- with(d1, sd(thr1))
  av_thr2 <- with(d1, mean(thr2))
  sd_thr2 <- with(d1, sd(thr2))
  av_t_thr <- with(d1, mean(t_thr))
  sd_t_thr <- with(d1, sd(t_thr))
  dr <- data.frame(av_alive = av_alive,
                   av_lambda = av_lambda,
                   av_thr1 = av_thr1,
                   av_thr2 = av_thr2,
                   av_t_thr = av_t_thr,
                   sd_lambda = sd_lambda,
                   sd_thr1 = sd_thr1,
                   sd_thr2 = sd_thr2,
                   sd_t_thr = sd_t_thr)

  dr1 <- read.delim("Run02_data.csv", sep = ",")
  dr1 <- rbind(dr1, dr)
  write.table(dr1, "Run02_data.csv", quote = FALSE,
              sep = ",", row.names = FALSE)
}
