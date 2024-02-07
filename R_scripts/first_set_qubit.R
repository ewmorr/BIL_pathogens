library(dplyr)

dat = read.csv("data/first_set_qubit_11082023.csv")
dat
str(dat)
dat = dat %>% filter(volume != "blank")
str(dat)
dat$volume = as.numeric(dat$volume)
str(dat)

mod1 = lm(concentration ~ site*volume, data = dat)
mod1
summary(mod1)

cor(dat[c("concentration", "volume")])
cor.test(dat$volume, dat$concentration)

plot(dat$volume ~ dat$concentration)
     