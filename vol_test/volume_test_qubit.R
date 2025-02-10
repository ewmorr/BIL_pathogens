library(dplyr)
library(ggplot2)
library(RColorBrewer)
source("R_scripts/ggplot_theme.txt")

conc = read.csv("data/volume_test_qubit_01022023.csv")
head(conc)

conc %>% filter(volume == "blank")
conc = conc %>% filter(volume != "blank")
conc$ID = paste(conc$site, conc$date, sep = " ")
conc$volume = as.numeric(conc$volume)

str(conc)

mod1 = lm(concentration ~ ID + volume, data = conc)
mod1
summary(mod1)
anova(mod1)

mod2 = lm(concentration ~ ID * volume, data = conc)
mod2
summary(mod2)
anova(mod2)
#no significant interaction
anova(mod1, mod2)

summary(mod1)
anova(mod1)

mod3 = lm(concentration ~ ID/volume, data = conc)
summary(mod3)
anova(mod3)

cor(conc[c("concentration", "volume")])
cor.test(conc$volume, conc$concentration)

plot(conc$volume ~ conc$concentration)


#add regression lines to the df based on the best model 
conc = cbind(conc, pred = predict(mod1))

#plot
p1 = ggplot(conc, aes(x = volume, y = concentration, color = ID)) +
    #geom_smooth(method = "lm", se = F, formula = ) +
    geom_line(mapping=aes(y=pred)) +
    geom_point(size = 2, position = position_jitter(width = 0.5, height = 0)) +
    scale_color_brewer(palette = "Paired") +
    scale_x_continuous(breaks = c(5,10,20,40,80)) +
    labs(
        y = expression(paste("DNA concentration (ng ",mu,"l"^-1,")")),
        x = "Volume trap fluid (ml)"     
    ) +
    my_gg_theme +
    theme(
        legend.text = element_text(size = 15)
    )
p1

pdf("figures/volume_test_DNA_conc.same_slope.pdf", width = 8, height = 6)
p1
dev.off()

mod4 = lm(concentration ~ ID, data = conc)
summary(mod4)
summary(mod1)

mod5 = lm(concentration ~ volume, data = conc)
summary(mod5)
summary(mod1)
