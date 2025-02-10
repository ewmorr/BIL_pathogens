library(dyplr)
library(ggplot2)

conc = read.csv("data/trap_wash_test_qubit_01022023.csv")
head(conc)

mod1 = aov(concentration ~ wash, data = conc)
summary(mod1)
TukeyHSD(mod1)


#summarize for plotting
conc.summary = conc %>% 
    group_by(wash) %>%
    summarize(
        avg = mean(concentration),
        se = sd(concentration)/sqrt(n()) 
    )

conc.summary$sig = c("b","b","a")
conc.summary$wash = factor(conc.summary$wash, levels = c("pre", "mid", "post"))

    
p1 = ggplot(conc.summary, aes(x = wash, y = avg)) +
    geom_point() +
    geom_errorbar(aes(ymin = avg-se, ymax = avg+se), width = 0.1 ) +
    scale_x_discrete(labels = c("pre" = "pre-wash", "mid" = "hose spray", "post" = "bleach + hose") ) +
    geom_text(aes(label = sig, y = avg+se+0.1)) +
    theme_bw() +
    labs(y = expression(paste("DNA concentration (ng ",mu,"l"^-1,")"))) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black")
    )
p1                     

ggplot(conc, aes(x = wash, y = concentration)) +
    geom_point() +
    #geom_errorbar(aes(ymin = avg-se, ymax = avg+se), width = 0.1 ) +
    scale_x_discrete(labels = c("pre" = "pre-wash", "mid" = "hose spray", "post" = "bleach + hose") ) +
    theme_bw() +
    labs(y = expression(paste("DNA concentration (ng ",mu,"l"^-1,")"))) +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10, color = "black")
    )


pdf("figures/trap_clean_test.pdf", height = 4, width = 6)
p1
dev.off()
