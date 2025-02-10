library(ggplot2)

vols = c(5,10,20,40,80)

saturating_fun = function(x, m, h){
    y = m * (x/(x+h))
    return(y)
}
saturating_fun(vols, 80, 20)

spp_num = saturating_fun(vols*10, 20, 5)
spp_num = c(50, 100, 150, 175, 180)

plot(spp_num ~ vols)

p1 = ggplot(data.frame(vols, spp_num), aes(vols, spp_num)) +

    geom_smooth(color = "black") +
    geom_point(size = 2) +
    theme_bw() +
    labs(x = "Volume fluid (ml)", y = "Number species")
p1

pdf("figures/accumulation_curve_example.pdf", width = 6, height = 4)
p1
dev.off()
