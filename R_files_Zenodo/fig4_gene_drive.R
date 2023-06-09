library(ggplot2)
library(ggpubr)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

c95mu7 = read.csv("gd_double_c95_mu7.csv", header = TRUE)
c95mu4 = read.csv("gd_double_c95_mu4.csv", header = TRUE)
c7mu7 = read.csv("gd_double_c7_mu7.csv", header = TRUE)
c7mu4 = read.csv("gd_double_c7_mu4.csv", header = TRUE)

c95mu7$sigma = as.factor(c95mu7$sigma)
c95mu4$sigma = as.factor(c95mu4$sigma)
c7mu7$sigma = as.factor(c7mu7$sigma)
c7mu4$sigma = as.factor(c7mu4$sigma)

c95mu7.plt <- ggplot(data=c95mu7) + theme_bw() +
  geom_line(aes(x = T, y = freq, linetype=type, color=sigma)) + 
  scale_linetype_manual(name="Variant",values=c(1,3)) +
  scale_color_manual(name="Inroduction rate",values=gg_color_hue(2)) +
  ylab("freq") +
  xlab("")

c95mu4.plt <- ggplot(data=c95mu4) + theme_bw() +
  geom_line(aes(x = T, y = freq, linetype=type, color=sigma)) + 
  scale_linetype_manual(name="Variant",values=c(1,3)) +
  scale_color_manual(name="Inroduction rate",values=gg_color_hue(2)) +
  ylab("") +
  xlab("")

c7mu7.plt <- ggplot(data=c7mu7) + theme_bw() +
  geom_line(aes(x = T, y = freq, linetype=type, color=sigma)) + 
  scale_linetype_manual(name="Variant",values=c(1,3)) +
  scale_color_manual(name="Inroduction rate",values=gg_color_hue(2)) +
  ylab("freq1") +
  xlab("time")

c7mu4.plt <- ggplot(data=c7mu4) + theme_bw() +
  geom_line(aes(x = T, y = freq, linetype=type, color=sigma)) + 
  scale_linetype_manual(name="Variant",values=c(1,3)) +
  scale_color_manual(name="Inroduction rate",values=gg_color_hue(2)) +
  ylab("") +
  xlab("time")

plt <- ggarrange(c95mu7.plt, c95mu4.plt, c7mu7.plt, c7mu4.plt, labels="auto",common.legend=TRUE, legend="right")

ggsave("fig5.pdf", plt, device = "pdf", )
plt
