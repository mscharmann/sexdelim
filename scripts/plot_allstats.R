library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

mydata <- read.table(args[1], sep = "\t", header = F, stringsAsFactors = FALSE)
colnames(mydata) <- c("chr","start","stop","log2_M_F_cov_ratio", "F_spec_kmers", "M_spec_kmers", "F_pi", "M_pi", "MvF_Fst", "MvF_dxy", "LD_r2", "GWAS_n_signif", "ZW_div", "XY_div", "het_m", "het_f", "NetDiv_F", "NetDiv_M")
mydata[mydata=="."] <- NA
mydata[, 2:ncol(mydata)] <- sapply(mydata[, 2:ncol(mydata)], as.numeric)
mydata$chr <- as.factor(mydata$chr)

top <- ggplot(data=mydata, aes(x=start, y=log2_M_F_cov_ratio, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x')

mid1 <- ggplot(data=mydata, aes(x=start, y=F_spec_kmers, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid2 <- ggplot(data=mydata, aes(x=start, y=M_spec_kmers, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid3 <- ggplot(data=mydata, aes(x=start, y=F_pi, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid4 <- ggplot(data=mydata, aes(x=start, y=M_pi, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid5 <- ggplot(data=mydata, aes(x=start, y=MvF_Fst, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid6 <- ggplot(data=mydata, aes(x=start, y=MvF_dxy, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid7 <- ggplot(data=mydata, aes(x=start, y=LD_r2, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid8 <- ggplot(data=mydata, aes(x=start, y=GWAS_n_signif, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid9 <- ggplot(data=mydata, aes(x=start, y=het_m, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid10 <- ggplot(data=mydata, aes(x=start, y=het_f, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid11 <- ggplot(data=mydata, aes(x=start, y=ZW_div, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid12 <- ggplot(data=mydata, aes(x=start, y=XY_div, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid13 <- ggplot(data=mydata, aes(x=start, y=NetDiv_F, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

bottom <- ggplot(data=mydata, aes(x=start, y=NetDiv_M, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.text.x=element_text(size=2, angle = 35)) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())


pdf(args[2], width = 20, height = 24)
plot_grid(top, mid1, mid2, mid3, mid4, mid5, mid6, mid7, mid8, mid9, mid10, mid11, mid12, mid13, bottom, ncol=1, align = "v")
dev.off()

