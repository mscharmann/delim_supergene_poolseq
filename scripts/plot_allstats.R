library(ggplot2)
library(cowplot)

args = commandArgs(trailingOnly=TRUE)

# catch error if table to be read is empty
my_func <- function(x){
		tryCatch({
        mydata <- read.table(x, sep = "\t", header = F, stringsAsFactors = FALSE)
        return(mydata)
        },
        error = function(error_condition) {
        pdf(args[2], width = 20, height = 24)
		
		plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')

		text(x = 0.5, y = 0.5, paste("no data selected"), cex = 1.6, col = "black")
		
		dev.off()
		q()
        }
    )
	}

mydata <- my_func(args[1])	

	
colnames(mydata) <- c("chr","start","stop","norm_coverage_ratio_log2", "p1_spec_kmers", "p2_spec_kmers", "p1_pi", "p2_pi", "dxy", "netDiv_p1", "netDiv_p2","Fst_Hudson","XY_like_div","ZW_like_div")
mydata[mydata=="."] <- NA
mydata[, 2:ncol(mydata)] <- sapply(mydata[, 2:ncol(mydata)], as.numeric)
mydata$chr <- as.factor(mydata$chr)

top <- ggplot(data=mydata, aes(x=start, y=norm_coverage_ratio_log2, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x')

mid1 <- ggplot(data=mydata, aes(x=start, y=p1_spec_kmers, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid2 <- ggplot(data=mydata, aes(x=start, y=p2_spec_kmers, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid3 <- ggplot(data=mydata, aes(x=start, y=p1_pi, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid4 <- ggplot(data=mydata, aes(x=start, y=p2_pi, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid5 <- ggplot(data=mydata, aes(x=start, y=dxy, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid6 <- ggplot(data=mydata, aes(x=start, y=netDiv_p1, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid7 <- ggplot(data=mydata, aes(x=start, y=netDiv_p2, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid8 <- ggplot(data=mydata, aes(x=start, y=Fst_Hudson, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

mid9 <- ggplot(data=mydata, aes(x=start, y=XY_like_div, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.title.x = element_blank(),
	axis.text.x = element_blank()) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())

bottom <- ggplot(data=mydata, aes(x=start, y=ZW_like_div, group = chr)) +
	theme_classic() +
	geom_line(size = 0.1) +
	theme(axis.text.x=element_text(size=2, angle = 35)) +
	facet_grid(~chr, scales = 'free_x', space = 'free_x') +
	theme(strip.background = element_blank(),
	strip.text.x = element_blank())


pdf(args[2], width = 20, height = 24)
plot_grid(top, mid1, mid2, mid3, mid4, mid5, mid6, mid7, mid8, mid9, bottom, ncol=1, align = "v")
dev.off()

