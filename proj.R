# NOTE: This script should only be used after results are generated with run.R

library(grid)
library(gridExtra)
library(ggplot2)
library(mp)

get_legend <- function(p) {
  tmp <- ggplot_gtable(ggplot_build(p))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

args <- commandArgs(T)

dataset <- args[1]
tech    <- args[2]
iter    <- args[3]
measure <- args[4]
out.file <- args[5]

X.file <- file.path("datasets", paste(dataset, "tbl", sep="."))
l.file <- file.path("datasets", paste(dataset, "labels", sep="."))
Y.file <- file.path("results", dataset, tech, paste("Y-", iter, ".tbl", sep=""))
sample.indices.file <- file.path("results", dataset, paste("sample-indices-", iter, ".tbl", sep=""))
Ys.file <- file.path("results", dataset, paste("Ys-", iter, ".tbl", sep=""))
Ysm.file <- file.path("results", dataset, paste("Ysm-", measure, "-", iter, ".tbl", sep=""))

X <- read.table(X.file)
Y <- read.table(Y.file)
l <- as.factor(read.table(l.file)$V1)
l <- l[!duplicated(X)]
X <- scale(unique(X))

sample.indices <- read.table(sample.indices.file)$V1
ls <- l[sample.indices]
Ys <- read.table(Ys.file)
Ysm <- read.table(Ysm.file)
Ym <- pekalska(dist(X), sample.indices, Ysm)

dfYs <- data.frame(x=Ys[,1], y=Ys[,2], Labels=ls)
pYs <- ggplot(dfYs) +
         theme_minimal() +
         theme(legend.position="none",
               plot.background=element_rect(fill="#ffffff", color="#000000"),
               axis.text=element_blank(),
               axis.title=element_blank(),
               panel.grid=element_blank()) +
         geom_point(aes(x=x, y=y, color=Labels), size=0.5, alpha=0.8) +
         scale_color_brewer(palette="Set1")
pYs <- ggplotGrob(pYs)

dfYsm <- data.frame(x=Ysm[,1], y=Ysm[,2], Labels=ls)
pYsm <- ggplot(dfYsm) +
          theme_minimal() +
          theme(legend.position="none",
                plot.background=element_rect(fill="#ffffff", color="#000000"),
                axis.text=element_blank(),
                axis.title=element_blank(),
                panel.grid=element_blank()) +
          geom_point(aes(x=x, y=y, color=Labels), size=0.5, alpha=0.8) +
          scale_color_brewer(palette="Set1")
pYsm <- ggplotGrob(pYsm)

dfY <- data.frame(x=Y[,1], y=Y[,2], Labels=l)
pY <- ggplot(dfY) +
        theme_minimal() +
        theme(legend.position="bottom",
              legend.title=element_blank(),
              legend.text=element_text(size=14),
              legend.background=element_rect(fill="#ffffff", color="#000000"),
              axis.text=element_blank(),
              axis.title=element_blank(),
              panel.grid=element_blank()) +
        geom_point(aes(x=x, y=y, color=Labels), alpha=0.8) +
        scale_color_brewer(palette="Set1") +
        annotation_custom(grob = pYs,
                          xmin=min(dfY$x), xmax=min(dfY$x) + (max(dfY$x) - min(dfY$x)) / 3,
                          ymin=max(dfY$y) - (max(dfY$y) - min(dfY$y)) / 3, ymax=max(dfY$y))

legend <- get_legend(pY)
pY <- pY + theme(legend.position="none")

dfYm <- data.frame(x=Ym[,1], y=Ym[,2], Labels=l)
pYm <- ggplot(dfYm) +
         theme_minimal() +
         theme(legend.position="none",
               axis.text=element_blank(),
               axis.title=element_blank(),
               panel.grid=element_blank()) +
         geom_point(aes(x=x, y=y, color=Labels), alpha=0.8) +
         scale_color_brewer(palette="Set1") +
         annotation_custom(grob = pYsm,
                           xmin=min(dfYm$x), xmax=min(dfYm$x) + (max(dfYm$x) - min(dfYm$x)) / 3,
                           ymin=min(dfYm$y), ymax=min(dfYm$y) + (max(dfYm$y) - min(dfYm$y)) / 3)

p <- grid.arrange(pY, pYm, legend,
                  ncol=2, nrow=2,
                  layout_matrix=rbind(c(3, 3), c(1, 2)),
                  widths=c(5, 5), heights=c(0.5, 4.5))

ggsave(out.file, p, width=10, height=5)
