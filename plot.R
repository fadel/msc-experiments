# Functions for plotting results from manipulation experiments.
# NOTE: This script should only be used after results are generated with run.R

library(cowplot)
library(gridExtra)
library(logging)
library(reshape2)
library(scales)

source("util.R")

# Plots results for experiments of all techniques for a given measure and
# dataset.
plot.measure <- function(ds, techniques, output.dir, measure, n.iter=30) {
  measure.df <- data.frame()
  scale.labels <- c()

  for (tech in techniques) {
    Y.measure  <- read.table(file.path(output.dir, ds$name, tech$name, paste(measure$name, "Y.tbl", sep="-")))
    Ym.measure <- read.table(file.path(output.dir, ds$name, tech$name, paste(measure$name, "Ym.tbl", sep="-")))

    measure.df <- rbind(measure.df,
                        data.frame(name=tech$name.pretty,
                                   type=paste("Y",  tech$name, sep="."),
                                   V1=Ym.measure - Y.measure)
                       )

    scale.labels <- c(scale.labels, tech$name.pretty)
  }

  p <- ggplot(measure.df) +
         background_grid(major="xy", minor="none") +
         theme(legend.position="none") +
         labs(x="", y=measure$name.pretty) +
         geom_boxplot(aes(type, V1, fill=name)) +
         scale_fill_brewer(palette="Set1", guide=guide_legend(title="")) +
         scale_y_continuous(limits=c(min(0, min(measure.df$V1)), max(0, max(measure.df$V1)))) +
         scale_x_discrete(labels=scale.labels)

  fname <- file.path(output.dir, "plots", ds$name, paste(measure$name, "pdf", sep="."))
  loginfo("Saving plot: %s", fname)
  save_plot(fname, p, base_aspect_ratio=1.3)
}

# Plot boxplots of techniques, one per measure diff per dataset.
plot.measures <- function(datasets, techniques, measures, output.dir, n.iter=30) {
  dir.create.safe(file.path(output.dir, "plots"))

  for (ds in datasets) {
    dir.create.safe(file.path(output.dir, "plots", ds$name))
    for (measure in measures) {
      if (is.null(ds$labels.file) && measure$name == "silhouette") {
        next
      }

      plot.measure(ds, techniques, output.dir, measure, n.iter)
    }
  }
}

# Same as above, but averages over all datasets.
plot.averages <- function(datasets, techniques, measures, output.dir, n.iter=30) {
  dir.create.safe(file.path(output.dir, "plots"))

  for (measure in measures) {
    measure.df <- data.frame()
    scale.labels <- c()
    for (tech in techniques) {
      measure.avg <- rep(0, n.iter)
      scale.labels <- c(scale.labels, tech$name.pretty)

      for (ds in datasets) {
        if (is.null(ds$labels.file) && measure$name == "silhouette") {
          next
        }

        Y.measure  <- read.table(file.path(output.dir, ds$name, tech$name, paste(measure$name, "Y.tbl", sep="-")))$V1
        Ym.measure <- read.table(file.path(output.dir, ds$name, tech$name, paste(measure$name, "Ym.tbl", sep="-")))$V1
        measure.avg <- measure.avg + (Ym.measure - Y.measure)
      }
      measure.avg <- measure.avg / length(datasets)
      measure.df <- rbind(measure.df, data.frame(tech=tech$name, V1=measure.avg))
    }

    p <- ggplot(measure.df) +
           background_grid(major="xy", minor="none") +
           theme(legend.position="none") +
           labs(x="", y=measure$name.pretty) +
           geom_boxplot(aes(tech, V1, fill=tech)) +
           scale_fill_brewer(palette="Set1", guide=guide_legend(title="")) +
           scale_y_continuous(limits=c(min(0, min(measure.df$V1)), max(0, max(measure.df$V1)))) +
           scale_x_discrete(labels=scale.labels)

    fname <- file.path(output.dir, "plots", paste(measure$name, "pdf", sep="."))
    loginfo("Saving plot: %s", fname)
    save_plot(fname, p, base_aspect_ratio=1.3)
  }
}

# Plot a single scatterplot of techniques and datasets, where x axis is the
# measure before manipulation and y axis is the measure after manipulation.
# Also adds a y=x line so that visual inspection is easier.
plot.scatter.measure <- function(measure, datasets, techniques, output.dir, n.iter=30) {
  measure.df <- data.frame()
  for (tech in techniques) {
    for (ds in datasets) {
      if (is.null(ds$labels.file) && measure$name == "silhouette") {
        next
      }

      base.path <- file.path(output.dir, ds$name, tech$name)
      fname <- file.path(base.path, paste(measure$name, "Y.tbl", sep="-"))
      Y.measure  <- read.table(fname)$V1
      fname <- file.path(base.path, paste(measure$name, "Ym.tbl", sep="-"))
      Ym.measure <- read.table(fname)$V1
      measure.df <- rbind(measure.df, data.frame(tech=tech$name.pretty,
                                                 dataset=ds$name.pretty,
                                                 x=mean(Y.measure),
                                                 y=mean(Ym.measure)))
    }
  }

  min.max <- min(max(measure.df$x), max(measure.df$y))
  p <- ggplot(measure.df) +
         background_grid(major="xy", minor="none") +
         theme(legend.position="right") +
         labs(x=paste(measure$name.pretty, "(before)", sep=" "),
              y=paste(measure$name.pretty, "(after)", sep=" ")) +
         geom_point(aes(x=x, y=y, color=tech, shape=dataset), alpha=0.8, size=3) +
         scale_color_brewer(palette="Set1", guide=guide_legend(title="Technique")) +
         scale_shape(guide=guide_legend(title="Dataset")) +
         geom_abline(intercept=0, slope=1)

  fname <- file.path(output.dir, "plots", paste(measure$name, "-scatter", ".pdf", sep=""))
  loginfo("Saving plot: %s", fname)
  save_plot(fname, p, base_aspect_ratio=1.5)

  p <- p +
    scale_x_log10(breaks=trans_breaks("log10", function(x) 10^x),
                  labels=trans_format("log10", math_format(10^ .x))) +
    scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),
                  labels=trans_format("log10", math_format(10^ .x))) +
    annotation_logticks()
  fname <- file.path(output.dir, "plots", paste(measure$name, "-scatter-log", ".pdf", sep=""))
  loginfo("Saving plot: %s", fname)
  save_plot(fname, p, base_aspect_ratio=1.5)

  p
}

# This function runs the scatterplot function above for all measures
plot.scatter <- function(datasets, techniques, measures, output.dir, n.iter=30) {
  dir.create.safe(file.path(output.dir, "plots"))

  for (measure in measures) {
    p <- plot.scatter.measure(measure, datasets, techniques, output.dir, n.iter)
  }
}

# Plot a single barplot of techniques and datasets, where the y axis shows the
# difference Ym-Y, with confidence intervals.
plot.ci.measure <- function(measure, datasets, techniques, output.dir, n.iter=30) {
  measure.df <- data.frame()
  for (tech in techniques) {
    for (ds in datasets) {
      if (is.null(ds$labels.file) && measure$name == "silhouette") {
        next
      }

      base.path <- file.path(output.dir, ds$name, tech$name)
      fname <- paste(measure$name, "Y.tbl", sep="-")
      Y.measure  <- read.table(file.path(base.path, fname))$V1
      fname <- paste(measure$name, "Ym.tbl", sep="-")
      Ym.measure <- read.table(file.path(base.path, fname))$V1
      measure.df <- rbind(measure.df, data.frame(tech=tech$name.pretty,
                                                 dataset=ds$name.pretty,
                                                 y=Ym.measure - Y.measure))
    }
  }

  p <- ggplot(measure.df) +
         background_grid(major="xy", minor="none") +
         theme(legend.position="right") +
         labs(x="", y=measure$name.pretty) +
         stat_summary(aes(x=tech, y=y, color=tech, shape=dataset), fun.data=ci.fun, position=position_dodge(width=0.75)) +
         scale_color_brewer(palette="Set1", guide=guide_legend(title="Technique")) +
         scale_shape(guide=guide_legend(title="Dataset")) +
         scale_x_discrete(expand=c(0, 0.01))

  fname <- file.path(output.dir, "plots", paste(measure$name, "-ci", ".pdf", sep=""))
  loginfo("Saving plot: %s", fname)
  save_plot(fname, p, base_aspect_ratio=1.65)

  p <- p + scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),
                         labels=trans_format("log10", math_format(10^ .x))) +
    annotation_logticks(sides="l")

  fname <- file.path(output.dir, "plots", paste(measure$name, "-ci-log", ".pdf", sep=""))
  loginfo("Saving plot: %s", fname)
  save_plot(fname, p, base_aspect_ratio=1.65)


  p
}

# This function runs the function above for all measures
plot.ci <- function(datasets, techniques, measures, output.dir, n.iter=30) {
  dir.create.safe(file.path(output.dir, "plots"))

  for (measure in measures) {
    p <- plot.ci.measure(measure, datasets, techniques, output.dir, n.iter)
  }
}

# Plot a single scatterplot of techniques and datasets, where x axis is the
# measure anipulation.
plot.evo.measure <- function(measure, datasets, techniques, output.dir) {
  measure.df <- data.frame()
  for (ds in datasets) {
    if (is.null(ds$labels.file) && measure$name == "silhouette") {
      next
    }
    for (tech in techniques) {
      base.path <- file.path(output.dir, ds$name, tech$name)
      fname <- paste(measure$name, "Ys-evo.tbl", sep="-")
      Ys.measure  <- read.table(file.path(base.path, fname))$V1
      fname <- paste(measure$name, "Y-evo.tbl", sep="-")
      Y.measure  <- read.table(file.path(base.path, fname))$V1
      measure.df <- rbind(measure.df, data.frame(tech=tech$name.pretty,
                                                 dataset=ds$name.pretty,
                                                 key=paste(tech$name, ds$name, sep="-"),
                                                 x=rev(Ys.measure),
                                                 y=rev(Y.measure)))
    }
  }

  p <- ggplot(measure.df) +
         background_grid(major="xy", minor="none") +
         theme(legend.position="right") +
         labs(x=paste(measure$name.pretty, "(Ys)", sep=" "),
              y=paste(measure$name.pretty, "(Y)", sep=" ")) +
         geom_point(aes(x=x, y=y, color=tech, shape=dataset), alpha=0.8, size=1.5) +
         geom_path(aes(x=x, y=y, color=tech, group=key)) +
         scale_color_brewer(palette="Set1", guide=guide_legend(title="Technique")) +
         scale_shape(guide=guide_legend(title="Dataset")) +
         geom_abline(intercept=0, slope=1)

  fname <- file.path(output.dir, "plots", paste(measure$name, "-evo", ".pdf", sep=""))
  loginfo("Saving plot: %s", fname)
  save_plot(fname, p, base_aspect_ratio=1.5)

  p
}

# This function runs the scatterplot function above for all measures
plot.evo <- function(datasets, techniques, measures, output.dir) {
  dir.create.safe(file.path(output.dir, "plots"))

  for (measure in measures) {
    p <- plot.evo.measure(measure, datasets, techniques, output.dir)
  }
}


# Experiment configuration
# Defines: datasets, techniques, measures, output.dir
source("config.R")

args <- commandArgs(T)

# logging setup
basicConfig()
addHandler(writeToFile,
           file=args[1],
           level="FINEST")

#plot.measures(datasets, techniques, measures, output.dir)
#plot.averages(datasets, techniques, measures, output.dir)
#plot.scatter(datasets, techniques, measures, output.dir)
#plot.ci(datasets, techniques, measures, output.dir)
plot.evo(datasets, techniques, measures, output.dir)
