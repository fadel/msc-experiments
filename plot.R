# Functions for plotting results from manipulation experiments.
# NOTE: This script should only be used after results are generated with run.R

require(cowplot)
require(gridExtra)
require(logging)

# Plots results for experiments of all techniques for a given measure and
# datasset.
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

  p
}

# This function runs the scatterplot function above for all measures
plot.scatter <- function(datasets, techniques, measures, output.dir, n.iter=30) {
  dir.create.safe(file.path(output.dir, "plots"))

  for (measure in measures) {
    p <- plot.scatter.measure(measure, datasets, techniques, output.dir, n.iter)
  }
}

# Experiment configuration
# Defines: datasets, techniques, measures, output.dir
source("config.R")

args <- commandArgs(T)

# Logging setup
basicConfig()
addHandler(writeToFile,
           file=args[1],
           level="FINEST")

plot.measures(datasets, techniques, measures, output.dir)
plot.averages(datasets, techniques, measures, output.dir)
plot.scatter(datasets, techniques, measures, output.dir)
