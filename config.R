# All files generated will be saved under this directory
output.dir <- "results"

datasets <- list(
                 list(name="wdbc",         name.pretty="WDBC",         data.file="datasets/wdbc.tbl",             labels.file="datasets/wdbc.labels",             scale=T),
                 list(name="mnist-2000",   name.pretty="MNIST",        data.file="datasets/mnist_2000.tbl",       labels.file="datasets/mnist_2000.labels",       scale=F),
                 list(name="segmentation", name.pretty="Segmentation", data.file="datasets/segmentation.tbl",     labels.file="datasets/segmentation.labels",     scale=T),
                 list(name="newsgroups",   name.pretty="Newsgroups",   data.file="datasets/newsgroups-500-3.tbl", labels.file="datasets/newsgroups-500-3.labels", scale=F),
                 list(name="faces",        name.pretty="Faces",        data.file="datasets/faces.tbl",            labels.file=NULL,                               scale=F)
                )

techniques <- list(
                   list(name="lamp-global",   name.pretty="LAMP-g",   fn="lamp", args=list()),
                   list(name="lamp-local-25", name.pretty="LAMP-25%", fn="lamp", args=list(cp=0.25)),
                   list(name="plmp",          name.pretty="PLMP",     fn="plmp", args=list()),
                   list(name="lsp",           name.pretty="LSP",      fn="lsp",  args=list(k=15)),
                   list(name="pekalska",      name.pretty="Pekalska", fn="pekalska.wrapper", args=list())
                  )

measures <- list(
                 list(name="silhouette", name.pretty="Silhouette"),
                 list(name="np",         name.pretty="Neighborhood preservation"),
                 list(name="stress",     name.pretty="Stress"),
                 list(name="precision",  name.pretty="Smoothed precision"),
                 list(name="recall",     name.pretty="Smoothed recall")
                )
