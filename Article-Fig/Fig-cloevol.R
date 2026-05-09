setwd("E:/#WWY_Project/TJUS_WES/TJUS_Analysis/TJUS_clonevol")
library(clonevol)
library(packcircles)
library(readr)
##输入数据(数据需要as.data.frame)
US041_gaussian_input <- read_delim("01.Data/US041.gaussian.input.tsv", 
                                      delim = "\t", escape_double = FALSE, 
                                    trim_ws = TRUE)
US041_gaussian_input <- as.data.frame(US041_gaussian_input)
vaf.col.names <- grep('.vaf', colnames(US041_gaussian_input), value=T)
sample.names <- gsub('.vaf', '', vaf.col.names)
US041_gaussian_input[, sample.names] <- US041_gaussian_input[, vaf.col.names]
vaf.col.names <- sample.names
sample.groups <- c('P', 'R')
names(sample.groups) <- vaf.col.names
clone.colors <- c('#2E86AB', '#A23B72', '#F18F01', '#C73E1D')



plot.variant.clusters(US041_gaussian_input,
                      cluster.col.name = 'cluster',
                      show.cluster.size = FALSE,
                      cluster.size.text.color = 'blue',
                      vaf.col.names = vaf.col.names,
                      vaf.limits = 200,
                      sample.title.size = 20,
                      violin = FALSE,
                      box = FALSE,
                      jitter = TRUE,
                      jitter.shape = 1,
                      jitter.color = clone.colors,
                      jitter.size = 3,
                      jitter.alpha = 1,
                      jitter.center.method = 'median',
                      jitter.center.size = 1,
                      jitter.center.color = 'darkgray',
                      jitter.center.display.value = 'none',
                      highlight = 'is.driver',
                      highlight.shape = 21,
                      highlight.color = 'blue',
                      highlight.fill.color = 'green',
                      highlight.note.col.name = 'gene',
                      highlight.note.size = 2,
                      order.by.total.vaf = FALSE)


str(US041_gaussian_input)
US041_gaussian_input$cluster <- as.factor(US041_gaussian_input$cluster)
US041_gaussian_input$is.driver <- as.integer(US041_gaussian_input$is.driver)
plot.pairwise(US041_gaussian_input, col.names = vaf.col.names,
              out.prefix = 'variants.pairwise.plot',
              colors = clone.colors)

plot.cluster.flow(US041_gaussian_input, vaf.col.names = vaf.col.names,
                  sample.names = c('Primary', 'Relapse'),
                  colors = clone.colors)






##推断克隆进化树
y = infer.clonal.models(variants = US041_gaussian_input,
                        cluster.col.name = 'cluster',
                        vaf.col.names = vaf.col.names,
                        sample.groups = sample.groups,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = 1,
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        clone.colors = clone.colors,
                        min.cluster.vaf = 0.01,
                        # min probability that CCF(clone) is non-negative
                        sum.p = 0.05,
                        # alpha level in confidence interval estimate for CCF(clone)
                        alpha = 0.05)
##将驱动基因事件映射到树中
y <- transfer.events.to.consensus.trees(y,
                                        US041_gaussian_input[US041_gaussian_input$is.driver,],
                                        cluster.col.name = 'cluster',
                                        event.col.name = 'gene')
##将基于节点的树转换为基于树枝的树
y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')
##把多个地块和树木画在一起
plot.clonal.models(y,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.highlight = 'is.driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'gene',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = 'output',
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 8,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,2))
