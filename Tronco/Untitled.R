#Install Tronco
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TRONCO")

required <- c("TRONCO")

#Vignettes
browseVignettes("TRONCO")

library(TRONCO)
data(aCML)
data(crc_maf)
data(crc_gistic)
data(crc_plain)

dataset = rename.gene(aCML, 'TET2', 'new name')
dataset = rename.type(dataset, 'Ins/Del', 'new type')
as.events(dataset, type = 'new type')

dataset = join.events(aCML, 
                      'gene 4',
                      'gene 88',
                      new.event='test',
                      new.type='banana',
                      event.color='yellow')

dataset = join.types(dataset, 'Nonsense point', 'Nonsense Ins/Del')
as.types(dataset)

dataset = delete.gene(aCML, gene = 'TET2')
dataset = delete.event(dataset, gene = 'ASXL1', type = 'Ins/Del')
dataset = delete.samples(dataset, samples = c('patient 5', 'patient 6'))
dataset = delete.type(dataset, type = 'Missense point')
view(dataset)

dataset = samples.selection(aCML, samples = as.samples(aCML)[1:3])
view(dataset)

dataset = events.selection(aCML,  filter.freq = .05, 
                           filter.in.names = c('EZH1','EZH2'), 
                           filter.out.names = 'SETBP1')

as.events(dataset)

## ----fig.width=7, fig.height=5.5, fig.cap="Multiple output from oncoprint can be captured as a gtable and composed via grid.arrange (package gridExtra). In this case we show  aCML data on top -- displayed after the as.alterations transformation -- versus a selected subdataset of events with a minimum frequency of 5%, force exclusion of SETBP1 (all events associated), and inclusion of EZH1 and EZH2.", results='hide'----
library(gridExtra)
grid.arrange(
  oncoprint(as.alterations(aCML, new.color = 'brown3'), 
            cellheight = 6, cellwidth = 4, gtable = TRUE,
            silent = TRUE, font.row = 6)$gtable,
  oncoprint(dataset, cellheight = 6, cellwidth = 4,
            gtable = TRUE, silent = TRUE, font.row = 6)$gtable, 
  ncol = 1)

## -----------------------------------------------------------------------------
view(aCML)

## -----------------------------------------------------------------------------
as.genotypes(aCML)[1:10,5:10]

## -----------------------------------------------------------------------------
as.events(aCML)[1:5, ]
as.events.in.sample(aCML, sample = 'patient 2')

## -----------------------------------------------------------------------------
as.genes(aCML)[1:8]

## -----------------------------------------------------------------------------
as.types(aCML)
as.colors(aCML)

## -----------------------------------------------------------------------------
head(as.gene(aCML, genes='SETBP1'))

## -----------------------------------------------------------------------------
as.samples(aCML)[1:10]

## -----------------------------------------------------------------------------
which.samples(aCML, gene='TET2', type='Nonsense point')

## -----------------------------------------------------------------------------
dataset = as.alterations(aCML)

## -----------------------------------------------------------------------------
view(dataset)

## -----------------------------------------------------------------------------
ngenes(aCML)
nevents(aCML)
nsamples(aCML)
ntypes(aCML)
npatterns(aCML)

## ----fig.width=6, fig.height=5, fig.cap="This plot gives a graphical visualization of the events that are in the dataset -- with a color per event type. It sorts samples to enhance exclusivity patterns among the events"----
oncoprint(aCML)

## ----fig.width=5, fig.height=5, fig.cap="This plot gives a graphical visualization of the events that are in the dataset -- with a color per event type. It it clusters samples/events"----
oncoprint(aCML, 
          legend = FALSE, 
          samples.cluster = TRUE, 
          gene.annot = list(one = list('NRAS', 'SETBP1'), two = list('EZH2', 'TET2')),
          gene.annot.color = 'Set2',
          genes.cluster = TRUE)

## -----------------------------------------------------------------------------
stages = c(rep('stage 1', 32), rep('stage 2', 32))
stages = as.matrix(stages)
rownames(stages) = as.samples(aCML)
dataset = annotate.stages(aCML, stages = stages)
has.stages(aCML)
head(as.stages(dataset))

## -----------------------------------------------------------------------------
head(as.stages(dataset))

## ----fig.width=6, fig.height=5------------------------------------------------
oncoprint(dataset, legend = FALSE)

## ----fig.width=6, fig.height=5, fig.cap="Example \texttt{oncoprint} output for aCML data with randomly annotated stages, in left, and samples clustered by group assignment in right -- for simplicity the group variable is again the stage annotation."----
oncoprint(dataset, group.samples = as.stages(dataset))

## -----------------------------------------------------------------------------
pathway = as.pathway(aCML,
                     pathway.genes = c('SETBP1', 'EZH2', 'WT1'),
                     pathway.name = 'MyPATHWAY',
                     pathway.color = 'red',
                     aggregate.pathway = FALSE)

## ----onco-pathway, fig.width=6.5, fig.height=2, fig.cap="Oncoprint output of a custom pathway called MyPATHWAY involving genes SETBP1, EZH2 and WT1; the genotype of each event is shown."----
oncoprint(pathway, title = 'Custom pathway',  font.row = 8, cellheight = 15, cellwidth = 4)

## ----fig.width=6.5, fig.height=1.8, fig.cap="Oncoprint output of a custom pair of pathways, with events shown"----
pathway.visualization(aCML, 
                      pathways=list(P1 = c('TET2', 'IRAK4'),  P2=c('SETBP1', 'KIT')),        
                      aggregate.pathways=FALSE,
                      font.row = 8)

## ----fig.width=6.5, fig.height=1, fig.cap="Oncoprint output of a custom pair of pathways, with events hidden"----
pathway.visualization(aCML, 
                      pathways=list(P1 = c('TET2', 'IRAK4'),  P2=c('SETBP1', 'KIT')),
                      aggregate.pathways = TRUE,
                      font.row = 8)

## ----eval=FALSE---------------------------------------------------------------
# library(rWikiPathways)
# # quotes inside query to require both terms
# my.pathways <- findPathwaysByText('SETBP1 EZH2 TET2 IRAK4 SETBP1 KIT')
# human.filter <- lapply(my.pathways, function(x) x$species == "Homo sapiens")
# my.hs.pathways <- my.pathways[unlist(human.filter)]
# # collect pathways idenifiers
# my.wpids <- sapply(my.hs.pathways, function(x) x$id)
# 
# pw.title<-my.hs.pathways[[1]]$name
# pw.genes<-getXrefList(my.wpids[1],"H")

## ----wikipathways, eval=FALSE-------------------------------------------------
# browseURL(getPathwayInfo(my.wpids[1])[2])
# browseURL(getPathwayInfo(my.wpids[2])[2])
# browseURL(getPathwayInfo(my.wpids[3])[2])

## ----include=FALSE------------------------------------------------------------
library(knitr)
opts_chunk$set(
  concordance = TRUE,
  background = "#f3f3ff"
)


## ----eval=FALSE---------------------------------------------------------------
# aCML = annotate.description(aCML, 'aCML data (Bioinf.)')

## -----------------------------------------------------------------------------
head(crc_maf[, 1:10])

## -----------------------------------------------------------------------------
dataset_maf = import.MAF(crc_maf)

## -----------------------------------------------------------------------------
dataset_maf = import.MAF(crc_maf, merge.mutation.types = FALSE)

## -----------------------------------------------------------------------------
dataset_maf = import.MAF(crc_maf, filter.fun = function(x){ x['Hugo_Symbol'] == 'APC'} )

## -----------------------------------------------------------------------------
dataset_maf = import.MAF(crc_maf, 
                         merge.mutation.types = FALSE, 
                         paste.to.Hugo_Symbol = c('MA.protein.change'))

## -----------------------------------------------------------------------------
crc_gistic

## -----------------------------------------------------------------------------
dataset_gistic = import.GISTIC(crc_gistic)

## -----------------------------------------------------------------------------
crc_plain

## -----------------------------------------------------------------------------
dataset_plain = import.genotypes(crc_plain, event.type='myVariant')

## ----results='hide', eval=FALSE-----------------------------------------------
# data = cbio.query(
#     genes=c('TP53', 'KRAS', 'PIK3CA'),
#     cbio.study = 'luad_tcga_pub',
#     cbio.dataset = 'luad_tcga_pub_cnaseq',
#     cbio.profile = 'luad_tcga_pub_mutations')


## -----------------------------------------------------------------------------
dataset = change.color(aCML, 'Ins/Del', 'dodgerblue4')
dataset = change.color(dataset, 'Missense point', '#7FC97F')
as.colors(dataset)

## -----------------------------------------------------------------------------
consolidate.data(dataset)

## -----------------------------------------------------------------------------
alterations = events.selection(as.alterations(aCML), filter.freq = .05)

## -----------------------------------------------------------------------------
gene.hypotheses = c('KRAS', 'NRAS', 'IDH1', 'IDH2', 'TET2', 'SF3B1', 'ASXL1')
aCML.clean = events.selection(aCML,
                              filter.in.names=c(as.genes(alterations), gene.hypotheses))
aCML.clean = annotate.description(aCML.clean, 
                                  'CAPRI - Bionformatics aCML data (selected events)')

## ----fig.width=8, fig.height=5.5, fig.cap="Data selected for aCML reconstruction annotated with the events which are part of a pattern that we will input to CAPRI."----
oncoprint(aCML.clean, gene.annot = list(priors = gene.hypotheses), sample.id = TRUE)

## -----------------------------------------------------------------------------
aCML.hypo = hypothesis.add(aCML.clean, 'NRAS xor KRAS', XOR('NRAS', 'KRAS'))

## ----eval=FALSE---------------------------------------------------------------
# aCML.hypo = hypothesis.add(aCML.hypo, 'NRAS or KRAS',  OR('NRAS', 'KRAS'))

## ----fig.width=6, fig.height=1, fig.cap="Oncoprint output to show the perfect (hard) exclusivity among NRAS/KRAS mutations in aCML"----
oncoprint(events.selection(aCML.hypo,
                           filter.in.names = c('KRAS', 'NRAS')),
          font.row = 8,
          ann.hits = FALSE)

## -----------------------------------------------------------------------------
aCML.hypo = hypothesis.add(aCML.hypo, 'SF3B1 xor ASXL1', XOR('SF3B1', XOR('ASXL1')),
                           '*')

## -----------------------------------------------------------------------------
as.events(aCML.hypo, genes = 'TET2') 
aCML.hypo = hypothesis.add(aCML.hypo,
                           'TET2 xor IDH2',
                           XOR('TET2', 'IDH2'),
                           '*')
aCML.hypo = hypothesis.add(aCML.hypo,
                           'TET2 or IDH2',
                           OR('TET2', 'IDH2'),
                           '*')

## ----fig.width=7, fig.height=2, fig.cap="{oncoprint} output to show the soft exclusivity among NRAS/KRAS mutations in aCML"----
oncoprint(events.selection(aCML.hypo,
                           filter.in.names = c('TET2', 'IDH2')),
          font.row = 8,
          ann.hits = FALSE)

## -----------------------------------------------------------------------------
aCML.hypo = hypothesis.add.homologous(aCML.hypo)

## -----------------------------------------------------------------------------
dataset = hypothesis.add.group(aCML.clean, OR, group = c('SETBP1', 'ASXL1', 'CBL'))

## ----fig.width=8, fig.height=6.5, fig.cap="oncoprint} output of the a dataset that has patterns that could be given as input to CAPRI to retrieve a progression model."----
oncoprint(aCML.hypo, gene.annot = list(priors = gene.hypotheses), sample.id = TRUE, 
          font.row=10, font.column=5, cellheight=15, cellwidth=4)

## -----------------------------------------------------------------------------
npatterns(dataset)
nhypotheses(dataset)

## -----------------------------------------------------------------------------
as.patterns(dataset)
as.events.in.patterns(dataset)
as.genes.in.patterns(dataset)
as.types.in.patterns(dataset)

## -----------------------------------------------------------------------------
head(as.hypotheses(dataset))
dataset = delete.hypothesis(dataset, event = 'TET2')
dataset = delete.pattern(dataset, pattern = 'OR_ASXL1_CBL')

## ----pattern-plot,fig.show='hide', fig.width=4, fig.height=2.2, fig.cap="Barplot to show an hypothesis: here we test genes SETBP1 and ASXL1 versus Missense point mutations of  CSF3R, which suggests that  that pattern does not 'capture' all the samples with  CSF3R mutations."----
tronco.pattern.plot(aCML,
                    group = as.events(aCML, genes=c('SETBP1', 'ASXL1')),
                    to = c('CSF3R', 'Missense point'),
                    legend.cex=0.8,
                    label.cex=1.0)

## ----pattern-plot-circos, fig.width=6, fig.height=6, fig.cap="Circos to show an hypothesis: here we test genes SETBP1 and ASXL1 versus Missense point mutations of  CSF3R. The combination of this and the previous  plots should allow to understand which pattern we shall write in an attempt to capture a potential causality relation between the pattern and the event."----
tronco.pattern.plot(aCML,
                    group = as.events(aCML, genes=c('TET2', 'ASXL1')),
                    to = c('CSF3R', 'Missense point'),
                    legend = 1.0,
                    label.cex = 0.8,
                    mode='circos')

## -----------------------------------------------------------------------------
model.capri = tronco.capri(aCML.hypo, boot.seed = 12345, nboot = 5)
model.capri = annotate.description(model.capri, 'CAPRI - aCML')

## -----------------------------------------------------------------------------
model.caprese = tronco.caprese(aCML.clean)
model.caprese = annotate.description(model.caprese, 'CAPRESE - aCML')

## -----------------------------------------------------------------------------
model.edmonds = tronco.edmonds(aCML.clean, nboot = 5, boot.seed = 12345)
model.edmonds = annotate.description(model.edmonds, 'MST Edmonds - aCML')

## -----------------------------------------------------------------------------
model.gabow = tronco.gabow(aCML.clean, nboot = 5, boot.seed = 12345)
model.gabow = annotate.description(model.gabow, 'MST Gabow - aCML')

## -----------------------------------------------------------------------------
model.chowliu = tronco.chowliu(aCML.clean, nboot = 5, boot.seed = 12345)
model.chowliu = annotate.description(model.chowliu, 'MST Chow Liu - aCML')

## -----------------------------------------------------------------------------
model.prim = tronco.prim(aCML.clean, nboot = 5, boot.seed = 12345)
model.prim = annotate.description(model.prim, 'MST Prim - aCML data')

gene.hypotheses = c('KRAS', 'NRAS', 'IDH1', 'IDH2', 'TET2', 'SF3B1', 'ASXL1')
alterations = events.selection(as.alterations(aCML), filter.freq = .05)
aCML.clean = events.selection(aCML,
    filter.in.names=c(as.genes(alterations), gene.hypotheses))
aCML.clean = annotate.description(aCML.clean, 
    'CAPRI - Bionformatics aCML data (selected events)')
aCML.hypo = hypothesis.add(aCML.clean, 'NRAS xor KRAS', XOR('NRAS', 'KRAS'))
aCML.hypo = hypothesis.add(aCML.hypo, 'SF3B1 xor ASXL1', XOR('SF3B1', XOR('ASXL1')),
    '*')
as.events(aCML.hypo, genes = 'TET2') 
aCML.hypo = hypothesis.add(aCML.hypo,
    'TET2 xor IDH2',
    XOR('TET2', 'IDH2'),
    '*')
aCML.hypo = hypothesis.add(aCML.hypo,
    'TET2 or IDH2',
    OR('TET2', 'IDH2'),
    '*')
aCML.hypo = hypothesis.add.homologous(aCML.hypo)
aCML.hypo = annotate.description(aCML.hypo, '')
aCML.clean = annotate.description(aCML.clean, '')
model.capri = tronco.capri(aCML.hypo, boot.seed = 12345, nboot = 5)
model.capri = annotate.description(model.capri, 'CAPRI - aCML')
model.caprese = tronco.caprese(aCML.clean)
model.caprese = annotate.description(model.caprese, 'CAPRESE - aCML')
model.edmonds = tronco.edmonds(aCML.clean, nboot = 5, boot.seed = 12345)
model.edmonds = annotate.description(model.edmonds, 'MST Edmonds - aCML')
model.gabow = tronco.gabow(aCML.clean, nboot = 5, boot.seed = 12345)
model.gabow = annotate.description(model.gabow, 'MST Gabow - aCML')
model.chowliu = tronco.chowliu(aCML.clean, nboot = 5, boot.seed = 12345)
model.chowliu = annotate.description(model.chowliu, 'MST Chow Liu - aCML')
model.prim = tronco.prim(aCML.clean, nboot = 5, boot.seed = 12345)
model.prim = annotate.description(model.prim, 'MST Prim - aCML data')

## -----------------------------------------------------------------------------
view(model.capri)

## ----fig.width=4,fig.height=4,warning=FALSE-----------------------------------
tronco.plot(model.capri, 
    fontsize = 12, 
    scale.nodes = 0.6, 
    confidence = c('tp', 'pr', 'hg'), 
    height.logic = 0.25, 
    legend.cex = 0.35, 
    pathways = list(priors = gene.hypotheses), 
    label.edge.size = 10)

## ----fig.width=7,fig.height=7,warning=FALSE, fig.cap="aCML data processed model by algorithms to extract models from individual patients, we show the otput of  CAPRESE,  and all algorithms based on Minimum Spanning Trees (Edmonds, Chow Liu and Prim). Only the  model retrieved by Chow Liu has two different edge colors as it was regularized with two different strategies: AIC and BIC."----
par(mfrow = c(2,2))
tronco.plot(model.caprese, fontsize = 22, scale.nodes = 0.6, legend = FALSE)
tronco.plot(model.edmonds, fontsize = 22, scale.nodes = 0.6, legend = FALSE)
tronco.plot(model.chowliu, fontsize = 22, scale.nodes = 0.6, legend.cex = .7)
tronco.plot(model.prim, fontsize = 22, scale.nodes = 0.6, legend = FALSE)

## -----------------------------------------------------------------------------
as.data.frame(as.parameters(model.capri))
has.model(model.capri)
dataset = delete.model(model.capri)

## -----------------------------------------------------------------------------
str(as.adj.matrix(model.capri))

## -----------------------------------------------------------------------------
marginal.prob = as.marginal.probs(model.capri)
head(marginal.prob$capri_bic)

## -----------------------------------------------------------------------------
joint.prob = as.joint.probs(model.capri, models='capri_bic')
joint.prob$capri_bic[1:3, 1:3]

## -----------------------------------------------------------------------------
conditional.prob = as.conditional.probs(model.capri, models='capri_bic')
head(conditional.prob$capri_bic)

## -----------------------------------------------------------------------------
str(as.confidence(model.capri, conf = c('tp', 'pr', 'hg')))

## ----selective-advantage------------------------------------------------------
as.selective.advantage.relations(model.capri)

## -----------------------------------------------------------------------------
model.boot = tronco.bootstrap(model.capri, nboot = 3, cores.ratio = 0)
model.boot = tronco.bootstrap(model.boot, nboot = 3, cores.ratio = 0, type = 'statistical')

## ----fig.width=4, fig.height=4, warning=FALSE, fig.cap="aCML model reconstructed by CAPRI with AIC / BIC as regolarizators and annotated with both non-parametric and statistical bootstrap scores. Edge thickness is proportional to the non-parametric scores."----

tronco.plot(model.boot, 
    fontsize = 12, 
    scale.nodes = .6,   
    confidence=c('sb', 'npb'), 
    height.logic = 0.25, 
    legend.cex = .35, 
    pathways = list(priors= gene.hypotheses), 
    label.edge.size=10)

## -----------------------------------------------------------------------------
as.bootstrap.scores(model.boot)
view(model.boot)

## ----fig.width=7, fig.height=7, fig.cap="Heatmap of the bootstrap scores for the CAPRI aCML model (via AIC regularization)."----

pheatmap(keysToNames(model.boot, as.confidence(model.boot, conf = 'sb')$sb$capri_aic) * 100, 
           main =  'Statistical bootstrap scores for AIC model',
           fontsize_row = 6,
           fontsize_col = 6,
           display_numbers = TRUE,
           number_format = "%d"
           )

## -----------------------------------------------------------------------------
model.boot = tronco.kfold.eloss(model.boot)
model.boot = tronco.kfold.prederr(model.boot, runs = 2, cores.ratio = 0)
model.boot = tronco.kfold.posterr(model.boot, runs = 2, cores.ratio = 0)

## -----------------------------------------------------------------------------
as.kfold.eloss(model.boot)
as.kfold.prederr(model.boot)
as.kfold.posterr(model.boot)

## -----------------------------------------------------------------------------
tabular = function(obj, M){
    tab = Reduce(
        function(...) merge(..., all = TRUE), 
            list(as.selective.advantage.relations(obj, models = M),
                as.bootstrap.scores(obj, models = M),
                as.kfold.prederr(obj, models = M),
                as.kfold.posterr(obj,models = M)))
  
    # merge reverses first with second column
    tab = tab[, c(2,1,3:ncol(tab))]
    tab = tab[order(tab[, paste(M, '.NONPAR.BOOT', sep='')], na.last = TRUE, decreasing = TRUE), ]
    return(tab)
}

head(tabular(model.boot, 'capri_bic'))

## ----fig.width=4,fig.height=4, warning=FALSE, fig.cap="aCML model reconstructed by CAPRI with  AIC/BIC as regolarizators and annotated with  non-parametric, as well as with entropy loss, prediction and posterior classification errors computed via cross-validation. Edge thickness is proportional to the non-parametric  scores."----
tronco.plot(model.boot, 
    fontsize = 12, 
    scale.nodes = .6, 
    confidence=c('npb', 'eloss', 'prederr', 'posterr'), 
    height.logic = 0.25, 
    legend.cex = .35, 
    pathways = list(priors= gene.hypotheses), 
    label.edge.size=10)

#Post reconstruction

gene.hypotheses = c('KRAS', 'NRAS', 'IDH1', 'IDH2', 'TET2', 'SF3B1', 'ASXL1')
alterations = events.selection(as.alterations(aCML), filter.freq = .05)
aCML.clean = events.selection(aCML,
                              filter.in.names=c(as.genes(alterations), gene.hypotheses))
aCML.clean = annotate.description(aCML.clean, 
                                  'CAPRI - Bionformatics aCML data (selected events)')
aCML.hypo = hypothesis.add(aCML.clean, 'NRAS xor KRAS', XOR('NRAS', 'KRAS'))
aCML.hypo = hypothesis.add(aCML.hypo, 'SF3B1 xor ASXL1', XOR('SF3B1', XOR('ASXL1')),
                           '*')
as.events(aCML.hypo, genes = 'TET2') 
aCML.hypo = hypothesis.add(aCML.hypo,
                           'TET2 xor IDH2',
                           XOR('TET2', 'IDH2'),
                           '*')
aCML.hypo = hypothesis.add(aCML.hypo,
                           'TET2 or IDH2',
                           OR('TET2', 'IDH2'),
                           '*')
aCML.hypo = hypothesis.add.homologous(aCML.hypo)
aCML.hypo = annotate.description(aCML.hypo, '')
aCML.clean = annotate.description(aCML.clean, '')
model.capri = tronco.capri(aCML.hypo, boot.seed = 12345, nboot = 5)
model.capri = annotate.description(model.capri, 'CAPRI - aCML')
model.caprese = tronco.caprese(aCML.clean)
model.caprese = annotate.description(model.caprese, 'CAPRESE - aCML')
model.edmonds = tronco.edmonds(aCML.clean, nboot = 5, boot.seed = 12345)
model.edmonds = annotate.description(model.edmonds, 'MST Edmonds - aCML')
model.gabow = tronco.gabow(aCML.clean, nboot = 5, boot.seed = 12345)
model.gabow = annotate.description(model.gabow, 'MST Gabow - aCML')
model.chowliu = tronco.chowliu(aCML.clean, nboot = 5, boot.seed = 12345)
model.chowliu = annotate.description(model.chowliu, 'MST Chow Liu - aCML')
model.prim = tronco.prim(aCML.clean, nboot = 5, boot.seed = 12345)
model.prim = annotate.description(model.prim, 'MST Prim - aCML data')

## -----------------------------------------------------------------------------
view(model.capri)

## ----fig.width=4,fig.height=4,warning=FALSE-----------------------------------
tronco.plot(model.capri, 
            fontsize = 12, 
            scale.nodes = 0.6, 
            confidence = c('tp', 'pr', 'hg'), 
            height.logic = 0.25, 
            legend.cex = 0.35, 
            pathways = list(priors = gene.hypotheses), 
            label.edge.size = 10)

## ----fig.width=7,fig.height=7,warning=FALSE, fig.cap="aCML data processed model by algorithms to extract models from individual patients, we show the otput of  CAPRESE,  and all algorithms based on Minimum Spanning Trees (Edmonds, Chow Liu and Prim). Only the  model retrieved by Chow Liu has two different edge colors as it was regularized with two different strategies: AIC and BIC."----
par(mfrow = c(2,2))
tronco.plot(model.caprese, fontsize = 22, scale.nodes = 0.6, legend = FALSE)
tronco.plot(model.edmonds, fontsize = 22, scale.nodes = 0.6, legend = FALSE)
tronco.plot(model.chowliu, fontsize = 22, scale.nodes = 0.6, legend.cex = .7)
tronco.plot(model.prim, fontsize = 22, scale.nodes = 0.6, legend = FALSE)

## -----------------------------------------------------------------------------
as.data.frame(as.parameters(model.capri))
has.model(model.capri)
dataset = delete.model(model.capri)

## -----------------------------------------------------------------------------
str(as.adj.matrix(model.capri))

## -----------------------------------------------------------------------------
marginal.prob = as.marginal.probs(model.capri)
head(marginal.prob$capri_bic)

## -----------------------------------------------------------------------------
joint.prob = as.joint.probs(model.capri, models='capri_bic')
joint.prob$capri_bic[1:3, 1:3]

## -----------------------------------------------------------------------------
conditional.prob = as.conditional.probs(model.capri, models='capri_bic')
head(conditional.prob$capri_bic)

## -----------------------------------------------------------------------------
str(as.confidence(model.capri, conf = c('tp', 'pr', 'hg')))

## ----selective-advantage------------------------------------------------------
as.selective.advantage.relations(model.capri)

## -----------------------------------------------------------------------------
model.boot = tronco.bootstrap(model.capri, nboot = 3, cores.ratio = 0)
model.boot = tronco.bootstrap(model.boot, nboot = 3, cores.ratio = 0, type = 'statistical')

## ----fig.width=4, fig.height=4, warning=FALSE, fig.cap="aCML model reconstructed by CAPRI with AIC / BIC as regolarizators and annotated with both non-parametric and statistical bootstrap scores. Edge thickness is proportional to the non-parametric scores."----

tronco.plot(model.boot, 
            fontsize = 12, 
            scale.nodes = .6,   
            confidence=c('sb', 'npb'), 
            height.logic = 0.25, 
            legend.cex = .35, 
            pathways = list(priors= gene.hypotheses), 
            label.edge.size=10)

## -----------------------------------------------------------------------------
as.bootstrap.scores(model.boot)
view(model.boot)

## ----fig.width=7, fig.height=7, fig.cap="Heatmap of the bootstrap scores for the CAPRI aCML model (via AIC regularization)."----

pheatmap(keysToNames(model.boot, as.confidence(model.boot, conf = 'sb')$sb$capri_aic) * 100, 
         main =  'Statistical bootstrap scores for AIC model',
         fontsize_row = 6,
         fontsize_col = 6,
         display_numbers = TRUE,
         number_format = "%d"
)

## -----------------------------------------------------------------------------
model.boot = tronco.kfold.eloss(model.boot)
model.boot = tronco.kfold.prederr(model.boot, runs = 2, cores.ratio = 0)
model.boot = tronco.kfold.posterr(model.boot, runs = 2, cores.ratio = 0)

## -----------------------------------------------------------------------------
as.kfold.eloss(model.boot)
as.kfold.prederr(model.boot)
as.kfold.posterr(model.boot)

## -----------------------------------------------------------------------------
tabular = function(obj, M){
  tab = Reduce(
    function(...) merge(..., all = TRUE), 
    list(as.selective.advantage.relations(obj, models = M),
         as.bootstrap.scores(obj, models = M),
         as.kfold.prederr(obj, models = M),
         as.kfold.posterr(obj,models = M)))
  
  # merge reverses first with second column
  tab = tab[, c(2,1,3:ncol(tab))]
  tab = tab[order(tab[, paste(M, '.NONPAR.BOOT', sep='')], na.last = TRUE, decreasing = TRUE), ]
  return(tab)
}

head(tabular(model.boot, 'capri_bic'))

## ----fig.width=4,fig.height=4, warning=FALSE, fig.cap="aCML model reconstructed by CAPRI with  AIC/BIC as regolarizators and annotated with  non-parametric, as well as with entropy loss, prediction and posterior classification errors computed via cross-validation. Edge thickness is proportional to the non-parametric  scores."----
tronco.plot(model.boot, 
            fontsize = 12, 
            scale.nodes = .6, 
            confidence=c('npb', 'eloss', 'prederr', 'posterr'), 
            height.logic = 0.25, 
            legend.cex = .35, 
            pathways = list(priors= gene.hypotheses), 
            label.edge.size=10)
