#!/usr/bin/env Rscript
# Three complementary, descriptive internal splicing figures.
# Local junction usage is NOT event-specific PSI; chains are partial, observed read segments.
suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(ggplot2)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) {
    stop('Usage: Rscript plot_junction_structure.R <relative_usage.tsv> <annotated_comparison.tsv> <reference_union_exons.tsv> <sample_A_chains.tsv> <sample_B_chains.tsv> <output_dir> <prefix>')
}
usage_file <- args[1]; anno_file <- args[2]; exon_file <- args[3]
case_chain_file <- args[4]; control_chain_file <- args[5]
outdir <- args[6]; prefix <- args[7]
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

u <- read_tsv(usage_file, show_col_types=FALSE, progress=FALSE)
j <- read_tsv(anno_file, show_col_types=FALSE, progress=FALSE)
e <- read_tsv(exon_file, show_col_types=FALSE, progress=FALSE)
check <- function(d, columns, label) {
    missing <- setdiff(columns, names(d))
    if (length(missing)) stop(label, ': missing columns: ', paste(missing, collapse=', '))
    if (!nrow(d)) stop(label, ': no rows')
}
check(u, c('junction_id', 'case_sample', 'control_sample',
    'case_read_support','control_read_support','reference_annotation_status',
    'case_donor_usage','control_donor_usage','case_donor_status','control_donor_status',
    'case_acceptor_usage','control_acceptor_usage','case_acceptor_status','control_acceptor_status',
    'donor_pooled_observed_junctions','acceptor_pooled_observed_junctions'), 'Usage')
check(j, c('junction_id','gene_label','chrom','intron_start0','intron_end0',
    'case_sample','control_sample','case_read_support','control_read_support',
    'reference_annotation_status'), 'Annotation')
check(e, c('chrom','start0','end0','strand'), 'Reference exons')
if (anyDuplicated(u$junction_id) || anyDuplicated(j$junction_id)) stop('Duplicate junction IDs')
if (!setequal(u$junction_id, j$junction_id)) stop('Usage and annotation contain different junction IDs')
if (length(unique(j$case_sample)) != 1L || length(unique(j$control_sample)) != 1L ||
    length(unique(j$chrom)) != 1L || length(unique(j$gene_label)) != 1L)
    stop('Expected one comparison, one gene, and one chromosome')
case_id <- unique(j$case_sample); control_id <- unique(j$control_sample)
chrom <- unique(j$chrom); gene <- unique(j$gene_label)
if (case_id == control_id) stop('Comparison samples must be different')
if (!identical(unique(u$case_sample), case_id) || !identical(unique(u$control_sample), control_id))
    stop('Usage and annotation sample IDs disagree')
if (any(e$chrom != chrom)) stop('Reference exon chromosome does not match input gene')
if (any(j$intron_start0 >= j$intron_end0) ||
    any(j$case_read_support < 0 | j$control_read_support < 0)) stop('Invalid junction coordinates or counts')

cols <- c('Ensembl-annotated'='#28537D', 'Absent from Ensembl 116'='#BE7C35')
annot_group <- function(x) ifelse(x == 'annotated_junction',
                                   'Ensembl-annotated', 'Absent from Ensembl 116')
common_theme <- theme_classic(base_size=12) + theme(
    plot.title=element_text(face='bold', size=15, margin=margin(b=12)),
    legend.position='bottom', legend.title=element_blank(),
    strip.background=element_blank(), strip.text=element_text(face='bold'),
    plot.margin=margin(15,20,15,15))
save_plot <- function(p, suffix, width=11, height=7) {
    stem <- file.path(outdir, paste0(prefix, '.', suffix))
    ggsave(paste0(stem, '.png'), p, width=width, height=height, dpi=300, bg='white')
    ggsave(paste0(stem, '.pdf'), p, width=width, height=height, device=cairo_pdf)
    cat('Wrote:', stem, '.png/.pdf\n')
}

# FIGURE 1: relative support for comparable donor/acceptor choices.
usage_rows <- list()
for (side in c('donor', 'acceptor')) {
    case_status <- u[[paste0('case_',side,'_status')]]
    control_status <- u[[paste0('control_',side,'_status')]]
    keep <- case_status == 'descriptive_comparable' & control_status == 'descriptive_comparable'
    keep[is.na(keep)] <- FALSE
    if (!any(keep)) next
    usage_rows[[side]] <- data.frame(
        site_type=if (side=='donor') 'Donor-site choices' else 'Acceptor-site choices',
        sample_A_usage=as.numeric(u[[paste0('case_',side,'_usage')]][keep]),
        sample_B_usage=as.numeric(u[[paste0('control_',side,'_usage')]][keep]),
        total_support=u$case_read_support[keep] + u$control_read_support[keep],
        annotation=annot_group(u$reference_annotation_status[keep]))
}
if (length(usage_rows)) {
    pdat <- bind_rows(usage_rows)
    pdat <- pdat[is.finite(pdat$sample_A_usage) & is.finite(pdat$sample_B_usage),]
    if (nrow(pdat)) {
        p_usage <- ggplot(pdat, aes(x=sample_A_usage,y=sample_B_usage)) +
            geom_abline(slope=1,intercept=0,linetype='dashed',colour='#9BA4AE',linewidth=.5) +
            geom_point(aes(colour=annotation,size=total_support),alpha=.7) +
            facet_wrap(~site_type,nrow=1) +
            scale_colour_manual(values=cols) +
            scale_size_area(max_size=5.5) +
            scale_x_continuous(limits=c(0,1),breaks=seq(0,1,.25),labels=function(v) paste0(v*100,'%')) +
            scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.25),labels=function(v) paste0(v*100,'%')) +
            coord_fixed(ratio=1) +
            labs(title='Relative internal splice-junction support',
                 x=paste0(case_id,' observed junction-share'),
                 y=paste0(control_id,' observed junction-share'),
                 size='Pooled supporting alignments') + common_theme
        save_plot(p_usage, 'relative_usage', width=12, height=6.7)
        cat('Junction points with comparable usage:',nrow(pdat),'(donor/acceptor points may repeat junctions)\n')
    }
} else {
    warning('No sites passed the comparable-site support filter; usage figure omitted')
}

# FIGURE 2: a local Sashimi-style junction-arc diagram, NOT a read-coverage Sashimi plot.
# Include arcs only where the individual sample supports the junction with >= N alignments.
# Filtering changes the visualisation, NOT the junction counts or annotated output tables.
min_arc_support_text <- Sys.getenv('JUNCTION_ARC_MIN_SUPPORT', unset='3')
if (!grepl('^[0-9]+$', min_arc_support_text))
    stop('JUNCTION_ARC_MIN_SUPPORT must be a positive integer')
min_arc_support <- suppressWarnings(as.integer(min_arc_support_text))
if (is.na(min_arc_support) || min_arc_support < 1L)
    stop('JUNCTION_ARC_MIN_SUPPORT must be a positive integer')
strong_arc_support <- 5L
# Prefer alternative splice-site groups and avoid choosing an extremely long intron.
j <- j %>% mutate(total_support=case_read_support + control_read_support,
                  annotation=annot_group(reference_annotation_status))
idx <- match(j$junction_id, u$junction_id)
alternative <- as.numeric(u$donor_pooled_observed_junctions[idx]) >= 2 |
               as.numeric(u$acceptor_pooled_observed_junctions[idx]) >= 2
span <- j$intron_end0 - j$intron_start0
eligible <- which(j$total_support > 0 & alternative & span <= 150000)
if (!length(eligible)) eligible <- which(j$total_support > 0 & span <= 150000)
if (!length(eligible)) eligible <- which(j$total_support > 0)
if (!length(eligible)) stop('No read-supported junctions for the arc plot')
seed <- j[eligible[which.max(j$total_support[eligible])],]
pad <- max(2000, min(12000, 0.1*(seed$intron_end0-seed$intron_start0)))
win_left <- max(0, floor(seed$intron_start0-pad))
win_right <- ceiling(seed$intron_end0+pad)
region <- j %>% filter(intron_start0>=win_left,intron_end0<=win_right,
                          pmax(case_read_support,control_read_support)>=min_arc_support) %>%
    arrange(desc(total_support),intron_start0) %>% slice_head(n=12)
if (!nrow(region)) stop('Cannot select a regional junction')
arc_parts <- list(); label_parts <- list()
for (k in seq_len(nrow(region))) {
    r <- region[k,]
    for (sample in c(case_id,control_id)) {
        n <- if (sample==case_id) r$case_read_support else r$control_read_support
        if (n < min_arc_support) next
        t <- seq(0,1,length.out=65)
        h <- .28 + .13*log1p(n)
        arc_parts[[length(arc_parts)+1L]] <- data.frame(
            sample=sample, junction_id=r$junction_id,
            x=r$intron_start0+t*(r$intron_end0-r$intron_start0),
            y=h*sin(pi*t),annotation=r$annotation,n=n,
            support_tier=ifelse(n>=strong_arc_support,'5+ reads','3-4 reads'))
        label_parts[[length(label_parts)+1L]] <- data.frame(
            sample=sample,x=(r$intron_start0+r$intron_end0)/2,
            y=h+.06,label=as.character(n))
    }
}
arc_dat <- bind_rows(arc_parts); label_dat <- bind_rows(label_parts)
arc_dat$sample <- factor(arc_dat$sample,levels=c(case_id,control_id))
label_dat$sample <- factor(label_dat$sample,levels=c(case_id,control_id))
ex_region <- e %>% filter(end0>win_left,start0<win_right)
p_arc <- ggplot() +
    geom_hline(yintercept=0,colour='#B9C3CA',linewidth=.5) +
    geom_rect(data=ex_region,
              aes(xmin=pmax(start0,win_left),xmax=pmin(end0,win_right),
                  ymin=-.12,ymax=-.025),inherit.aes=FALSE,fill='#626E79') +
    geom_path(data=arc_dat,
              aes(x=x,y=y,group=interaction(sample,junction_id),
                  colour=annotation,linewidth=support_tier),
              alpha=.90,lineend='round') +
    geom_text(data=label_dat,aes(x=x,y=y,label=label),size=3,check_overlap=TRUE) +
    facet_grid(sample ~ .) +
    scale_colour_manual(values=cols,drop=FALSE) +
    scale_linewidth_manual(values=c('3-4 reads'=.55,'5+ reads'=1.25),guide='none') +
    scale_x_continuous(labels=function(v) format(round(v/1e6,3),nsmall=3,trim=TRUE),
                       limits=c(win_left,win_right),expand=expansion(mult=c(.005,.005))) +
    labs(title=paste(gene,'internal splice-junction arcs'),
         x=paste0(chrom,' genomic coordinate (Mb)'),y=NULL) +
    common_theme +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),axis.ticks.y=element_blank())
save_plot(p_arc, 'junction_arcs', width=11.5, height=7.7)
cat('Arc region (0-based half-open):',chrom,win_left,win_right,
    '; eligible junctions:',nrow(region),
    '; visible sample-specific arcs:',length(label_parts),
    '; minimum support per sample:',min_arc_support,'\n')

# FIGURE 3: individual read-supported consecutive internal chains, not full isoforms.
load_chains <- function(path, sample) {
    z <- read_tsv(path, show_col_types=FALSE, progress=FALSE)
    check(z,c('sample_id','chrom','junction_chain_start0_end0','n_junctions','read_support'),path)
    z <- z[z$sample_id == sample & z$chrom == chrom & z$n_junctions >= 2 & z$read_support >= 1,]
    if (!nrow(z)) stop('No multi-junction chains in ',path,' for ',sample)
    z$starts <- vector('list',nrow(z)); z$ends <- vector('list',nrow(z))
    for (i in seq_len(nrow(z))) {
        pieces <- strsplit(z$junction_chain_start0_end0[i],',',fixed=TRUE)[[1]]
        pairs <- strsplit(pieces,'-',fixed=TRUE)
        if (any(lengths(pairs)!=2L)) stop('Malformed junction chain: ',z$junction_chain_start0_end0[i])
        coords <- as.integer(unlist(pairs))
        if (anyNA(coords)) stop('Non-integer junction coordinate')
        starts <- coords[seq(1,length(coords),2)]
        ends <- coords[seq(2,length(coords),2)]
        if (any(ends <= starts) || is.unsorted(starts) || length(starts)!=z$n_junctions[i])
            stop('Invalid junction chain coordinates')
        z$starts[[i]] <- starts; z$ends[[i]] <- ends
    }
    z$left <- vapply(z$starts,function(v) as.numeric(min(v)),numeric(1))
    z$right <- vapply(z$ends,function(v) as.numeric(max(v)),numeric(1))
    z
}
a <- load_chains(case_chain_file,case_id)
b <- load_chains(control_chain_file,control_id)
chains <- bind_rows(a,b) %>% mutate(span=right-left) %>% arrange(desc(read_support),span)
seed_chains <- chains %>% filter(span<=250000)
if (!nrow(seed_chains)) seed_chains <- chains
chain_seed <- seed_chains[1,]
chain_left <- max(0,chain_seed$left - 20000)
chain_right <- chain_seed$right + 20000
selected <- chains %>% filter(left >= chain_left,right <= chain_right) %>%
    group_by(sample_id) %>% arrange(desc(read_support),span,.by_group=TRUE) %>%
    slice_head(n=6) %>% ungroup()
if (!nrow(selected)) stop('No chains overlap selected structure region')
selected <- selected %>% arrange(sample_id,desc(read_support),span)
selected$track <- seq_len(nrow(selected))
selected$track_label <- paste0(selected$sample_id,' · chain ',
                              ave(selected$track,selected$sample_id,FUN=seq_along),
                              ' (n=',selected$read_support,')')
bridges <- list(); chain_arcs <- list()
for (i in seq_len(nrow(selected))) {
    starts <- selected$starts[[i]]; ends <- selected$ends[[i]]
    for (k in seq_along(starts)) {
        chain_arcs[[length(chain_arcs)+1L]] <- data.frame(
            x=starts[k],xend=ends[k],y=selected$track[i],sample=selected$sample_id[i])
        if (k < length(starts) && ends[k] < starts[k+1]) {
            bridges[[length(bridges)+1L]] <- data.frame(
                xmin=ends[k],xmax=starts[k+1],ymin=selected$track[i]-.08,
                ymax=selected$track[i]+.08)
        }
    }
}
chain_arc_dat <- bind_rows(chain_arcs)
bridge_dat <- if (length(bridges)) bind_rows(bridges) else data.frame(xmin=numeric(),xmax=numeric(),ymin=numeric(),ymax=numeric())
ref_chain <- e %>% filter(end0>chain_left,start0<chain_right)
chain_cols <- setNames(c('#28537D','#15958E'),c(case_id,control_id))
p_chain <- ggplot() +
    geom_rect(data=ref_chain,
              aes(xmin=pmax(start0,chain_left),xmax=pmin(end0,chain_right),
                  ymin=.15,ymax=.38),fill='#737E87') +
    geom_segment(data=selected,aes(x=left,xend=right,y=track,yend=track),
                 colour='#B6C1CA',linewidth=.55) +
    geom_rect(data=bridge_dat,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
              fill='#63717D',inherit.aes=FALSE) +
    geom_curve(data=chain_arc_dat,
               aes(x=x,xend=xend,y=y,yend=y,colour=sample),
               curvature=-.28,linewidth=.8,inherit.aes=FALSE) +
    scale_colour_manual(values=chain_cols) +
    scale_y_continuous(breaks=c(.25,selected$track),
                       labels=c('Reference exon union',selected$track_label),
                       limits=c(0,max(selected$track)+.9),expand=expansion(mult=c(0,.02))) +
    scale_x_continuous(labels=function(v) format(round(v/1e6,3),nsmall=3,trim=TRUE)) +
    coord_cartesian(xlim=c(chain_left,chain_right)) +
    labs(title=paste(gene,'read-supported internal junction chains'),
         x=paste0(chrom,' genomic coordinate (Mb)'),y=NULL) + common_theme +
    theme(legend.position='none',axis.text.y=element_text(size=9))
save_plot(p_chain, 'junction_chains', width=12.5, height=max(6,3.5+.48*nrow(selected)))
cat('Chain region (0-based half-open):',chrom,chain_left,chain_right,'; tracks:',nrow(selected),'\n')
cat('FIGURE VALIDATION COMPLETED\n')
