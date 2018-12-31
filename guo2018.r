source("pipe.r")

rl(scdb_dir="scrna_db_Li2018/", scfigs_dir="figs_Li2018/")

def_vars()

tgconfig::set_param("mcp_heatmap_alt_side_text", F,  "metacell")
tgconfig::set_param("mcell_mc2d_width", 1200,  "metacell")
tgconfig::set_param("mcell_mc2d_height", 1200,  "metacell")

tot_reads_per_cell = 2e4

# build all objects 
md = read.table("data/Guo2018_cell_clust_md.txt", header=T, sep="\t", stringsAsFactors = F)

tpm = fread("data/GSE99254_NSCLC.TCell.S12346.TPM.txt", sep="\t", header=T, stringsAsFactors=F, check.names=F, data.table=F)
tpm = tpm[!is.na(tpm$symbol), ]

full_mat = round(tpm[, c(-1, -2)] * (tot_reads_per_cell / 1e6))
rownames(full_mat) = tpm$symbol

types2set = list(N="Normal", P="Blood", "T"="Tumor")
md$batch_set_id = unlist(types2set[substr(md$sampleType, 1, 1)])

sample2cd25 = list(C="all", H="low", R="high", Y="int")

gate_type = substr(md$sampleType, 3, 3)
md$CD8_gate = ifelse(gate_type == "C", "+", "-")
md$CD4_gate = ifelse(gate_type == "C", "-", "+")
md$CD25_gate = unlist(sample2cd25[substr(md$sampleType, 3, 3)])

md[md$Patient == "P1118" & md$CD25_gate == "low" | md$Patient == "P1208" & md$sampleType == "TTH", "CD25_gate"] = "all"

md$amp_batch_id = paste(md$Patient, md$sampleType, sep="_")

pats_info = read.table("data/Guo2018_patients_info.txt", header=T, sep="\t", check.names=F)
fmd = merge(md, pats_info, by.x="Patient", by.y="Patient ID")
rownames(fmd) = fmd$UniqueCell_ID

full_mat = full_mat[, rownames(fmd)]

# load mat, create a metacell object by Guo 2018 original clusters
scdb_add_mat(guo_id, scm_new_matrix(Matrix(as.matrix(full_mat)), fmd))

mat_g_full = scdb_mat(guo_id)

mcell_mat_ignore_cells(guo_filt_id, guo_id, fmd[fmd$majorCluster == "filtered", "UniqueCell_ID"])
mat_g = scdb_mat(guo_filt_id)

mc_factor = as.factor(mat_g@cell_metadata[mat_g@cells, 'majorCluster'])
guo_mc = as.numeric(mc_factor)
names(guo_mc) = mat_g@cells

mcell_new_mc(guo_filt_id, guo_mc, character(0), mat_g)
mc_g = scdb_mc(guo_filt_id)

color_key = read.table("data/Guo2018_cl_colors.txt", header=T, sep="\t")
color_key$gene = ""
mc_g@color_key = color_key
g_col2group = get_mc_col2group(mc_g)
g_group2col = get_mc_group2col(mc_g)
mc_g@colors = as.vector(g_group2col[levels(mc_factor)])
scdb_add_mc(guo_filt_id, mc_g)

# Create lateral gene sets to blacklist 
#
# Wrapper for creating cell cycle, stress and IFN gene sets
create_lateral_gset = function(gset_nm = "guo2018_lateral", suff="filtered", cor_thresh=0.1, sz_cor_thresh=0.2, cl_to_keep=NULL) 
{
	genes_anchors = c('MKI67', 'PCNA', 'TOP2A', 'SMC4', 'MCM3', 'TXN', 'HSP90AB1', 'FOS', 'JUN', 'ISG15', 'OAS1', 'WARS', 'IFIT1')
	tab_fn = paste(.scfigs_base, "guo2018_lateral_gmods.txt", sep="/")
	
	
	if (is.null(cl_to_keep)) {
		mcell_mat_rpt_cor_anchors(mat_id=guo_filt_id, gene_anchors = genes_anchors, cor_thresh = cor_thresh, gene_anti = c(), tab_fn = tab_fn, sz_cor_thresh = sz_cor_thresh)
		gcor_mat = read.table(tab_fn, header=T)
		foc_genes = apply(gcor_mat[, genes_anchors], 1, which.max)
		
		gset = gset_new_gset(sets = foc_genes, desc = "Guo2018 Cell cycle, IFN1 and stress correlated genes")
		
		scdb_add_gset(gset_nm, gset)
		
		sub_mat_id = paste(guo_filt_id, gset_nm, sep="_")
		mcell_mat_ignore_genes(new_mat_id = sub_mat_id, mat_id = guo_filt_id, ig_genes = names(foc_genes), reverse = T)
		
		mcell_gset_split_by_dsmat(gset_id = gset_nm, mat_id = sub_mat_id, K = 24)
		mcell_plot_gset_cor_mats(gset_id = gset_nm, scmat_id = sub_mat_id)
	} else {
		mcell_gset_remove_clusts(gset_id = gset_nm, filt_clusts = cl_to_keep, new_id = paste(gset_nm, suff, sep="_"), reverse=T)
	}
	
}

# Creating here a lateral gene set, yet will later use the one used in the paper for consistency (differences are because of different initial random seeding)
create_lateral_gset(gset_nm = "guo2018_lateral_new")

# manually selected the clusters of genes to keep (here from the gene set used in the paper)
create_lateral_gset(cl_to_keep = c(2, 5, 6, 7, 10, 11, 12, 13, 15, 24))

# create metacells from their cells	(creating here an object but will use the one used in the paper)
mel_basic_mat2mc(guo_filt_id, c('Tumor', 'Normal', 'Blood'), lateral_gset_id=guo_lat_gset_id, force_new=T, name="_new")

guo_conf = colorize_by_confusion_mat(mc_id=guo_our_mc_clean_id, graph_id=guo_our_mc_id)
guo_conf = colorize_by_confusion_mat(mc_id=guo_our_mc_clean_id, graph_id=guo_our_mc_id, res=guo_conf, supmc_file="config/guo2018_supmc.txt", marks_file="config/guo2018_marks.txt")

# manually assign the 2 CD4-GZMK mcs
mc_guo_ours = scdb_mc(guo_our_mc_clean_id)
cd4_gzmk_col = "#448ec1"	
mc_guo_ours@colors[c(20, 64)] = cd4_gzmk_col
mc_guo_ours@color_key = rbind(mc_guo_ours@color_key, data.frame(gene="", group="CD4-GZMK", color=cd4_gzmk_col))
scdb_add_mc(guo_our_mc_clean_id, mc_guo_ours)

# basic plots
md_fields = c('Patient', 'CD8_gate', 'CD4_gate', 'CD25_gate', 'Sex', 'Stage', 'Histology', 'Smoker')
mc_guo_ours = scdb_mc(guo_our_mc_clean_id)

go_col2group = get_mc_col2group(mc_guo_ours)

guo_mc_ord = order(ord_by_id[[guo_our_mc_clean_id]][go_col2group[mc_guo_ours@colors]] + 1e-5 * seq_along(mc_guo_ours@colors))

mel_basic_mc_mc2d_plots(mc_id=guo_our_mc_clean_id, mat_id=guo_our_mc_id, graph_id=guo_our_mc_id, lateral_gset_id = guo_lat_gset_id, metadata_fields_to_export=md_fields, mc_ord=guo_mc_ord)

tgconfig::set_param("mcell_mc2d_width", 500,  "metacell")
tgconfig::set_param("mcell_mc2d_cex", 0.6,  "metacell")
for (mdf in c(md_fields, 'batch_set_id')) {
	mcell_mc2d_plot_by_factor(guo_our_mc_clean_id, guo_our_mc_id, mdf, T)
}

mat_guo = scdb_mat(guo_our_mc_id)
mc_guo = scdb_mc(guo_filt_id)
g_col2group = get_mc_col2group(mc_guo)
g_group2col = get_mc_group2col(mc_guo)

mc_guo_ours = scdb_mc(guo_our_mc_clean_id)
go_col2group = get_mc_col2group(mc_guo_ours)
go_group2col = get_mc_group2col(mc_guo_ours)

# Melanoma T/NK cells
mc_t_nk = scdb_mc(tumor_t_nk_id)
t_col2group = get_mc_col2group(mc_t_nk)
t_group2col = get_mc_group2col(mc_t_nk)

# guo mc cor to our melanoma T/NK mcs
lfp_gu = log2(mc_guo_ours@mc_fp)
lfp_t = log2(mc_t_nk@mc_fp)

f_cd8_mc = mc_t_nk@colors %in% t_group2col[cd8_nms]
f_treg_mc = mc_t_nk@colors %in% t_group2col['Treg']

gg = intersect(rownames(lfp_t), rownames(lfp_gu))

mc_t_nk_ann = data.frame(row.names=1:ncol(lfp_t), mel_group=t_col2group[mc_t_nk@colors])
mc_ann = data.frame(row.names=1:ncol(mc_guo_ours@mc_fp), mc_group=go_col2group[mc_guo_ours@colors])

t_nk_mc_ord = order(ord_by_id[[tumor_t_nk_id]][t_col2group[mc_t_nk@colors]] + 1e-5 * seq_along(mc_t_nk@colors))
guo_mc_ord = order(ord_by_id[[guo_our_mc_clean_id]][go_col2group[mc_guo_ours@colors]] + 1e-5 * seq_along(mc_guo_ours@colors))

# breakdown of mc groups to sample source
src_colors = list(Blood="#ffd966", Normal="#6dc983", Tumor="#ed7d31")
grp_src = table(go_col2group[mc_guo_ours@colors[mc_guo_ours@mc]], mat_guo@cell_metadata[names(mc_guo_ours@mc), 'batch_set_id'])
grp_src_n = grp_src / rowSums(grp_src)
grp_src_n = grp_src_n[names(ord_by_id[[guo_our_mc_clean_id]]), c('Blood', 'Normal', 'Tumor')]

.plot_start(scfigs_fn(guo_our_mc_clean_id, "mc_grp_by_source"), 700, 700)
layout(matrix(1:2, nrow=1), width=c(1,3))
par(mar=c(4,10,4,0))
barplot(rep(1, nrow(grp_src_n)), col=go_group2col[rownames(grp_src_n)], horiz=T, xaxt='n', yaxt='n', ylim=c(0, nrow(grp_src_n) * 1.2 + 0.5))
mtext(rownames(grp_src_n), 2, at=seq(0.7, by=1.2, length=nrow(grp_src_n)), las=2, line=0.5, cex=2)
par(mar=c(4,0.5,4, 2))
barplot(t(grp_src_n), horiz=T, col=unlist(src_colors[colnames(grp_src_n)]), yaxt='n', ylim=c(0, nrow(grp_src_n) * 1.2 + 0.5))
legend("topright", legend=colnames(grp_src_n), fill=unlist(src_colors[colnames(grp_src_n)]), ncol=ncol(grp_src_n), bty='n', cex=1.5)
dev.off()

ug = unique(g_col2group)
ugo = unique(go_col2group)

# clust - mc-group conf heatmap
clust_g_ann = data.frame(row.names=ug, cl_group=ug)
mc_g_ann = data.frame(row.names=ugo, mc_group=ugo)

grp_comp = table(go_col2group[mc_guo_ours@colors[mc_guo_ours@mc]], g_col2group[mc_guo@colors[mc_guo@mc[names(mc_guo_ours@mc)]]])
grp_comp = grp_comp[names(ord_by_id[[guo_our_mc_clean_id]]), ]
grp_comp_n = grp_comp / rowSums(grp_comp)

.plot_start(scfigs_fn(guo_our_mc_clean_id, "clust_mcGroup_cell_sharing"), max(800, 200 + 20*nrow(mc_g_ann)), max(800, 200 + 20*nrow(clust_g_ann)))
cord = order(apply(grp_comp_n, 2, which.max))
pheatmap(t(grp_comp_n[, cord]), cluster_rows=F, cluster_cols=F, cellwidth=20, cellheight=20, color=colorRampPalette(c('white', 'red'))(101), display_numbers=t(grp_comp[, cord]), annotation_row=clust_g_ann, annotation_col=mc_g_ann, annotation_color=list(cl_group=g_group2col[colnames(grp_comp_n)[cord]], mc_group=go_group2col[rownames(grp_comp_n)]))
dev.off()

# clust - mc conf heatmap
mc_comp = table(mc_guo_ours@mc, g_col2group[mc_guo@colors[mc_guo@mc[names(mc_guo_ours@mc)]]])
mc_comp = mc_comp[guo_mc_ord, ]
mc_comp_n = mc_comp / rowSums(mc_comp)

.plot_start(scfigs_fn(guo_our_mc_clean_id, "clust_mc_cell_sharing"), max(800, 500 + 5*nrow(mc_ann)), max(800, 200 + 20*nrow(clust_g_ann)))

cord = order(ifelse(grepl("other|DN|DP|diverse|iNKT", colnames(mc_comp_n), perl=T), 100, 1) * apply(mc_comp_n, 2, which.max))
pheatmap(t(mc_comp_n[, cord]), cluster_rows=F, cluster_cols=F, cellwidth=5, cellheight=20, color=colorRampPalette(c('white', 'red'))(101), annotation_row=clust_g_ann, annotation_col=mc_ann, annotation_color=list(cl_group=g_group2col[colnames(mc_comp_n)[cord]], mc_group=go_group2col[names(ord_by_id[[guo_our_mc_clean_id]])]))
dev.off()

# diff expression of cells in a specific Guo cluster: breakdown by those who are in mcs unique to this cluster and those that are in mcs that also contain substantial number of cells from another cluster
cluster_diff_expr_by_mc_groups = function(mc_comp, cl1, cl2, min_cells_in_mc=10, re_to_ignore="^MIR[0-9]|^SNOR[A-Z][0-9]", lat_gset_id=guo_lat_gset_id, min_abs_enr=1.5, max_genes=70) 
{
	cl_ind = mat_guo@cell_metadata[names(mc_guo_ours@mc), 'majorCluster'] == cl1
	
	mcs1 = names(which(mc_comp[, cl1] >= min_cells_in_mc))
	mcs2 = names(which(mc_comp[, cl2] >= min_cells_in_mc))
	
	cells1 = names(mc_guo_ours@mc)[mc_guo_ours@mc %in% setdiff(mcs1, mcs2) & cl_ind]
	cells2 = names(mc_guo_ours@mc)[mc_guo_ours@mc %in% intersect(mcs1, mcs2) & cl_ind]

	message(sprintf("comparing %d cells unique to %s (%d mcs) to %d cells (%d mcs) sharing mcs with %s, out of total %d %s cells", length(cells1), cl1, length(setdiff(mcs1, mcs2)), length(cells2), length(intersect(mcs1, mcs2)), cl2, sum(cl_ind), cl1))
	
	de = diff_expr(mc_guo_ours, mat_guo@mat, mcs1=NULL, mcs2=NULL, nms1=cells1, nms2=cells2, geo_mean=T, min_max_umi=50, reg=5)
	if (!is.null(re_to_ignore)) {
		de = de[grep(re_to_ignore, de$gene, perl=T, invert=T), ]
	}
	if (!is.null(lat_gset_id)) {
		lat_gset = scdb_gset(lat_gset_id)
		de = de[!(de$gene %in% names(lat_gset@gene_set)), ]
	}
	de_f = de %>% filter(abs(enr) >= min_abs_enr) 
	if (!is.null(de_f) & nrow(de_f) > 0) {
		message(sprintf("%d out of %d passed enr filter", nrow(de_f), nrow(de)))
	
		stopifnot(nrow(de_f) < max_genes)
		.plot_start(scfigs_fn(guo_our_mc_clean_id, sprintf("diff_expr_by_mcs_%s_unique_vs_%s_shared", cl1, cl2)), 300, max(300, 100 + max_genes * 12))
		par(mar=c(4,12,4,1))
		barplot(de_f$enr, names.arg=de_f$gene, las=2, horiz=T, col=ifelse(de_f$enr > 0, g_group2col[cl1], g_group2col[cl2]), main=sprintf("%s unique vs\n%s shared", cl1, cl2), ylim=c(0, max_genes))
		dev.off()
	}
}

clusts1 = c("CD8_C5-ZNF683", "CD4_C4-CD69" , "CD4_C4-CD69")
clusts2 = c("CD8_C6-LAYN"  , "CD4_C2-ANXA1", "CD4_C6-GZMA")
for (i in seq_along(clusts1)) {
	cluster_diff_expr_by_mc_groups(mc_comp, clusts1[i], clusts2[i])
}

# compare mean umi (geom-mean) per group - melanoma vs lung (norm by diff in naives)
mat_t_nk = scdb_mat(tumor_t_nk_id)
mat_t_nk_ds = scm_downsamp(mat_t_nk@mat, 500)

mel_lat_gset = scdb_gset(lateral_gset_id)
mel_lat_genes = names(mel_lat_gset@gene_set)
all_genes_ann = annotate_genes(rownames(lfp_t))
gf = apply(lfp_t[, f_cd8_mc], 1, max) > 1 & !grepl("TF", all_genes_ann$type) & (! (rownames(lfp_t) %in% mel_lat_genes))
dysf_genes =  tail(sort(cor(lfp_t['LAG3', f_cd8_mc], t(lfp_t[gf, f_cd8_mc]))[1,]), 30)
cyto_genes =  tail(sort(cor(lfp_t['FGFBP2', f_cd8_mc], t(lfp_t[gf, f_cd8_mc]))[1,]), 30)

treg_gf = apply(lfp_t[, f_treg_mc], 1, max) > 1 & !grepl("TF", all_genes_ann$type)& (! (rownames(lfp_t) %in% mel_lat_genes))
treg_genes =  tail(sort(cor(lfp_t["IL2RA", f_treg_mc], t(lfp_t[treg_gf, f_treg_mc]))[1,]), 30)

# mel vs corresponding lung mc grp gene comparison
compare_mel_to_lung_corresponding_group = function(mc_t_nk, mat_t_nk_ds, mel_grp='dysfunctional', lung_grp='CD8-dysf', mel_genes=names(dysf_genes), lung_ref_grp='CD8-GZMK', de_min_enr=2, de_min_max_umi_per_cell=0.005, text_cex=1, size=600, re_to_filt="^MIR[0-9]|^SNOR[A-Z][0-9]", umi_reg=64, min_xlim=3, mel_naive='naive', lung_naive='naive-CXCR6', n_bp_genes=450, gene_labels_to_show=NA, mel_reg=10, lung_reg=10, mat_guo_as_ds=NULL, points_col=NULL, name="")
{
	message(sprintf("compare %s (by %s) mel to %s (by %s) lung", mel_grp, mel_naive, lung_grp, lung_naive))
	gg = intersect(rownames(mc_t_nk@mc_fp), rownames(mc_guo_ours@mc_fp))
	if (!is.null(re_to_filt)) {
		gg = grep(re_to_filt, gg, invert=T, perl=T, v=T)
	}
		
	if (is.null(lung_ref_grp)) {
		lung_genes = NULL
	} else {
		lung_de = diff_expr(mc_guo_ours, mat_guo@mat, which(go_col2group[mc_guo_ours@colors] == lung_grp), which(mc_guo_ours@colors == go_group2col[lung_ref_grp]), filter_outlier_genes=T, geo_mean=T, min_max_umi=0.2, reg=0.01)
	
		lung_genes = intersect(gg, lung_de[lung_de$enr >= de_min_enr, 'gene'])
	}
	
	mel_stat = rowSums(mat_t_nk@mat[gg, mc_t_nk@colors[mc_t_nk@mc] == t_group2col[mel_grp]])
	lung_stat = rowSums(mat_guo@mat[gg, mc_guo_ours@colors[mc_guo_ours@mc] == go_group2col[lung_grp]])
	
	norm_to = min(sum(mel_stat), sum(lung_stat))
	mel_n = mel_stat * (norm_to / sum(mel_stat))
	lung_n = lung_stat * (norm_to / sum(lung_stat))
	
	x = log2((mel_n + lung_n)/2)
	y = log2( (mel_n + umi_reg) / (lung_n + umi_reg) )
	
	mel_genes = intersect(mel_genes, gg)
	lung_genes = intersect(lung_genes, gg)
	joined_genes = intersect(mel_genes, lung_genes)
	rest_genes = setdiff(gg, union(mel_genes, lung_genes))
	gene_cols = rep('grey', length(gg))
	names(gene_cols) = gg
	gene_cols[mel_genes] = 'darkblue'
	gene_cols[lung_genes] = 'darkred'
	gene_cols[joined_genes] = 'darkgreen'
	
	.plot_start(scfigs_fn(guo_our_mc_clean_id, sprintf("mel_%s_vs_lung_%s_genes_log_tot", mel_grp, lung_grp)), size, size)
	plot(x, y, pch=19, cex=0.5, col=gene_cols, main=sprintf("Mel %s vs Lung %s", mel_grp, lung_grp), xlab="mean mel + lung", ylab="log2(mel/lung)", xlim=c(min_xlim, max(x)))
	text(x[rest_genes], y[rest_genes], rest_genes, pos=2, cex=text_cex, col="grey")
	if (length(mel_genes) > 0) {
		text(x[mel_genes], y[mel_genes], mel_genes, pos=2, cex=text_cex, col='darkblue')
	}
	if (length(lung_genes) > 0) {
		text(x[lung_genes], y[lung_genes], lung_genes, pos=2, cex=text_cex, col='darkred')
	}
	if (length(joined_genes) > 0) {
		text(x[joined_genes], y[joined_genes], joined_genes, pos=2, cex=text_cex, col='darkgreen')
	}
	abline(h=0, lty=2)
	if (length(mel_genes) > 0 | length(lung_genes) > 0) {
		legend("topleft", legend=c(paste(mel_grp, "enr"), paste(lung_grp, "enr"), "common"), fill=c('darkblue', 'darkred', 'darkgreen'), bty='n')
	}
	dev.off()
	
	
	mel_naive_stat = rowSums(mat_t_nk@mat[gg, mc_t_nk@colors[mc_t_nk@mc] == t_group2col[mel_naive]])
	lung_naive_stat = rowSums(mat_guo@mat[gg, mc_guo_ours@colors[mc_guo_ours@mc] == go_group2col[lung_naive]])
	norm_to = min(sum(mel_naive_stat), sum(mel_stat), sum(lung_naive_stat), sum(lung_stat))
	mel_stat = mel_stat * (norm_to / sum(mel_stat))
	mel_naive_stat = mel_naive_stat * (norm_to / sum(mel_naive_stat))
	lung_stat = lung_stat * (norm_to / sum(lung_stat))
	lung_naive_stat = lung_naive_stat * (norm_to / sum(lung_naive_stat))
	reg = 50
	x1 = log2( (mel_stat + reg) / (mel_naive_stat + reg))
	y1 = log2((lung_stat + reg) / (lung_naive_stat + reg))
	
	mel_de = diff_expr(mc_t_nk, mat_t_nk_ds, which(mc_t_nk@colors == t_group2col[mel_grp]), which(mc_t_nk@colors == t_group2col[mel_naive]), geo_mean=T, min_max_umi=0, reg=mel_reg)
	mel_gs  = intersect(gg, mel_de[mel_de$enr >= de_min_enr & pmax(mel_de$tot1, mel_de$tot2) >= sum(mc_t_nk@colors[mc_t_nk@mc] == t_group2col[mel_grp]) * de_min_max_umi_per_cell, 'gene'])
	
	if (is.null(mat_guo_as_ds)) {
		mat_guo_as_ds = mat_guo@mat
	}
	lung_de = diff_expr(mc_guo_ours, mat_guo_as_ds, which(mc_guo_ours@colors == go_group2col[lung_grp]), which(mc_guo_ours@colors == go_group2col[lung_naive]), geo_mean=T, min_max_umi=0, reg=lung_reg)
	lung_gs  = intersect(gg, lung_de[lung_de$enr >= de_min_enr & pmax(lung_de$tot1, lung_de$tot2) >= sum(mc_guo_ours@colors[mc_guo_ours@mc] == go_group2col[lung_grp]) * de_min_max_umi_per_cell, 'gene'])
	joined_gs = intersect(mel_gs, lung_gs)
	all_gs = union(lung_gs, mel_gs)
	x1 = mel_de[gg, 'enr']
	names(x1) = gg
	y1 = lung_de[gg, 'enr']
	names(y1) = gg
	write.table(data.frame(row.names=all_gs, mel_enr=mel_de[all_gs, 'enr'], lung_enr=lung_de[all_gs, 'enr']), scfigs_fn(guo_our_mc_clean_id, sprintf("mel_%s_by_%s_vs_lung_%s_by_%s_genes", mel_grp, mel_naive, lung_grp, lung_naive), ext="txt"), quote=F, sep="\t")
	
	.plot_start(scfigs_fn(guo_our_mc_clean_id, sprintf("mel_%s_by_%s_vs_lung_%s_by_%s_genes%s", mel_grp, mel_naive, lung_grp, lung_naive, name)), size, size)
	if (is.null(points_col)) {
		gene_cols = ifelse(gg %in% joined_gs, 'darkgreen', ifelse(gg %in% lung_gs, 'darkred', ifelse(gg %in% mel_gs, 'darkblue', 'lightgrey')))
	} else {
		gene_cols = rep(points_col, length(gg))
	}
	names(gene_cols) = gg
	plot(x1, y1, pch=19, col=gene_cols, main="", xlab=paste(mel_grp, "by", mel_naive), ylab=paste(lung_grp, "by", lung_naive), xlim=c(min(x1), max(x1) + ifelse(is.null(gene_labels_to_show), 0, 0.5)))
	points(x1[mel_gs], y1[mel_gs], col=gene_cols[mel_gs], pch=19)
	points(x1[lung_gs], y1[lung_gs], col=gene_cols[lung_gs], pch=19)
	points(x1[joined_gs], y1[joined_gs], col=gene_cols[joined_gs], pch=19)
	
	if (!is.null(gene_labels_to_show)) {
		if (is.na(gene_labels_to_show)){
			gene_labels_to_show = all_gs
		}
		text(x1[gene_labels_to_show], y1[gene_labels_to_show], gene_labels_to_show, cex=0.9, pos=ifelse(x1[gene_labels_to_show] > y1[gene_labels_to_show], 4, 2))
	}
	grid(col='black')
	legend("topleft", legend=c(paste(mel_grp, "enr"), paste(lung_grp, "enr"), "common"), pch=19, col=c('darkblue', 'darkred', 'darkgreen'), bty='n')
	dev.off()
	
	x1b = (x1 + y1)/2
	y1b = y1 - x1
	x1 = x1b
	y1 = y1b
	.plot_start(scfigs_fn(guo_our_mc_clean_id, sprintf("mel_%s_by_%s_vs_lung_%s_by_%s_genes%s_v2", mel_grp, mel_naive, lung_grp, lung_naive, name)), size, size)
	gene_cols = ifelse(gg %in% joined_gs, 'darkgreen', ifelse(gg %in% lung_gs, 'darkred', ifelse(gg %in% mel_gs, 'darkblue', 'lightgrey')))
	names(gene_cols) = gg
	plot(x1, y1, pch=19, col=gene_cols, main="", xlab=paste(mel_grp, "by", mel_naive), ylab=paste(lung_grp, "by", lung_naive), xlim=c(0, max(x1) + ifelse(is.null(gene_labels_to_show), 0, 0.5)))
	points(x1[mel_gs], y1[mel_gs], col=gene_cols[mel_gs], pch=19)
	points(x1[lung_gs], y1[lung_gs], col=gene_cols[lung_gs], pch=19)
	points(x1[joined_gs], y1[joined_gs], col=gene_cols[joined_gs], pch=19)
	
	if (!is.null(gene_labels_to_show)) {
		if (is.na(gene_labels_to_show)){
			gene_labels_to_show = all_gs
		}
		text(x1[gene_labels_to_show], y1[gene_labels_to_show], gene_labels_to_show, cex=0.7, pos=4)
	}
	grid(col='black')
	legend("topleft", legend=c(paste(mel_grp, "enr"), paste(lung_grp, "enr"), "common"), pch=19, col=c('darkblue', 'darkred', 'darkgreen'), bty='n')
	dev.off()
	
	list(x=x, y=y, x1=x1, y1=y1, mel_gs=mel_gs, lung_gs=lung_gs, mel_de=mel_de, lung_de=lung_de)
}

ml_dysf = compare_mel_to_lung_corresponding_group(mc_t_nk, mat_t_nk_ds, umi_reg = 64)
ml_cyto = compare_mel_to_lung_corresponding_group(mc_t_nk, mat_t_nk_ds, mel_grp='effector2', lung_grp='CD8-cyto', mel_genes=names(cyto_genes), umi_reg = 64)
ml_dysf_by_cyto = compare_mel_to_lung_corresponding_group(mc_t_nk, mat_t_nk_ds, mel_grp = 'dysfunctional', lung_grp = 'CD8-dysf', mel_genes = names(dysf_genes), lung_ref_grp='CD8-cyto', mel_naive='effector2', lung_naive = 'CD8-cyto', umi_reg = 64)
ml_treg = compare_mel_to_lung_corresponding_group(mc_t_nk, mat_t_nk_ds, mel_grp='Treg', lung_grp='tumor-Treg', lung_ref_grp='blood-Treg', mel_genes=names(treg_genes), umi_reg = 64)

y = read.table(scfigs_fn(guo_our_mc_clean_id, "mel_dysfunctional_by_naive_vs_lung_CD8-dysf_by_naive-CXCR6_genes", ext="txt"), header=T)
y = y[order(y$lung_enr - y$mel_enr), ]
.plot_start(scfigs_fn(guo_our_mc_clean_id, "mel_dysfunctional_by_naive_vs_lung_CD8-dysf_by_naive-CXCR6_genes_mel_enr_barplot"), 300, 400)
par(mar=c(4,8,1,1))
barplot(t(y[1:20, ]), beside=T, horiz=T, las=2, col=c(t_group2col['dysfunctional'], 'grey'), ylim=c(0, 22 * 3))
legend("topleft", legend=c('Melanoma', 'Lung'), fill=c(t_group2col['dysfunctional'], 'grey'), bty='n', ncol=2)
dev.off()

dysf_genes_to_show = unique(c('KLRK1', 'KLRC4-KLRK1', 'CD8A', 'KLRD1', 'KLRC4', 'KLRC1', 'CD8B', 'KLRC3', 'NKG7', 'KLRC2', 'GZMB', 'CCL3', 'ANAPC1P1', 'ZNF683', 'KRT86', 'TNFRSF9', 'RRM2', 'VCAM1', 'KRT81', 'CCL4', 'KIR2DL4', 'CCL3', 'VCAM1', 'TNFRSF9', 'TNS3', 'NKG7', 'GZMB', 'SPRY2', 'CCL4', 'CCL4L1', 'CXCL13', 'AKAP5', 'LAG3', 'CD8A', 'IFNG', 'ITGA2', 'PHLDA1', 'SLC2A8', 'ID3', 'ZBED2', 'EGR2', 'CSF1', 'TNFSF9', 'EOMES', 'TNFSF4', 'SPRY2', 'PDCD1', 'CTLA4'))
ml_dysf = compare_mel_to_lung_corresponding_group(mc_t_nk, mat_t_nk_ds, umi_reg = 64, gene_labels_to_show = dysf_genes_to_show, points_col='grey', name="_manual")

treg_genes_to_show = c('BATF', 'CCL22', 'CCR8', 'CD177', 'CSF1', 'CSF2RB', 'CTLA4', 'CX3CR1', 'EBI3', 'ENTPD1', 'FANK1', 'FOXP3', 'GNG8', 'IKZF2', 'IL12RB2', 'IL1R1', 'IL1R2', 'IL1RN', 'IL2RA', 'LAIR2', 'LAYN', 'LTA', 'PHLDA1', 'TIGIT', 'TNFRSF13B', 'TNFRSF18', 'TNFRSF1B', 'TNFRSF4', 'TNFRSF8', 'TNFRSF9', 'VDR', 'ZBED2', 'ZBTB32', 'PDCD1')
ml_treg = compare_mel_to_lung_corresponding_group(mc_t_nk, mat_t_nk_ds, mel_grp='Treg', lung_grp='tumor-Treg', lung_ref_grp='blood-Treg', mel_genes=names(treg_genes), umi_reg = 64, gene_labels_to_show = treg_genes_to_show, points_col='grey', name="_manual")

y = read.table(scfigs_fn(guo_our_mc_clean_id, "mel_Treg_by_naive_vs_lung_tumor-Treg_by_naive-CXCR6_genes", ext="txt"), header=T)
y = y[order(y$lung_enr - y$mel_enr), ]
.plot_start(scfigs_fn(guo_our_mc_clean_id, "mel_Treg_by_naive_vs_lung_tumor-Treg_by_naive-CXCR6_genes_mel_enr_barplot"), 300, 400)
par(mar=c(4,8,1,1))
barplot(t(y[1:3, ]), beside=T, horiz=T, las=2, col=c(t_group2col['Treg'], 'grey'), ylim=c(0, 22 * 3))
legend("topleft", legend=c('Melanoma', 'Lung'), fill=c(t_group2col['Treg'], 'grey'), bty='n', ncol=2)
dev.off()

# prolifiration analysis
create_lateral_gset(suff="cc_filt", cl_to_keep=5:7)
guo_cc_gset = scdb_gset("guo2018_lateral_cc_filt")
guo_cc_genes = names(guo_cc_gset@gene_set)

mel_cc_gset = scdb_gset("mel_cc_filt")
mel_cc_genes = intersect(gg, names(mel_cc_gset@gene_set))

f_cc = colSums(mat_guo@mat[guo_cc_genes, names(mc_guo_ours@mc)]) /  colSums(mat_guo@mat[, names(mc_guo_ours@mc)])
f_by_mel_cc = colSums(mat_guo@mat[mel_cc_genes, names(mc_guo_ours@mc)]) /  colSums(mat_guo@mat[, names(mc_guo_ours@mc)])

cc_reg = 1e-3
lung_cc_cutoff = 0.015

lf_cc = log2(cc_reg + f_cc)
lf_by_mel_cc = log2(cc_reg + f_by_mel_cc)
.plot_start(scfigs_fn(guo_our_mc_clean_id, "f_cc_lung_vs_mel_gsets"), 400, 400)
plot(lf_cc, lf_by_mel_cc, pch=19, col=mc_guo_ours@colors[mc_guo_ours@mc[names(lf_cc)]], cex=0.5, xlab="% cc UMIs (log2) by lung gset", ylab="% cc UMIs (log2) by mel gset")
abline(v=log2(lung_cc_cutoff + cc_reg), lty=2)
abline(v=log2(f_cc_cutoff + cc_reg), lty=2, col='red')
dev.off()

# % prolif per group
p_cc_by_grp = sort(tapply(f_cc >= lung_cc_cutoff, go_col2group[mc_guo_ours@colors[mc_guo_ours@mc[names(f_cc)]]], mean))
n_cc_by_grp = tapply(f_cc >= lung_cc_cutoff, go_col2group[mc_guo_ours@colors[mc_guo_ours@mc[names(f_cc)]]], sum)

.plot_start(scfigs_fn(guo_our_mc_clean_id, "p_cc_per_grp"), 300, 600)
par(mar=c(5, 8, 1, 8))
barplot(p_cc_by_grp, col=go_group2col[names(p_cc_by_grp)], las=2, horiz=T)
mtext(n_cc_by_grp[names(p_cc_by_grp)], 4, line=0.5, at = seq(0.7, by=1.2, length=length(n_cc_by_grp)), las=2)
dev.off()

# %prolif on CD8-dysf strats
f_cd8_dysf_mc = mc_guo_ours@colors == go_group2col['CD8-dysf']
f_cd8_dysf_c = f_cd8_dysf_mc[mc_guo_ours@mc]

lung_lat_gset = scdb_gset(guo_lat_gset_id)
lung_lat_genes = names(lung_lat_gset@gene_set)

layn_de = diff_expr(mc_guo_ours, mat_guo@mat, which(go_col2group[mc_guo_ours@colors] == 'CD8-dysf'), which(mc_guo_ours@colors == go_group2col['CD8-GZMK']), filter_outlier_genes=T, geo_mean=T, min_max_umi=0.2, reg=0.1)
layn_genes = layn_de %>% filter(pmin(tot1, tot2) >= 175 & enr >= 1.5 & !(gene %in% lung_lat_genes))
layn_gs = layn_genes$enr
names(layn_gs) = layn_genes$gene

.plot_start(scfigs_fn(guo_our_mc_clean_id, "lung_dysf_genes"), 300, 100+nrow(layn_genes)*10)
par(mar=c(4,10, 1, 1))
barplot(sort(layn_gs), horiz=T, las=2, xlab='CD8 dysf over GZMK')
dev.off()

f_layn = colSums(mat_guo@mat[names(layn_gs), names(mc_guo_ours@mc)]) /  colSums(mat_guo@mat[, names(mc_guo_ours@mc)])

layn_cut = cut(f_layn[f_cd8_dysf_c], breaks=quantile(f_layn[f_cd8_dysf_c], seq(0, 1, by=0.2)), include.lowest=T, labels=paste0("Q", seq(0.2, 1, by=0.2)))
layn_strats = as.character(layn_cut)
names(layn_strats) = names(f_layn[f_cd8_dysf_c])

yy_dysf = tapply(f_cc[f_cd8_dysf_c] >= lung_cc_cutoff, layn_strats[names(mc_guo_ours@mc)][f_cd8_dysf_c], mean)
.plot_start(scfigs_fn(guo_our_mc_clean_id, "p_cc_by_CD8-dysf_strats"), 400, 400)
barplot(yy_dysf, col=go_group2col['CD8-dysf'], ylab="% prolif", main="CD8-dysf cells by %dysf UMIs")
dev.off()

# group patient composition
.plot_start(scfigs_fn(guo_our_mc_clean_id, "patient_grp_comp"), 900, 500)
layout(matrix(c(4, 1:3), nrow=1))
par(mar=c(4,1,4,4))
cell_by_src = split(names(mc_guo_ours@mc), mat_guo@cell_metadata[names(mc_guo_ours@mc), 'batch_set_id'])
patient_ord = NULL
for (src in c('Tumor', 'Normal', 'Blood')) {
	nms = cell_by_src[[src]]
	pat_grp = table(factor(mat_guo@cell_metadata[nms, 'Patient'], levels=unique(mat_guo@cell_metadata$Patient)), factor(go_col2group[mc_guo_ours@colors[mc_guo_ours@mc[nms]]], levels=names(ord_by_id[[guo_our_mc_clean_id]])))
	pat_grp_n = pat_grp / rowSums(pat_grp)
	if (is.null(patient_ord)) {
		patient_ord = rownames(pat_grp_n)[order(pat_grp_n[, 'CD8-dysf'])]
	}
	pat_grp_n = pat_grp_n[patient_ord, names(ord_by_id[[guo_our_mc_clean_id]])]
	barplot(t(pat_grp_n), col=go_group2col[colnames(pat_grp_n)], las=2, horiz=T, yaxt=ifelse(src == 'Tumor', 's', 'n'), cex.names=2, cex.main=2)
	title(main=src, cex.main=2)
	mtext(text = rowSums(pat_grp[rownames(pat_grp_n), ]), side = 4, line=0.5, at=seq(0.7, by=1.2, len=nrow(pat_grp)), las=2)
}
dev.off()
