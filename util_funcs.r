# 
#' Load library and additional packages, init metacell working directories
#'
#' @param scdb_dir 
#' @param scfigs_dir 
#'
#' @export
#'
rl = function(scdb_dir="scrna_db", scfigs_dir="figs_paper")
{
	library(Matrix)
	library(MASS)
	library(pheatmap)
	library(flowCore)
	library(metacell)
  library(dplyr)
	library(glmnet)
	library(plotrix)
	library(data.table)
	library(tgconfig)
	
  tgconfig::override_params("config/mel_all.yaml",package="metacell")

  dir.create(scdb_dir, showWarnings = F)
  scdb_init(scdb_dir, force_reinit=T)
  
  dir.create(scfigs_dir, showWarnings = F)
  scfigs_init(scfigs_dir)
}

#' Define global variables
#'
#' @return
#' @export
#'
def_vars = function() {
	
	# TODO remove redundant vars
	all_id <<- "mel_all"
	all_facs_id <<- paste0(all_id, "_facs")
	filt_id <<- "mel_filt"
	
	tumor_id <<- "mel_filt_Tumor"
	tumor_clean_id <<- "mel_filt_Tumor_outClean"
	tumor_t_nk_id <<- "mel_filt_Tumor_outClean_submc_T_NK"
	tumor_non_t_nk_id <<- "mel_filt_Tumor_outClean_submc_non_T_NK"
	tumor_non_t_nk_merged_id <<- "mel_filt_Tumor_outClean_submc_non_T_NK_merged"
	
	
	pbmc_tumors_comp_id <<- "mel_filt_tumors_pbmc_comp"
	pbmc_tumors_comp_clean_id <<- "mel_filt_tumors_pbmc_comp_outClean"
	
	mito_gset_id <<- "human_mito"
	ig_gset_id <<- "human_IG"
	ncrna_gset_id <<- "human_ncRNA"
	ncrp_gset_id <<- "human_ncRP"
	
	lateral_gset_id  <<- "mel_lateral"
	gset_nm <<- "Human_stress"
	
	# Guo 2018 vars
	guo_id <<- "guo2018_tpm_scaled"
	guo_filt_id <<- paste0(guo_id, "_filt")
	
	guo_our_mc_id <<- "guo2018_tpm_scaled_filt_Tumor_Normal_Blood"
	guo_our_mc_clean_id <<- "guo2018_tpm_scaled_filt_Tumor_Normal_Blood_outClean"
	
	guo_lat_gset_id <<- "guo2018_lateral_filtered"
	
	# Tirosh 2016 vars
	tirosh_clean_id <<- "tirosh2016_tpm_scaled_unresolved_malignant_non-malignant_outClean"

	pbmc_pats <<- c("h1---N", "h1---IFN100", "h1---IFN1000", "h2---N", "h2---IFN100", "h2---IFN1000", "p13_PBMC-3-(S)C-N", "p13-3-(S)C-N")
	
	all_tumor_pats <<- c("p1-4-LN-N",
	"p4-4-LN-1IT",
	"p2-4-LN-2IT",
	"p5-4-LN-2IT",
	"p3-4-LN-mIT",
	"p6-4-(S)C-1IT",
	"p9-4-(S)C-2IT",
	"p8-4-(S)C-mIT",
	"p10-4-MSC-2IT",
	"p11-3-(S)C-N",
	"p12-3-(S)C-N-1", "p12-3-(S)C-N-2",
	"p13-3-(S)C-N",
	"p15-3-(S)C-N",
	"p21-3-(S)C-N",
	"p17-3-(S)C-N-1", "p17-3-(S)C-N-2",
	"p16-3-(S)C-T",
	"p18-3-(S)C-T",
	"p19-3-(S)C-T",
	"p20-3-(S)C-1IT",
	"p23-3-LN-N",
	"p24-3-LN-N-2",
	"p25-3-p(S)C-N",
	"p26-2-p(S)C-N",
	"p27-4-(S)C-1IT",
	"p28-4-(S)C-2IT") 

	ggroup_cols <<- c('grey', 'red', 'royalblue', 'mediumorchid4', '#B3EE3A', 'gold')
	names(ggroup_cols) <<- c('None', 'ligand', 'receptor', 'TF', 'cytotoxic', 'exhaustion')

	ggroup_ord <<- 1:length(ggroup_cols)
	names(ggroup_ord) <<- names(ggroup_cols)

	all_nms <<- c("T", "NK", "B", "plasma", "myeloid", "pDC", "osteoclast", "melanocyte", "erythrocyte")

	t_nms <<- c("effector2", "dysfunctional", "effector1", "naive", "Treg", "Tfh", "dysf-cd4", "em-cd4")
	
	non_t_nms <<- c("macrophage", "osteoclast", "monocyte", "non-classic-monocyte", "mature-DC", "immature-DC", "DC", "pDC", "plasma", "B", "NK")
		
	non_t_merged_nms <<- c("macrophage", "osteoclast", "monocyte", "DC", "pDC", "plasma", "B", "NK")

	guo_our_nms <<- c("CD8-MAIT", "CD8-cyto", "CD8-dysf", "CD8-dysf-KLRC4", "CD8-GZMK", "naive-CXCR6", "naive-TCF7", "CD4-GZMK", "blood-Treg", "tumor-Treg", "Tfh", "CD4-cyto")
	
	tirosh_nms <<- c("Tumor", "Tumor/CAF", "Endothelial" , "B", "myeloid", "NK", "CD8-cyto", "naive", "CD8-transitional", "CD8-dysf", "CD4-Treg", "Tfh")
	
	ord_non_t_nms <<- 1:(length(non_t_nms)-1)
	names(ord_non_t_nms) <<- grep("NK", non_t_nms, invert=T, v=T)
	
	ord_non_t_merged_nms <<- 1:(length(non_t_merged_nms)-1)
	names(ord_non_t_merged_nms) <<- grep("NK", non_t_merged_nms, invert=T, v=T)
	
	ord_guo_our_nms <<- seq_along(guo_our_nms)
	names(ord_guo_our_nms) <<- guo_our_nms
	
	ord_tirosh_nms <<- seq_along(tirosh_nms)
	names(ord_tirosh_nms) <<- tirosh_nms
	
	cd8_nms <<- c("effector1", "effector2", "dysfunctional")
	
	ord_t_nk_nms <<- 1:(length(t_nms)+1)
	names(ord_t_nk_nms) <<- c("NK", t_nms)
	
	ord_all_nms <<- 1:length(all_nms)
	names(ord_all_nms) <<- all_nms
	ord_by_id <<- list()
	ord_by_id[[tumor_t_nk_id]] <<- ord_t_nk_nms
	ord_by_id[[tumor_non_t_nk_id]] <<- ord_non_t_nms
	ord_by_id[[tumor_non_t_nk_merged_id]] <<- ord_non_t_merged_nms
	ord_by_id[[tumor_clean_id]] <<-ord_all_nms
	ord_by_id[[guo_our_mc_clean_id]] <<- ord_guo_our_nms
	ord_by_id[[tirosh_clean_id]] <<- ord_tirosh_nms
	
	
	myl_nms <<- c("osteoclast", "pDC", "mature-DC", "DC", "immature-DC", "macrophage", "monocyte", "non-classic-monocyte")
	
	# classify patients by reactivity assay
	reac_pats <<- c("p13-3-(S)C-N","p23-3-LN-N","p27-4-(S)C-1IT","p4-4-LN-1IT","p9-4-(S)C-2IT", "p28-4-(S)C-2IT")
	non_reac_pats <<- c("p5-4-LN-2IT",  "p10-4-MSC-2IT", "p8-4-(S)C-mIT", "p3-4-LN-mIT")
	na_reac_pats <<- setdiff(all_tumor_pats, c(reac_pats, non_reac_pats))
	reac_pat_cols <<- c('darkgreen', 'grey', 'red')
	names(reac_pat_cols) <<- c('reactive', 'unknown', 'non-reactive')
	
	reac_col2type <<- c('reactive', 'unknown', 'non-reactive')
	names(reac_col2type) <<- c('darkgreen', 'grey', 'red')
	
	pat_reac_col <<- c(rep(reac_pat_cols['reactive'], length(reac_pats)), rep(reac_pat_cols['non-reactive'], length(non_reac_pats)), rep(reac_pat_cols['unknown'], length(na_reac_pats)))
	names(pat_reac_col) <<- c(reac_pats, non_reac_pats, na_reac_pats)
	
	min_t_nk_cells_per_patient <<- 100
	
	f_cc_cutoff <<- 0.0325
	
	t_nk_genes <<-  c('CD3D', 'CD4', 'CD8A', 'CD8B', 'TRAC', 'TRBC2', 'CTLA4', 'FOXP3', 'CCL5', 'NKG7', 'GZMA', 'TNFRSF4', 'IL7R', 'TCF7', 'IFNG', 'TIGIT', 'LAG3', 'PDCD1', 'KIR2DL4', 'ISG15', 'ITM2A', 'GZMB', 'CCL3', 'CCL4', 'CXCL13', 'GZMK', 'GZMH',  'KLRB1', 'CCR5', 'PRF1', 'IL32', 'FGFBP2', 'FCGR3A', 'CX3CR1', 'CD5', 'GNLY')
}



#' Compute differential gene expression between two groups of cells on a downsampled umi matrix
#'
#' @param mc metacell object 
#' @param mat_ds downsampled umi matrix object to use
#' @param mcs1 metacell ids of first group
#' @param mcs2 metacell ids of second group 
#' @param reg regulation value (added to total umi expression before calculating log-ratio)
#' @param min_max_umi min scaled total umis to filter genes by (taking max of the 2 groups)
#' @param nms1 cell names in first group (overriding mcs1)
#' @param nms2 cell names in second group (overriding mcs2)
#' @param filter_outlier_genes filter out genes that are enriched in small number of mcs
#' @param compare_top_mc_to_n_highest rank of mc to compare top expressing mc to, if filtering genes
#' @param max_top_to_n_highest_ratio if filtering genes, threshold on ratio of top to n'th mcs
#' @param verbose 
#' @param geo_mean apply geometric mean instead of sum of total umis
#' @param geo_mean_per_cell geometric mean per cell. False by default, and then multiply gene geometric mean by number of cells to get "total umi"
#'
#' @return data frame with gene name, number of scaled umis in each group, enrichment value
#' @export
#'
diff_expr = function(mc, mat_ds, mcs1=NULL, mcs2=NULL, reg=5, min_max_umi=50, nms1=NULL, nms2=NULL, filter_outlier_genes=F, compare_top_mc_to_n_highest=3, max_top_to_n_highest_ratio=3, verbose=T, geo_mean=F, geo_mean_per_cell=F)
{
	if (is.null(nms1)) {
		nms1 = names(mc@mc)[mc@mc %in% mcs1]
	}
	if (is.null(nms2)) {
		nms2 = names(mc@mc)[mc@mc %in% mcs2]
	}
	nms1 = intersect(colnames(mat_ds), nms1)
	nms2 = intersect(colnames(mat_ds), nms2)
	if (verbose) {
		message(sprintf("comparing %d vs %d cells", length(nms1), length(nms2)))
	}
	
	if (geo_mean) {
		df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), mu1=apply(mat_ds[, nms1], 1, function(y) {exp(mean(log(1+y)))-1}), mu2=apply(mat_ds[, nms2], 1, function(y) {exp(mean(log(1+y)))-1}), stringsAsFactors = F)
		df$tot1 = df$mu1 * length(nms1)
		df$tot2 = df$mu2 * length(nms2)
	} else {
		df = data.frame(row.names=rownames(mat_ds), gene=rownames(mat_ds), tot1=Matrix::rowSums(mat_ds[, nms1]), tot2=Matrix::rowSums(mat_ds[, nms2]), stringsAsFactors = F)
	}
	
	norm_by = min(sum(df$tot1), sum(df$tot2))
	df$tot1 = df$tot1 / sum(df$tot1) * norm_by
	df$tot2 = df$tot2 / sum(df$tot2) * norm_by
	df = df[pmax(df$tot1, df$tot2) >= min_max_umi, ]
	
	if (geo_mean && geo_mean_per_cell) {
		df$enr = log2( (df$mu1 + reg) / (df$mu2 + reg))
	} else {
		df$enr = log2( (df$tot1 + reg) / (df$tot2 + reg))
	}
	df = df[order(df$enr, decreasing=T), ]
	
	if (filter_outlier_genes) {
		fp = mc@mc_fp[intersect(rownames(mc@mc_fp), df$gene), ]
		if (!is.null(mcs1)) {
			fp = fp[, mcs1]
		}
		gmax = apply(fp, 1, max)
		gnext = apply(fp, 1, function(v) head(tail(sort(v), n=compare_top_mc_to_n_highest), n=1) )
		df[rownames(fp), 'out_r'] =  gmax/gnext
		to_filt = !is.na(df$out_r) & df$out_r > max_top_to_n_highest_ratio & df$enr > 0
		if (sum(to_filt) > 0) {
			if (verbose) {
				message(sprintf("filtering %d outlier genes:", sum(to_filt)))
			}
			print(df[to_filt, ])
			df = df[!to_filt, ]
		}
	}
	df	        
}


#' Scatter plot of 2 genes log footprint (metacell enrichment) values
#'
#' @param nm1 first gene name (or just x-lab title if x is supplied)
#' @param nm2 second gene name (or just x-lab title if y is supplied)
#' @param lfp gene-metacell enrichment matrix. Would usually be log2(mc@mc_fp) if mc is the metacell object
#' @param cols metacell colors. Vector length equal to ncol(lfp), Would usually be mc@colors
#' @param ofn output file name (plot to screen if NULL)
#' @param x x-values (optional, take lfp[nm1, ] if null)
#' @param y y-values (optional, take lfp[nm2, ] if null)
#' @param show_mc_ids show ids of metacells
#' @param cex 
#' @param cex.lab 
#' @param add_grid 
#' @param main 
#' @param xlim 
#' @param ylim 
#'
#' @export
#'
plt = function(nm1, nm2, lfp, cols, ofn=NULL, x=NULL, y=NULL, show_mc_ids=T, cex=3, cex.lab=1, add_grid=F, main="", xlim=NULL, ylim=NULL) {
	if (is.null(x)) {
		x = lfp[nm1, ]
	}
	if (is.null(y)) {
		y = lfp[nm2, ]
	}
	if (!is.null(ofn)) {
		.plot_start(ofn, 450, 450)
	}
	if (is.null(xlim)) {
		xlim = range(x)
	}
	if (is.null(ylim)) {
		ylim = range(y)
	}
	plot(x, y, pch=21, cex=cex, bg=cols, xlab=nm1, ylab=nm2, cex.lab=cex.lab, main=main, xlim=xlim, ylim=ylim)
	if (show_mc_ids) {
		text(x, y, colnames(lfp), cex=cex/4)
	}
	if (add_grid) {
		grid(col='black', lwd=0.5)
	}
	if (!is.null(ofn)) {
		dev.off()
	}
}


#' Perform logicle transform on columns of input data frame using flowCore::logicleTransform
#'
#' @param df 
#'
#' @return transformed data frame
#' @export
#'
.logicle_transform = function(df) {
	logicle  <- flowCore::logicleTransform("myTransform")
	
	if (is.null(ncol(df))) {
		df = as.data.frame(df)
		colnames(df) = 'v'
		rownames(df) = 1:nrow(df)
	}
	# add dummy column to keep it a df
	#ff = .df_to_ff(cbind(df, data.frame(dummy=df[,1])))
	ind = rowSums(is.na(df)) == 0
	m = as.matrix(df[ind, ])
	colnames(m) = colnames(df)
	ff = flowCore::flowFrame(m)
	
	ff_t = flowCore::transform(ff, flowCore::transformList(colnames(df), logicle))
	
	v = ff_t@exprs[, colnames(df)]
	if (ncol(df) == 1) {
		names(v) = rownames(df)[ind]
	}
	v
}


#' metacell group composition barplot by patient
#'
#' @param mc_id id of input metacell object
#' @param mat_id id of input mat objecy
#' @param groups to filter
#' @param name added to output figure
#' @param col2group optional named vector mapping colors to group names 
#' @param group2col optional named vector mapping group names to colors
#' @param pat_grp cell count table, patients on rows, groups on columns (optional, computed if null)
#' @param pat_grp_n normalized pat_grp (fraction of cells per patient, optional, computed if null)
#' @param pats_ord patient order (optional)
#' @param pat_groups seperate figures per patient group (optional)
#' @param clust_pats h-clust patients? (vector length equal pat_groups)
#' @param min_pat_cells min number of cells per patient (vector length equal pat_groups)
#'
#' @return order of patients
#' @export
#'
mc_group_composition_barplots = function(mc_id, mat_id, groups=c('Treg', 'Tfh', 'Naive', 'GZMK', 'GZMH', 'Naive-CT', 'NK', 'AKAP5'), name=NULL, col2group=NULL, group2col=NULL, pat_grp=NULL, pat_grp_n=NULL, pats_ord=NULL, pat_groups = list(Tumors=tumor_pats, PBMC=pbmc_pats), clust_pats = c(T,  F), min_pat_cells = c(100, 1))
{
	mc = scdb_mc(mc_id)
	mat = scdb_mat(mat_id)
	stopifnot(!is.null(mc) & !is.null(mat))
	
	if (is.null(col2group)) {
		col2group = as.character(mc@color_key$group)
		names(col2group) = as.character(mc@color_key$color)
	}
	
	if (is.null(group2col)) {
		group2col = as.character(mc@color_key$color)
		names(group2col) = as.character(mc@color_key$group)
	}
	
	if (is.null(pat_grp)) {
		pat_grp = table(mat@cell_metadata[names(mc@mc), 'PatientID'], col2group[mc@colors[mc@mc]])
		pat_grp = pat_grp[, groups]
	}
	
	if (is.null(pat_grp_n)) {
		pat_grp_n = pat_grp / rowSums(pat_grp)
	}
	
	pat_grp_n_ord = list()
	bp_spaces = list()
	valid_ind = list()

	for (i in seq_along(pat_groups)) {
		g_nm = names(pat_groups)[i]
		pat_grp_n2 = pat_grp_n[pat_groups[[i]], ]

		if (clust_pats[i]) {
			hc = hclust(dist(cor(t(pat_grp_n2))), method='ward.D2')
			ord = hc$order
		} 
		else if (is.null(pats_ord[[g_nm]])) {
			ord = 1:nrow(pat_grp_n2)
		}
		else {
			ord = pats_ord[[g_nm]]
		}
		pats_ord[[g_nm]] = rownames(pat_grp_n2[ord, ])
		
		pat_grp_n_ord [[ g_nm ]]  = pat_grp_n2[ord, ]
		bp_spaces[[g_nm]] = NULL
		valid_ind[[g_nm]] = rowSums(pat_grp[rownames(pat_grp_n2[ord,]), ]) >= min_pat_cells[i]
	}
	if (length(intersect(c('PBMC', 'Tumors'), names(pat_groups))) == 2) {
		pat_grp_n_ord[["Tumors_and_PBMC"]] = rbind(pat_grp_n_ord[["PBMC"]][1:7, ], pat_grp_n_ord[["Tumors"]])
		bp_spaces[["Tumors_and_PBMC"]] = c(rep(0.2, 7), 0.8, rep(0.2, nrow(pat_grp_n_ord[["Tumors"]])-1))
		valid_ind[["Tumors_and_PBMC"]] = c(valid_ind[["PBMC"]][1:7], valid_ind[["Tumors"]])
	}
	for (i in seq_along(pat_grp_n_ord)) {
		g_nm = names(pat_grp_n_ord)[i]
		pat_grp_n2 = pat_grp_n_ord[[g_nm]]
		
		pat_grp_n2[!valid_ind[[g_nm]] , ] = NA
		rownames(pat_grp_n2) = sprintf("%s (%d)", rownames(pat_grp_n2), rowSums(pat_grp[rownames(pat_grp_n2),]))
		
		.plot_start(scfigs_fn(mc_id, sprintf('mc_group_%s%s', g_nm, ifelse(is.null(name), "", paste0("_", name)))), 500, max(30 * nrow(pat_grp_n2), 600))
		par(mar=c(4,20, 4, 2))

		barplot(t(pat_grp_n2), col=group2col[colnames(pat_grp_n2)], las=2, horiz=T, main=g_nm, cex.names=1.5, space=bp_spaces[[g_nm]])
		dev.off()
		
		subtype_cor = cor(log2(1e-3+pat_grp_n2[rowSums(is.na(pat_grp_n2)) == 0, ]))
		.plot_start(scfigs_fn(mc_id, sprintf('mc_group_cor_%s%s', g_nm, ifelse(is.null(name), "", paste0("_", name)))), 100 + ncol(subtype_cor) * 50, 100 + ncol(subtype_cor) * 50)
		subtype_cor[is.na(subtype_cor)] = 0
		pheatmap(subtype_cor, display_numbers=T, breaks=seq(-1, 1, length=101), treeheight_row=20, treeheight_col=0, fontsize=16)
		dev.off()
	}
	
	pats_ord
}

#' Wrapper function for assigning groups (and colors) to metacells based on their confusion matrix clustering. 
#' Run first to create the clustered confusion matrix, use it and the mc_sup object to classify sup_ids and create the supmc_file (and optionally the marks_file) and then rerun it to assign colors to the metacell object
#' 
#' @param mc_id id of the metacell object (both used as the input and as the output object)
#' @param graph_id id of the cell knn graph object
#' @param supmc_file file name assigning name and colors to clusters of metacells (optional, if null then no coloring is done. tab delimited with these columns: supid, color, name)
#' @param marks_file file name assigning name and colors to metacells by thresholding expression of genes (optional, will override colors assigned by the supmc_file, tab-delimited with these columns: name, gene, color, T_fold. T_fold is the threshold on the gene log2 fp value).
#' @param res return value of this function (can be supplied after the first run to shorten run-time)
#' @param show_mc_ids in heatmap plot
#'
#' @return list containining mc_hc (hierrarchical clustering of the metacells) and mc_sup (info on derived clusters of metacells) 
#' 
colorize_by_confusion_mat = function(mc_id = t_nk_id, graph_id=filt_id, supmc_file=NULL, marks_file=NULL, res=NULL, show_mc_ids=F) 
{
	# Cluster metacells by confusion matrix
	if (is.null(res)) {	
		mc_hc = mcell_mc_hclust_confu(mc_id=mc_id,
																	graph_id=graph_id)
	}
	else {
		mc_hc = res$mc_hc
	}
	
	# Annotate metacell clusters (globally and locally enriched genes)
	if (is.null(res)) {
		mc_sup = mcell_mc_hierarchy(mc_id=mc_id,
																mc_hc=mc_hc, T_gap=0.04)
	}
	else {
		mc_sup = res$mc_sup
	}

	# Colorize metacells based on the manually created supmc_file (and optionally by the marks_file)
	if (!is.null(supmc_file)) {
		mc_colorize_sup_hierarchy(mc_id=mc_id,
															supmc = mc_sup,
															supmc_key = supmc_file,
															gene_key= marks_file)
	}

	# generate metacell clusters heatmap
	mcell_mc_plot_hierarchy(mc_id=mc_id,
													graph_id=graph_id,
													mc_order=mc_hc$order,
													sup_mc = mc_sup,
													width=1600, heigh=3000, min_nmc=2, show_mc_ids=show_mc_ids)
	
	list(mc_hc=mc_hc, mc_sup=mc_sup)
}

#' Wrapper for mcell_mc_plot_marks
#'
#' @param mc_id 
#' @param mat_id 
#' @param lateral_gset_id 
#' @param min_gene_fold 
#' @param k_per_clust 
#' @param text_cex 
#'
#' @export
#'
mc_gen_marks_heatmaps = function(mc_id, mat_id, lateral_gset_id=NULL, min_gene_fold=tgconfig::get_param("scm_mc_mark_min_gene_fold", "metacell"), k_per_clust=tgconfig::get_param("scm_mc_mark_k_per_clust", "metacell"), text_cex=tgconfig::get_param("mcp_heatmap_text_cex", "metacell")) {
	
	orig_min_gene_fold = tgconfig::get_param("scm_mc_mark_min_gene_fold", "metacell")
	orig_k_per_clust = tgconfig::get_param("scm_mc_mark_k_per_clust", "metacell")
	orig_text_cex = tgconfig::get_param("mcp_heatmap_text_cex", "metacell")
	
	tgconfig::set_param("scm_mc_mark_min_gene_fold", min_gene_fold, "metacell")
	tgconfig::set_param("scm_mc_mark_k_per_clust", k_per_clust, "metacell")
	tgconfig::set_param("mcp_heatmap_text_cex", text_cex, "metacell")
	
	if (is.null(lateral_gset_id)) {
		mcell_gset_from_mc_markers(gset_id=paste0(mc_id, "_marks"), mc_id=mc_id)
		for (plot_cells in c(T,F)) {
			mcell_mc_plot_marks(mc_id, paste0(mc_id, "_marks"), mat_id, plot_cells=plot_cells)
		}
	}
	else {
		mcell_gset_from_mc_markers(gset_id=paste0(mc_id, "_marks"), mc_id=mc_id, blacklist_gset_id=lateral_gset_id)
		mcell_gset_from_mc_markers(gset_id=paste0(mc_id, "_marks_lateral"), mc_id=mc_id, filt_gset_id=lateral_gset_id)
		
		for (plot_cells in c(T,F)) {
			mcell_mc_plot_marks(mc_id, paste0(mc_id, "_marks"), mat_id, lateral_gset_id=paste0(mc_id, "_marks_lateral"), plot_cells=plot_cells)
		}
	}
	
	tgconfig::set_param("scm_mc_mark_min_gene_fold", orig_min_gene_fold, "metacell")
	tgconfig::set_param("scm_mc_mark_k_per_clust", orig_k_per_clust, "metacell")
	tgconfig::set_param("mcp_heatmap_text_cex", orig_text_cex, "metacell")
	
}

#' Assign input genes to groups (ligand, receptor, cytotoxic, TF, exhausted)
#'
#' @param genes input genes
#' @param tfs_file file name of TFs (single column, no header)
#' @param add_tfs vector of TFs to add
#' @param remove_tfs vector of TFs to remove
#'
#' @return data frame with genes, assigned group and color
#' @export
#'
annotate_genes = function(genes, tfs_file="data/Lambert2018_TF_names_v_1.01.txt", add_tfs=c('TOX', paste0('ID', 1:4)), remove_tfs=paste0("TET", 1:3)) 
{
	tfs = read.table(tfs_file, stringsAsFactors = F)$V1
	if (!is.null(add_tfs)) {
		tfs = union(tfs, add_tfs)
	}
	if (!is.null(remove_tfs)) {
		tfs = setdiff(tfs, remove_tfs)
	}
	ggroups = list(ligand=grep("^IL[0-9]+$|^CCL[0-9]|^CXCL[0-9]|^CX3CL[0-9]|^TNFSF|^XCL", genes, v=T, perl=T),
									 receptor=grep("^IL[0-9]+R|^CCR[0-9]|^CXCR[0-9]|^CX3CR[0-9]|^TNFRSF|^XCR", genes, v=T, perl=T),
									 #putativeTF=intersect(genes, rownames(tfs)[!tfs$Validated]),
									 TF=intersect(genes, tfs),
									 cytotoxic=grep("^GZM|^KLR|^CTS|GNLY|PRF1", genes, v=T, perl=T),
									 exhaustion=c('ENTPD1', 'CTLA4', 'TIGIT', 'ICOS', 'BTLA', 'HAVCR2', 'PDCD1', 'LAG3', 'CD244', 'CD96', 'TNFRSF9'))
	
	df = data.frame(row.names=genes, color=rep(ggroup_cols['None'], length(genes)), type=rep('None', length(genes)), stringsAsFactors = F)

	for (gg in names(ggroups)) {
		cgs = intersect(genes, ggroups[[gg]])
		if (length(cgs) > 0) {
			color_crash = names(which(df[cgs, 'color'] != ggroup_cols['None']))
			if (length(color_crash) > 0) {
				message(sprintf("Coloring by %s, overwriting %d previously colored (%s)", gg, length(color_crash), paste0(color_crash, collapse=", ")))
			}
			df[cgs, 'color'] = ggroup_cols[gg]
			df[cgs, 'type'] = gg
		}
	}
	df
}

#' Differential gene expression between groups of mcs (one group against the rest). Also computes total umi count per differentially expressed genes.
#'
#' @param mc_id 
#' @param mat_id 
#' @param lateral_gset_id 
#' @param n_ds 
#' @param mat_ds 
#' @param min_log2_enr 
#' @param min_mean_umi 
#' @param groups 
#' @param filter_outlier_genes 
#' @param compare_top_mc_to_n_highest 
#' @param max_top_to_n_highest_ratio 
#' @param vgel_group_ord 
#' @param blist_genes 
#' @param gen_vgels 
#' @param plot_contours 
#'
#' @return
#' @export
#'
gene_diff_expr_by_groups = function(mc_id, mat_id, lateral_gset_id, n_ds=500, mat_ds=NULL, min_log2_enr=1, min_mean_umi=0.1, groups=c('effector1', 'effector2', 'dysfunctional'), filter_outlier_genes=T,  compare_top_mc_to_n_highest=3, max_top_to_n_highest_ratio=3, vgel_group_ord=ord_t_nk_nms, blist_genes=NULL, gen_vgels=F, plot_contours=F) 
{
	mc = scdb_mc(mc_id)
	stopifnot(!is.null(mc))
	
	mat = scdb_mat(mat_id)
	stopifnot(!is.null(mat))
	
	if (is.null(mat_ds)) {
		mat_ds = scm_downsamp(mat@mat, n_ds)
	}
	
	col2group = as.character(mc@color_key$group)
	names(col2group) = as.character(mc@color_key$color)
	
	group2col = as.character(mc@color_key$color)
	names(group2col) = as.character(mc@color_key$group)
	
	lateral_gset = scdb_gset(lateral_gset_id)
	
	mat_ds = mat_ds[setdiff(rownames(mat_ds), unique(c(mat@ignore_genes, names(lateral_gset@gene_set), blist_genes))), ]
	
	de = list()
	mcsz = table(mc@mc)
	cells = c()
	for (g in groups) {
		message(g)		
		
		n_cells = sum(mcsz[mc@colors == group2col[g]])
		cells = c(cells, names(mc@mc)[mc@colors[mc@mc] == group2col[[g]]])
		
		de[[g]] = diff_expr(mc, mat_ds, which(mc@colors == group2col[g]), which(mc@colors %in% group2col[setdiff(groups, g)]), filter_outlier_genes = filter_outlier_genes,  compare_top_mc_to_n_highest=compare_top_mc_to_n_highest, max_top_to_n_highest_ratio=max_top_to_n_highest_ratio) %>%
				filter(tot1 >= n_cells * min_mean_umi & enr >= min_log2_enr)
	}
	if (gen_vgels) {
		mcell_mc_plot_vgels(mc_id, unique(unlist(lapply(de, function(x) x$gene))), reorder_preset=group2col[names(vgel_group_ord)], lane_w=30, height=200, plot_title=T)
	}
	
	bad_genes = names(which(table(unlist(lapply(de, function(x) x$gene))) > 1))
	if (length(bad_genes) > 0) {
		message(sprintf("removing %d multiple group genes (%s)", length(bad_genes), paste0(bad_genes, collapse=", ")))

		for (g in groups) {
			deg = de[[g]]
			de[[g]] = deg[!deg$gene %in% bad_genes, ]
		}
	}
	genes_ann = data.frame(row.names=unlist(lapply(de, function(x) x$gene)), stringsAsFactors=F)
	
	all_cells = cells
	totu = matrix(NA, nrow=length(groups), ncol=length(all_cells), dimnames=list(groups, all_cells))
	totu_ds = totu
	ann_cols = list()
	maxg = max(unlist(lapply(de, nrow)))
	
	for (g in groups) {
		curr_g = de[[g]]$gene
		if (length(curr_g) > 1) {
			totu[g, ] = Matrix::colSums(mat@mat[curr_g, all_cells])
			totu_ds[g, ] = Matrix::colSums(mat_ds[curr_g, all_cells])
			
			.plot_start(scfigs_fn(mc_id, sprintf("gr_enr_genes_%s", g)), 300, 100+ 12 * maxg)
			par(mar=c(4,8,4,1))
			barplot(rev(de[[g]]$enr), horiz=T, names.arg=rev(de[[g]]$gene), xlab='enr (log2)', las=2, ylim=c(0, maxg+2), width=1, col=group2col[g])
			dev.off()
			
		}
		else if (length(curr_g) == 1) {
			totu[g, ] = mat@mat[curr_g, all_cells]
			totu_ds[g, ] = mat_ds[curr_g, all_cells]
		}
		genes_ann[, g] = ifelse(rownames(genes_ann) %in% de[[g]]$gene, 1, 0)
		ann_cols[[g]] = c('white', group2col[g])
		names(ann_cols[[g]]) = 0:1
	}
	if (colSums(genes_ann) > 0) {
		genes_ann = genes_ann[, colSums(genes_ann) > 0]
		gs_cor_c = cor(t(log2(1+7*as.matrix(mat_ds[rownames(genes_ann), cells]))))
		ghc = hclust(dist(gs_cor_c), method='ward.D2')
		
		
		.plot_start(scfigs_fn(mc_id, sprintf("gr_enr_genes_%s", gsub(":", "_", paste0(groups, collapse="_")))), max(700, 400+nrow(gs_cor_c)* 8), max(700, 200+nrow(gs_cor_c)* 8))
		z_gs_cor_c = quantile(gs_cor_c, 0.98)
		pheatmap(pmin(pmax(gs_cor_c[ghc$order, ghc$order], -z_gs_cor_c), z_gs_cor_c), cluster_rows=F, cluster_cols=F, cellwidth=8, cellheight=8, annotation_col=genes_ann, annotation_row=genes_ann, annotation_colors = ann_cols)
		dev.off()
	}
	
	totu_n = log2(1 + 7000 * t(t(totu) / colSums(mat@mat[, colnames(totu)])))
	all_gc = split(all_cells, col2group[mc@colors[mc@mc[all_cells]]])
	if (plot_contours) {
		for (i in seq(1, length(groups)-1)) {
			for (j in seq(i+1, length(groups))) {
				nm1 = groups[i]
				nm2 = groups[j]
				if (sum(totu[nm1,]) > 0 & sum(totu[nm2, ]) > 0) {
					.plot_start(scfigs_fn(mc_id, paste0(gsub(":.*", "", nm1, perl=T), "_vs_", gsub(":.*", "", nm2, perl=T), "_totg")), length(groups) * 300, 300)
					layout(matrix(1:length(groups), nrow=1, byrow = T))
					par(mar=c(4,4,3,1))
					for(g in groups) {
					#for (g in c(t_nms, 'NK')) {
						plot(totu_n[nm1,], totu_n[nm2, ], main=sprintf("%s (%d)", g, length(all_gc[[g]])), xlim=range(totu_n[nm1,]), ylim=range(totu_n[nm2,]), pch=19, cex=0.3, col='grey', xlab=nm1, ylab=nm2)
						points(totu_n[nm1, all_gc[[g]]], totu_n[nm2, all_gc[[g]]], col=group2col[g], pch=19, cex=0.3); 
						nbins = max(25, round(length(all_gc[[g]])/100))
						d = MASS::kde2d(totu_n[nm1, all_gc[[g]]], totu_n[nm2, all_gc[[g]]], n=nbins)
						contour(d$x, d$y, d$z, drawlabels=F, add=T)
					}
					dev.off()
				}
			}
		}
	}
	
	list(de=de, totu=totu, totu_ds=totu_ds)
}

#' Plot metacell information (Fig 3b)
#'
#' @param mc_pat cell per metacell (column) and patient (row)
#' @param mc_colors 
#' @param mc_groups 
#' @param is_mc_t_nk 
#' @param mat_id 
#' @param ofn 
#' @param order_mcs_by 
#' @param pats_order 
#' @param groups_order 
#' @param width 
#' @param height 
#' @param min_mc_f_to_show 
#' @param miss_col 
#' @param mc_pat_zlim 
#' @param min_cells_in_mc_per_pat 
#' @param max_log2_p 
#' @param show_is_t 
#' @param cex_f 
#' @param ncells_ylim 
#' @param pat_comp_colspec 
#'
#' @export
#'
plot_mc_metadata = function(mc_pat, mc_colors, mc_groups, is_mc_t_nk, mat_id, ofn, order_mcs_by="group", pats_order=NULL, groups_order=NULL, width=NULL, height=720, min_mc_f_to_show=0.5, miss_col='grey90', mc_pat_zlim=0.8, min_cells_in_mc_per_pat=2, max_log2_p=10, show_is_t=F, cex_f=0.8, ncells_ylim=NULL, pat_comp_colspec=c('white', RColorBrewer::brewer.pal(n=9, 'YlOrRd')))
{

	mat = scdb_mat(mat_id)
	stopifnot(!is.null(mat))
	
	mc_pat_n = mc_pat / rowSums(mc_pat)

	if (order_mcs_by == 'group') {
		mc_order = order(as.numeric(ordered(mc_groups, levels=groups_order)) + seq_along(mc_groups) * 1e-6)
	}
	else if (mc_order_by == 'hclust') {
		mc_hc = hclust(dist(cor(t(mc_pat_n))), method='ward.D2')
		mc_order = mc_hc$order
	}
	
	ann_col = tgconfig::get_param("mcp_metadata_annot_colors", "metacell")
	if (is.null(width)) {
		width = 150 + length(mc_colors) * 4
	}
	
	.plot_start(ofn, w=width, h=height)
	heights = c(2,3,10)
	if (show_is_t) {
		heights = c(heights, 1)
	}
	heights = c(heights, 0.7)
	
	layout(matrix(1:(4 + ifelse(show_is_t, 1, 0)), ncol=1), heights=heights)
	par(mar=c(0.5,20,2,1))
	
	# 1 mc size "barplot"
	mcsz = rowSums(mc_pat)
	if (is.null(ncells_ylim)) {
		ncells_ylim = c(0, max(mcsz))
	}
	image(seq_along(mcsz), seq(ncells_ylim[1], ncells_ylim[2], length=100), matrix(0, length(mcsz), 100), xaxt='n', xlab="", ylab="", col=NA, cex.axis=1.5*cex_f)
	title(ylab="#cells", line=5, cex.lab=2*cex_f)
	rect(seq_along(mcsz)-0.5, 0, seq_along(mcsz)+0.5, mcsz[mc_order], border=NA, col='lightblue')
	abline(h=median(mcsz), lty=2)
	
	# 2 mc patients contribution
	image(1:nrow(mc_pat_n), seq(0, 1, length=ncol(mc_pat_n)), mc_pat_n, col=NA, xaxt='n', xlab="", ylab="", frame.plot=F, cex.axis=1.5*cex_f, las=2)
	title(ylab="patient contribution", line=5, cex.lab=2*cex_f)
	for (i in 1:nrow(mc_pat_n)) {
		v = sort(mc_pat_n[mc_order[i], ], decreasing=T)
		cumv = cumsum(v)
		rect(xleft=i-0.5, xright=i+0.5, ybottom=c(0, cumv[-length(cumv)]), ytop=cumv, border='lightgrey', col=c('darkgreen', 'lightgreen', rep('white', ncol(mc_pat_n)-2)))
	}
	
	par(mar=c(0.5,20,0.5,1))
	
	# 3 mc on patient tab
	if (is.null(pats_order)) {
		pats_order = colnames(mc_pat_n)
	}
	mc_pat_n = mc_pat_n[, pats_order]

	reg = min(mc_pat_n[mc_pat_n > 0])
	
	pat_comp_colors = colorRampPalette(pat_comp_colspec)(200)
	message(sprintf("trimming %d vals to %.2f", sum(mc_pat_n > mc_pat_zlim), mc_pat_zlim))
	image(log2(reg+pmin(mc_pat_n[mc_order, ], mc_pat_zlim)), zlim=log2(reg + c(0, mc_pat_zlim)), xaxt='n', xlab="", ylab="", yaxt="n", col=pat_comp_colors)
	mtext(pats_order, side=2, line=0.5, at=seq(0, 1, length=length(pats_order)), las=2, cex=1.5*cex_f)	
	abline(v=-0.5 / (length(mc_groups)-1) + cumsum(rle(mc_groups[mc_order])$lengths)/(length(mc_groups)-1))
	abline(h=-0.5 / (ncol(mc_pat_n)-1) + seq(0, 1, length=ncol(mc_pat_n)), lwd=0.5)
	
	# 4 is mc T/NK?
	if (show_is_t) {
		image(as.matrix(1:length(mc_colors)), col=ifelse(is_mc_t_nk[mc_order], 'orange', 'white'), xaxt='n', yaxt='n', xlab="", ylab="T/NK", cex.lab=2*cex_f)
	}
	
	# 5 mc colors
	image(as.matrix(1:length(mc_colors)), col=mc_colors[mc_order], xaxt='n', yaxt='n', xlab="", ylab="")
	
	md = mat@cell_metadata %>% 
		group_by(PatientID, stage, location) %>%
		summarize(ncells=length(PatientID)) %>%
		data.frame
	rownames(md) = md$PatientID
	dev.off()

	.plot_start(gsub(".png", "_log_legend.png", ofn), 200, 800)
	image(t(log2(reg + seq(0, mc_pat_zlim, length=length(pat_comp_colors)))), col=pat_comp_colors, xaxt='n', yaxt='n')
	box(lty=1, lwd=3)
	axis(4, at=seq(0, 1, length=3), labels=round(seq(0, mc_pat_zlim, length=3), 2), las=2)
	dev.off()

}


#' Wrapper function for creating gene sets from the supplied anchor genes
#'
#' @param mat_id 
#' @param gene_anchors 
#' @param gset_nm 
#' @param cor_thresh 
#' @param gene_anti 
#' @param sz_cor_thresh 
#' @param nclusts 
#'
#' @export
#'
select_gene_modules_by_anchor_genes = function(mat_id, gene_anchors, gset_nm, cor_thresh = 0.1, gene_anti=c(), sz_cor_thresh=0.7, nclusts=20) 
{
		tab_fn = sprintf("%s/%s.txt", scfigs_dir(mat_id, "gmods_by_anchors"), gset_nm)
		message("mcell_mat_rpt_cor_anchors")
		mcell_mat_rpt_cor_anchors(mat_id=mat_id,
															gene_anchors = gene_anchors,
															cor_thresh = cor_thresh,
															gene_anti = gene_anti,
															tab_fn = tab_fn,
															sz_cor_thresh=sz_cor_thresh)
	
	
		foc_gcor = read.table(tab_fn, sep="\t", h=T, stringsAsFactors=F)
		foc_genes = apply(foc_gcor[,setdiff(colnames(foc_gcor),c("sz_cor","max","neg_max"))], 1, which.max)
		gset = tgGeneSets(foc_genes, gset_nm)
		
		scdb_add_gset(gset_nm, gset)
		
		sub_mat_id = paste(mat_id, gset_nm, sep="_")
		
		mcell_mat_ignore_genes(sub_mat_id, mat_id, names(foc_genes), reverse=T)
		
		message("mcell_gset_split_by_dsmat")
		mcell_gset_split_by_dsmat(gset_nm, sub_mat_id, nclusts)
		mcell_plot_gset_cor_mats(gset_nm, sub_mat_id)
	
}


#' Control for clone sharing by mc group. Sample cells from each patient according to the probablity to recover a tcr by %tcr_genes, re-shuffling clonal assignment, yet maintaining number and size of clones. repeat n_iters times
#'
#' @param mc_id 
#' @param t_tcr produced by tcr_clones.r
#' @param mat_id 
#' @param n_iters 
#' @param tcr_genes 
#' @param sym_mat symmetrize matrix
#' @param patients 
#' @param ignore_groups 
#'
#' @return summed counts from all patients and iterations
#' @export
#'
shuffle_clones_mc_pairs_by_patient = function(mc_id, t_tcr, mat_id=mc_id, n_iters=100, tcr_genes=c('TRAC', 'TRBC2'), sym_mat=T, patients=NULL, ignore_groups='NK') 
{
	mc = scdb_mc(mc_id)
	mat = scdb_mat(mat_id)
	
	group2col = get_mc_group2col(mc)
	
	if (is.null(ignore_groups)) {
	 	cells = names(mc@mc) 
	} else {
		cells = names(mc@mc)[!(mc@colors[mc@mc] %in% ignore_groups)]
		t_tcr = t_tcr[!(t_tcr$mc_grp %in% ignore_groups), ]
	}
	tcr_sc = colSums(mat@mat[tcr_genes, cells]) / colSums(mat@mat[, cells])
	tcr_sc_bin = cut(tcr_sc, breaks=c(0, seq(min(tcr_sc[tcr_sc > 0]), quantile(tcr_sc, 0.98), length=8), max(tcr_sc)), include.lowest=T)
	names(tcr_sc_bin) = names(tcr_sc)
	
	p_tcr_by_bin = tapply(names(tcr_sc) %in% rownames(t_tcr), tcr_sc_bin, mean)
	
	valid_clones = t_tcr %>% group_by(clone_id) %>% summarise(clone_size = length(Well_ID)) %>% filter(clone_size > 1) %>% data.frame 
	t_tcr_gt2 = t_tcr %>% inner_join(valid_clones, by="clone_id")
	
	t_tcr_gt2_by_pat = split(t_tcr_gt2, t_tcr_gt2$PatientID)
	mc_by_pat = split(mc@mc[cells], mat@cell_metadata[cells, 'PatientID'])
	
	mc_mc = matrix(0, length(mc@colors), length(mc@colors), dimnames=list(seq_along(mc@colors), seq_along(mc@colors)))
	
	if (is.null(patients)) {
		patients = names(t_tcr_gt2_by_pat)
	}
	n_pairs = 0
	
	for (p in patients) {
		
		pat_mcs = mc_by_pat[[p]]
		pat_p_mcs = p_tcr_by_bin[tcr_sc_bin[names(pat_mcs)]]
		
		pat_cln_ids = t_tcr_gt2_by_pat[[p]]$clone_id
		
		if (!is.null(pat_cln_ids)) {
			inds = combn(seq_along(pat_cln_ids), 2)
			inds = inds[,  pat_cln_ids[inds[1,]] == pat_cln_ids[inds[2,]]]		
			
			message(sprintf("%s, %d cells, %d clones, %d edges", p, length(pat_mcs), length(unique(pat_cln_ids)), ncol(inds)))
	
			for (i in 1:n_iters) {
				s_mcs = sample(pat_mcs, length(pat_cln_ids), replace=F, prob=pat_p_mcs)	
				n_pairs = n_pairs + ifelse(is.null(nrow(inds)), 1, ncol(inds))
				if (is.null(nrow(inds))) {
					mc_mc[ s_mcs[inds[1]], s_mcs[inds[2]]] = mc_mc[ s_mcs[inds[1]], s_mcs[inds[2]]] + 1
				} else {
					mc_pairs = data.frame(a=s_mcs[inds[1,]], b=s_mcs[inds[2,]]) %>% group_by(a, b) %>% summarize(n=length(a)) %>% data.frame
					mc_mc[ cbind(mc_pairs$a, mc_pairs$b)] = mc_mc[ cbind(mc_pairs$a, mc_pairs$b)] + mc_pairs$n
				}
			}
		}
	}
	
	message(sprintf("n_pairs = %d used %d (%f)", n_pairs, sum(mc_mc), sum(mc_mc)/n_pairs))
	if (sym_mat) {
		mc_mc = mc_mc + t(mc_mc)
	}
	
	return(mc_mc)
}


#' Return mc color to group "dictionary"
#'
#' @param mc 
#'
#' @export
#'
get_mc_col2group = function(mc) {
	col2group = as.character(mc@color_key$group)
	names(col2group) = as.character(mc@color_key$color)
	col2group
}

#' Return mc group to color "dictionary"
#'
#' @param mc 
#'
#' @export
#'
get_mc_group2col = function(mc) {
	group2col = as.character(mc@color_key$color)
	names(group2col) = as.character(mc@color_key$group)
	group2col
}


#' Generate a pseudo-metacell object by using input metacells group names as new metacell assignments
#'
#' @param mc_id 
#' @param mat_id 
#' @param suffix 
#' @param ignore_mcs 
#' @param min_gene_fold 
#' @param k_per_clust 
#' @param text_cex 
#' @param group_ord 
#' @param rebuild 
#'
#' @export
#'
mc_by_color_group = function(mc_id=tumor_t_nk_id, mat_id=mc_id, suffix="_by_color_groups", ignore_mcs=NULL, min_gene_fold=1, k_per_clust=5, text_cex=4, group_ord=ord_t_nk_nms, rebuild=F)
{
	mc_g_id = paste0(mc_id, suffix)
	mc_g = scdb_mc(mc_g_id)
	if (is.null(mc_g) | rebuild) {
		mc = scdb_mc(mc_id)
		if (!is.null(ignore_mcs)) {
			mc = mc_set_outlier_mc(mc, ignore_mcs)
		}
		col2group = get_mc_col2group(mc)
		group2col = get_mc_group2col(mc)
		
		mat = scdb_mat(mat_id)
		mcg = mc
		mcg@mc = group_ord[col2group[mc@colors[mc@mc]]]
		names(mcg@mc) = names(mc@mc)
		
		mc_g = tgMCCov(mcg@mc, mc@outliers, mat)
		mc_g@colors = group2col[names(group_ord)]
		scdb_add_mc(mc_g_id, mc_g)
	}
	
	mc_gen_marks_heatmaps(mc_g_id, mat_id, lateral_gset_id=NULL, min_gene_fold=min_gene_fold, k_per_clust=k_per_clust, text_cex=text_cex) 
	
}

#' Plot %umis of a list of genes per patient
#'
#' @param mc_ids can supply several mc_ids, will sum up data per patient for all of them
#' @param genes 
#' @param genes_nm 
#' @param genes_col 
#' @param mat_ids 
#' @param reg 
#' @param min_cells 
#' @param summary_func 
#'
#' @export
#'
plot_f_tot_genes_by_group = function(mc_ids, genes, genes_nm, genes_col, mat_ids=mc_ids, reg=1e-3, min_cells=20, summary_func=median)
{
	cells = c()
	f_genes = c()
	patients = c()
	groups = c()
	
	group2col = c()
	
	for (i in seq_along(mc_ids)) {
		mc_id = mc_ids[i]
		mat_id = mat_ids[i]
		mc = scdb_mc(mc_id)
		c_col2group = get_mc_col2group(mc)
		c_group2col = get_mc_group2col(mc)
	
		mat = scdb_mat(mat_id)

		cells = c(cells, names(mc@mc))
		f_genes = c(f_genes, colSums(mat@mat[genes, names(mc@mc)]) / colSums(mat@mat[, names(mc@mc)]))
		patients = c(patients, mat@cell_metadata[names(mc@mc), 'PatientID'])
		groups = c(groups, c_col2group[mc@colors[mc@mc]])
		
		group2col = c(group2col, c_group2col)
	}
	df = data.frame(row.names=cells, f_genes=f_genes, patient=patients, group=groups, l_f_genes=log2(reg + f_genes))
	
	odir = sprintf("%s/f_umis_on_%s_genes", .scfigs_base, genes_nm)
	dir.create(odir, showWarnings = F)
	pats_ord = names(sort(tapply(df$f_genes, df$patient, summary_func)))
	
	m = matrix(NA, length(unique(df$patient)), length(unique(df$group)), dimnames=list(unique(df$patient), unique(df$group)))
	m_counts = m
	m_counts[is.na(m_counts)] = 0
	
	for (g in unique(df$group)) {
		c_df = df %>% filter(group == g)
		
		c_med = tapply(c_df$f_genes, c_df$patient, median)
		pat_n = table(c_df$patient)
		
		m_counts[names(pat_n), g] = pat_n
		m[names(c_med), g] = c_med
		
		o_pat = ordered(c_df$patient, levels=pats_ord)
		
		.plot_start(sprintf("%s/%s_%s_min_%d_cells.png", odir, genes_nm, g, min_cells), 700, 400)
		par(mar=c(10, 4, 3, 1))
		boxplot(l_f_genes ~ o_pat, c_df, notch=T, boxfill=ifelse(pat_n[pats_ord] >= min_cells, group2col[g], 'white'), border=ifelse(pat_n[pats_ord] >= min_cells, 'black', 'white'), main=paste(genes_nm, g), las=2, ylab='%umi (log2)', pch=19, cex=0.5)
		dev.off()
	}
	
	o_pat = ordered(df$patient, levels=pats_ord)
	pat_n = table(df$patient)
	.plot_start(sprintf("%s/%s_all_min_%d_cells.png", odir, genes_nm, min_cells), 700, 400)
	par(mar=c(10, 4, 3, 1))
	boxplot(l_f_genes ~ o_pat, df, notch=T, boxfill=ifelse(pat_n[pats_ord] >= min_cells, genes_col, 'white'), border=ifelse(pat_n[pats_ord] >= min_cells, 'black', 'white'), main=genes_nm, las=2, ylab='%umi (log2)', pch=19, cex=0.5)
	dev.off()
	
	.plot_start(sprintf("%s/%s_all_by_reac_min_%d_cells.png", odir, genes_nm, min_cells), 700, 400)
	par(mar=c(10, 4, 3, 1))
	boxplot(l_f_genes ~ o_pat, df, notch=T, boxfill=ifelse(pat_n[pats_ord] >= min_cells, pat_reac_col[pats_ord], 'white'), border=ifelse(pat_n[pats_ord] >= min_cells, 'black', 'white'), main=genes_nm, las=2, ylab='%umi (log2)', pch=19, cex=0.5)
	dev.off()
	
	m_n = log2( (reg + m) / (reg + mean(m, na.rm = T)))
	z_m_n = quantile(abs(m_n), 0.99, na.rm=T)
	.plot_start(sprintf("%s/%s_group_enr.png", odir, genes_nm), 700, 700)
	pheatmap(pmin(pmax(m_n[pats_ord, ], -z_m_n), z_m_n), breaks=seq(-z_m_n, z_m_n, len=100), cluster_rows=F, cellwidth=15, cellheight=15)
	dev.off()

	m_pat_n = log2( (reg + m) / (reg + rowMeans(m, na.rm = T)))
	z_m_pat_n = quantile(abs(m_pat_n), 0.99, na.rm=T)
	.plot_start(sprintf("%s/%s_group_pat_enr.png", odir, genes_nm), 700, 700)
	pheatmap(pmin(pmax(m_pat_n[pats_ord, ], -z_m_pat_n), z_m_pat_n), breaks=seq(-z_m_pat_n, z_m_pat_n, len=100), cluster_rows=F, cellwidth=15, cellheight=15)
	dev.off()
	
	tapply(df$l_f_genes, df$patient, summary_func)
}


#' set metacell params
#'
#' @param params named list (param-name --> value)
#'
#' @return list of previous params
#' @export
#'
override_metacell_params = function(params)
{
	prev_params = list()
	for (nm in names(params)) {
		prev_params[[nm]] = tgconfig::get_param(nm, "metacell")
		tgconfig::set_param(nm, params[[nm]], "metacell")
	}
	prev_params
}

#' Restore metacell params after calling to override_metacell_params
#'
#' @param prev_params the return value of override_metacell_params
#'
#' @export
#'
restore_metacell_params = function(prev_params)
{
	for (nm in names(prev_params)) {
		tgconfig::set_param(nm, prev_params[[nm]], "metacell")
	}
}


#' Compute Pielou's evenness 
#'
#' @param v vector of clone sizes
#'
#' @return
#' @export
#'
clonality_from_cln_counts_vec = function(v)
{
	v = v [ v > 0]
	v_n = v / sum(v)
	1 + sum(v_n * log(v_n))/log(length(v)) 
}


#' Boxplot-like plot of value per patient stratified by %cells in a reference group
#'
#' @param mc_id 
#' @param x 
#' @param nm 
#' @param col 
#' @param mat_id 
#' @param grp 
#' @param f_grp_strats 
#' @param odir 
#' @param ylab 
#' @param width 
#' @param height 
#'
#' @export
#'
#' @examples
.plot_stat_by_f_grp_helper = function(mc_id, x, nm, col, mat_id=mc_id, grp='dysfunctional', f_grp_strats=c(0.2, 0.4), odir=NULL, ylab='fraction', width=180, height=400) 
{
	mc = scdb_mc(mc_id)
	group2col = get_mc_group2col(mc)
	col2group = get_mc_col2group(mc)
	mat = scdb_mat(mat_id)
	
	pat_grp = table(mat@cell_metadata[names(mc@mc), 'PatientID'], col2group[mc@colors[mc@mc]])
	nms = intersect(rownames(pat_grp), names(x))
	x = x[nms]
	pat_grp = pat_grp[nms, ]
	pat_grp_n = pat_grp / rowSums(pat_grp)
	
	pat_f_grp = pat_grp_n[, grp]
	grp_pat_strat = ifelse(pat_f_grp < f_grp_strats[1], 'low', ifelse(pat_f_grp < f_grp_strats[2], 'mid', 'high'))
	o_grp_strat = ordered(grp_pat_strat, levels=c('low', 'mid', 'high'))
	
	.plot_start(scfigs_fn(mc_id, sprintf("f_%s_strat_by_f_%s", nm, grp), dir=odir), width, height)
	
	stripchart(x ~ o_grp_strat, method='jitter', vertical=T, pch=21, cex=1.5, las=2, main=sprintf("%s\n(%.2f)", nm, cor(x, pat_grp_n[, grp], method='spearman')), bg=col, ylab=ylab)
	segments(seq(0.75, by=1, len=3), tapply(x,  o_grp_strat, median), seq(1.25, by=1, len=3), tapply(x,  o_grp_strat, median), lwd=2)
	mw = wilcox.test(x[grp_pat_strat == 'low'], x[grp_pat_strat == 'high'])
	pv = mw$p.value
	mtext(ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", ""))), side=3, line=0, cex=1.5, at=2)
	dev.off()
}

#' Barplot like plot, bar per clone, showing clone group composition
#'
#' @param t_tcr produced by tcr_clones.r
#' @param cln_grp 
#' @param cln_size_range 
#' @param sort_by 
#' @param grp_ord 
#' @param group2col 
#' @param assert_on_dup_patients 
#' @param odir 
#'
#' @export
#'
plot_cln_group_composition = function(mc_id, t_tcr, cln_grp, cln_size_range, sort_by='hclust', grp_ord=names(ord_t_nk_nms), group2col, assert_on_dup_patients=T, odir=.scfigs_base) 
{
	cln_sz = rowSums(cln_grp)
	ind = cln_sz >= cln_size_range[1] & cln_sz <= cln_size_range[2]
	cln_grp = cln_grp[ind,]
	cln_sz = cln_sz[ind]
	
	cln2pat = unique(t_tcr[, c('clone_id', 'UnifPatientID')])
	if (max(table(cln2pat$clone_id)) > 1) {
		if (assert_on_dup_patients) {
			stop("found cross-patient clones!")
		} else {
			valid_clns = names(which(table(cln2pat$clone_id) == 1))
			cln2pat	= cln2pat[cln2pat$clone_id %in% valid_clns, ]
		}
	}
	rownames(cln2pat) = cln2pat$clone_id
	pats = cln2pat[as.character(rownames(cln_grp)), 'UnifPatientID']
	
	cln_grp_n = cln_grp / rowSums(cln_grp)
	if (sort_by == 'hclust') {
		hc = hclust(dist(cor(t(cln_grp_n))), method='ward.D2')
		cln_ord = hc$order
	} else if (sort_by == 'patient_type') {
		cln_ord = order(as.numeric(ordered(pats, levels=unique(pats))) + 1e-1 * cln_grp_n[, 'dysfunctional'] + 1e-2 * cln_grp_n[, 'effector1'] + 1e-3 * cln_grp_n[, 'effector2'] + 1e-4 * cln_grp_n[, 'Tfh'] + 1e-5 * cln_grp_n[, 'Treg'])
		
	}
	grp_ord = intersect(grp_ord, colnames(cln_grp))
	cln_grp = cln_grp[cln_ord, grp_ord]
	cln_sz = cln_sz[cln_ord]
	pats = pats[cln_ord]
	
	
	.plot_start(scfigs_fn(mc_id, sprintf("clone_group_comp_%d_to_%d_%s", cln_size_range[1], cln_size_range[2], sort_by), odir), 100 + 10 * nrow(cln_grp), 400)
	barplot(t(cln_grp), col=group2col[grp_ord], las=2, names.arg = pats)
	mtext(rownames(cln_grp), 3, at=seq(0.6, by=1.2, length=nrow(cln_grp)), las=2, line=1)
	dev.off()
	rownames(cln_grp)
}

#' Linear regression model with lasso regularization on input metacell based score from TFs. Uses glmnet package.
#'
#' @param mc_id 
#' @param dysf_sc First response vector 
#' @param cyto_sc Second response vector (predicting first - second)
#' @param all_genes_ann 
#' @param min_max_abs_lfp
#' @param lambda 
#' @param f 
#' @param min_cor_for_barplot 
#' @param nm1 
#' @param nm2 
#' @param grp1 
#' @param grp2 
#' @param nfolds 
#' @param min_vars 
#' @param max_genes_for_plot_size 
#' @param lateral_gset_id 
#' @param gen_plots 
#'
#' @return
#' @export
#'
#' @examples
model_by_tfs = function(mc_id, dysf_sc, cyto_sc, all_genes_ann, min_max_abs_lfp=1, lambda=NULL, f=NULL, min_cor_for_barplot=0.3, nm1='dysf', nm2='cyto', grp1='dysfunctional', grp2='effector2', nfolds=10, min_vars=3, max_genes_for_plot_size=-1, lateral_gset_id='mel_lateral', gen_plots=T) 
{
	mc = scdb_mc(mc_id)
	
	group2col = as.character(mc@color_key$color)
	names(group2col) = as.character(mc@color_key$group)
	
	lat_gset = scdb_gset(lateral_gset_id)
	lfp_t = log2(mc@mc_fp)
	
	if (is.null(f)) {
		f = rep(TRUE, ncol(lfp_t))
	}
	
	lfp_tfs = lfp_t[setdiff(rownames(all_genes_ann)[grepl('TF', all_genes_ann$type)], names(lat_gset@gene_set)), f]		
	lfp_tfs = lfp_tfs[apply(abs(lfp_tfs), 1, max) > min_max_abs_lfp, ]
	
	if (!is.null(cyto_sc)) {	
		cyto_tfs = sort(cor(cyto_sc[f], t(lfp_tfs))[1,])
		cyto_tfs = cyto_tfs[cyto_tfs > min_cor_for_barplot]
		if (gen_plots) {
			nbars = max(length(cyto_tfs), max_genes_for_plot_size)
			.plot_start(scfigs_fn(mc_id, paste0("fig2_", nm2, "_tfs_barplot")), 200, 100+nbars*12)
			par(mar=c(4,8,1,1))
			barplot(cyto_tfs, horiz=T, las=2, col=all_genes_ann[names(cyto_tfs), 'color'], ylim=c(0, nbars+2), xlab=paste('cor to', nm2))
			dev.off()
		}
	}
	
	dysf_tfs = sort(cor(dysf_sc[f], t(lfp_tfs))[1,])
	dysf_tfs = dysf_tfs[dysf_tfs > min_cor_for_barplot]
	if (gen_plots) {
		nbars = max(length(dysf_tfs), max_genes_for_plot_size)
		.plot_start(scfigs_fn(mc_id, paste0("fig2_", nm1, "_tfs_barplot")), 200, 100+nbars*12)
		par(mar=c(4,8,1,1))
		barplot(dysf_tfs, horiz=T, las=2, col=all_genes_ann[names(dysf_tfs), 'color'], ylim=c(0, nbars+2), xlab=paste('cor to', nm1))
		dev.off()
	}
	if (is.null(cyto_sc)) {	
		obs = dysf_sc[f]
	}
	else {
		obs = dysf_sc[f] - cyto_sc[f]
	}
	
	if (is.null(nfolds)) {
		nfolds = length(obs) # leave one out
		#message(paste("nfolds", nfolds))
	}
	cvfit = cv.glmnet(t(lfp_tfs), obs, nfolds=nfolds)
	
	if (is.null(lambda)) {
		lambda = cvfit$lambda[max(which(cvfit$lambda == cvfit$lambda.1se), min(which(cvfit$nzero >= min_vars)))]
		
		#message(sprintf("%.2f lambda", lambda))
	}
	
	if (gen_plots) {
		.plot_start(scfigs_fn(mc_id, paste0("fig2_", nm1, "_vs_", nm2, "_TF_model_CV")), 400, 400)
		plot(cvfit)
		abline(v=log(lambda))
		dev.off()
	}
	
	cvcoef = coef(cvfit, s=lambda)
	
	pos_feats = neg_feats = r2 = NULL
	if (sum(cvcoef != 0) > 1) {
		pred = predict(cvfit, newx=t(lfp_tfs), s=lambda)
		r2 = cor(obs, pred)**2
		if (gen_plots) {
			plt(sprintf("%s - %s (obs)", nm1, nm2), sprintf("%s - %s (model)", nm1, nm2), lfp_tfs, col=mc@colors[f], x=obs, y=pred, add_grid=T, cex.lab=1.5, cex=2, show_mc_ids=F, ofn=scfigs_fn(mc_id, paste0("fig2_", nm1, "_vs_", nm2, "_obs_vs_model")), main=sprintf("R^2 = %.4f", r2))
			
		}
		cvcoef = cvcoef[order(cvcoef[,1]),]
		cvcoef = cvcoef[cvcoef != 0]
		cvcoef = cvcoef[intersect(names(cvcoef), rownames(lfp_tfs))]
		if (gen_plots) {
			nbars = max(length(cvcoef), max_genes_for_plot_size)
			.plot_start(scfigs_fn(mc_id, paste0("fig2_", nm1, "_vs_", nm2, "_model_coef")), 100+nbars*12, 200)
			par(mar=c(8,4,1,1))
			barplot(cvcoef, las=2, col=ifelse(cvcoef > 0, group2col[grp1], group2col[ifelse(is.null(cyto_sc), 'white', grp2)]), ylab='coeffecient', xlim=c(0, nbars+2))
			dev.off()
		}
		#print(cvcoef)
		
		neg_feats = names(which(cvcoef < 0))
		if (length(neg_feats) > 0 & gen_plots) {
			for (i in seq_along(neg_feats)) { 
				plt(paste(nm2, "score"), neg_feats[i], lfp_t[, f], cols=mc@colors[f], x=if (is.null(cyto_sc)) { dysf_sc[f] } else { cyto_sc[f] }, show_mc_ids=F, add_grid=T, ofn=scfigs_fn(mc_id, sprintf("tfs_%s_grad_lfp_%d_%s", nm2, i, neg_feats[i]), scfigs_dir(mc_id, "fig2_grad_lfps")), cex.lab=1.5)
			}
		}
		
		pos_feats = names(which(cvcoef > 0))
		if (length(pos_feats) > 0 & gen_plots) {
			for (i in seq_along(pos_feats)) { 
				plt(paste(nm1, "score"), pos_feats[i], lfp_t[, f], cols=mc@colors[f], x=dysf_sc[f], show_mc_ids=F, add_grid=T, ofn=scfigs_fn(mc_id, sprintf("tfs_%s_grad_lfp_%d_%s", nm1, i, pos_feats[i]), scfigs_dir(mc_id, "fig2_grad_lfps")), cex.lab=1.5)
			}
		}
	}
	list(pos_tfs=pos_feats, neg_tfs=neg_feats, coef=cvcoef, r2=r2)
}

#' Wrapper func for plotting FACS index per patient and mc group
#'
#' @param fmd data frame with FACS indes (probably the @cell_metadata slot of the mat object)
#' @param mc_id 
#' @param patients 
#' @param ab_x 
#' @param ab_y 
#' @param groups 
#' @param fig_pref 
#'
#' @export
#'
plot_facs_by_patient_wrap = function(fmd, mc_id=tumor_clean_id, patients=c("p2-4-LN-2IT", "p4-4-LN-1IT", "p9-4-(S)C-2IT", "p11-3-(S)C-N", "p13-3-(S)C-N", "p27-4-(S)C-1IT"), ab_x="CD45", ab_y="CD3", groups=list("T"="T", NK="NK", B_plasma=c("B", "plasma"), myeloid=c("myeloid", "pDC", "osteoclast")), fig_pref="S1A") {
	mc = scdb_mc(mc_id)
	group2col = get_mc_group2col(mc)
	col2group = get_mc_col2group(mc)
	
	.plot_start(scfigs_fn(mc_id, sprintf("%s_FACS_%s_vs_%s", fig_pref, ab_y, ab_x)), 200 * length(groups), 200 * length(patients))
	layout(matrix(1:(length(groups)*length(patients)), nrow=length(patients), ncol=length(groups), byrow=T))
	par(mar=c(4,4,3,1))
	
	ab_x_c = paste0(ab_x, "_Ab")
	ab_y_c = paste0(ab_y, "_Ab")
	
	for (pat in patients) {
		c_fmd = fmd[fmd$PatientID == pat, ]
		for (i in seq_along(groups)) {
			g = groups[[i]]
			g_ind = rownames(c_fmd)[col2group[mc@colors[mc@mc[rownames(c_fmd)]]] %in% g]
			
			c_vals = .logicle_transform(c_fmd[, c(ab_x_c, ab_y_c)])
			g_ind = intersect(g_ind, rownames(c_vals))
			
			plot(c_vals[, ab_x_c], c_vals[, ab_y_c], pch=19, cex=0.5, col="grey90", xlab=paste(ab_x,'(logicle)'), ylab=paste(ab_y, '(logicle)'), main=paste(pat, names(groups)[i], 'sp', c_fmd[1, 'sp']))
			points(c_vals[g_ind, ab_x_c], c_vals[g_ind, ab_y_c], pch=19, cex=0.5, col=mc@colors[mc@mc[g_ind]])
			if (length(g_ind) >= 200) {
				d = MASS::kde2d(c_vals[g_ind, ab_x_c], c_vals[g_ind, ab_y_c], n=100)
				contour(d$x, d$y, d$z, drawlabels=F, add=T, col='grey30')
			}
		}
	}
	dev.off()
}

#' Compute genes geometric mean (and total umis) per patient for cells belong to a given mc group
#'
#' @param mc_id 
#' @param group 
#' @param min_cells 
#' @param mat_id 
#' @param min_pat_tot_umi 
#' @param min_pats_passing_tot_cutoff 
#'
#' @return
#' @export
#'
genes_geomean_by_patients_for_group = function(mc_id, group, min_cells=50, mat_id=mc_id, min_pat_tot_umi=10, min_pats_passing_tot_cutoff=4) 
{
	mc = scdb_mc(mc_id)
	mat = scdb_mat(mat_id)
	
	col2group = get_mc_col2group(mc)
	
	cells = names(mc@mc)[which(col2group[mc@colors[mc@mc]] == group)]
	pats = mat@cell_metadata[cells, 'PatientID']
	valid_pats = names(which(table(pats) >= min_cells))
	
	ind = pats %in% valid_pats
	
	cells = cells[ind]
	pats = pats[ind]
	
	tot = metacell:::.row_stats_by_factor(mat@mat[, cells], pats, Matrix::rowSums)
	
	f = rowSums(tot >= min_pat_tot_umi) >= min_pats_passing_tot_cutoff
	
	geom = metacell:::.row_stats_by_factor(mat@mat[f, cells], pats, function(y) {exp(Matrix::rowMeans(log(1+y)))-1})
	
	list(geom=geom, tot=tot)
}

#' Plot patient mc group composition with metadata (supp fig 4C-D)
#'
#' @param mc_id 
#' @param name 
#' @param filter_groups 
#' @param patients 
#' @param sort_pats_by 
#' @param min_cells 
#' @param mat_id 
#' @param md_nms 
#'
#' @return
#' @export
#'
#' @examples
plot_patients_group_comp = function(mc_id, name, filter_groups=NULL, patients=NULL, sort_pats_by=NULL, min_cells=0, mat_id=mc_id, md_nms = c('stage', 'location', 'treatment')) 
{
	mc = scdb_mc(mc_id)
	mat = scdb_mat(mat_id)
	col2group = get_mc_col2group(mc)
	group2col = get_mc_group2col(mc)
	
	pat_grp = table(mat@cell_metadata[names(mc@mc), 'PatientID'], col2group[mc@colors[mc@mc]])
	if (!is.null(filter_groups)) {
		pat_grp = pat_grp[, filter_groups]
	}
	pat_grp = pat_grp[, intersect(names(ord_by_id[[mc_id]]), colnames(pat_grp))]
	pat_grp = pat_grp[rowSums(pat_grp) >= min_cells, ]
	if (!is.null(patients)) {
		pat_grp = pat_grp[intersect(patients, rownames(pat_grp)), ]
	}
	pat_grp_n = pat_grp / rowSums(pat_grp)
	if (!is.null(sort_pats_by)) {
		if (sort_pats_by == 'hclust') {
			hc = hclust(dist(cor(t(pat_grp_n))), method='ward.D2')
			pat_grp_n = pat_grp_n[hc$order, ]
		} else if (sort_pats_by %in% colnames(pat_grp_n)) {
			pat_grp_n = pat_grp_n[order(pat_grp_n[, sort_pats_by]), ]
		} else {
			stop(paste("unknown sort by", sort_pats_by))
		}
	}
	
	f_inf_tab = read.table("data/processed_f_infiltrate.txt", header=T)
	
	pats_md = unique(mat@cell_metadata[names(mc@mc), c('PatientID', md_nms, 'Patient CD3+ % from CD45+')])
	rownames(pats_md) = pats_md$PatientID
	pats_md$f_inflitrate = f_inf_tab[ rownames(pats_md), 'f_inflitrate']
	pats_md[, 5] = as.numeric(pats_md[, 5])
	# manual correction, following update from the pathologist !
	pats_md['p3-4-LN-mIT', 'location'] = '(S)C'
	pats_md['p4-4-LN-1IT', 'location'] = '(S)C'
	.plot_start(scfigs_fn(mc_id, paste0("mc_comp_", name)), 300 + 25 * nrow(pat_grp_n), 560)
	layout(matrix(1:(length(md_nms)+4), ncol=1), heights=c(280, rep(20, length(md_nms)), 60, 60, 100))
	
	par(mar=c(0.5,10,3,1))
	barplot(t(pat_grp_n), col= group2col[colnames(pat_grp_n)], ylab='fraction', las=2, cex.names=1.5, cex.axis=1.5, cex.lab=1.5, xaxt='n')
	
	par(mar=c(0.2, 10, 0.2, 1))
	md_dict =  tgconfig::get_param("mcp_metadata_annot_colors", "metacell")
	
	for (nm in md_nms) {
		barplot(rep(1, nrow(pat_grp_n)), col=unlist(md_dict[[nm]][pats_md[rownames(pat_grp_n), nm]]), yaxt='n', border=NA)
		mtext(nm, 2, las=2)
	}
	par(mar=c(0.5, 10,0.5,1))
	barplot(pats_md[rownames(pat_grp_n), 5], col='royalblue', las=2, cex.axis=1.5, cex.lab=1.5, ylab='%CD3', border = NA)
	abline(h=0)
	barplot(ifelse(is.na(pats_md[rownames(pat_grp_n), 6]), max(pats_md[rownames(pat_grp_n), 6], na.rm=T), pats_md[rownames(pat_grp_n), 6]), col=ifelse(is.na(pats_md[rownames(pat_grp_n), 6]), 'grey90', 'royalblue4'), las=2, cex.axis=1.5, cex.lab=1.5, ylab='%Inflitrate', names.arg = rownames(pat_grp_n), cex.names=1.5, border=NA)
	abline(h=0)
	dev.off()

}

#' Compute umis fraction of a set of genes per given strata of cells, with chisq based p-values
#'
#' @param mat_id 
#' @param genes 
#' @param strats 
#' @param name 
#' @param col 
#'
#' @export
#'
f_genes_by_strats_with_chisq_pval = function(mat_id, genes, strats, name, col)
{
	mat = scdb_mat(mat_id)
	
	if (length(genes) > 1) {
		gu = tapply(Matrix::colSums(mat@mat[genes, names(strats)]), strats, sum)
	} else {
		gu = tapply(mat@mat[genes, names(strats)], strats, sum)
	}
	totu = tapply(Matrix::colSums(mat@mat[, names(strats)]), strats, sum)
	
	f_u = gu / totu
	
	.plot_start(scfigs_fn(mat_id, paste0(name, "f_genes_by_strats")), length(unique(strats))*30 + 100, 600)
	layout(matrix(1:3, ncol=1), heights=c(400, 70, 70))
	par(mar=c(4, 4, 4, 1))
	barplot(f_u, col=col, ylab='%UMIs', las=2, main=name)
	par(mar=c(0, 4, 0, 1))
	pv = sapply(2:length(gu), function(i) { x = chisq.test(matrix(c(gu[i-1], totu[i-1], gu[i],
																																	totu[i]), ncol=2)); x$p.value })
	stars = ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", "NS"))) 
	
	barplot(rep(1, length(f_u)), col=NA, border=NA, yaxt='n')
	text(seq(1.2, by=1.2, length=length(f_u)-1), rep(0.3, length(pv)-1), stars, cex=1.5, srt=0)
	barplot(rep(1, length(f_u)), col=NA, border=NA, yaxt='n')
	text(seq(1.2, by=1.2, length=length(f_u)-1), rep(0.5, length(pv)-1), sprintf("%.2e", pv), cex=1.5, srt=45)
	
	dev.off()
	
}

#' Boxplot of gene lfp value (per mc) over given mc groups
#'
#' @param mc_id 
#' @param gene 
#' @param groups 
#'
#' @export
#'
compare_lfp_on_gene = function(mc_id, gene, groups=NULL)
{
	mc = scdb_mc(mc_id)
	lfp = log2(mc@mc_fp[gene, ])
	col2group = get_mc_col2group(mc)
	group2col = get_mc_group2col(mc)
	
	if (is.null(groups)) {
		groups = ord_by_id[[mc_id]]		
	}
	groups = intersect(groups, unique(col2group[mc@colors]))
	df = data.frame(lfp=lfp, group = col2group[mc@colors])
	grp_o = ordered(df$group, levels=groups)
	.plot_start(scfigs_fn(mc_id, sprintf("comp_%s_on_%s", gene, paste0(groups, collapse="_"))),
							200 + length(groups) * 40, 300)
	par(mar=c(8, 4, 4, 2))
	boxplot(lfp ~ grp_o, df, las=2, main=gene, ylab="fp (log2)", pch=30)
	stripchart(lfp ~ grp_o, df, vertical=T, method="jitter", pch=19, cex=0.7, add=T, col=group2col[groups], jitter=0.3)
	dev.off()
}

# wrap for opening a plot (png, ps or pdf)
.plot_start = function(fn, w, h)
	
{
	device = get_param("mc_plot_device", "metacell")
	res = get_param("mc_plot_ppi", "metacell")
	pointsize = get_param("mc_plot_pointsize", "metacell")
	
	if (device == "png") {
		png(filename=sub("ps$", "png", fn), width=w, height=h, res=res, pointsize = pointsize)
	}
	else if (device == "ps") {
		postscript(file=sub("png$", "ps", fn), width=w/res, height=h/res)
	}
	else if (device == "pdf") {
		pdf(file=sub("png$", "pdf", fn), width=w/res, height=h/res)
	}
	else if (device == "svg") {
		RSvgDevice::devSVG(file=sub("png$", "svg", fn), width=w/res, height=h/res)
	}
	
	else {
		stop(sprintf("unknown output device type: %s", device))
	}
}
