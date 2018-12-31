source("util_funcs.r")

# download data files
if (!file.exists("data")) {
	download.file("http://www.wisdom.weizmann.ac.il/~lubling/Li2018/data_Li2018.tar.gz", "data_Li2018.tar.gz")
	
	system("tar xfz data_Li2018.tar.gz")
}

if (!file.exists("data/GSE99254_NSCLC.TCell.S12346.TPM.txt")) {
	download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE99254&format=file&file=GSE99254%5FNSCLC%2ETCell%2ES12346%2ETPM%2Etxt%2Egz", "data/GSE99254_NSCLC.TCell.S12346.TPM.txt.gz")
	system("gunzip data/GSE99254_NSCLC.TCell.S12346.TPM.txt.gz")
}

if (!file.exists("scrna_db_Li2018")) {
	download.file("http://www.wisdom.weizmann.ac.il/~lubling/Li2018/scrna_db_Li2018.tar.gz", "scrna_db_Li2018.tar.gz")
	
	system("tar xfz scrna_db_Li2018.tar.gz")
}

# init - load the metacell package and define global variables
dir.create("scrna_db_Li2018", showWarnings = F)
dir.create("figs_Li2018", showWarnings = F)
rl(scdb_dir="scrna_db_Li2018", scfigs_dir="figs_Li2018")

def_vars()

#
# Load all MARS tables into a master mat. Fix mixed plates metadata info.
mel_build_master_mat = function(correct_mixed_plates=T, force_new=F)
{
	# load all melanoma MARS plates (tumor, pbmc, tumor_blocks)
	mcell_import_multi_mars(mat_nm=all_id, dataset_table_fn="data/MelanomaSampleIndex_all_valid_210718.txt", base_dir="data/umi.tab", force=force_new)
	message("import mat done")
	
	full_m = scdb_mat(all_id)
	
	mcell_plot_batch_stats(all_id)
	min_umis_cutoff = mcell_plot_umis_per_cell(all_id)
	message("batch stats done")
	
	# some plates have cells with different treatment. Assign the correct treatment on the cell level, and remove empty wells
	message("correcting treatment, PatientID and Short Sample Name for subsets of cells")

	w2c = read.table("data/wells_cells.txt", header=T, stringsAsFactors=T, sep="\t")
	transl_tab = read.table("data/mixed_plates_treatment_info.txt", header=T)
	
	cell2treat = reshape2::melt(as.matrix(transl_tab), value.name="treatment") %>% 
		left_join(w2c, by=c("Var1" = "well_coordinates", "Var2" = "Amp_batch_ID")) %>% 
		select(Well_ID, treatment)
	
	empty_cells = as.character(cell2treat[cell2treat$treatment == 'Empty', 'Well_ID'])
	
	cells = as.character(cell2treat[cell2treat$treatment != 'Empty', 'Well_ID'])
	treat = as.character(cell2treat[cell2treat$treatment != 'Empty', 'treatment'])
	
	mdc = full_m@cell_metadata[cells, ]
	mdc[, 'treatment'] = treat
	mdc[, 'PatientID'] = paste(mdc[, 'new PID'], mdc$stage, mdc$location, mdc$treatment, sep="-")
	mdc$Patient = ifelse(is.na(mdc$n_tumor), mdc$PatientID, paste(mdc$PatientID, mdc$n_tumor, sep="-"))
	
	mdc[, 'Sample Short Name']  = paste0(mdc$PatientID, ifelse(mdc$CD3 == 1, "_CD3", ""), ifelse(mdc$CD45 == 1, "_CD45", ""), ifelse(mdc$live == 1, '_live', ""), "_", mdc$Processing)
	
	full_m@cell_metadata[cells, ] = mdc
	
	scdb_add_mat(all_id, full_m)
	
	mcell_mat_ignore_cells(all_id, all_id, empty_cells)
	
	mcell_add_gene_stat(all_id, all_id, force=force_new)
	}

#
# Build blaclisted gene sets
mel_build_blacklist_gsets_by_gene_nms = function() 
{
	full_m = scdb_mat(all_id)
	
	# mitochondrial gene set
	mt_cands = grep("^MT-|^MTRN|^MTAT|^MTND|^MRP", full_m@genes, v=T, perl=T)
	mito = data.table::fread("data/MitoCarta2_human.txt", header=T, sep="\t", stringsAsFactors=F)
	mt_both = intersect(mt_cands, mito$Symbol)
	mt_cands = setdiff(mt_cands, mt_both)
	mitocarta = setdiff(mito$Symbol, mt_both)
	mt_genes = c(rep('regexp MT', length(mt_cands)), rep('MitoCarta2', length(mitocarta)), rep('MitoCarta2 and regexp MT', length(mt_both)))
	names(mt_genes) = c(mt_cands, mitocarta, mt_both)
	scdb_add_gset(mito_gset_id, gset_new_gset(mt_genes, 'mitochondrial genes'))
	
	# IG gene set
	ig_nms = grep("^IGK|^IGL|^IGJ|^IGH|^IGBP|^IGSF", full_m@genes, v=T, perl=T)
	ig_genes = rep("IG", length(ig_nms))
	names(ig_genes) = ig_nms
	scdb_add_gset(ig_gset_id, gset_new_gset(ig_genes, 'IG genes'))
	
	# ncRNA gene set
	ncrna_nms = c('MALAT1', 'XIST', 'NEAT1', 'hsa-mir-6723')
	ncrna_genes = rep("ncRNA", length(ncrna_nms))
	names(ncrna_genes) = ncrna_nms
	scdb_add_gset(ncrna_gset_id, gset_new_gset(ncrna_genes, "ncRNA genes"))
	
	# RP pseudo genes set
	ncrp_nms = grep("^RP[0-9]+-", full_m@genes, v=T, perl=T)
	ncrp_genes = rep("ncRP", length(ncrp_nms))
	names(ncrp_genes) = ncrp_nms
	scdb_add_gset(ncrp_gset_id, gset_new_gset(ncrp_genes, "RP##- genes"))
	
}

#
# build clean mat: remove blaclisted genes (mitochondrial, ncRNA, RP[0-9], IG) and small cells
mel_build_blist_filtered_master_mat = function(min_umis_post_gene_ignore = 500, max_mito_f = 0.6, force_new=T)
{
	full_m = scdb_mat(all_id)
	
	blist_gsets = sapply(c(mito_gset_id, ig_gset_id, ncrna_gset_id, ncrp_gset_id), scdb_gset)
	blist_genes = unlist(lapply(blist_gsets, function(gs) names(gs@gene_set)))

	uc = Matrix::colSums(full_m@mat)
	
	mt_genes = intersect(names(blist_gsets[[mito_gset_id]]@gene_set), full_m@genes)
	mito_f = Matrix::colSums(full_m@mat[mt_genes, ]) / uc
	
	# filter cells with low counts, large MT fraction, ignore MT genes, RP##- and mega-strong RNA genes
	mcell_mat_ignore_genes(filt_id, all_id, blist_genes)
	
	filt_mat = scdb_mat(filt_id)
	
	mcell_mat_ignore_cells(filt_id, filt_id, union(filt_mat@ignore_cells, names(uc)[mito_f >= max_mito_f | Matrix::colSums(filt_mat@mat) <= min_umis_post_gene_ignore]))
	
	# working on the filtered mat
	mcell_add_gene_stat(filt_id, filt_id, force=force_new)
}

#
# metacell generation wrapper on selected batch set values
mel_basic_mat2mc = function(mat_id, batch_sets, lateral_gset_id, force_new=T, 
														T_vm = 0.2, T_tot = 200, T_top3 = 3,
														cgraph_knn = 200, cgraph_downsamp=T, 
														bootstrap_n_resamp=500, bootstrap_p_resamp=0.75,
														mc_K=30, min_mc_size=30, mc_alpha=2, cells=NULL, name="") 
{
	mat = scdb_mat(mat_id)
	stopifnot(!is.null(mat))
	
	# create new mat
	if (is.null(cells)) {
		bs_id = paste0(batch_sets, collapse="_")
		new_mat_id = sprintf("%s_%s%s", mat_id, bs_id, name)
		
		message(sprintf("will submat %s into %s", mat_id, new_mat_id))
		
		s_mat = scdb_mat(new_mat_id)
		if (is.null(s_mat) | force_new) {
			mcell_mat_ignore_cells(new_mat_id, mat_id, mat@cells[ mat@cell_metadata[mat@cells, 'batch_set_id'] %in% batch_sets], reverse=T)
		
			mcell_add_gene_stat(new_mat_id, new_mat_id, force=force_new)
		}
	} else {
		new_mat_id = sprintf("%s_%s", mat_id, name)
		s_mat = scdb_mat(new_mat_id)
		if (is.null(s_mat) | force_new) {
			mcell_mat_ignore_cells(new_mat_id, mat_id, cells, reverse=T)
			mcell_add_gene_stat(new_mat_id, new_mat_id, force=force_new)
		}
	}
	# select genes to affect graph creation
	mcell_gset_filter_varmean(new_mat_id, new_mat_id, T_vm=T_vm, force_new=T)
	mcell_gset_filter_cov(new_mat_id, new_mat_id, T_tot=T_tot, T_top3=T_top3)
	
	if (!is.null(lateral_gset_id)) {
		marker_gset = scdb_gset(new_mat_id)	
		lateral_gset = scdb_gset(lateral_gset_id)
		marker_gset = gset_new_restrict_gset(marker_gset, lateral_gset, inverse=T, "cgraph markers w/o lat genes")
		scdb_add_gset(new_mat_id, marker_gset)
	}
	
	# create cgraph
	mcell_add_cgraph_from_mat_bknn(new_mat_id, new_mat_id, new_mat_id, K=cgraph_knn, dsamp=cgraph_downsamp)
	message(new_mat_id, " cgraph done")
	
	# bootstrap
	mcell_coclust_from_graph_resamp(new_mat_id, new_mat_id,
																	min_mc_size=round(cgraph_knn/5),
																	p_resamp=bootstrap_p_resamp,
																	n_resamp=bootstrap_n_resamp)
	message(new_mat_id, " bootstrap done")
	
	# mc from coclust matrix
	mcell_mc_from_coclust_balanced(new_mat_id, new_mat_id, new_mat_id, K=mc_K, min_mc_size=min_mc_size, alpha=mc_alpha)
	
	# find and remove outliers
	T_lfc = 3000 # impossible cutoff - to make it only split heterogenous mcs
	
	# breaking down heterogenous metacells and removing outliers by extreme gene expression
	mcell_mc_split_filt(new_mc_id=sprintf("%s_outClean", new_mat_id), new_mat_id, new_mat_id, T_lfc=T_lfc, plot_mats=F)
	
}

#
#
mel_basic_mc_mc2d_plots = function(mc_id, mc2d_id=mc_id, graph_id=mc_id, mat_id=mc_id, lateral_gset_id=NULL, metadata_fields_to_export=c('PatientID', 'stage', 'location', 'treatment'), mc_ord=NULL, plot_2d=T)
{

	if (is.null(lateral_gset_id)) {
		mcell_gset_from_mc_markers(gset_id=mc_id, mc_id=mc_id)
		
		mcell_mc_plot_marks(mc_id, mc_id, mat_id, plot_cells=T)
		mcell_mc_plot_marks(mc_id, mc_id, mat_id, plot_cells=F)
	} else {
		mcell_gset_from_mc_markers(gset_id=mc_id, mc_id=mc_id, blacklist_gset_id=lateral_gset_id)
		mcell_gset_from_mc_markers(gset_id=paste0(mc_id, "_lateral"), mc_id=mc_id, filt_gset_id=lateral_gset_id)
	
	
		mcell_mc_plot_marks(mc_id, mc_id, mat_id, lateral_gset_id=paste0(mc_id, "_lateral"), plot_cells=T, mc_ord=mc_ord)
		mcell_mc_plot_marks(mc_id, mc_id, mat_id, lateral_gset_id=paste0(mc_id, "_lateral"), plot_cells=F, mc_ord=mc_ord)
	}
	
	if (plot_2d) {
		mcell_mc2d_force_knn(mc2d_id, mc_id, graph_id)
		mcell_mc2d_plot(mc2d_id=mc2d_id)	
	}
	
	mcell_mc_export_tab(mc_id=mc_id, gstat_id=mat_id, mat_id=mat_id, T_fold=2, metadata_fields=metadata_fields_to_export)
}

#
# Proliferation analysis for the given metacell object (mc_id)
mel_mc_prolif_score = function(mc_id=tumor_t_nk_id, mat_id=mc_id, mc2d_id=mc_id, f_cc_cutoff=0.0325, plot_mode='circle_size', cc_gset_id="mel_cc_filt", add_genes_by_names=F, pats_order=NULL, strat_by=NULL, patient_fname='PatientID')
{
	mc = scdb_mc(mc_id)
	mat = scdb_mat(mat_id)
	mc2d = scdb_mc2d(mc2d_id)
	
	col2group = get_mc_col2group(mc)
	group2col = get_mc_group2col(mc)
	
	height = tgconfig::get_param("mcell_mc2d_gene_height", "metacell")
	width = tgconfig::get_param("mcell_mc2d_gene_width", "metacell")
	mc_cex = tgconfig::get_param("mcell_mc2d_gene_mc_cex", "metacell")
	sc_cex = tgconfig::get_param("mcell_mc2d_gene_cell_cex", "metacell")

	cells = names(mc@mc)
	
	cc_gset = scdb_gset(cc_gset_id)
	cc_genes = names(cc_gset@gene_set)
	if (add_genes_by_names) {
		cc_genes = union(cc_genes, grep('^HIST|^CENP|^SMC[0-9]', mat@genes, v=T, perl=T))
	}
	diff_genes = setdiff(cc_genes, mat@genes)
	if (length(diff_genes) > 0) {
		message(sprintf("warning: missing %d genes (%s)", length(diff_genes), paste0(diff_genes, collapse=", ")))
		cc_genes = setdiff(cc_genes, diff_genes)
	}
	tot_cc_c = Matrix::colSums(mat@mat[cc_genes, cells])
	f_cc_c =  tot_cc_c / Matrix::colSums(mat@mat[, cells])
	tot_cc = tapply(f_cc_c > f_cc_cutoff, mc@mc[cells], sum)
	tara = tapply(f_cc_c, mc@mc[cells], length) - tot_cc
	
	df = data.frame(row.names=names(tot_cc), tot_cc=tot_cc, tara=tara, f_cc=tot_cc/(tot_cc+tara))
	df = df[order(df$tot_cc / (df$tot_cc + df$tara)), ]
	
	.plot_start(scfigs_fn(mc2d_id, sprintf("prolif_frac_by_%s%s_%s", cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""), plot_mode)), width, height)
	plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=sc_cex, col='lightgrey', xlab="", ylab="")
	if (plot_mode == 'pie') {
		for (i in 1:nrow(df)) {
			c_mc = as.numeric(rownames(df)[i])
			plotrix::floating.pie(mc2d@mc_x[c_mc], mc2d@mc_y[c_mc], c(df[i, 'tot_cc'], df[i, 'tara'])+0.1, radius=diff(range(mc2d@sc_x, na.rm=T)) / 50, col=c(mc@colors[c_mc], 'white'))
		}
	}
	else if (plot_mode == "circle_size") {
		mcs = as.numeric(rownames(df))
		points(mc2d@mc_x[mcs], mc2d@mc_y[mcs], pch=21, bg=mc@colors[mcs], cex=mc_cex * log10(1 + 100*df$f_cc)/2)
		fs = seq(0.05, length=1 + ceiling(max(df$f_cc)/0.15), by=0.15)

		legend("topleft", legend=paste0(fs * 100, "%"), pch=21, pt.bg='darkgrey', col='black', pt.cex=mc_cex * log10(1 + 100*fs)/2, bty='n', y.intersp=2, inset=0.04)
	}
	#text(mc2d@mc_x, mc2d@mc_y, 1:length(mc2d@mc_x), cex=0.7)
	dev.off()
	
	df_cc = data.frame(row.names=names(f_cc_c), f_cc=f_cc_c, tot_cc=tot_cc_c, patient=mat@cell_metadata[names(f_cc_c), patient_fname], group=col2group[mc@colors[mc@mc[names(f_cc_c)]]], stringsAsFactors = F)
	
	# %prolif per group
	.plot_start(scfigs_fn(mc_id, sprintf("prolif_frac_by_%s%s_group", cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 450, 250)
	f_cc_by_grp = tapply(df_cc$f_cc >= f_cc_cutoff, df_cc$group, mean)
	tot_cc_by_grp = tapply(df_cc$f_cc >= f_cc_cutoff, df_cc$group, length)
	f_cc_by_grp = f_cc_by_grp[names(ord_by_id[[mc_id]])]
	tot_cc_by_grp = tot_cc_by_grp[names(f_cc_by_grp)]

	par(mar=c(8,4,4,1))
	barplot(f_cc_by_grp, col=group2col[names(f_cc_by_grp)], las=2, ylab='% prolif cells')
	mtext(tot_cc_by_grp, 3, at=seq(0.6, by=1.2, length=length(tot_cc_by_grp)), las=2, line=0.5)
	dev.off()
	
	# %prolif per patient
	.plot_start(scfigs_fn(mc_id, sprintf("prolif_frac_by_%s%s_patient", cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 450, 250)
	f_cc_by_pat = tapply(df_cc$f_cc >= f_cc_cutoff, df_cc$patient, mean)
	tot_cc_by_pat = tapply(df_cc$f_cc >= f_cc_cutoff, df_cc$patient, length)
	if (is.null(pats_order)) {
		f_cc_by_pat = sort(f_cc_by_pat)
		pats_order = names(f_cc_by_pat)
	} else {
		f_cc_by_pat = f_cc_by_pat[pats_order]
	}
	tot_cc_by_pat = tot_cc_by_pat[names(f_cc_by_pat)]
	par(mar=c(8,4,4,1))
	barplot(f_cc_by_pat, col="lightblue", las=2, ylab='% prolif cells')
	mtext(tot_cc_by_pat, 3, at=seq(0.6, by=1.2, length=length(tot_cc_by_pat)), las=2, line=0.5)
	dev.off()
	
	# qc - cc gene contrib by grp on prolif cells
	df_cc_pos = df_cc[df_cc$f_cc >= f_cc_cutoff, ]
	pos_cells = rownames(df_cc_pos)
	cc_g_contrib = metacell:::.row_stats_by_factor(mat@mat[cc_genes, pos_cells], df_cc_pos$group, Matrix::rowSums)
	cc_tot_umi = tapply(Matrix::colSums(mat@mat[cc_genes, pos_cells]), df_cc_pos$group, sum)
	cc_g_contrib_n = t(t(cc_g_contrib) / as.numeric(cc_tot_umi))
	cc_g_contrib_n_ord = cc_g_contrib_n[order(rowMeans(cc_g_contrib_n)), ]
	
	.plot_start(scfigs_fn(mc_id, sprintf("cc_gene_contrib_by_group_%s%s", cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 900, 400)
	par(mar=c(8,4,1,1))
	boxplot(t(log2(1e-3+cc_g_contrib_n_ord)), las=2, pch=NA, cex=0.7, ylab="mean umi over cc cells (log2)")
	for (i in 1:ncol(cc_g_contrib_n_ord)) {
		points(1:nrow(cc_g_contrib_n_ord) + (i - ncol(cc_g_contrib_n_ord)/2)*0.04, log2(1e-3+cc_g_contrib_n_ord[, i]), pch=19, cex=0.7, col=group2col[colnames(cc_g_contrib_n_ord)][i])
	}
	dev.off()
	
	# qc - per gene, umi lfp per group/patient
	cc_cells = rownames(df_cc_pos)
	all_m = mat@mat[cc_genes, cells]
	cc_m = mat@mat[cc_genes, cc_cells]
	
	.plot_m_fp = function(m, grps, nm, grps_ord, fp_reg=0.1) {
		grp_geomean = metacell:::.row_stats_by_factor(m, grps, function(y) {exp(Matrix::rowMeans(log(1+y)))-1})	
		grp_meansize = tapply(Matrix::colSums(m), grps, mean)
		g_fp = t(1000*t(grp_geomean)/as.vector(grp_meansize))
		lfp = log2((fp_reg+g_fp)/apply(fp_reg+g_fp, 1, median))
		ghc = hclust(dist(cor(t(lfp))), method='ward.D2')
		.plot_start(scfigs_fn(mc_id, sprintf("cc_gene_lfp_by_%s_%s%s", nm, cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 500, 500)
		pheatmap(pmin(pmax(lfp[ghc$order, intersect(grps_ord, colnames(lfp))], -3), 3), breaks=seq(-3, 3, len=100), cluster_rows=F, cluster_cols=F, main=nm)
		dev.off()
	}
	.plot_m_fp(all_m, df_cc$group, "group_all_cells", grps_ord=names(ord_by_id[[mc_id]]))
	.plot_m_fp(cc_m, df_cc_pos$group, "group_prolif_cells", grps_ord = names(ord_by_id[[mc_id]]))
	.plot_m_fp(all_m, df_cc$patient, "patient_all_cells", grps_ord = pats_order)
	.plot_m_fp(cc_m, df_cc_pos$patient, "patient_prolif_cells", grps_ord = pats_order)
	
																			 
	# prolif strat on %dysf %treg umis
	if (!is.null(strat_by)) {
		for (nm in names(strat_by)) {
			strats = strat_by[[nm]]
			.plot_start(scfigs_fn(mc_id, sprintf("prolif_f_umis_strat_by_%s_%s%s", nm, cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 300, 300)
			par(mar=c(8,4,4,1))
			boxplot(log2(1e-3+df_cc[names(strats), 'f_cc']) ~ strats, notch=T, main=paste("cc umis by", nm), las=2, pch=19)
			dev.off()
			
			.plot_start(scfigs_fn(mc_id, sprintf("f_prolif_strat_by_%s_%s%s", nm, cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 300, 300)
			par(mar=c(8,4,4,1))
			y = tapply(df_cc[names(strats), 'f_cc'] >= f_cc_cutoff, strats, mean)
			barplot(y, main=paste("%prolif by", nm), las=2, col=group2col[nm])
			dev.off()

			.plot_start(scfigs_fn(mc_id, sprintf("f_prolif_strat_by_%s_per_type_%s%s", nm, cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 300, 300)
			par(mar=c(8,4,4,1))
			grps = df_cc[names(strats), 'group']
			c_f_cc = df_cc[names(strats), 'f_cc']
			y = matrix(NA, length(unique(grps)), length(unique(strats)), dimnames=list(unique(grps), sort(unique(strats))))
			for (grp in rownames(y)) {
				tt = tapply(c_f_cc[grps == grp] >= f_cc_cutoff, strats[grps == grp], mean)
				y[grp, names(tt)] = tt
			}
			n_y = table(grps, strats)
			y[ n_y[rownames(y), colnames(y)] < 250] = 0
			barplot(y, main=paste("%prolif by", nm), las=2, col=group2col[rownames(y)], beside=T)
			dev.off()

			.plot_start(scfigs_fn(mc_id, sprintf("tot_prolif_strat_by_%s_%s%s", nm, cc_gset_id, ifelse(add_genes_by_names, "_add_by_name", ""))), 300, 300)
			par(mar=c(8,4,4,1))
			y = tapply(df_cc[names(strats), 'tot_cc'], strats, mean)
			barplot(y, main=paste("Mean prolif UMIs by", nm), las=2, col=group2col[nm])
			dev.off()
		}
	}
	
	list(prolif_df=df, f_cc=f_cc_c, df_cc=df_cc, f_cc_by_grp=f_cc_by_grp, tot_cc_by_grp=tot_cc_by_grp)
}

#
# gsets
build_lateral_gene_sets = function() 
{

	# The 3 gene sets for cell-cycle, IFN and stress genes were build using the following commands. Since the process involves random downsampling of the umi matrix the results can be slightly different, so for consistency we're supplying the gene sets used in the paper under the scrna_db_Li2018 folder.
	#	select_gene_modules_by_anchor_genes(filt_id, c('MKI67', 'HIST1H1D', 'PCNA', 'SMC4', 'MCM3'), gset_nm="mel_cc", cor_thresh=0.1, sz_cor_thresh=0.1)
	#select_gene_modules_by_anchor_genes(filt_id, c('ISG15', 'OAS1', 'WARS', 'IFIT1'), gset_nm="mel_ifn", cor_thresh=0.1, sz_cor_thresh=0.1)
	#select_gene_modules_by_anchor_genes(filt_id, c("TXN", "HSP90AB1", "HSPA1A", "FOS", "HIF1A"), gset_nm="mel_stress", cor_thresh=0.1, sz_cor_thresh=0.1, nclusts=24)
	
	# Stopped at this point to manually select the following gene clusters to filter
	mcell_gset_remove_clusts("mel_cc", filt_clusts=c(4, 7, 11, 12, 14, 16:18), new_id = "mel_cc_filt", reverse=T)
	mcell_gset_remove_clusts("mel_ifn", filt_clusts=c(5, 14:17, 19), new_id = "mel_ifn_filt", reverse=T)
	mcell_gset_remove_clusts("mel_stress", filt_clusts=c(9, 10, 16, 17, 18, 22), new_id = "mel_stress_filt", reverse=T)
	
	lat_gsets = lapply(c('mel_cc_filt', 'mel_ifn_filt', 'mel_stress_filt'), scdb_gset)
	new_lat_genes = unlist(lapply(1:length(lat_gsets), function(i) { gs = lat_gsets[[i]]@gene_set; v = rep(i, length(gs)); names(v) = names(gs); return(v) }))
	
	# manually add some genes to make sure they're in
	mat_f = scdb_mat(filt_id)
	add_cc_genes = grep('^HIST|^CENP|^SMC[0-9]', mat_f@genes, v=T, perl=T)
	add_ifn_genes = grep('^IFI', mat_f@genes, v=T, perl=T)
	add_genes = rep(1:2, times=c(length(add_cc_genes), length(add_ifn_genes)))
	names(add_genes) = c(add_cc_genes, add_ifn_genes)
	add_genes = add_genes[setdiff(names(add_genes), names(new_lat_genes))]
	scdb_add_gset(lateral_gset_id, gset_new_gset(c(new_lat_genes, add_genes), 'lateral: CC, IFN, stress'))
	
}

#
# Process the umi table (fiter cells and genes) and generate the global metacell object for all Tumor samples. Due to varying random seeding, the exact partitioning of cells to metacells may differ upon each run. For consistency, the downstream analysis is using the metacell objects used by us to generate the figures. The current function plots a metacell-to-metacell cell count matrix to demonstrate that the differences between the metacell objects are marginal. 
build_metacells = function(rebuild_data=F) 
{
	metacell:::.set_seed(42)
	mel_build_master_mat(force_new=rebuild_data)
	mel_build_blacklist_gsets_by_gene_nms()
	mel_build_blist_filtered_master_mat(force_new = rebuild_data)
	build_lateral_gene_sets()
	
	# build mc from all Tumor samples
	mel_basic_mat2mc(filt_id, 'Tumor', lateral_gset_id, name="_new")
	
	
	# Tumor samples with matching PBMC samples
	# 
	# Note!!! The published version used all cells of patients p13, p17 and p27 with their matching PBMC samples. Current version filters out cells that didn't pass QC (mainly small cells). Results are similar yet cleaner.
	mat_filt = scdb_mat(filt_id)
	pbmc_pats = grep("p13-|p17-|p27-|PBMC", unique(mat_filt@cell_metadata[mat_filt@cells, 'PatientID']), v=T, perl=T)

	pbmc_cells = rownames(mat_filt@cell_metadata[mat_filt@cells, ])[mat_filt@cell_metadata[mat_filt@cells, 'PatientID'] %in% pbmc_pats]
	
	mel_basic_mat2mc(filt_id, cells=pbmc_cells, name="tumors_pbmc_comp", cgraph_knn=100, lateral_gset_id = lateral_gset_id)
	
}

#
# generate a heatmap (and a table) cell membership between the 2 given metacells, showing annotations of each mc if they exist
compare_mc_partitioning = function(mc1_id, mc2_id, do_log=F, ord1='by_col', ord2='by_whichMax') 
{
	mc1 = scdb_mc(mc1_id)
	mc2 = scdb_mc(mc2_id)
	
	ann_cols = NA
	mc1_ann = NA
	g1_ann = NA
	
	if (length(mc1@colors) == ncol(mc1@mc_fp)) {
		col2group1 = get_mc_col2group(mc1)
		group2col1 = get_mc_group2col(mc1)
		mc1_ann = data.frame(row.names=seq_along(mc1@colors), group1=col2group1[mc1@colors], stringsAsFactors = F)
		g1_ann = data.frame(row.names=unique(mc1_ann$group1), group1=unique(mc1_ann$group1), stringsAsFactors = F)
		ann1 = group2col1[unique(mc1_ann$group1)]
		ann_cols = list(group1=ann1)
	}
	
	mc2_ann = NA
	g2_ann = NA
	if (length(mc2@colors) == ncol(mc2@mc_fp)) {
		col2group2 = get_mc_col2group(mc2)
		group2col2 = get_mc_group2col(mc2)
		mc2_ann = data.frame(row.names=seq_along(mc2@colors), group2=col2group2[mc2@colors], stringsAsFactors = F)
		g2_ann = data.frame(row.names=unique(mc2_ann$group2), group2=unique(mc2_ann$group2), stringsAsFactors = F)
		ann2 = group2col2[unique(mc2_ann$group2)]
		ann_cols[['group2']] = ann2
	}
	
	cells = intersect(names(mc1@mc), names(mc2@mc))
	message(sprintf("Comparing %s (%d cells) with %s (%d cells), %d common cells", mc1_id, length(mc1@mc), mc2_id, length(mc2@mc), length(cells)))
	
	tab = table(mc1@mc[cells], mc2@mc[cells])

	tab_n = tab / rowSums(tab)
	tab_n2 = t(t(tab) / colSums(tab))

	rord = 1:nrow(tab_n)
	cord = 1:ncol(tab_n)

	gaps1 = gaps2 = NULL
	if (ord1 == "by_col" && !is.null(nrow(mc1_ann))) {
		fac1 = factor(mc1_ann$group1, levels=names(sort(tapply(apply(tab_n, 1, which.max), mc1_ann$group1, mean))))
		rord = order(as.numeric(fac1) + rord * 1e-6)
		ann_cols[['group1']] = ann_cols[['group1']][levels(fac1)]
		gaps1 = which(diff(as.numeric(fac1[rord])) > 0)
	} else if (ord1 == "by_whichMax") {
		rord = order(apply(tab_n, 1, which.max))
	}
	tab_n = tab_n[rord, ]	
	tab_n2 = tab_n2[rord, ]
	
	if (ord2 == "by_col" && !is.null(nrow(mc2_ann))) {
		fac2 = factor(mc2_ann$group2, levels=names(sort(tapply(apply(tab_n2, 2, which.max), mc2_ann$group2, mean))))
		cord = order(as.numeric(fac2) + cord * 1e-6)
		ann_cols[['group2']] = ann_cols[['group2']][levels(fac2)]
		gaps2 = which(diff(as.numeric(fac2[cord])) > 0)
	} else if (ord2 == "by_whichMax") {
		cord = order(apply(t(t(tab) / colSums(tab)), 2, which.max))
	}
	tab_n = tab_n[, cord]
	tab_n2 = tab_n2[, cord]
	
	ofn_base = sprintf("memb_comp_with_%s_%s_%s_%s", mc2_id, ord1, ord2, ifelse(do_log, "log", "f"))
	csize = 2
	png(scfigs_fn(mc1_id, ofn_base), width=400 + ncol(tab_n) * csize, height=200 + nrow(tab_n) * csize)
	pheatmap(tab_n, color=colorRampPalette(c('white', 'darkred'))(101), cluster_rows=F, cluster_cols=F, cellwidth=csize, cellheight=csize, annotation_row=mc1_ann, annotation_col=mc2_ann, annotation_colors=ann_cols, show_colnames=F, show_rownames=F, border_col=NA, gaps_row=gaps1, gaps_col=gaps2)
	dev.off()
	
	tt = NULL
	row_size = col_size = csize
	if (!is.null(nrow(mc1_ann))) {
		tt = apply(tab, 2, function(v) tapply(v, as.character(mc1_ann$group1), sum) )
		row_size=10
	}
	if (!is.na(mc2_ann)) {
		if (is.null(tt)) {
			tt = tab
		}
		tt = apply(tt, 1, function(v) tapply(v, as.character(mc2_ann$group2), sum) )
		col_size = 10
	}
	
	if (!is.null(tt)) {
		png(scfigs_fn(mc1_id, paste0(ofn_base, "_collaped")), width=max(800, 400 + ncol(tt) * col_size), height=max(800, 200 + nrow(tt) * row_size))
		pheatmap(log2(1+tt), color=colorRampPalette(c('white', 'darkred'))(101), cluster_rows=T, cluster_cols=T, cellwidth=col_size, cellheight=row_size, annotation_row=g1_ann, annotation_col=g2_ann, annotation_colors=ann_cols, show_colnames=!is.null(nrow(g2_ann)), show_rownames=!is.null(nrow(g1_ann)), border_col=NA, treeheight_row=2, treeheight_col=2)
		dev.off()
		
	}
}

#
#
fig_s_pbmc_plots = function(paper_version=F) 
{
	if (paper_version) {
		mc_id = paste0(pbmc_tumors_comp_id, "_paper")
		graph_id = mc_id
	} else {
		mc_id = paste0(pbmc_tumors_comp_clean_id, "_nosmall")
		graph_id = paste0(pbmc_tumors_comp_id, "_nosmall")
	}
	mat_id = graph_id 
	lateral_gset_id = paste0(graph_id, "_lateral")
	
	# annotate the metacells
	conf_pbmc_res <<- colorize_by_confusion_mat(mc_id = mc_id, graph_id = graph_id)
	
	# Manual stop here to generate the supmc and marks files
	if (paper_version) {
		conf_pbmc_res <<- colorize_by_confusion_mat(mc_id= mc_id, graph_id=graph_id, supmc_file="config/pbmc_tumor_supmc_paper.txt", marks_file="config/pbmc_tumor_marks_paper.txt", res=conf_pbmc_res)
	} else {
		conf_pbmc_res <<- colorize_by_confusion_mat(mc_id= mc_id, graph_id= graph_id, supmc_file="config/pbmc_tumor_supmc_nosmall.txt", marks_file="config/pbmc_tumor_marks_nosmall.txt", res=conf_pbmc_res)
	}
	
	mel_basic_mc_mc2d_plots(mc_id=mc_id, mat_id=mat_id, graph_id=graph_id, lateral_gset_id = lateral_gset_id)
	
	#mcell_mc2d_plot(mc_id)
	
	mel_plot_e_gc_barplots(mc_id, "all_hm_genes", genes=t_nk_genes, ncol=2, panel_height=48, panel_width=280, ord_first_by_color=T)

	mat_pbmc = scdb_mat(mat_id)
	mc_pbmc = scdb_mc(mc_id)
	p_col2group = get_mc_col2group(mc_pbmc)
	p_group2col = get_mc_group2col(mc_pbmc)
	md = mat_pbmc@cell_metadata[names(mc_pbmc@mc), ]
	md$group = p_col2group[mc_pbmc@colors[mc_pbmc@mc]]
	md[, 'Cell type'] = gsub(" ", "", md[, 'Cell type'])
	md = md[grep("CCR7", md[, 'Cell type'], invert=T), ]
	cell_by_type = split(rownames(md), md[, 'Cell type'])
	
	if (paper_version) {
		pbmc_grp_nms = c("NK", "effector2" , "dysfunctional", "effector1", "naive", "Treg", "Tfh", "em", "macrophage", "monocyte", "DC", "plasma",	"B")
	} else {
		pbmc_grp_nms = c("NK", "effector2", "dysfunctional", "effector1", "naive", "naive-like", "Treg", "Tfh", "macrophage", "monocyte", "DC", "pDC", "plasma",	"B")
	}
	
	all_pbmc_pat_nms = c("p13-3-(S)C-N", "p17-3-(S)C-N-1", "p17-3-(S)C-N-2", "p27-4-(S)C-1IT", "p13_PBMC-3-(S)C-N", "p17_PBMC-3-(S)C-N", "p27_PBMC-4-(S)C-1IT")
	
	mcell_mc2d_plot_by_factor(mc2d_id=mc_id, mat_id=mat_id, meta_field="PatientID", single_plot=F, neto_points=T, filter_values=all_pbmc_pat_nms)
	
	pbmc_pat_nms = grep("^p17", all_pbmc_pat_nms, v=T, invert=T)
	
	.plot_start(scfigs_fn(mc_id, "pat_grp_comp"), 600, 300)
	layout(matrix(c(3,1:2), nrow=1))
	for (ctype in c('CD3', 'CD45')) {
		dum = expand.grid(pbmc_pat_nms, pbmc_grp_nms)
		p_pat_grp = table(c(as.character(md[cell_by_type[[ctype]], 'PatientID']), as.character(dum[,1])), c(as.character(md[cell_by_type[[ctype]], 'group']), as.character(dum[,2])))
		p_pat_grp = p_pat_grp -1
		p_pat_grp_n = p_pat_grp / rowSums(p_pat_grp)
		par(mar=c(4,1,4,1))
		barplot(t(p_pat_grp_n[pbmc_pat_nms, intersect(pbmc_grp_nms, colnames(p_pat_grp_n))]), horiz=T, col=p_group2col[pbmc_grp_nms], las=2, yaxt=ifelse(ctype == 'CD3', 's', 'n'), main=ctype)
	}
	dev.off()

	metacell:::plot_color_bar(vals=pbmc_grp_nms, cols=p_group2col[pbmc_grp_nms], title="", fig_fn=scfigs_fn(mc_id, "group_colors_legend"))

	
}

#
# Wrapper for mcell_mc2d_plot_gene
mel_mc2d_plot_genes_of_interest = function(mc_id, mc2d_id=mc_id, min_max_lfp=1)
{
	mc = scdb_mc(mc_id)
	lfp = log2(mc@mc_fp)
	
	genes_of_interest = unique(unlist(c('CCR6', 'ENTPD1', 'MAF', 'SSR4', 'LCP1', 'CTLA4', 'TIGIT', 'ICOS', 'ICOSLG', 'BTLA', 'HAVCR2', 'PDCD1', 'BATF', 'RBPJ', 'LAG3', 'FOXP3', 'PRF1', 'GNLY', 'FCER1G', 'TCF7', 'FGFBP2', 'MKI67', sapply(c('^CTS[A-Z]$', '^CD[0-9]', '^GZM', '^KLR', '^CCL', '^CCR', '^CX3', '^CXC', '^IL[0-9]', '^ID[0-9]', '^TGFB', '^IFN', '^IKZ', '^TNFSF', '^TNFRSF', '^XCL', '^XCR'), grep, x=rownames(lfp), perl=T, v=T)))) 
	genes_of_interest = genes_of_interest[ apply(lfp[genes_of_interest, ], 1, max) >= min_max_lfp ]
	
	for (gene in genes_of_interest) {
		mcell_mc2d_plot_gene(mc2d_id, gene, show_legend=T, neto_points=T)
	}
	
}


#
# Fig 1D metacell barplot per gene
mel_plot_e_gc_barplots = function(mc_id, name, genes=NULL, ncolumns=2, panel_height=50, panel_width=300, ord_first_by_color=T, n_ideal_umi=1000) 
{
	mc = scdb_mc(mc_id)
	col2group = get_mc_col2group(mc)
	
	if (is.null(genes)) {
		marks_gset = scdb_gset(mc_id)
		genes = names(marks_gset@gene_set)
	}
	e_gc = mc@e_gc[genes, ] * n_ideal_umi
	
	if (ord_first_by_color) {
		e_gc = e_gc[, order(as.numeric(ordered(col2group[mc@colors], levels=names(ord_by_id[[mc_id]]))) + as.numeric(colnames(e_gc)) * 1e-6)]
	}
	
	.plot_start(scfigs_fn(mc_id, sprintf("mc_geom_mean_%s", name)), w=ncolumns * panel_width, h=panel_height * ceiling(length(genes) / ncolumns))
	layout(matrix(1:(length(genes) + length(genes) %% 2), ncol=ncolumns))
	par(mar=c(0.5, 12, 0.5, 1))
	
	for (g in genes) {
		barplot(e_gc[g, ], border=NA, col=mc@colors[as.numeric(colnames(e_gc))], xaxt='n', yaxt='n', space=0)
		yaxp = par("yaxp")
		axis(2, yaxp=c(yaxp[1], yaxp[2], 1), las=2, cex=1)
		mtext(g, 2, line=1.5, cex=1.2, las=2)
	}
	dev.off()
	
}

#
# Export breakdown of cells to groups and patients in the generated metacell objects
mel_collect_table_stats = function() 
{
	mc_ids = c(NA, NA, NA, tumor_clean_id, tumor_t_nk_id, tumor_non_t_nk_id)
	mat_ids = c(all_id, filt_id, tumor_id, tumor_id, tumor_t_nk_id, tumor_non_t_nk_id)
	
	for (i in seq_along(mc_ids)) {
		message(mat_ids[i])
		mat = scdb_mat(mat_ids[i])
		if (!is.null(mat)) {
			if (is.na(mc_ids[i])) {
				pat_count = table(mat@cell_metadata[mat@cells, 'PatientID'])
				write.table(pat_count, sprintf("%s/%s.cells_per_patient.txt", .scfigs_base, mat_ids[i]), sep="\t", quote=F)
			}
			else {
				mc = scdb_mc(mc_ids[i])
				col2group = get_mc_col2group(mc)
				pat_grp = table(mat@cell_metadata[names(mc@mc), 'PatientID'], col2group[mc@colors[mc@mc]])
				write.table(pat_grp, sprintf("%s/%s.cells_per_patient_and_group.txt", .scfigs_base, mc_ids[i]), sep="\t", quote=F)
				
				md = mat@cell_metadata[names(mc@mc), ]
				md$mc = mc@mc
				md$mc_group = col2group[mc@colors[mc@mc]]
				
				write.table(md, sprintf("%s/%s.metadata_mc_membership.txt", .scfigs_base, mc_ids[i]), sep="\t", quote=F)
			}
		}
	}
	
	for (gset_id in c("human_IG", "human_mito", "human_ncRNA", "human_ncRP", "mel_cc_filt", "mel_ifn_filt", "mel_stress_filt")) {
		gset = scdb_gset(gset_id)
		write.table(sort(names(gset@gene_set)), sprintf("%s/gene_set.%s.txt", .scfigs_base, gset_id), sep="\t", quote=F)
	}          
	
}

#
#
generate_figs = function() 
{
	conf_res = colorize_by_confusion_mat(mc_id=tumor_clean_id, graph_id=tumor_id)
	
	# after manual build of supmc/marks file and generation of the supmc and marks files 
	conf_res = colorize_by_confusion_mat(mc_id=tumor_clean_id, graph_id=tumor_id, supmc_file="config/all_tumor_supmc_paper.txt", marks_file="config/all_tumor_marks_paper.txt", res=conf_res)
	
	# remove erytrocytes mc
	mc_a = scdb_mc(tumor_clean_id)
	scdb_add_mc(paste0(tumor_clean_id, "_all"), mc_a)
	a_group2col = get_mc_group2col(mc_a)
	mc_a = mc_set_outlier_mc(mc_a, which(mc_a@colors == a_group2col['erythrocyte']))
	
	# split tumor mc's to T/NK and B/plasma/myeloid/osteoclasts/macrophages/pDC
	mc_a = scdb_mc(tumor_clean_id)
	a_group2col = get_mc_group2col(mc_a)
	
	col2grp = c(rep('T_NK', 2), rep('non_T_NK', 5), rep('rest', 2))
	names(col2grp) = a_group2col[c('T', 'NK', 'B', 'plasma', 'myeloid', 'osteoclast', 'pDC', 'melanocyte', 'erythrocyte')]
	mcell_mc_split_by_color_group(tumor_clean_id, tumor_id, col2grp=col2grp)
	
	# T/NK mcs
	t_sup <<- colorize_by_confusion_mat(mc_id=tumor_t_nk_id, graph_id=tumor_id)
	
	# after manual build of supmc/marks file
	t_sup <<- colorize_by_confusion_mat(mc_id=tumor_t_nk_id, graph_id=tumor_id, supmc_file="config/t_nk_tumor_supmc_paper.txt", marks_file="config/t_nk_tumor_marks_paper.txt", res=t_sup)
	
	mcell_mc_plot_hierarchy(mc_id=tumor_t_nk_id, graph_id=tumor_id, mc_order=t_sup$mc_hc$order, sup_mc = t_sup$mc_sup, width=2000, height=3000, min_nmc=2, shades=colorRampPalette((c('white', 'darkred', 'black'))), plot_grid=F)
	
	mel_basic_mc_mc2d_plots(mc_id=tumor_t_nk_id, mat_id=tumor_id, graph_id=tumor_id, lateral_gset_id = lateral_gset_id)
	
	# non T/NK mcs: B, plasma, myeloid, osteoclast, pDC
	nt_sup <<- colorize_by_confusion_mat(mc_id=tumor_non_t_nk_id, graph_id=tumor_id)
	
	# after manual build of supmc/marks file
	nt_sup <<- colorize_by_confusion_mat(mc_id=tumor_non_t_nk_id, graph_id=tumor_id, supmc_file="config/non_t_nk_tumor_supmc_paper.txt", marks_file="config/non_t_nk_tumor_marks_paper.txt", res=nt_sup)
	
	mcell_mc_plot_hierarchy(mc_id=tumor_non_t_nk_id, graph_id=tumor_id, mc_order=nt_sup$mc_hc$order, sup_mc = nt_sup$mc_sup, width=2100, height=3000, min_nmc=2, shades=colorRampPalette((c('white', 'darkred', 'black'))), plot_grid=F)
	
	mel_basic_mc_mc2d_plots(mc_id=tumor_non_t_nk_id, mat_id=tumor_id, graph_id=tumor_id, lateral_gset_id = lateral_gset_id)
	
	
	# load objects
	mc_t_nk = scdb_mc(tumor_t_nk_id)
	lfp_t = log2(mc_t_nk@mc_fp)
	
	mat_t_nk = scdb_mat(tumor_t_nk_id)
	mat_t_nk_ds = scm_downsamp(mat_t_nk@mat, 500)
	
	t_group2col = get_mc_group2col(mc_t_nk)
	t_col2group = get_mc_col2group(mc_t_nk)
	
	mc_non_t_nk = scdb_mc(tumor_non_t_nk_id)
	lfp_nt = log2(mc_non_t_nk@mc_fp)
	
	mat_non_t_nk = scdb_mat(tumor_non_t_nk_id)
	mat_non_t_nk_ds = scm_downsamp(mat_non_t_nk@mat, 500)
	
	nt_group2col = get_mc_group2col(mc_non_t_nk)
	nt_col2group = get_mc_col2group(mc_non_t_nk)
	
	mat_all = scdb_mat(tumor_id)
	
	global_col2group = c(t_col2group, nt_col2group)
	global_group2col = c(t_group2col, nt_group2col)

	# Fig 6D boxplots
	g_6d_foc = c('effector2', 'naive', 'effector1', 'dysfunctional')
	f_6d_foc = t_col2group[mc_t_nk@colors] %in% g_6d_foc
	mc_6d_o = ordered(t_col2group[mc_t_nk@colors[f_6d_foc]], levels=g_6d_foc)
	for (gene in c('ENTPD1', 'ITGAE', 'KLRG1')) {
		.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("%s_on_cd8_naive_mcs", gene)), 400, 500)
		par(mar=c(10, 4, 4, 1))
		boxplot(lfp_t[gene, f_6d_foc] ~ mc_6d_o, las=2, main=gene, pch=26, ylab="lfp")
		stripchart(lfp_t[gene, f_6d_foc] ~ mc_6d_o, pch=19, cex=1, method="jitter", jitter=0.3, vertical=T, col=t_group2col[g_6d_foc], add=T)
		
		dev.off()
	}

	####
	
	prev_params = override_metacell_params(list(mcell_mc2d_plot_key=F, mcell_mc2d_height=1500, mcell_mc2d_width=1500, mcell_mc2d_cex=0.5))
	mcell_mc2d_plot(tumor_clean_id, plot_edges=T)
	mcell_mc2d_plot(tumor_t_nk_id, plot_edges=F)
	mcell_mc2d_plot(tumor_t_nk_id, plot_edges=T)
	restore_metacell_params(prev_params)
	
	mc_t_nk_pats = table(mc_t_nk@mc, mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'])
	tumor_pats = intersect(all_tumor_pats, names(which(colSums(mc_t_nk_pats) >= min_t_nk_cells_per_patient)))
	
	metacell:::plot_color_bar(vals=names(ord_t_nk_nms), cols=t_group2col[names(ord_t_nk_nms)], title="", fig_fn=scfigs_fn(tumor_t_nk_id, "group_colors_legend"))
	
	# group composition barplots (for Fig 3B)
	pats_grps_t_nk = table( mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'], t_col2group[mc_t_nk@colors[mc_t_nk@mc]])
	pats_grps_t_nk = pats_grps_t_nk[tumor_pats, ]
	pats_by_dysf = tumor_pats[order(pats_grps_t_nk[, 'dysfunctional']/ rowSums(pats_grps_t_nk[, t_nms]))]
	t_pats_ord = mc_group_composition_barplots(tumor_t_nk_id, tumor_t_nk_id, t_nms, 'T_clust', t_col2group, t_group2col, pat_groups=list(Tumors=tumor_pats), clust_pats=T, min_pat_cells=80)
	t_pats_dysf_ord = mc_group_composition_barplots(tumor_t_nk_id, tumor_t_nk_id, t_nms, 'T_by_dysf', t_col2group, t_group2col, pat_groups=list(Tumors=pats_by_dysf), clust_pats=F, min_pat_cells=80)
	
	# plot specific genes lfp on the 2D projection (Fig 1E)
	mel_mc2d_plot_genes_of_interest(tumor_t_nk_id)
	
	# run proliferation analysis on T/NK (Fig 5A-B)
	prolif_df = mel_mc_prolif_score()
	
	# bargraphs of key genes (Fig 1D)
	mel_plot_e_gc_barplots(tumor_t_nk_id, "all_hm_genes", genes=t_nk_genes, ncol=2, panel_height=48, panel_width=280, ord_first_by_color=T)
	
	# Fig 1F
	mc_pat = table(mc_t_nk@mc, mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'])
	v = rowSums(mc_pat >= 2)
	.plot_start(scfigs_fn(tumor_t_nk_id, "patients_per_mc_ge2_cells"), 400, 240)
	hist(v, length(unique(v)), col='royalblue', main="", xlab="#patients")
	dev.off()

	mel_collect_table_stats()
	
	#
	# Fig 2 and related sup figs
	#
	tot_t_nk_per_pat = table(mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'])

	lat_gset = scdb_gset(lateral_gset_id)
	lat_genes = names(lat_gset@gene_set)
	
	# genes cor heatmap on CD8 cells or metacells
	f_cd8_mc = mc_t_nk@colors %in% t_group2col[cd8_nms]
	f_cd8_c = f_cd8_mc[mc_t_nk@mc]
	
	cd8_genes = setdiff(names(which(apply(abs(lfp_t[, f_cd8_mc]), 1, max) > log2(3))), lat_genes)
	
	cd8_genes = grep("^RP", cd8_genes, invert=T, v=T)
	all_genes_ann = annotate_genes(rownames(lfp_t))
	g_ann = data.frame(row.names=cd8_genes, type=all_genes_ann[cd8_genes, 'type'])
	
	cd8_cor_c = tgs_cor(t(as.matrix(mat_t_nk_ds[cd8_genes, f_cd8_c])), spearman=T)
	cd8_cor_mc = cor(t(lfp_t[cd8_genes, f_cd8_mc]))
	diag(cd8_cor_c) = NA
	diag(cd8_cor_mc) = NA
	
	# Fig 2a: CD8 gene-gene cor heatmaps
	blwtrd_cols = colorRampPalette(c('blue', 'white', 'red'))(101)
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_cd8_c_genes_cor"), max(700, 300 + length(cd8_genes) * 12), max(700, 300 + length(cd8_genes) * 12))
	pheatmap(pmin(pmax(cd8_cor_c, -0.1), 0.1), clustering_method="ward.D2", cutree_rows=15, cutree_cols=15, treeheight_col=0, treeheight_row=0, cellwidth=10, cellheight=10, fontsize=12, col=blwtrd_cols, show_colnames=F) #annotation_row=g_ann, annotation_colors=list(type=ggroup_cols), 
	dev.off()
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_cd8_mc_genes_cor"), max(700, 300 + length(cd8_genes) * 12), max(700, 300 + length(cd8_genes) * 12))
	pheatmap(pmin(pmax(cd8_cor_mc, -0.7), 0.7), clustering_method="ward.D2", cutree_rows=15, cutree_cols=15, treeheight_col=0, treeheight_row=0, cellwidth=10, cellheight=10, fontsize=12, col=blwtrd_cols, show_colnames=F) #, annotation_row=g_ann, annotation_colors=list(type=ggroup_cols)
	dev.off()
	
	# Find high variance non-TF and non blacklisted genes
	gf = apply(lfp_t[, f_cd8_mc], 1, max) > 1 & !grepl("TF", all_genes_ann$type) & (! (rownames(lfp_t) %in% lat_genes))
	
	dysf_anchor = 'LAG3'
	dysf_genes =  tail(sort(cor(lfp_t[dysf_anchor, f_cd8_mc], t(lfp_t[gf, f_cd8_mc]))[1,]), 30)
	dysf_sc = colMeans(lfp_t[names(dysf_genes), ])
	
	# S2B - cor to lag3 - mc vs sc
	dysf_anchor_mc_cor = cor(lfp_t[dysf_anchor, f_cd8_mc], t(lfp_t[setdiff(rownames(lfp_t)[gf], dysf_anchor), f_cd8_mc]))[1,]
	dysf_anchor_c_cor = cor(as.matrix(mat_t_nk_ds[dysf_anchor, f_cd8_c]), as.matrix(t(mat_t_nk_ds[setdiff(rownames(lfp_t)[gf], dysf_anchor), f_cd8_c])), method='spearman')[1,]

	.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("fig_s2b_genes_cor_to_dysf_anchor_mc_vs_sc_labels")), 400, 400)
	lab_f = dysf_anchor_c_cor >= min(dysf_anchor_c_cor[names(dysf_genes)], na.rm=T)
	plot(dysf_anchor_mc_cor, dysf_anchor_c_cor, pch=19, cex=0.5, col=ifelse(names(dysf_anchor_c_cor) %in% names(dysf_genes), 'red', 'black'), xlim=c(0, max(dysf_anchor_mc_cor)), ylim=c(0, max(dysf_anchor_c_cor)), xlab=sprintf("cor to %s on mc", dysf_anchor), ylab=sprintf("cor to %s on cells", dysf_anchor), main=sprintf("R^2 = %.2f", cor(dysf_anchor_mc_cor, dysf_anchor_c_cor, use='pair')**2))
	text(dysf_anchor_mc_cor[lab_f], dysf_anchor_c_cor[lab_f], names(dysf_anchor_c_cor[lab_f]), cex=0.7, col=ifelse(names(dysf_anchor_c_cor[lab_f]) %in% names(dysf_genes), 'red', 'black'), pos=2)
	dev.off()
	
	.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("fig2_genes_cor_to_dysf_anchor_mc_vs_sc")), 400, 400)
	plot(dysf_anchor_mc_cor, dysf_anchor_c_cor, pch=19, cex=0.5, col=ifelse(names(dysf_anchor_c_cor) %in% names(dysf_genes), 'red', 'black'), xlab=sprintf("cor to %s on mc", dysf_anchor), ylab=sprintf("cor to %s on cells", dysf_anchor))
	dev.off()
	
	# helper func
	plot_genes_by_anchors = function(anchor_gene, n=30, name='dysfunctional')
	{
		gf = apply(lfp_t[, f_cd8_mc], 1, max) > 1 & !grepl("TF", all_genes_ann$type) & (! (rownames(lfp_t) %in% lat_genes))
		
		genes =  tail(sort(cor(lfp_t[anchor_gene, f_cd8_mc], t(lfp_t[gf, f_cd8_mc]))[1,]), n)
		
		.plot_start(scfigs_fn(tumor_t_nk_id, paste0("fig2_", name, "_genes_by_", anchor_gene, "_barplot")), 200, 500)
		par(mar=c(4,8,1,1))
		#barplot(dysf_genes, horiz=T, las=2, col=all_genes_ann[names(dysf_genes), 'color'])
		barplot(genes, horiz=T, col=t_group2col[name], yaxt='n', xlab=paste('cor to', anchor_gene))
		mtext(names(genes), 2, at=seq(0.6, len=length(genes), by=1.2), las=2, line=0.5)
		dev.off()
	}
	
	# Alternative anchors as a naive control
	for (anchor_gene in c("LAG3", "HAVCR2")) {
		plot_genes_by_anchors(anchor_gene)
	}
	
	for (i in seq_along(dysf_genes)) { 
		plt("Dysf-score", names(dysf_genes)[i], lfp_t[, f_cd8_mc], cols=mc_t_nk@colors[f_cd8_mc], x=dysf_sc[f_cd8_mc], show_mc_ids=F, add_grid=T, ofn=scfigs_fn(tumor_t_nk_id, sprintf("dysf_grad_lfp_%d_%s", i, names(dysf_genes)[i]), scfigs_dir(tumor_t_nk_id, "fig2_grad_lfps")), cex.lab=1.5)
	}
	
	cyto_anchor = 'FGFBP2'
	cyto_genes =  tail(sort(cor(lfp_t[cyto_anchor, f_cd8_mc], t(lfp_t[gf, f_cd8_mc]))[1,]), 30)
	cyto_sc = colMeans(lfp_t[names(cyto_genes), ])
	
	# S2C-D - robustness
	g_cor_to_dysf = sapply(rownames(lfp_t)[gf], function(nm) cor(lfp_t[nm, f_cd8_mc], colMeans(lfp_t[setdiff(names(dysf_genes), nm), f_cd8_mc])))
	g_cor_to_cyto = sapply(rownames(lfp_t)[gf], function(nm) cor(lfp_t[nm, f_cd8_mc], colMeans(lfp_t[setdiff(names(cyto_genes), nm), f_cd8_mc])))
	
	dysf_anchor_full_mc_cor = cor(lfp_t[dysf_anchor, f_cd8_mc], t(lfp_t[gf, f_cd8_mc]))[1,]
	dysf_anchor_ranks = length(dysf_anchor_full_mc_cor) + 1 - rank(dysf_anchor_full_mc_cor)
	yy = tail(sort(g_cor_to_dysf), 60)
	.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("S2C_dysf_robustness")), 300, 700)
	par(mar=c(4,12,1,1))
	barplot(yy, horiz=T, col=ifelse(names(yy) %in% names(dysf_genes), t_group2col['dysfunctional'], 'darkgrey'), las=2, border=NA, cex.names=0.7, names.arg=sprintf("%s (%d)", names(yy), dysf_anchor_ranks[names(yy)]))
	dev.off()
	
	cyto_anchor_full_mc_cor = cor(lfp_t[cyto_anchor, f_cd8_mc], t(lfp_t[gf, f_cd8_mc]))[1,]
	cyto_anchor_ranks = length(cyto_anchor_full_mc_cor) + 1 - rank(cyto_anchor_full_mc_cor)
	yy = tail(sort(g_cor_to_cyto), 60)
	.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("S2D_cyto_robustness")), 300, 700)
	par(mar=c(4,12,1,1))
	barplot(yy, horiz=T, col=ifelse(names(yy) %in% names(cyto_genes), t_group2col['effector2'], 'darkgrey'), las=2, border=NA, cex.names=0.7, names.arg=sprintf("%s (%d)", names(yy), cyto_anchor_ranks[names(yy)]))
	dev.off()
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_cyto_genes_barplot"), 200, 500)
	par(mar=c(4,8,1,1))
	#barplot(cyto_genes, horiz=T, las=2, col=all_genes_ann[names(cyto_genes), 'color'])
	barplot(cyto_genes, horiz=T, col=t_group2col['effector2'], yaxt='n', xlab=paste('cor to', cyto_anchor))
	mtext(names(cyto_genes), 2, at=seq(0.6, len=length(cyto_genes), by=1.2), col=ifelse(all_genes_ann[names(cyto_genes), 'color'] == 'grey', 'black', 'black'), las=2, line=0.5)
	dev.off()
	
	for (anchor_gene in c("FGFBP2", "CX3CR1")) {
		plot_genes_by_anchors(anchor_gene, name="effector2")
	}
	
	for (i in seq_along(cyto_genes)) { 
		plt("Eff2-score", names(cyto_genes)[i], lfp_t[, f_cd8_mc], cols=mc_t_nk@colors[f_cd8_mc], x=cyto_sc[f_cd8_mc], show_mc_ids=F, add_grid=T, ofn=scfigs_fn(tumor_t_nk_id, sprintf("cyto_grad_lfp_%d_%s", i, names(cyto_genes)[i]), scfigs_dir(tumor_t_nk_id, "fig2_grad_lfps")), cex.lab=1.5)
	}
	
	plt("Dysf-score", "Eff2-score", lfp_t[, f_cd8_mc], mc_t_nk@colors[f_cd8_mc], y=cyto_sc[f_cd8_mc], x=dysf_sc[f_cd8_mc], ofn=scfigs_fn(tumor_t_nk_id, "fig2_plt_cyto_vs_dysf_cd8"), show_mc_ids=F, cex=3, add_grid=T)
	
	tot_mc_umis = tapply(Matrix::colSums(mat_t_nk@mat[, names(mc_t_nk@mc)]), mc_t_nk@mc, sum)
	plt("Dysf-score", "total mc umis (log2)", lfp_t[, f_cd8_mc], mc_t_nk@colors[f_cd8_mc], x=dysf_sc[f_cd8_mc], y=log2(tot_mc_umis[f_cd8_mc]), ofn=scfigs_fn(tumor_t_nk_id, "fig2_plt_tot_mc_umi_vs_dysf_cd8"), show_mc_ids=F, cex=3, add_grid=T)
	plt("Eff2-score", "total mc umis (log2)", lfp_t[, f_cd8_mc], mc_t_nk@colors[f_cd8_mc], x=cyto_sc[f_cd8_mc], y=log2(tot_mc_umis[f_cd8_mc]), ofn=scfigs_fn(tumor_t_nk_id, "fig2_plt_tot_mc_umi_vs_cyto_cd8"), show_mc_ids=F, cex=3, add_grid=T)
	
	# %cyto UMIs vs %dysf UMIs on single cells (+ %cell-cycle)
	cc_res = mel_mc_prolif_score()
	
	cd8_cells_f = f_cd8_c
	
	tot = Matrix::colSums(mat_t_nk@mat[, cd8_cells_f])
	tot_dysf = Matrix::colSums(mat_t_nk@mat[names(dysf_genes), cd8_cells_f])
	tot_cyto = Matrix::colSums(mat_t_nk@mat[names(cyto_genes), cd8_cells_f])
	max_dysf_cyto_f_umi = 0.05
	sc_tot = data.frame(row.names=names(tot), f_cyto=tot_cyto/tot, f_dysf=tot_dysf/tot, group= t_col2group[mc_t_nk@colors[mc_t_nk@mc[cd8_cells_f]]], patient=mat_t_nk@cell_metadata[names(tot), 'PatientID'], f_cc=cc_res$f_cc[names(tot)])
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_cyto_vs_dysf_cells_cd8"), 900, 300)
	layout(matrix(1:3, 1, 3))
	cd8_nms_to_plot = cd8_nms
	par(mar=c(5,5,1,1))
	
	for (g in cd8_nms_to_plot) {
		plot(pmin(sc_tot$f_dysf, max_dysf_cyto_f_umi), pmin(sc_tot$f_cyto, max_dysf_cyto_f_umi), pch=19, cex=0.4, col='darkgrey', xlab='% Dysf UMIs', ylab='% Eff2 UMIs')
		points(pmin(sc_tot[sc_tot$group == g, 'f_dysf'], max_dysf_cyto_f_umi), pmin(sc_tot[sc_tot$group == g, 'f_cyto'], max_dysf_cyto_f_umi), pch=19, cex=0.4, col=t_group2col[g]) 
	}
	dev.off()
	
	# f_dysf on CD8s by patient + color by reactivity
	sc_tot = sc_tot[sc_tot$patient %in% tumor_pats, ]
	pats_dysf_load = tapply(sc_tot$f_dysf, sc_tot$patient, median)
	pats_by_cd8_f_dysf = names(sort(pats_dysf_load))
	o_dysf = ordered(sc_tot$patient, levels=pats_by_cd8_f_dysf)
	
	f_dysf_eff2_reg = 1e-3
	sc_tot$l_f_dysf = log2(f_dysf_eff2_reg + sc_tot$f_dysf)
	sc_tot$l_f_eff2 = log2(f_dysf_eff2_reg + sc_tot$f_cyto)
	sc_tot$l_f_cc = log2(f_dysf_eff2_reg + sc_tot$f_cc)
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig6_f_cyto_vs_dysf_on_cd8_by_patient"), 900, 400)
	par(mar=c(10,5,1,1))
	boxplot(l_f_dysf  ~ o_dysf, sc_tot, las=2, boxfill='white', border='white', xlim=c(-5.5, length(pats_by_cd8_f_dysf) + 0.5), ylim=c(log2(f_dysf_eff2_reg) - 1, max(sc_tot[, c('l_f_dysf', 'l_f_eff2')])), ylab='%UMIs (log2)')
	
	boxplot(l_f_dysf ~ o_dysf, sc_tot, varwidth=F, notch=F, pch=19, cex=0.3, col=t_group2col['dysfunctional'], boxwex=0.25, at=seq_along(pats_by_cd8_f_dysf)-0.15, add=T, xaxt='n', yaxt='n', las=2)
	boxplot(l_f_eff2 ~ o_dysf, sc_tot, varwidth=F, notch=F, pch=19, cex=0.3, col=t_group2col['effector2'], boxwex=0.25, at=seq_along(pats_by_cd8_f_dysf)+0.15, add=T, xaxt='n', yaxt='n', las=2)
	points(seq_along(pats_by_cd8_f_dysf), rep(log2(f_dysf_eff2_reg) - 0.5, length(pats_by_cd8_f_dysf)), pch=21, cex=3, bg=pat_reac_col[pats_by_cd8_f_dysf], col=ifelse(pats_by_cd8_f_dysf %in% na_reac_pats, 'white', 'black'))
	legend("topleft", legend=c('dysfunctional', 'effector2'), fill=t_group2col[c('dysfunctional', 'effector2')], bty='n', cex=1.5)
	leg_labs = reac_col2type[intersect(reac_pat_cols, pat_reac_col)]
	
	legend("bottomleft", legend=leg_labs, pch=21, pt.bg=names(leg_labs), pt.cex=2, bty='n', cex=1.5)
	dev.off()
	
	pat_cd8_comp = table(mat_t_nk@cell_metadata[names(tot), 'PatientID'], t_col2group[mc_t_nk@colors[mc_t_nk@mc[names(tot)]]])
	pat_cd8_comp_n = pat_cd8_comp/ rowSums(pat_cd8_comp)
	pat_cd8_comp_n = pat_cd8_comp_n[pats_by_cd8_f_dysf, intersect(colnames(pat_cd8_comp_n), t_nms)]
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig6_cd8_comp_by_patient"), 900, 250)
	par(mar=c(10,5,1,1))
	barplot(t(pat_cd8_comp_n), las=2, co=t_group2col[colnames(pat_cd8_comp_n)], xlim=c(-5.5, length(pats_by_cd8_f_dysf)*1.2 + 0.5), ylab='fraction')
	dev.off()
	
	dysf_vs_cyto_model = model_by_tfs(tumor_t_nk_id, dysf_sc, cyto_sc, all_genes_ann, max_genes_for_plot_size = 30, f=f_cd8_mc)
	
	# Fig S4B: top cor genes to dysf load on groups per patient
	group_per_patient_top_cor_to_dysf_helper = function(group, n_genes=50) 
	{
		pat_geom = genes_geomean_by_patients_for_group(tumor_non_t_nk_id, group=group)
		curr_pats = intersect(colnames(pat_geom$geom), names(which(!is.na(pats_dysf_load))))
		cg = tail(sort(cor(pats_dysf_load[curr_pats], t(pat_geom$geom[, curr_pats]), method='spearman')[1,]), n_genes)
		tot_cg = rowSums(pat_geom$tot[names(cg), ])
		
		.plot_start(scfigs_fn(tumor_non_t_nk_id, sprintf("S4B_%s_by_pat_top_g_dysf_load_cor_barplot", group)), 300, 700)
		par(mar=c(4,15,1,1))
		barplot(cg, horiz=T, col=nt_group2col[group], yaxt='n', xlab='cor to patient dysf load')
		mtext(sprintf("%s (%d)", names(cg), tot_cg), 2, at=seq(0.6, len=length(cg), by=1.2), las=2, line=0.5)
		dev.off()
		
	}
	for (gg in c('monocyte', 'B')) {
		group_per_patient_top_cor_to_dysf_helper(gg)
	}
	
	### CD4 (Treg)
	
	# genes cor heatmap
	f_treg_mc = mc_t_nk@colors %in% t_group2col['Treg']
	f_treg_c = f_treg_mc[mc_t_nk@mc]
	treg_genes = setdiff(names(which(apply(lfp_t[, f_treg_mc], 1, max) > log2(3))), lat_genes)
	
	treg_genes = grep("^RP", treg_genes, invert=T, v=T)
	g_ann = data.frame(row.names=treg_genes, type=all_genes_ann[treg_genes, 'type'])
	
	treg_cor_c = tgs_cor(t(as.matrix(mat_t_nk_ds[treg_genes, f_treg_c])), spearman=T)
	treg_cor_mc = cor(t(lfp_t[treg_genes, f_treg_mc]))
	
	# Treg gene-gene cor heatmaps
	treg_cutree = 5
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_treg_c_genes_cor"), max(700, 300 + length(treg_genes) * 12), max(700, 300 + length(treg_genes) * 12))
	pheatmap(pmin(pmax(treg_cor_c, -0.1), 0.1), clustering_method="ward.D2", cutree_rows=treg_cutree, cutree_cols=treg_cutree, treeheight_col=0, annotation_row=g_ann, annotation_colors=list(type=ggroup_cols), cellwidth=10, cellheight=10, fontsize=12)
	dev.off()
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_treg_mc_genes_cor"), max(700, 300 + length(treg_genes) * 12), max(700, 300 + length(treg_genes) * 12))
	pheatmap(pmin(pmax(treg_cor_mc, -0.7), 0.7), clustering_method="ward.D2", cutree_rows=treg_cutree, cutree_cols=treg_cutree, treeheight_col=0, annotation_row=g_ann, annotation_colors=list(type=ggroup_cols), cellwidth=10, cellheight=10, fontsize=12)
	dev.off()
	
	# remove TFs for model 
	gf = apply(lfp_t[, f_treg_mc], 1, max) > 1 & !grepl("TF", all_genes_ann$type)& (! (rownames(lfp_t) %in% lat_genes))
	#gf = apply(lfp_t[, f_treg_mc], 1, max) > 1
	
	treg_anchor = "IL2RA"
	treg_genes =  tail(sort(cor(lfp_t[treg_anchor, f_treg_mc], t(lfp_t[gf, f_treg_mc]))[1,]), 30)
	treg_sc = colMeans(lfp_t[names(treg_genes), ])
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_Treg_genes_barplot"), 200, 500)
	par(mar=c(4,8,1,1))
	barplot(treg_genes, horiz=T, col=t_group2col['Treg'], yaxt='n', xlab=paste('cor to', treg_anchor))
	mtext(names(treg_genes), 2, at=seq(0.6, len=length(treg_genes), by=1.2), col=ifelse(all_genes_ann[names(treg_genes), 'color'] == 'grey', 'black', 'black'), las=2, line=0.5)
	dev.off()
	
	for (i in seq_along(treg_genes)) { 
		plt("Treg-score", names(treg_genes)[i], lfp_t[, f_treg_mc], cols=mc_t_nk@colors[f_treg_mc], x=treg_sc[f_treg_mc], show_mc_ids=F, add_grid=T, ofn=scfigs_fn(tumor_t_nk_id, sprintf("Treg_grad_lfp_%d_%s", i, names(treg_genes)[i]), scfigs_dir(tumor_t_nk_id, "fig2_grad_lfps")), cex.lab=1.5)
	}
	
	plt("Treg-score", "total mc umis (log2)", lfp_t[, f_treg_mc], mc_t_nk@colors[f_treg_mc], x=treg_sc[f_treg_mc], y=log2(tot_mc_umis[f_treg_mc]), ofn=scfigs_fn(tumor_t_nk_id, "fig2_plt_tot_mc_umi_vs_treg_Treg"), show_mc_ids=F, cex=3, add_grid=T)
	
	treg_model = model_by_tfs(tumor_t_nk_id, treg_sc, NULL, all_genes_ann, f=f_treg_mc, nm1='Treg', nm2='None', grp1 = 'Treg', grp2=NULL)
	
	# Compare known exhaustion genes cor to Dysf gradient and to Treg gradient
	fig2_hand_selected_genes_dysf_and_treg_cor = function(size=4, label_cex=0.7, xlim=NULL, ylim=NULL) 
	{
		ex_genes = c("AKAP5", "BTLA", "CD96", "CXCL13","FOXP3", "ICOS", "IFNG", "IKZF2", "IL1R1", "IL1R2", "IL2RA", "IL2RB", "IL32", "PHLDA1", "SNAP47", "TBC1D4", "TFRC", "TNFRSF18","TRAF5", "BATF", "CSF1", "CTLA4", "ENTPD1", "HAVCR2", "LAG3", "PDCD1", "TNFRSF1B", "TNFRSF9", "TIGIT") 
		#ex_genes = c( "BIRC3", "CARD16", "CCL3", "CCND2", "CD244", "CD27", "CD7",  "CTSC",  "F5", "FAM3C", "FASLG", "GADD45A", "GALNT1", "GBP2", "HNRNPA1P21",  "KAT2B", "KIR2DL4", "LAYN", "LYST", "MIR155HG", "NAMPT", "PIM2", "PTMS", "RDH10", "RTKN2", "S100A4", "SIRPG",  "TNIP3", "TNS3", "TOX", "TPP1", "UGP2")
		cd8_exhaust_c = cor(dysf_sc[f_cd8_mc], t(lfp_t[ex_genes, f_cd8_mc]))
		treg_exhaust_c = cor(treg_sc[f_treg_mc], t(lfp_t[ex_genes, f_treg_mc]))
		svg(sprintf("%s/%s.%s.svg", .scfigs_base, tumor_t_nk_id, "exhaust_genes_Treg_vs_Dysf_grads"),5, 5)
		#.plot_start(scfigs_fn(tumor_t_nk_id, "exhaust_genes_Treg_vs_Dysf_grads"), size, size)
		par(mar=c(6,6,6,6))
		plot(cd8_exhaust_c, treg_exhaust_c, pch=21, asp=1, xlim=xlim, ylim=ylim, col='black', bg=all_genes_ann[ex_genes, 'color'], xlab="cor to dysf grad", ylab="cor to Treg grad")
		text(cd8_exhaust_c, treg_exhaust_c, ex_genes, pos=ifelse(cd8_exhaust_c > 0.5, 2, 4), cex=label_cex)
		#grid(col="black", lwd=0.5)
		
		dev.off()
		list(dysf=cd8_exhaust_c, treg=treg_exhaust_c)
	}
	dysf_vs_treg_selected_genes = fig2_hand_selected_genes_dysf_and_treg_cor()
	
	# Single genes TF models
	fig2_tf_model_per_gene = function(min_max_coef=0.1, min_r2=0.5, z=0.4, target_g_coef=NULL) {
		tfs = rownames(all_genes_ann)[grepl("TF", all_genes_ann$type)]
		target_genes = unique(c(names(dysf_genes), names(cyto_genes), names(treg_genes)))
		target_genes_ann = data.frame(row.names=target_genes, 
																	dysf=ifelse(target_genes %in% names(dysf_genes), 1, 0),
																	cyto=ifelse(target_genes %in% names(cyto_genes), 1, 0),
																	Treg=ifelse(target_genes %in% names(treg_genes), 1, 0),
																	type=all_genes_ann[target_genes, 'type'], stringsAsFactors = F)
		tfs_ann = data.frame(row.names=tfs, type=all_genes_ann[tfs, 'type'])
		
		dysf_cols = c('white', t_group2col["dysfunctional"])
		cyto_cols = c('white', t_group2col["effector2"])
		treg_cols = c('white', t_group2col["Treg"])
		names(dysf_cols) = names(cyto_cols) = names(treg_cols) = 0:1
		
		if (is.null(target_g_coef)) {
			target_g_coef = matrix(0, length(target_genes), length(tfs), dimnames=list(target_genes, tfs))

			for (g in target_genes) {
				g_model = model_by_tfs(tumor_t_nk_id, lfp_t[g, ], NULL, all_genes_ann, nm1=g, nm2='None', min_cor_for_barplot = 0, gen_plots = F, f=f_cd8_mc)
				#message(sprintf("%s R^2= %.2f", g, g_model$r2))
				if (length(g_model$coef) > 0 & g_model$r2 >= min_r2) {
					target_g_coef[g, names(g_model$coef)] = g_model$coef
				}
			}
		}
		target_g_coef_f = target_g_coef[apply(target_g_coef, 1, max) > min_max_coef, apply(target_g_coef, 2, max) > min_max_coef]
		
		message(sprintf("Remained with %d (out of %d) target genes and %d (out of %d) TFs", nrow(target_g_coef_f), ifelse(is.null(target_g_coef), -1, nrow(target_g_coef)), ncol(target_g_coef_f), ifelse(is.null(target_g_coef), -1, ncol(target_g_coef))))
		chc = hclust(dist(cor(target_g_coef_f)), method='ward.D2')
		rhc = hclust(dist(cor(t(target_g_coef_f))), method='ward.D2')
		
		.plot_start(scfigs_fn(tumor_t_nk_id, "fig2_tf_per_gene"), max(700, 300 + ncol(target_g_coef_f)*10), max(700, 300 + nrow(target_g_coef_f)*10))
		
		pheatmap(pmin(pmax(t(target_g_coef_f[rhc$order, chc$order]), -z), z), breaks=seq(-z, z, length=100), color=colorRampPalette(c('blue', 'white', 'red'))(101), cluster_rows=F, cluster_cols=F, annotation_col=target_genes_ann, annotation_row=tfs_ann, annotation_colors=list(type=ggroup_cols, dysf=dysf_cols, cyto=cyto_cols, Treg=treg_cols), cellwidth=9, cellheight=9, fontsize=9)
		dev.off()
		
		target_g_coef
	}
	
	target_g_coef_f = fig2_tf_model_per_gene()
	
	plt("dysf-score", "Treg-score", lfp_t, cols=mc_t_nk@colors, x=dysf_sc, y=treg_sc, ofn=scfigs_fn(tumor_t_nk_id, "fig2_plt_treg_vs_dysf_all_mc"), show_mc_ids=F, add_grid=T, cex=2)
	
	# global Tfh vs Treg and Dysf genes: Create mc by groups, splitting Tregs and cd8-dysf to low and high expressed groups
	mc_t_nk_sep = mc_t_nk
	mc_t_nk_sep@colors[dysf_sc > 0.7] = 'gold4'
	mc_t_nk_sep@colors[114] = t_group2col['dysf-cd4'] # keep the dysf-cd4 mc
	mc_t_nk_sep@colors[lfp_t['FOXP3', ] > 2.25] = 'red4'
	mc_t_nk_sep@color_key = rbind(mc_t_nk_sep@color_key, data.frame(gene=c('', ''), group=c('dysf-hi', 'Treg-hi'), color=c('gold4', 'red4')))
	
	treg_dysf_hi_suff = "_Treg_dysf_hi"
	treg_dysf_hi_nms = c("Treg", "Treg-hi", "Tfh", "dysf-cd4", "naive", "effector1", "effector2", "em-cd4", "dysfunctional", "dysf-hi", "NK")
	treg_dysf_hi_ord = seq_along(treg_dysf_hi_nms)
	names(treg_dysf_hi_ord) = treg_dysf_hi_nms
	tumor_t_nk_sep_id = paste0(tumor_t_nk_id, treg_dysf_hi_suff)
	scdb_add_mc(tumor_t_nk_sep_id, mc_t_nk_sep)
	
	tumor_t_nk_sep_group_id = paste0(tumor_t_nk_sep_id, "_by_color_groups")
	
	mc_by_color_group(mc_id=tumor_t_nk_sep_id, mat_id=tumor_t_nk_id, group_ord=treg_dysf_hi_ord, rebuild=T)
	
	mc_t_nk_g = scdb_mc(tumor_t_nk_sep_group_id)
	lfp_t_g = log2(mc_t_nk_g@mc_fp)
	colnames(lfp_t_g) = treg_dysf_hi_nms
	n_genes_per_group = 20
	
	g1s = c('Treg-hi', 'dysf-hi', 'dysf-hi', 'Treg-hi', 'dysf-hi', 'NK')
	g2s = c('Tfh', 'Tfh', 'Treg-hi', 'Treg', 'dysfunctional', 'effector2')
	for (i in 1:length(g1s)) {
		g1 = g1s[i]
		g2 = g2s[i]
		g12_g = union(names(tail(sort(lfp_t_g[, g1]), n_genes_per_group)), names(tail(sort(lfp_t_g[, g2]), n_genes_per_group)))
		.plot_start(scfigs_fn(tumor_t_nk_sep_group_id, sprintf("%s_vs_%s_top_%d_genes_test", g1, g2, n_genes_per_group)),400, 400)
		
		plot(lfp_t_g[g12_g, g1], lfp_t_g[g12_g, g2], pch=21, col='black', bg=all_genes_ann[g12_g, 'color'], xlab=paste(g1, "fp (log2)"), ylab=paste(g2, "fp (log2)"), cex=0.8)
		abline(a=0, b=1, lty=3, col='grey')
		abline(h=0, lty=2)
		abline(v=0, lty=2)
		text(lfp_t_g[g12_g, g1], lfp_t_g[g12_g, g2], g12_g, cex=0.8, pos=ifelse(lfp_t_g[g12_g, g1] > 0, 2, 4))
		dev.off()
	}
	
	# diff expr Treg (enr over naive) vs Dysf (enr over naive)
	t_nk_c_by_grp = split(names(mc_t_nk@mc), t_col2group[mc_t_nk@colors[mc_t_nk@mc]])
	group_enr_over_naive_comparison_helper = function(grp1='Treg', grp2='dysfunctional', zlim=c(-2, 4.5)) 
	{
		grp1_vs_naive = diff_expr(mc_t_nk, mat_t_nk_ds, NULL, NULL, nms1=t_nk_c_by_grp[[grp1]], nms2=t_nk_c_by_grp[['naive']], min_max_umi=0)
		grp2_vs_naive = diff_expr(mc_t_nk, mat_t_nk_ds, NULL, NULL, nms1=t_nk_c_by_grp[[grp2]], nms2=t_nk_c_by_grp[['naive']], min_max_umi=0)
		
		t1 = grp1_vs_naive %>% filter(tot1 > length(t_nk_c_by_grp[[grp1]]) * 0.05) 
		t2 = grp2_vs_naive %>% filter(tot1 > length(t_nk_c_by_grp[[grp2]]) * 0.05) 
		g12b = union(t1$gene, t2$gene)
		
		g12b_s = union(g12b[pmax(grp1_vs_naive[g12b, 'enr'], grp2_vs_naive[g12b, 'enr']) >= 2],   intersect(g12b, c(
			rownames(all_genes_ann)[all_genes_ann$type != 'None'],
			grep("^IL[0-9]|^IFN", g12b, v = T, perl = T)
		)))
		
		for (add_labels in c(T, F)) {
			.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("fig2K_genes_enr_over_naive_%s_vs_%s%s", grp2, grp1, ifelse(add_labels, "_labels", ""))), 600, 600)
			plot(pmin(pmax(grp1_vs_naive[g12b, 'enr'], zlim[1]), zlim[2]), pmin(pmax(grp2_vs_naive[g12b, 'enr'], zlim[1]), zlim[2]), pch=19, cex=0.5, xlab=paste(grp1, 'over Naive enr (log2)'), ylab=paste(grp2, 'over Naive enr (log2)'))
			if (add_labels) {
				text(pmin(pmax(grp1_vs_naive[g12b_s, 'enr'], zlim[1]), zlim[2]), pmin(pmax(grp2_vs_naive[g12b_s, 'enr'], zlim[1]), zlim[2]), g12b_s, cex=0.7, pos=3, col=all_genes_ann[g12b_s, 'color'])
			}
			dev.off()
		}
	}
	group_enr_over_naive_comparison_helper()
	group_enr_over_naive_comparison_helper(grp1='Treg', grp2='Tfh')
	group_enr_over_naive_comparison_helper(grp1='Tfh', grp2='dysfunctional')
	
	############################
	#
	# Myeloid and patient related analysis (S3, Fig 3)
	#
	############################
	message("Starting myeloid related analysis...")
	
	# create non-T mc for plots - merging colors of im/mature-DC --> DC and non-classic-monocyte --> monocyte
	mc_non_t_nk_merged = mc_non_t_nk
	mc_non_t_nk_merged@colors = rep("white", n=length(mc_non_t_nk@colors))
	mc_non_t_nk_merged@color_key = data.frame(group=c(), gene=c(), color=c(), priority=c(), T_fold=c())
	scdb_add_mc(tumor_non_t_nk_merged_id, mc_non_t_nk_merged)
	nt_sup = colorize_by_confusion_mat(mc_id = tumor_non_t_nk_merged_id, graph_id=tumor_id, 
																		 supmc_file=paste0("config/non_t_nk_tumor_supmc_merged_grps.txt"), 
																		 marks_file=paste0("config/non_t_nk_tumor_marks_merged_grps.txt"), res=nt_sup) 
	mc_non_t_nk_merged = scdb_mc(tumor_non_t_nk_merged_id)
	mc2d_non_t_nk_merged = scdb_mc2d(tumor_non_t_nk_id)
	mc2d_non_t_nk_merged@mc_id = tumor_non_t_nk_merged_id
	scdb_add_mc2d(tumor_non_t_nk_merged_id, mc2d_non_t_nk_merged)
	
	prev_params = override_metacell_params(list(mcell_mc2d_plot_key=F, mcell_mc2d_height=1500, mcell_mc2d_width=1500, mcell_mc2d_cex=0.5))
	mcell_mc2d_plot(tumor_non_t_nk_merged_id, plot_edges=T)
	mcell_mc2d_plot(tumor_non_t_nk_merged_id, plot_edges=F)
	restore_metacell_params(prev_params)
	
	
	f_myl_mc = nt_col2group[mc_non_t_nk@colors] %in% myl_nms
	f_myl_c = f_myl_mc[mc_non_t_nk@mc]
	
	# bargraphs of key genes
	non_t_nk_genes =  c('CD74', 'CD14', 'MMP9', 'HLA-DRA', 'MS4A1', 'TYROBP', 'CD79A', 'MZB1', 'XBP1', 'IRF7', 'IRF8', 'TCF4', 'GZMB', 'TCL1A', 'LILRA4', 'LYZ', 'CDKN1C', 'IDO1', 'ISG15', 'CLEC10A', 'CD1C', 'CD1A', 'APOE', 'C1QA', 'FCGR3A', 'S100A8', 'S100A12', 'CD274')
	mel_plot_e_gc_barplots(tumor_non_t_nk_merged_id, "selected_hm_genes", genes=non_t_nk_genes, ncol=2, panel_height=48, panel_width=280, ord_first_by_color=T)
	
	metacell:::plot_color_bar(vals=rev(non_t_merged_nms), cols=global_group2col[rev(non_t_merged_nms)], title="", fig_fn=scfigs_fn(tumor_non_t_nk_merged_id, "group_colors_legend"))
	
	# patient composition barplots
	global_pat_grp = table(mat_all@cell_metadata[c(names(mc_non_t_nk@mc), names(mc_t_nk@mc)), 'PatientID'], global_col2group[c(mc_non_t_nk@colors[mc_non_t_nk@mc], mc_t_nk@colors[mc_t_nk@mc])])
	global_pat_grp = global_pat_grp[tumor_pats, ]
	
	stopifnot(length(t_nms) + length(non_t_nms) == ncol(global_pat_grp) & length(intersect(t_nms, non_t_nms)) == 0)
	
	# normalize #cells separetaly for T/NK and non-T/NK
	global_pat_grp_norm =  global_pat_grp
	global_pat_grp_norm[, t_nms] = global_pat_grp_norm[, t_nms] / rowSums(global_pat_grp_norm[, t_nms])
	global_pat_grp_norm[, non_t_nms] = global_pat_grp_norm[, non_t_nms] / rowSums(global_pat_grp_norm[, non_t_nms])
	
	non_t_pats_ord = mc_group_composition_barplots(tumor_non_t_nk_id, tumor_id, groups=non_t_nms, name="non_T", col2group=global_col2group, group2col=global_group2col, pat_groups=list(Tumors=tumor_pats), pat_grp=global_pat_grp[, non_t_nms], pat_grp_n = global_pat_grp_norm[, non_t_nms], clust_pats=T, min_pat_cells=50)
	
	pats_by_dysf = tumor_pats[order(global_pat_grp_norm[, 'dysfunctional']/ rowSums(global_pat_grp_norm[, t_nms]))]
	
	all_nms = c(t_nms, non_t_nms)
	t_pats_dysf_ord = mc_group_composition_barplots(tumor_non_t_nk_id, tumor_id, groups=all_nms, name='all_by_dysf', global_col2group, global_group2col, pat_groups=list(Tumors=pats_by_dysf), pat_grp = global_pat_grp[, all_nms], pat_grp_n=global_pat_grp_norm[, all_nms], clust_pats=F, min_pat_cells=80)
	
	mc_gen_marks_heatmaps(tumor_non_t_nk_id, tumor_non_t_nk_id, lateral_gset_id)
	
	# stratify patients by %dysf
	pat_f_dysf = global_pat_grp_norm[, 'dysfunctional']
	dysf_pat_strat = ifelse(pat_f_dysf < 0.2, 'low', ifelse(pat_f_dysf < 0.4, 'mid', 'high'))
	o_dysf_strat = ordered(dysf_pat_strat, levels=c('low', 'mid', 'high'))
	
	.plot_f_grp_by_dysf_helper = function(mc_id, x, nm, col) 
	{
		.plot_start(scfigs_fn(mc_id, sprintf("f_%s_strat_by_f_dysf", nm)), 200, 300)
		
		stripchart(x ~ o_dysf_strat, method='jitter', vertical=T, pch=21, cex=2, las=2, main=sprintf("%s (%.2f)", nm, cor(x, global_pat_grp_norm[, 'dysfunctional'], method='spearman')), bg=col, ylab='fraction')
		segments(seq(0.75, by=1, len=3), tapply(x,  o_dysf_strat, median), seq(1.25, by=1, len=3), tapply(x,  o_dysf_strat, median), lwd=2)
		mw = wilcox.test(x[dysf_pat_strat == 'low'], x[dysf_pat_strat == 'high'])
		pv = mw$p.value
		mtext(ifelse(pv < 0.001, "***", ifelse(pv < 0.01, "**", ifelse(pv < 0.05, "*", ""))), side=3, line=0, cex=2, at=2)
		dev.off()
	}
	for (g in unique(names(t_group2col))) {
		.plot_f_grp_by_dysf_helper(tumor_t_nk_id, global_pat_grp_norm[, g], g, t_group2col[g])
	}
	for (g in unique(names(nt_group2col))) {
		.plot_f_grp_by_dysf_helper(tumor_non_t_nk_id, global_pat_grp_norm[, g], g, nt_group2col[g])
	}
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "f_dysf_pat_strats"), 900, 400)
	par(mar=c(12,6,3,1))
	barplot(sort(global_pat_grp_norm[, 'dysfunctional']), col=global_group2col['dysfunctional'], ylab='fraction', las=2, cex.names=1.5, cex.axis=1.5, cex.lab=1.5)
	dysf_strats_sz = table(o_dysf_strat)
	ds_sz_cs = cumsum(dysf_strats_sz)
	axis(3, labels=c('low', 'mid', 'high'), las=1, at=1.2 * ((ds_sz_cs + c(0, ds_sz_cs[-3])) / 2), tick=F, cex.axis=1.5)
	abline(v=0.1+ds_sz_cs[1:2]*1.2)
	dev.off()
	
	# Fig 4C
	f_inf_tab = read.table("data/processed_f_infiltrate.txt", header=T)
	md_nms = c('stage', 'location', 'treatment')
	pats_md = unique(mat_t_nk@cell_metadata[names(mc_t_nk@mc), c('PatientID', md_nms, 'Patient CD3+ % from CD45+')])
	rownames(pats_md) = pats_md$PatientID
	pats_md$f_inflitrate = f_inf_tab[ rownames(pats_md), 'f_inflitrate']
	pats_md[, 5] = as.numeric(pats_md[, 5])
	
	# manual correction following new input from the pathologist!
	pats_md['p3-4-LN-mIT', 'location'] = '(S)C'
	pats_md['p4-4-LN-1IT', 'location'] = '(S)C'
	.plot_start(scfigs_fn(tumor_t_nk_id, "f_dysf_pat_strats_with_metadata"), 700, 460)
	layout(matrix(1:7, ncol=1), heights=c(180, 20, 20, 20, 60, 60, 100))
	
	par(mar=c(0.5,10,3,1))
	barplot(global_pat_grp_norm[pats_by_dysf, 'dysfunctional'], col=global_group2col['dysfunctional'], ylab='fraction', las=2, cex.names=1.5, cex.axis=1.5, cex.lab=1.5, xaxt='n')
	dysf_strats_sz = table(o_dysf_strat)
	ds_sz_cs = cumsum(dysf_strats_sz)
	axis(3, labels=c('low', 'medium', 'high'), las=1, at=1.2 * ((ds_sz_cs + c(0, ds_sz_cs[-3])) / 2), tick=F, cex.axis=1.5)
	abline(v=0.1+ds_sz_cs[1:2]*1.2)
	
	par(mar=c(0.2, 10, 0.2, 1))
	md_dict =  tgconfig::get_param("mcp_metadata_annot_colors", "metacell")
	
	for (nm in md_nms) {
		barplot(rep(1, length(pats_by_dysf)), col=unlist(md_dict[[nm]][pats_md[pats_by_dysf, nm]]), yaxt='n', border=NA)
		mtext(nm, 2, las=2)
	}
	par(mar=c(0.5, 10,0.5,1))
	barplot(pats_md[pats_by_dysf, 5], col='royalblue', las=2, cex.axis=1.5, cex.lab=1.5, ylab='%CD3', border = NA)
	abline(h=0)
	barplot(ifelse(is.na(pats_md[pats_by_dysf, 6]), max(pats_md[pats_by_dysf, 6], na.rm=T), pats_md[pats_by_dysf, 6]), col=ifelse(is.na(pats_md[pats_by_dysf, 6]), 'grey90', 'royalblue4'), las=2, cex.axis=1.5, cex.lab=1.5, ylab='%Inflitrate', names.arg = pats_by_dysf, cex.names=1.5, border=NA)
	abline(h=0)
	dev.off()
	
	# mc composition by patient
	mc_pat_t = table(mc_t_nk@mc, mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'])
	mc_pat_nt = table(mc_non_t_nk@mc, mat_non_t_nk@cell_metadata[names(mc_non_t_nk@mc), 'PatientID'])
	max_mcsz = max(c(table(mc_t_nk@mc), table(mc_non_t_nk@mc)))
	
	plot_mc_metadata(rbind(mc_pat_t, mc_pat_nt), mc_colors=c(mc_t_nk@colors, mc_non_t_nk_merged@colors), mc_groups=global_col2group[c(mc_t_nk@colors, mc_non_t_nk_merged@colors)], is_mc_t_nk=NULL, mat_id=tumor_id, scfigs_fn(tumor_id, "mc_metadata"), pats_order=rownames(global_pat_grp_norm)[order(global_pat_grp_norm[, 'dysfunctional'])], groups_order=names(c(ord_by_id[[tumor_t_nk_id]], ord_by_id[[tumor_non_t_nk_merged_id]])), show_is_t=F)
	
	plot_mc_metadata(mc_pat_t, mc_colors=mc_t_nk@colors, mc_groups=t_col2group[mc_t_nk@colors], is_mc_t_nk=NULL, mat_id=tumor_t_nk_id, scfigs_fn(tumor_t_nk_id, "mc_metadata"), pats_order=rownames(global_pat_grp_norm)[order(global_pat_grp_norm[, 'dysfunctional'])], groups_order=names(ord_by_id[[tumor_t_nk_id]]), show_is_t=F, ncells_ylim=c(0, max_mcsz))
	plot_mc_metadata(mc_pat_nt, mc_colors=mc_non_t_nk@colors, mc_groups=nt_col2group[mc_non_t_nk@colors], is_mc_t_nk=NULL, mat_id=tumor_non_t_nk_id, scfigs_fn(tumor_non_t_nk_id, "mc_metadata"), pats_order=rownames(global_pat_grp_norm)[order(global_pat_grp_norm[, 'dysfunctional'])], groups_order=names(ord_by_id[[tumor_non_t_nk_id]]), show_is_t=F, ncells_ylim=c(0, max_mcsz))
	
	# Fig S4C - T composition of treatment naive
	treatment_naive_pats = mat_t_nk@cell_metadata[names(mc_t_nk@mc), ] %>% filter(treatment=='N') %>% dplyr::select(PatientID) %>% unique
	plot_patients_group_comp(tumor_t_nk_id, "treatment_naive", filter_groups=setdiff(t_nms, 'NK'), patients=treatment_naive_pats$PatientID, sort_pats_by='dysfunctional', min_cells=100) 
	
	# Fig S4D - pairs of metastases side by side
	plot_patients_group_comp(tumor_t_nk_id, "metastases_pairs", filter_groups=setdiff(t_nms, 'NK'), patients=c('p12-3-(S)C-N-1', 'p12-3-(S)C-N-2', 'p17-3-(S)C-N-1', 'p17-3-(S)C-N-2')) 
	
	nt_genes_ann = annotate_genes(rownames(lfp_nt))
	tfs = rownames(nt_genes_ann)[grepl("TF", nt_genes_ann$type)]
	
	# genes cor heatmap
	message("gen myeloid genes heatmap..")
	
	mono_mac_dc_genes = setdiff(names(which(apply(abs(lfp_nt[, f_myl_mc]), 1, max) > 3)), lat_genes)
	mono_mac_dc_genes = grep("^RP", mono_mac_dc_genes, invert=T, v=T)
	
	nt_g_ann = data.frame(row.names=mono_mac_dc_genes, type=nt_genes_ann[mono_mac_dc_genes, 'type'])
	
	myl_cor_c = tgs_cor(t(as.matrix(mat_non_t_nk_ds[mono_mac_dc_genes, f_myl_c])), spearman=T)
	myl_cor_mc = cor(t(lfp_nt[mono_mac_dc_genes, f_myl_mc]))
	
	# myl gene-gene cor heatmaps
	blwtrd_cols = colorRampPalette(c('blue', 'white', 'red'))(101)
	.plot_start(scfigs_fn(tumor_non_t_nk_id, "fig3_myl_c_genes_cor"), max(700, 300 + length(mono_mac_dc_genes) * 12), max(700, 300 + length(mono_mac_dc_genes) * 12))
	pheatmap(pmin(pmax(myl_cor_c, -0.3), 0.3), clustering_method="ward.D2", cutree_rows=10, cutree_cols=10, treeheight_col=0, treeheight_row=0, cellwidth=10, cellheight=10, fontsize=12, col=blwtrd_cols, show_colnames=F) #annotation_row=g_ann, annotation_colors=list(type=ggroup_cols), 
	dev.off()
	
	.plot_start(scfigs_fn(tumor_non_t_nk_id, "fig3_myl_mc_genes_cor"), max(700, 300 + length(mono_mac_dc_genes) * 12), max(700, 300 + length(mono_mac_dc_genes) * 12))
	pheatmap(pmin(pmax(myl_cor_mc, -1), 1), clustering_method="ward.D2", cutree_rows=10, cutree_cols=10, treeheight_col=0, treeheight_row=0, cellwidth=10, cellheight=10, fontsize=12, col=blwtrd_cols, show_colnames=F) #, annotation_row=g_ann, annotation_colors=list(type=ggroup_cols)
	dev.off()
	
	##### 
	# monocyte - DC - macrophage gradients
	myl_groups = c('macrophage', 'monocyte', 'DC')
	myl_de = gene_diff_expr_by_groups(tumor_non_t_nk_id, tumor_non_t_nk_id, lateral_gset_id, n_ds=500, mat_ds=mat_non_t_nk_ds, groups=myl_groups, filter_outlier_genes=F, vgel_group_ord=ord_non_t_nms, min_log2_enr=2, blist_genes=c(lat_genes, tfs))
	
	myl_tot = colSums(mat_non_t_nk@mat[, colnames(myl_de$totu)])
	
	.plot_start(scfigs_fn(tumor_non_t_nk_id, "fig_S3_myeloid_signatures_f_umis_on_cells"), 1200, 400)
	par(mfrow=c(1,3))
	par(mar=c(5,5,5,5))
	cg = split(colnames(myl_de$totu), nt_col2group[mc_non_t_nk@colors[mc_non_t_nk@mc[colnames(myl_de$totu)]]])
	for (g in myl_groups) {
		plot(myl_de$totu['monocyte', ]/myl_tot- myl_de$totu['DC',]/myl_tot, myl_de$totu['macrophage',]/myl_tot - myl_de$totu['DC',]/myl_tot, pch=19, cex=0.7, cex.main=2, cex.lab=2, cex.axis=2, col='lightgrey', xlab='monocyte - DC', ylab='macrophage - DC', main=g)
		ctotu = myl_de$totu[, cg[[g]]]
		cmyl_tot = myl_tot[cg[[g]]]
		points(ctotu['monocyte',] / cmyl_tot - ctotu['DC',] / cmyl_tot, ctotu['macrophage',]/ cmyl_tot - ctotu['DC',] / cmyl_tot, col=nt_group2col[g], pch=19, cex=0.7) 
	}
	dev.off()
	
	myl_genes = list()
	myl_scores = list()
	f_myl_g_mc = list()
	
	for (g in myl_groups) {
		myl_genes[[g]] = myl_de$de[[g]]$gene
		myl_scores[[g]] = colMeans(lfp_nt[myl_genes[[g]], ])
		f_myl_g_mc[[g]] =  mc_non_t_nk@colors == nt_group2col[g]
		
		for (i in seq_along(myl_genes[[g]])) { 
			gene = myl_genes[[g]][i]
			plt(paste(g, "score"), gene, lfp_nt[, f_myl_mc], cols=mc_non_t_nk@colors[f_myl_mc], x=myl_scores[[g]][f_myl_mc], show_mc_ids=F, add_grid=T, ofn=scfigs_fn(tumor_non_t_nk_id, sprintf("%s_grad_lfp_%d_%s", g, i, gene), scfigs_dir(tumor_non_t_nk_id, "fig_s3_myloid_grads")), cex.lab=1.5)
		}
	}
	myl_x = myl_scores[['monocyte']] - myl_scores[['DC']]
	myl_y = myl_scores[['macrophage']] - myl_scores[['DC']]
	plt("monocyte - DC", "macrophage - DC", lfp_nt[, f_myl_mc], cols=mc_non_t_nk@colors[f_myl_mc], x=myl_x[f_myl_mc], y=myl_y[f_myl_mc], show_mc_ids=T, add_grid=T, ofn=scfigs_fn(tumor_non_t_nk_id, "fig_s3_myeloid_signatures_f_umis_on_myl_mcs"), cex.lab=1.5)
	plt("monocyte - DC", "macrophage - DC", lfp_nt, cols=mc_non_t_nk@colors, x=myl_x, y=myl_y, show_mc_ids=T, add_grid=T, ofn=scfigs_fn(tumor_non_t_nk_id, "fig_s3_myeloid_signatures_f_umis_on_all_mcs"), cex.lab=1.5)
	
	# model monoc-mac-DC gradients with TFs
	nt_tf_models = list()
	for (g in myl_groups) {
		nt_tf_models[[g]] = model_by_tfs(tumor_non_t_nk_id, myl_scores[[g]], NULL, nt_genes_ann[rownames(lfp_nt), ], f=f_myl_g_mc[[g]], nm1=g, nm2="None", grp1=g, grp2=NULL, nfolds=min(sum(f_myl_g_mc[[g]]), 10), max_genes_for_plot_size=30)
	}
	mono_vs_mac_model = model_by_tfs(tumor_non_t_nk_id, myl_scores[['monocyte']], myl_scores[['macrophage']], nt_genes_ann[rownames(lfp_nt), ], f = f_myl_g_mc[['macrophage']] | f_myl_g_mc[['monocyte']], nm1='monocyte', nm2='macrophage', grp1='monocyte', grp2='macrophage', max_genes_for_plot_size = 30)
	
	# lateral modules signatures on patients
	ifn_gset = scdb_gset("mel_ifn_filt")
	ifn_genes = names(ifn_gset@gene_set)
	
	stress_gset = scdb_gset("mel_stress_filt")
	hs_sets = unique(stress_gset@gene_set [ grepl("HSP", names(stress_gset@gene_set))])
	er_sets = unique(stress_gset@gene_set [ grepl("FOS|JUN", names(stress_gset@gene_set))])
	stopifnot(length(intersect(hs_sets, er_sets)) == 0)
	hs_genes = names(stress_gset@gene_set)[stress_gset@gene_set %in% hs_sets]
	er_genes = names(stress_gset@gene_set)[stress_gset@gene_set %in% er_sets]
	
	med_ifn = plot_f_tot_genes_by_group(c(tumor_t_nk_id, tumor_non_t_nk_id), ifn_genes, "IFN", "lightblue")
	med_hs = plot_f_tot_genes_by_group(c(tumor_t_nk_id, tumor_non_t_nk_id), hs_genes, "heat_shock", "orange")
	med_er = plot_f_tot_genes_by_group(c(tumor_t_nk_id, tumor_non_t_nk_id), er_genes, "JUN_FOS", "darkgreen")
	
	# Fig 3A: 2d proj per patient T and non-T, only for patients with enough cells
	mc_t_nk_pats = table(mc_t_nk@mc, mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'])
	mc_non_t_nk_pats = table(mc_non_t_nk@mc, mat_non_t_nk@cell_metadata[names(mc_non_t_nk@mc), 'PatientID'])
	min_cells_for_pat_spec_2d_plot = 500
	stopifnot(all(colnames(mc_t_nk_pats) == colnames(mc_non_t_nk_pats)))
	big_pats = names(which(colSums(mc_t_nk_pats) >= min_cells_for_pat_spec_2d_plot & colSums(mc_non_t_nk_pats) >= min_cells_for_pat_spec_2d_plot))
	
	prev_params = override_metacell_params(list(mcell_mc2d_plot_key=F, mcell_mc2d_height=500, mcell_mc2d_width=500, mcell_mc2d_cex=1))
	mcell_mc2d_plot_by_factor(tumor_t_nk_id, tumor_t_nk_id, 'PatientID', single_plot=F, filter_values=big_pats, neto_points=T)
	mcell_mc2d_plot_by_factor(tumor_non_t_nk_id, tumor_non_t_nk_id, 'PatientID', single_plot=F, filter_values=big_pats, neto_points=T)
	restore_metacell_params(prev_params)
	
	# monocyte strats
	monoc_genes = myl_genes[['monocyte']]
	f_monoc_mc = f_myl_g_mc[['monocyte']]
	f_monoc_cc = f_monoc_mc[mc_non_t_nk@mc]
	f_monoc_neto = Matrix::colSums(mat_non_t_nk@mat[monoc_genes, f_monoc_cc]) / Matrix::colSums(mat_non_t_nk@mat[, f_monoc_cc])
	monoc_neto_cut = cut(f_monoc_neto, breaks=quantile(f_monoc_neto, seq(0, 1, by=0.2)), include.lowest=T, labels=paste0("Q", seq(0.1, 1, by=0.2)))
	monoc_neto_strats = as.character(monoc_neto_cut)
	names(monoc_neto_strats) = names(f_monoc_neto)
	
	for (g in c('LYZ', 'MXD1', 'MAF', 'CEBPB', 'KLF4', 'PTPRE', 'NFKBIA')) {
		f_genes_by_strats_with_chisq_pval(tumor_non_t_nk_id, g, monoc_neto_strats, paste0(g, "_on_monocytes"), nt_group2col['monocyte'])
	}
	f_genes_by_strats_with_chisq_pval(tumor_non_t_nk_id, monoc_genes, monoc_neto_strats, "MonocGenes_on_monocytes", nt_group2col['monocyte'])
	
	###############
	#
	# TCR analysis
	#
	###############
	tcr = read.csv("data/TCR_pipeline_result.csv", header=T, stringsAsFactors=F)

	# update clone id by CDR3 nucleotide instead of AA sequence	
	tcr$prev_clone_id = tcr$clone_id
	fac_clones = ordered(tcr$CDR3_first, levels=unique(tcr$CDR3_first))
	tcr$clone_id = as.numeric(fac_clones)

	tcr_cells = intersect(tcr$Well_ID, names(mc_t_nk@mc))
	rownames(tcr) = tcr$Well_ID
	
	t_tcr = tcr[tcr_cells, ]
	t_tcr$mc = mc_t_nk@mc[tcr_cells]
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "clone_sizes_barplot"), 300, 300)
	barplot(log2(table(pmin(table(t_tcr$clone_id), 12))), xlab='clone size', ylab='count (log2)')
	dev.off()
	
	t_tcr$PatientID = mat_t_nk@cell_metadata[t_tcr$Well_ID, 'PatientID']
	t_tcr$UnifPatientID = gsub("-.*", "", t_tcr$PatientID)
	t_tcr$mc_grp = t_col2group[mc_t_nk@colors[mc_t_nk@mc[t_tcr$Well_ID]]]
	t_tcr$seq_batch = mat_t_nk@cell_metadata[t_tcr$Well_ID, 'seq_batch_id']
	
	# pr TCR detection by expression of TCR genes (TRAC + TRBC2)
	tcr_sc = Matrix::colSums(mat_t_nk@mat[c('TRAC', 'TRBC2'), ]) / Matrix::colSums(mat_t_nk@mat)
	tcr_sc_bin = cut(tcr_sc, breaks=c(0, seq(min(tcr_sc[tcr_sc > 0]), quantile(tcr_sc, 0.98), length=8), max(tcr_sc)), include.lowest=T)
	
	f_t = mc_t_nk@colors[mc_t_nk@mc] != t_group2col['NK']
	
	p_tcr_t = tapply(names(tcr_sc)[f_t] %in% rownames(t_tcr), tcr_sc_bin[f_t], mean)
	p_tcr_a = tapply(names(tcr_sc) %in% rownames(t_tcr), tcr_sc_bin, mean)
	
	p_tcr_grp_t = table(tcr_sc_bin[f_t], t_col2group[mc_t_nk@colors[mc_t_nk@mc[f_t]]])
	p_tcr_grp_a = table(tcr_sc_bin, t_col2group[mc_t_nk@colors[mc_t_nk@mc]])
	
	p_tcr_grp_t_n = p_tcr_grp_t / rowSums((p_tcr_grp_t))
	p_tcr_grp_a_n = p_tcr_grp_a / rowSums((p_tcr_grp_a))
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "p_tcr_detect"), 600, 500)
	layout(matrix(1:2, nrow=1), widths=c(3,1))
	par(mar=c(4,4,4,3))
	t_neto_nms = setdiff(t_nms, 'NK')
	barplot(t(p_tcr_grp_t_n[, t_neto_nms]), col=t_group2col[t_neto_nms], las=2, main='T cells: % TRAC + TRBC2', ylab='fraction', names.arg=100*as.numeric(gsub("]", "", gsub(".*,", "", rownames(p_tcr_grp_t_n)))))
	points(seq(0.7, length=length(p_tcr_t), by=1.2), p_tcr_t, pch=19, cex=2)
	axis(4, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2), las=2)
	par(mar=c(4,1,4,1))
	plot(1:5, col=NA, xlab="", ylab="", xaxt='n', yaxt='n', bty='n')
	legend("topleft", legend="Pr(detect)", pch=19, pt.cex=2, bty='n', col='black')
	dev.off()
	
	if ('NK' %in% t_tcr$mc_grp) {
		.plot_start(scfigs_fn(tumor_t_nk_id, "p_tcr_detect_with_NK"), 600, 500)
		layout(matrix(1:2, nrow=1), widths=c(3,1))
		par(mar=c(4,4,4,3))
		barplot(t(p_tcr_grp_a_n[, names(ord_t_nk_nms)]), col=t_group2col[names(ord_t_nk_nms)], las=2, main='T/NK cells: % TRAC + TRBC2', ylab='fraction', names.arg=100*as.numeric(gsub("]", "", gsub(".*,", "", rownames(p_tcr_grp_a_n)))))
		points(seq(0.7, length=length(p_tcr_a), by=1.2), p_tcr_a, pch=19, cex=2)
		axis(4, at=seq(0, 1, by=0.2), labels=seq(0, 1, by=0.2), las=2)
		par(mar=c(4,1,4,1))
		plot(1:5, col=NA, xlab="", ylab="", xaxt='n', yaxt='n', bty='n')
		legend("topleft", legend="Pr(detect)", pch=19, pt.cex=2, bty='n', col='black')
		dev.off()
	}
	
	# clones by patients
	# detecting jumper cells (same well, seq batch, clone, different amp_batch) which are likely cross-plate contamination - and remove them
	jumper_cells = t_tcr %>% 
		group_by(clone_id, cell_name) %>% 
		summarize(n_batches=length(unique(Amp.Batch)), n_seq_batches=length(unique(seq_batch)), n_patients=length(unique(PatientID)), min_reads=min(reads_freq), max_reads=max(reads_freq)) %>% filter(n_batches > 1) %>% data.frame
	
	write.table(jumper_cells, sprintf("%s/jumper_cells.txt", .scfigs_base), quote=F, sep="\t")
	t_tcr = t_tcr %>% anti_join(jumper_cells, by=c('clone_id', 'cell_name')) %>% data.frame
	rownames(t_tcr) = t_tcr$Well_ID
	
	valid_clones = t_tcr %>% group_by(clone_id) %>% summarise(clone_size = length(Well_ID)) %>% filter(clone_size > 1) %>% data.frame 
	t_tcr_ge2 = t_tcr %>% inner_join(valid_clones, by="clone_id")
	t_tcr_singles = t_tcr[setdiff(t_tcr$Well_ID, t_tcr_ge2$Well_ID), ]
	
	# mc mc enrichment - looking at all pairs of cells sharing a clone
	mc_pairs = do.call('rbind', lapply(split(t_tcr_ge2$mc, t_tcr_ge2$clone_id), function(x) { t(combn(x, 2)) }))
	mc_mc = table(c(seq_along(mc_t_nk@colors), mc_pairs[,1]), c(seq_along(mc_t_nk@colors), mc_pairs[,2]))
	diag(mc_mc) = diag(mc_mc) - 1
	mc_mc = mc_mc + t(mc_mc)
	
	cln_shr_ignore_grps = c('NK', 'dysf-cd4', 'em-cd4')
	e_mc_mc = shuffle_clones_mc_pairs_by_patient(tumor_t_nk_id, t_tcr, n_iters=200, ignore_groups=cln_shr_ignore_grps)
	stopifnot(rownames(mc_mc)== rownames(e_mc_mc) & colnames(mc_mc) == colnames(e_mc_mc))
	
	e_mc_mc = (sum(mc_mc) / sum(e_mc_mc)) * e_mc_mc 
	mc_mc_enr = log2( (mc_mc + 1e-3) / (e_mc_mc + 1e-3))
	
	mc_ann = data.frame(row.names=seq_along(mc_t_nk@colors), group=t_col2group[mc_t_nk@colors])
	hc = hclust(dist(cor(mc_mc_enr)), method='ward.D2')
	.plot_start(scfigs_fn(tumor_t_nk_id, "mc_mc_clone_sharing_norm_by_patient_shuffle"), 200+nrow(mc_mc) * 6, 100+nrow(mc_mc)*6)
	pheatmap(pmin(pmax(mc_mc_enr[hc$order, hc$order], -4), 4),  annotation_row=mc_ann, annotation_col=mc_ann, annotation_colors=list(group=t_group2col), cluster_rows=F, cluster_cols=F, cellwidth=6, cellheight=6)
	dev.off()
	
	# collapse mc-mc pairs counts on groups
	grp_f = t_col2group[mc_t_nk@colors]
	o_mc_mc_grp = metacell:::.row_stats_by_factor(t(metacell:::.row_stats_by_factor(mc_mc, grp_f, rowSums)), grp_f, rowSums)
	
	e_mc_mc_grp = metacell:::.row_stats_by_factor(t(metacell:::.row_stats_by_factor(e_mc_mc, grp_f, rowSums)), grp_f, rowSums)
	
	grp_reg=5
	full_o_mc_mc_grp = o_mc_mc_grp
	o_mc_mc_grp = o_mc_mc_grp[rowSums(o_mc_mc_grp) > 0 & !(rownames(o_mc_mc_grp) %in% cln_shr_ignore_grps), colSums(o_mc_mc_grp) > 0 & !(colnames(o_mc_mc_grp) %in% cln_shr_ignore_grps)]
	e_mc_mc_grp = e_mc_mc_grp[rownames(o_mc_mc_grp), colnames(o_mc_mc_grp)]
	enr_mc_mc_grp = log2((o_mc_mc_grp+grp_reg) / (e_mc_mc_grp+grp_reg))
	grp_mc_mc_hc = hclust(dist(cor(enr_mc_mc_grp)), method='ward.D2')
	grp_ann = data.frame(row.names=rownames(o_mc_mc_grp), group = rownames(o_mc_mc_grp))
	
	.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("mc_mc_clone_sharing_norm_group_agg_by_patient_shuffle_reg%d", grp_reg)), 200+nrow(o_mc_mc_grp)*25, 100+nrow(o_mc_mc_grp)*25)
	enr_mc_mc_grp_o = enr_mc_mc_grp[grp_mc_mc_hc$order, grp_mc_mc_hc$order]
	pheatmap(pmin(pmax(enr_mc_mc_grp_o, -3), 3), breaks=seq(-3, 3, length=100), display_numbers=round(enr_mc_mc_grp_o, 1), annotation_row=grp_ann, annotation_col=grp_ann, annotation_colors=list(group=t_group2col), cluster_rows=F, cluster_cols=F, cellwidth=25, cellheight=25)
	dev.off()
	
	for (min_grp_size_in_cln in 1:2) {
		aa = do.call("rbind", lapply(split(t_tcr_ge2$mc_grp, t_tcr_ge2$clone_id), function(grps) { ugrps = names(which(table(grps)>=min_grp_size_in_cln)); ugrpsi = expand.grid(ugrps, ugrps); cbind(as.character(ugrpsi[,1]), as.character(ugrpsi[,2])) }))
		aa_t = table(aa[,1], aa[,2])
		aa_full = matrix(0, nrow(full_o_mc_mc_grp), ncol(full_o_mc_mc_grp), dimnames=list(rownames(full_o_mc_mc_grp), colnames(full_o_mc_mc_grp)))
		aa_full[rownames(aa_t), colnames(aa_t)] = aa_t
		write.table(aa_full, sprintf("%s/cln_grp_sharing_min%d_grp_sz_in_cln.txt", .scfigs_base, min_grp_size_in_cln), quote=F, sep="\t")
	}
	
	.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("mc_mc_clone_sharing_norm_group_agg_by_patient_shuffle_reg%d_tree", grp_reg)), 300, 200)
	par(mar=c(8,4,1,1))
	hcd = as.dendrogram(grp_mc_mc_hc)
	plot(hcd)
	dev.off()
	
	# clones mc group composition barplots
	cln_grp = table(t_tcr$clone_id, t_tcr$mc_grp)
	mid_clns_ord = plot_cln_group_composition(tumor_t_nk_id, t_tcr, cln_grp, c(8,20), group2col=t_group2col)
	top_clns_ord = plot_cln_group_composition(tumor_t_nk_id, t_tcr, cln_grp, c(21,max(rowSums(cln_grp))), group2col=t_group2col)
	
	# clone size by patients colored by dominant group 
	pat_grp = table(mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'], t_col2group[mc_t_nk@colors[mc_t_nk@mc]])
	pat_grp_n = pat_grp / rowSums(pat_grp)
	pats_by_dysf = rownames(pat_grp_n)[order(pat_grp_n[, 'dysfunctional'])]
	clones = t_tcr %>% group_by(clone_id) %>% summarize(patient=names(tail(sort(table(PatientID)), 1)), clone_size=length(clone_id), max_group=names(tail(sort(table(mc_grp)), 1))) %>% data.frame
	o = ordered(clones$patient, levels=intersect(pats_by_dysf, unique(clones$patient)))
	gnx = rnorm(nrow(clones), 0, 0.1)
	gny = rnorm(nrow(clones), 0, 0.1)
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "clone_size_by_patient_ord_by_f_dysf"), 600, 300)
	par(mar=c(10,4,1,1))
	plot(as.numeric(o)+gnx, log2(clones$clone_size)+gny, pch=19, col=t_group2col[clones$max_group], cex=0.9, xlab="", xaxt='n', ylab='clone size (log2) + noise')
	axis(1, labels=levels(o), at=seq_along(levels(o)), las=2)
	dev.off()
	
	# clones shared across patients
	cln_pat = table(t_tcr$PatientID, t_tcr$clone_id)
	
	dup_cln_pat = cln_pat[ , colSums(cln_pat > 0) > 1]
	dup_cln_pat = dup_cln_pat[rowSums(dup_cln_pat) > 0, ]
	
	l_dup_cln_pat = log2(1+dup_cln_pat)
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "clones_across_patients_hc"), 1600, 800)
	hc = hclust(dist(l_dup_cln_pat), method='ward.D2')
	pheatmap(l_dup_cln_pat[hc$order, ], color=colorRampPalette(c('white', colorRampPalette(c('blue', 'red', 'black'))(6)))(101), number_format="%d", treeheight_col=0, treeheight_row=0, fontsize_row=11, fontsize_col=5)
	dev.off()
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "clones_across_patients"), 1600, 800)
	pheatmap(l_dup_cln_pat, cluster_rows=F, color=colorRampPalette(c('white', colorRampPalette(c('blue', 'red', 'black'))(6)))(101), number_format="%d", treeheight_col=0, treeheight_row=0, fontsize_row=11, fontsize_col=5)
	dev.off()

	clns_pat = table(t_tcr$clone_id, t_tcr$PatientID)
	t_pat_grp = table(mat_t_nk@cell_metadata[names(mc_t_nk@mc), 'PatientID'], t_col2group[mc_t_nk@colors[mc_t_nk@mc]])
	tumor_pats = intersect(all_tumor_pats, names(which(rowSums(t_pat_grp) >= min_t_nk_cells_per_patient)))
	
	t_pat_grp_n = t_pat_grp / rowSums(t_pat_grp)
	all_pats_by_dysf = rownames(t_pat_grp_n)[order(t_pat_grp_n[, 'dysfunctional'])]
	pats_by_dysf = intersect(all_pats_by_dysf, tumor_pats)
	
	
	# Patient clone-size composition and type compostion by clone size
	fig4_cd = function(max_clnsz_for_comp = 3, min_clones_per_pat = 0, min_cells_per_pat_clnsz=10, ds_tcrs=F, ds_n=NULL ) {
		
		if (ds_tcrs) {
			tot_pats = table(t_tcr$PatientID)
			if (is.null(ds_n)) {
				ds_n = median(tot_pats)
			}
			pat_inds = split(1:nrow(t_tcr), t_tcr$PatientID)
			t_tcr_ds = do.call("rbind", lapply(names(which(tot_pats >= ds_n)), function(pat) t_tcr[ sample(pat_inds[[pat]], ds_n, replace=F), ] ))
			clns_pat = table(t_tcr_ds$clone_id, t_tcr_ds$PatientID)
			
		} else {
			clns_pat = table(t_tcr$clone_id, t_tcr$PatientID)
		}
		
		pat_cln_comp = matrix(0, ncol(clns_pat), max_clnsz_for_comp+1, dimnames=list(colnames(clns_pat), 0:max_clnsz_for_comp))
		for (pat in colnames(clns_pat)) {
			tab = table(pmin(clns_pat[, pat], max_clnsz_for_comp))
			pat_cln_comp[pat, as.character(names(tab))] = tab
		}
		pat_cln_comp = pat_cln_comp[, -1]
		pat_cln_comp = pat_cln_comp[rowSums(pat_cln_comp) >= min_clones_per_pat, ]
		pat_cln_comp_n = pat_cln_comp / rowSums(pat_cln_comp)
		
		pats_ord = rownames(pat_cln_comp_n)[order(pat_cln_comp_n[, "1"], decreasing = T)]
		
		pat_cln_comp = pat_cln_comp[pats_ord, ]
		pat_cln_comp_n = pat_cln_comp_n[pats_ord, ]
		clns_pat = clns_pat[, pats_ord]
		
		.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("fig4a_pats_comp%s_maxClSz%d_minClperPat%d_minCellsperPatClnSz%d.png", ifelse(ds_tcrs, sprintf("_ds_to_%d", ds_n), ""), max_clnsz_for_comp, min_clones_per_pat, min_cells_per_pat_clnsz)), 200 + nrow(pat_cln_comp_n)*30, 330 + 50 * max_clnsz_for_comp)
		layout(matrix(seq(max_clnsz_for_comp+3, 1), ncol=1), heights=c(rep(50, max_clnsz_for_comp), 40, 40, 250))
		
		# breakdown of clone sizes by patient
		par(mar=c(15,8,0.5,1))
		clnsz_cols = RColorBrewer::brewer.pal(n=max_clnsz_for_comp, 'Greens')
		barplot(t(pat_cln_comp_n), ylab="", col=clnsz_cols, las=2, cex.names=1.5, cex.axis = 1.5)
		
		# N tcrs per patient
		par(mar=c(0.5,8,1,1))
		tot_tcrs = colSums(clns_pat)
		barplot(tot_tcrs[pats_ord], xaxt='n', yaxt='n', ylab="#TCRs", col='lightblue4', cex.lab=1.5)
		yaxp = par("yaxp")
		axis(2, yaxp=c(yaxp[1], yaxp[2], 2), las=2, cex.axis=1.2)
		
		# N clones per patient
		tot_clns = rowSums(pat_cln_comp)
		barplot(tot_clns[pats_ord], xaxt='n', yaxt='n', ylab="#clones", col='lightblue', cex.lab=1.5)
		yaxp = par("yaxp")
		axis(2, yaxp=c(yaxp[1], yaxp[2], 2), las=2, cex.axis=1.2)
		
		# pie plots: breakdown to groups by clone size
		par(mar=c(0.5,8,0.5,1))
		for (clnsz in 1:max_clnsz_for_comp) {
			barplot(rowSums(pat_cln_comp), xaxt='n', yaxt='n', ylab="", xlab="", col='white', border=NA, ylim=c(0, 1), cex.lab=1.5)
			for (i in 1:nrow(pat_cln_comp)) {
				pat = rownames(pat_cln_comp)[i]
				c_clns = rownames(clns_pat)[pmin(clns_pat[, pat] , max_clnsz_for_comp) == clnsz]
				c_mc_grp = t_tcr[t_tcr$PatientID == pat & t_tcr$clone_id %in% c_clns, 'mc_grp']
				if (length(c_mc_grp) >= min_cells_per_pat_clnsz) {
					gtab = table(c_mc_grp)
					plotrix::floating.pie(0.6 + (i-1) * 1.2, 0.5, gtab, radius=0.4, col=t_group2col[names(gtab)])
				}
			}
			rect(-0.7, 0.25, -0.2, 0.75, col=clnsz_cols[clnsz])
			mtext(paste("|cln|", ifelse(clnsz == max_clnsz_for_comp, ">=", "="), clnsz), 2, line=0.5, at=0.5, las=2, cex=1)
		}
		dev.off()
		
		# stats by f_dysf strats
		for (i in 1:max_clnsz_for_comp) {
			.plot_stat_by_f_grp_helper(tumor_t_nk_id, pat_cln_comp_n[,i], paste0(ifelse(ds_tcrs, sprintf("ds_to_%d_", ds_n), ""), "f_clns_size_", i), clnsz_cols[i], odir=.scfigs_base)
		}
		.plot_stat_by_f_grp_helper(tumor_t_nk_id, apply(clns_pat, 2, clonality_from_cln_counts_vec), paste0(ifelse(ds_tcrs, sprintf("ds_to_%d_", ds_n), ""), 'clonality'), 'grey', odir=.scfigs_base, ylab='clonality')

		pat_cln_comp_n
	}
	pat_cln_comp_n = fig4_cd()
	
	# mc fraction size 1 clones vs >1 clones
	f_tcr_eq1_mc = tapply(names(mc_t_nk@mc) %in% t_tcr_singles$Well_ID, mc_t_nk@mc, mean)
	f_tcr_ge2_mc = tapply(names(mc_t_nk@mc) %in% t_tcr_ge2$Well_ID, mc_t_nk@mc, mean)
	.plot_start(scfigs_fn(tumor_t_nk_id, "fig4e_f_mc_by_clone_size"), 300, 300)
	par(mar=c(4,4,1,1))
	plot(f_tcr_eq1_mc, f_tcr_ge2_mc, pch=21, bg=mc_t_nk@colors, xlab='% cells from n=1 clones', ylab='% cells from n>=2 clones', cex=1.5)
	abline(a=0, b=1)
	dev.off()
	
	# clones composition
	for (sort_by in c('hclust', 'patient_type')) {
		mid_clns_ord = plot_cln_group_composition(tumor_t_nk_id, t_tcr, cln_grp, c(8,20), group2col=t_group2col, sort_by=sort_by, odir=.scfigs_base)
		top_clns_ord = plot_cln_group_composition(tumor_t_nk_id, t_tcr, cln_grp, c(21,max(rowSums(cln_grp))), group2col=t_group2col, sort_by=sort_by, odir=.scfigs_base)
	}
	
	# proliferation
	f_dysf_mc = mc_t_nk@colors == t_group2col['dysfunctional']
	f_dysf_cc = f_dysf_mc[mc_t_nk@mc]
	
	f_dysf_neto = Matrix::colSums(mat_t_nk@mat[names(dysf_genes), f_dysf_cc]) / Matrix::colSums(mat_t_nk@mat[, f_dysf_cc])
	dysf_neto_cut = cut(f_dysf_neto, breaks=quantile(f_dysf_neto, seq(0, 1, by=0.1)), include.lowest=T, labels=paste0("Q", seq(0.1, 1, by=0.1)))
	dysf_neto_strats = as.character(dysf_neto_cut)
	names(dysf_neto_strats) = names(f_dysf_neto)
	
	f_eff1_mc = mc_t_nk@colors == t_group2col['effector1']
	f_eff1_cc = f_eff1_mc[mc_t_nk@mc]
	f_eff1_neto = Matrix::colSums(mat_t_nk@mat[names(dysf_genes), f_eff1_cc]) / Matrix::colSums(mat_t_nk@mat[, f_eff1_cc])
	eff1_neto_cut = cut(f_eff1_neto, breaks=quantile(f_eff1_neto, seq(0, 1, by=0.2)), include.lowest=T, labels=paste0("Q", seq(0.2, 1, by=0.2)))
	eff1_neto_strats = as.character(eff1_neto_cut)
	names(eff1_neto_strats) = names(f_eff1_neto)
	
	f_eff1_dysf_mc = mc_t_nk@colors == t_group2col['dysfunctional'] | mc_t_nk@colors == t_group2col['effector1']
	f_eff1_dysf_cc = f_eff1_dysf_mc[mc_t_nk@mc]
	f_eff1_dysf = Matrix::colSums(mat_t_nk@mat[names(dysf_genes), f_eff1_dysf_cc]) / Matrix::colSums(mat_t_nk@mat[, f_eff1_dysf_cc])
	eff1_dysf_cut = cut(f_eff1_dysf, breaks=quantile(f_eff1_dysf, seq(0, 1, by=0.1)), include.lowest=T, labels=paste0("Q", seq(0.1, 1, by=0.1)))
	eff1_dysf_strats = as.character(eff1_dysf_cut)
	names(eff1_dysf_strats) = names(f_eff1_dysf)
	
	tt = table(t_col2group[mc_t_nk@colors[mc_t_nk@mc[f_eff1_dysf_cc]]], eff1_dysf_strats)
	tt_n = t(t(tt) / colSums(tt))
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "transitional_dysf_stratas_comp"), 300, 300)
	par(mar=c(8,4,4,1))
	barplot(tt_n, las=2, col=t_group2col[rownames(tt_n)])
	dev.off()
	
	f_treg = Matrix::colSums(mat_t_nk@mat[names(treg_genes), f_treg_c]) / Matrix::colSums(mat_t_nk@mat[, f_treg_c])
	treg_cut = cut(f_treg, breaks=quantile(f_treg, seq(0, 1, by=0.2)), include.lowest=T, labels=paste0("Q", seq(0.2, 1, by=0.2)))
	treg_strats = as.character(treg_cut)
	names(treg_strats) = names(f_treg)
	
	cc_gset = scdb_gset("mel_cc_filt")
	cc_genes = names(cc_gset@gene_set)
	
	t_cc_res = mel_mc_prolif_score(mc_id=tumor_t_nk_id, f_cc_cutoff = f_cc_cutoff, pats_order=rownames(pat_cln_comp_n), strat_by=list(dysfunctional=dysf_neto_strats, Treg=treg_strats, effector1=eff1_neto_strats))
	nt_cc_res = mel_mc_prolif_score(mc_id=tumor_non_t_nk_id, f_cc_cutoff = f_cc_cutoff, pats_order=rownames(pat_cln_comp_n))
	
	# unified T/ non-T %prolif
	min_cells_per_group_for_f_prolif=200
	f_cc_by_grp = c(t_cc_res$f_cc_by_grp, nt_cc_res$f_cc_by_grp)
	tot_cc_by_grp = c(t_cc_res$tot_cc_by_grp, nt_cc_res$tot_cc_by_grp)
	grp_ind = tot_cc_by_grp >= min_cells_per_group_for_f_prolif
	f_cc_by_grp = sort(f_cc_by_grp[grp_ind], decreasing = T)
	
	.plot_start(scfigs_fn(tumor_id, sprintf("f_prolif_by_grp_ge%d_cells", min_cells_per_group_for_f_prolif)), 700, 300)
	par(mar=c(10,4,4,1))
	barplot(f_cc_by_grp, col=global_group2col[names(f_cc_by_grp)], las=2, ylab='% prolif cells')
	mtext(tot_cc_by_grp[names(f_cc_by_grp)], 3, at=seq(0.6, by=1.2, length=length(f_cc_by_grp)), las=2, line=0.5)
	dev.off()
	
	# heatmap of strats lfp high var genes
	strats_lfp_heatmap = function(strats, group='dysfunctional', fp_reg=0.1, min_lfp_diff=0.7, add_genes='TCF7') 
	{
		mc_t_nk = scdb_mc(tumor_t_nk_id)
		lfp_t = log2(mc_t_nk@mc_fp)
		
		geom = metacell:::.row_stats_by_factor(mat_t_nk@mat[rownames(lfp_t), names(strats)], strats, function(y) {exp(Matrix::rowMeans(log(1+y)))-1})
		geom_mean = tapply(Matrix::colSums(mat_t_nk@mat[, names(strats)]), strats, mean)
		geom_n = 1000 * t(t(geom) / as.vector(geom_mean))
		strats_lfp = log2((fp_reg+geom_n)/apply(fp_reg+geom_n, 1, median))
		
		f = !grepl("^RPL|^RPS", rownames(strats_lfp), perl=T) & apply(strats_lfp, 1, function(v) { max(v) - min(v) } ) >= min_lfp_diff
		if (!is.null(add_genes)) {
			f[add_genes] = T
		}
		strats_lfp_f = strats_lfp[f, ]
		o = order(apply(strats_lfp_f, 1, which.max) + 1e-3 * apply(strats_lfp_f, 1, max))
		
		.plot_start(scfigs_fn(tumor_t_nk_id, sprintf("genes_var_on_%s_strats", group)), 350, 100 + sum(f)*10)
		pheatmap(pmin(pmax(strats_lfp_f[o, ], -1), 1), cluster_rows=F, cluster_cols=F, cellwidth=8, cellheight=6, fontsize_row=7, fontsize_col=9, gaps_row=cumsum(rle(apply(strats_lfp_f[o, ], 1, which.max))$lengths))
		dev.off()
	}
	
	strats_lfp_heatmap(dysf_neto_strats)
	strats_lfp_heatmap(treg_strats, group='Treg')
	
	# %genes on dysf strats with chisq-pval
	for (g in c('TIGIT', 'TOX2', 'ID3', 'CD8A')) {
		f_genes_by_strats_with_chisq_pval(tumor_t_nk_id, g, dysf_neto_strats, paste0(g, "_on_dysf"), t_group2col['dysfunctional'])
	}
	f_genes_by_strats_with_chisq_pval(tumor_t_nk_id, names(dysf_genes), dysf_neto_strats, "DysfGenes_on_dysf", t_group2col['dysfunctional'])

	# compare prolif vs not (diff expresion) per dysf strata (to rule out that the prolif weakly dysf cells are actually not dysf)
	t_f_cc = t_cc_res$f_cc
	for (qs in sort(unique(dysf_neto_strats))) {
		dysf_curr_q_cells = names(which(dysf_neto_strats == qs))
		pos_c = names(which(t_f_cc[dysf_curr_q_cells] >= f_cc_cutoff))
		neg_c = names(which(t_f_cc[dysf_curr_q_cells] < f_cc_cutoff))
		dysf_qs_de = diff_expr(mc_t_nk, mat_t_nk_ds, NULL, NULL, nms1=pos_c, nms2=neg_c)
		de_f = dysf_qs_de %>% filter(abs(enr) > 1)
		stopifnot(nrow(de_f) <= 40)
		png(scfigs_fn(tumor_t_nk_id, paste0("dysf_strat_de_by_prolif_", qs)), 200, 100 + 40 * 15)
		par(mar=c(4,8,4,1))
		barplot(de_f$enr, horiz=T, names.arg=de_f$gene, las=2, col=ifelse(de_f$gene %in% cc_genes, 'grey', ifelse(de_f$enr > 0, 'blue', 'red')), main=sprintf("dysf %s prolif vs not", qs), ylim=c(0, 40))
		dev.off()
		message(sprintf("%s: %.2f prolif cells", qs, mean(t_f_cc[dysf_curr_q_cells] >= f_cc_cutoff)))
	}

	is_prolif = t_f_cc[names(dysf_neto_strats)] >= f_cc_cutoff
	
	tot_dysf_cyto = data.frame(row.names=names(mc_t_nk@mc), tot_cyto=Matrix::colSums(mat_t_nk@mat[names(cyto_genes), names(mc_t_nk@mc)]), tot_dysf=Matrix::colSums(mat_t_nk@mat[names(dysf_genes), names(mc_t_nk@mc)]), group=t_col2group[mc_t_nk@colors[mc_t_nk@mc]], stringsAsFactors=F)
	tot_dysf_cyto[names(which(is_prolif)), 'group'] = 'dysf-prolif'
	tot_dysf_cyto[names(which(!is_prolif)), 'group'] = 'dysf-non-prolif'
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "sum_dysf_by_grp_breaking_dysf_by_prolif"), 400, 300)
	par(mar=c(10,4,1,1))
	boxplot(log2(1+tot_dysf) ~ group, tot_dysf_cyto, las=2, notch=T, ylab="sum DYSF umis (log2")
	dev.off()
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "sum_cyto_by_grp_breaking_dysf_by_prolif"), 400, 300)
	par(mar=c(10,4,1,1))
	boxplot(log2(1+tot_cyto) ~ group, tot_dysf_cyto, las=2, notch=T, ylab="sum cyto umis (log2")
	dev.off()
	
	# clonality per patient by reactivity
	clonality = apply(clns_pat, 2, clonality_from_cln_counts_vec)
	.plot_start(scfigs_fn(tumor_t_nk_id, "pat_clonality_by_reac"), w=180, h=400)
	par(mar=c(8,4,1,1))
	stripchart(clonality ~ reac_col2type[ pat_reac_col[names(clonality)]], cex=1.5, las=2, pch=21, bg='grey', vertical=T, method='jitter', ylab='clonality')
	dev.off()

	#################################
	#
	# FACS antibody index sort plots
	#
	#################################
	plot_ab_pair_by_mc_group = function(mat_id, mc_id, panel, ab1, ab2, group_ord=NULL, cex=0.7)
	{
		mat = scdb_mat(mat_id)
		stopifnot(!is.null(mat))
		
		mc = scdb_mc(mc_id)
		stopifnot(!is.null(mc))
		
		fdf =  mat@cell_metadata[names(mc@mc), ]
		
		c_by_p = split(rownames(fdf), fdf$PatientID)
		c_by_p[['all']] = rownames(fdf)
		
		for (nm in names(c_by_p)) {
			df = fdf[c_by_p[[nm]], ]
			if (nrow(df) > 0) {
				df = df[!is.na(df[,'Staining panel'] ) & df[, 'Staining panel'] == panel, paste0(c(ab1, ab2), "_Ab")]
				df = df[rowSums(is.na(df)) == 0, ]
				if (nrow(df) > 0) {
					df = .logicle_transform(df[rowSums(is.na(df)) == 0, ])
					
					col2group = as.character(mc@color_key$group)
					names(col2group) = as.character(mc@color_key$color)
					
					group2col = as.character(mc@color_key$color)
					names(group2col) = as.character(mc@color_key$group)
					
					gcells = split(names(mc@mc), col2group[mc@colors[mc@mc]])
					
					if (is.null(group_ord)) {
						group_ord = sort(names(gcells))
					}
					
					.plot_start(scfigs_fn(mc_id, sprintf("panel_%d_%s_%s_vs_%s", panel, nm, ab1, ab2)), 250 * length(group_ord), 300)
					layout(matrix(seq_along(group_ord), nrow=1))
					
					smooth =F 
					message(sprintf("panel %d, %s vs %s, total %d cells", panel, ab1, ab2, nrow(df)))

					for (group in group_ord) {
						
						curr_cells = intersect(gcells[[group]], rownames(df))
						if (smooth) {
							smoothScatter(df[curr_cells, 1], df[curr_cells, 2], xlab=ab1, ylab=ab2, main=sprintf("%d: %s (%d)", panel, group, length(gcells[[group]])), xlim=range(df[,1]), ylim=range(df[,2]))
						}
						else {
							plot(df[, 1], df[, 2], xlab=ab1, ylab=ab2, pch=19, cex=cex, cex.main=4*cex, cex.lab=3*cex, col=ifelse(length(curr_cells) > 0, 'darkgray', NA), main=sprintf("%d: %s (%d)", panel, group, length(curr_cells)))
							if (length(curr_cells) > 0) {
								points(df[curr_cells, 1], df[curr_cells, 2], pch=19, cex=cex, col=group2col[group])
								if (length(curr_cells) > 50) {
									d = MASS::kde2d(df[curr_cells, 1], df[curr_cells, 2], n=100)
									contour(d$x, d$y, d$z, drawlabels=F, add=T)
								}
							}
						}
					}
					dev.off()	
				}
			}
		}
	}
	
	fmat_id = all_facs_id
	mc_all_id = tumor_clean_id #"mel_filt_clean_all_sep_nk_broad"
	mc_t_nk_id = tumor_t_nk_id #"mel_filt_clean_all_broad_submc_T_NK" 
	
	mcell_add_mars_facs_data(fmat_id, all_id, "data/facs_index")
	
	fmat = scdb_mat(fmat_id)
	t_nk_nms = names(ord_by_id[[mc_t_nk_id]])
	
	
	for (spanel in c(2, 6, 9, 14, 15)) {
		plot_ab_pair_by_mc_group(fmat_id, mc_all_id, spanel, "CD3", "CD45", group_ord=all_nms)
		plot_ab_pair_by_mc_group(fmat_id, mc_t_nk_id, spanel, "CD3", "CD45", group_ord=t_nk_nms)
	}
	for (spanel in c(6,9,14,15)) {
		for(ab2 in c('CD103', 'PD1')) {
			plot_ab_pair_by_mc_group(fmat_id, mc_t_nk_id, spanel, "CD8", ab2, group_ord=t_nk_nms)
		}
	}
	for (ab2 in c('CD56', 'CD11b', 'PDL1')) {
		plot_ab_pair_by_mc_group(fmat_id, mc_t_nk_id, 2, "CD8", ab2, group_ord=t_nk_nms)
	}
	plot_ab_pair_by_mc_group(fmat_id, mc_t_nk_id, 6, "CD8", 'CD137', group_ord=t_nk_nms)
	
	for (spanel in c(2, 6, 9, 14)) {
		plot_ab_pair_by_mc_group(fmat_id, mc_t_nk_id, spanel, "CD4", "CD8", group_ord=t_nk_nms)
	}
	
	plot_ab_pair_by_mc_group(fmat_id, mc_all_id, 15, "CD45", "CD45RA", group_ord=all_nms)
	for (ab2 in c('PD1', 'CD103', 'DAPI', 'CCR7')) {
		plot_ab_pair_by_mc_group(fmat_id, mc_t_nk_id, 15, "CD8", ab2, group_ord=t_nk_nms)
	}
	
	valid_spanels = c(2,6,9)
	md = fmat@cell_metadata[names(mc_t_nk@mc), ]
	md$sp = md[, 'Staining panel']
	md = md[md$sp %in% valid_spanels & !is.na(md$CD4_Ab) & !is.na(md$CD8_Ab), ]
	
	vn = .logicle_transform(md[, c('CD4_Ab', 'CD8_Ab')])
	colnames(vn) = c('CD4_n', 'CD8_n')
	md = cbind(md, vn)
	cd4_cutoffs = list("2"=1.8, "6"=1.4, "9"=1.6)
	cd8_cutoffs = list("2"=2, "6"=0.8, "9"=0.8)
	md$facs_type = NA
	for (sp in valid_spanels) {
		md[md$sp == sp & md$CD4_n >= cd4_cutoffs[[as.character(sp)]] & md$CD8_n < cd8_cutoffs[[as.character(sp)]], 'facs_type'] = 'CD4'
		md[md$sp == sp & md$CD4_n < cd4_cutoffs[[as.character(sp)]] & md$CD8_n >= cd8_cutoffs[[as.character(sp)]], 'facs_type'] = 'CD8'
		md[md$sp == sp & md$CD4_n < cd4_cutoffs[[as.character(sp)]] & md$CD8_n < cd8_cutoffs[[as.character(sp)]], 'facs_type'] = 'DN'
	}
	
	facs_cols = c('red', 'orange', 'lightseagreen')
	names(facs_cols) = c('CD4', 'CD8', 'DN')
	
	png(scfigs_fn(tumor_t_nk_id, "CD4_vs_CD8_by_panel", scfigs_dir(tumor_t_nk_id, "facs_idx")), w=900, h=300)
	par(mfrow=c(1,3))
	par(mar=c(3,3,3,3))
	ccols = facs_cols[md$facs_type]
	ccols[is.na(ccols)] = 'lightgray'
	for (sp in valid_spanels) {
		plot(md[md$sp == sp, 'CD4_n'], md[md$sp == sp, 'CD8_n'], pch=19, cex=0.5, col=ccols[md$sp == sp], main=paste("panel", sp))
		abline(v=cd4_cutoffs[[as.character(sp)]], lty=2)
		abline(h=cd8_cutoffs[[as.character(sp)]], lty=2) 
	}
	dev.off()
	
	cells = intersect(rownames(md), names(mc_t_nk@mc))
	facs_cells = split(cells, md[cells, 'facs_type'])
	
	mc2d = scdb_mc2d(tumor_t_nk_id)
	
	png(scfigs_fn(tumor_t_nk_id, "2d_proj_by_facs_type", scfigs_dir(tumor_t_nk_id, "facs_idx")), w=900, h=300)
	par(mfrow=c(1,3))
	par(mar=c(1,1,3,1))
	for (type in names(facs_cells)) {
		plot(mc2d@sc_x, mc2d@sc_y, pch=19, cex=0.5, col='lightgray', main=type)
		points(mc2d@sc_x[facs_cells[[type]]], mc2d@sc_y[facs_cells[[type]]], pch=19, cex=0.5, col=facs_cols[type]) }
	dev.off()
	
	#######################
	#
	# Additional supp figs
	#
	#######################
	fmd = fmat@cell_metadata[names(mc_a@mc), ]
	fmd$sp = fmd[, 'Staining panel']
		
	# S1 - FACS on selected patients
	plot_facs_by_patient_wrap(fmd, patients=c("p2-4-LN-2IT", "p11-3-(S)C-N", "p27-4-(S)C-1IT"))
	
	# S1c - MC size distributions
	mcsz = table(mc_a@mc)
	.plot_start(scfigs_fn(tumor_clean_id, "S1C_mc_size_distr"), 450, 450)
	hist(mcsz, 30, col='navyblue', xlab="#cells in metacell", ylab="#metacells")
	dev.off()
	
	# S1d - %ribo UMIs vs tot UMIs
	mc_f_ribo = tapply(Matrix::colSums(mat_all@mat[grep("^RPS|^RPL", mat_all@genes, perl=T, v=T), names(mc_a@mc)]), mc_a@mc, sum) / tapply(Matrix::colSums(mat_all@mat[, names(mc_a@mc)]), mc_a@mc, sum)
	mc_mu_tot_umi = tapply(Matrix::colSums(mat_all@mat[, names(mc_a@mc)]), mc_a@mc, mean)
	.plot_start(scfigs_fn(tumor_clean_id, "S1D_mc_f_ribo_vs_tot_umi"), 450, 450)
	plot(log2(mc_mu_tot_umi), mc_f_ribo, pch=21, bg=mc_a@colors, xlab="mc mean tot UMIs (log2)", ylab="mc %ribo UMIs", cex=2)
	dev.off()
	
	# S1f FACS CD4 vs CD8 on T subsets
	nms = setdiff(names(ord_by_id[[tumor_t_nk_id]]), 'em-cd4')
	names(nms) = nms
	plot_facs_by_patient_wrap(fmd, mc_id=tumor_t_nk_id, patients=c("p4-4-LN-1IT", "p9-4-(S)C-2IT", "p13-3-(S)C-N"), ab_x="CD4", ab_y="CD8", groups=as.list(nms), fig_pref="S1F")
	
	# S1g - compare naive CD4 to naive CD8 (by FACS)
	naive_cells = names(mc_t_nk@mc)[which(mc_t_nk@colors[mc_t_nk@mc] == t_group2col['naive'])]
	naive_cd4 = intersect(naive_cells, facs_cells[['CD4']])
	naive_cd4_m = rowMeans(mat_t_nk@mat[, naive_cd4])
	
	naive_cd8 = intersect(naive_cells, facs_cells[['CD8']])
	naive_cd8_m = rowMeans(mat_t_nk@mat[, naive_cd8])
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "S1G_naive_cd4_vs_cd8_mu_UMIs"), 450, 450)
	plot(log2(naive_cd4_m), log2(naive_cd8_m), pch=19, cex=0.5, col='navyblue', main=sprintf("R^2 = %.2f", cor(naive_cd4_m, naive_cd8_m) **2), xlab=sprintf('mean naive CD4 (%d cells) UMIs (log2)', length(naive_cd4)), ylab=sprintf('mean naive CD8 (%d cells) UMIs (log2)', length(naive_cd8)))
	yy = abs(log2((1e-2+naive_cd4_m) / (1e-2+naive_cd8_m)))
	f = yy > 3 | (yy > 1 & log2(pmin(naive_cd4_m, naive_cd8_m)) > -4)
	text(log2(naive_cd4_m[f]), log2(naive_cd8_m[f]), names(naive_cd8_m)[f], cex=0.7, pos=ifelse(naive_cd4_m[f] > naive_cd8_m[f], 4, 2))
	dev.off()
	
	# S1h - diff expr naive vs rest
	non_naive_cells = setdiff(names(mc_t_nk@mc), naive_cells)
	non_naive_cells = setdiff(non_naive_cells, names(mc_t_nk@mc)[which(mc_t_nk@colors[mc_t_nk@mc] == t_group2col['NK'])])
	naive_vs_rest = diff_expr(mc_t_nk, mat_t_nk_ds, NULL, NULL, nms1=naive_cells, nms2=non_naive_cells, min_max_umi=200)
	yy = naive_vs_rest %>% filter(enr > 1)
	yy$lf_tot1 = log2(yy$tot1 / length(naive_cells))
	.plot_start(scfigs_fn(tumor_t_nk_id, "S1H_naive_vs_rest_T"), 450, 450)
	plot(yy$lf_tot1, yy$enr, pch=19, cex=0.5, col='navyblue', xlab='mean naive UMIs (log2)', ylab='naive enr (log2)', main='Naive vs rest of T')
	text(yy$lf_tot1, yy$enr, yy$gene, pos=ifelse(yy$lf_tot1 > -2, 2, 4), cex=0.7)
	dev.off()
	
	#S1I: lfp of dysf/CD4 metacells (scatter of log total vs. lfp?)
	.plot_start(scfigs_fn(tumor_t_nk_id, "S1I_dysf_cd4_top_lfp_genes_barplot"), 200, 900)
	par(mar=c(4,8,1,1))
	top_dysf_cd4 = tail(sort(lfp_t[grep("^RPL|^RPS", rownames(lfp_t), perl=T, invert=T), mc_t_nk@colors == t_group2col['dysf-cd4']]), 70)
	barplot(top_dysf_cd4, horiz=T, col=t_group2col['dysf-cd4'], yaxt='n', xlab='enr (log2)')
	mtext(names(top_dysf_cd4), 2, at=seq(0.6, len=length(top_dysf_cd4), by=1.2), las=2, line=0.5)
	dev.off()
	
	#S1J: lfp of memory CD4 metacells (scatter of log total vs. lfp?)
	.plot_start(scfigs_fn(tumor_t_nk_id, "S1I_em_cd4_top_mean_lfp_genes_barplot"), 200, 900)
	par(mar=c(4,8,1,1))
	top_em_cd4 = tail(sort(rowMeans(lfp_t[grep("^RPL|^RPS", rownames(lfp_t), perl=T, invert=T), mc_t_nk@colors == t_group2col['em-cd4']])), 70)
	barplot(top_em_cd4, horiz=T, col=t_group2col['em-cd4'], yaxt='n', xlab='enr (log2)')
	mtext(names(top_em_cd4), 2, at=seq(0.6, len=length(top_em_cd4), by=1.2), las=2, line=0.5)
	dev.off()
		
	# S3B - marker heat map for non-T
	mel_basic_mc_mc2d_plots(mc_id=tumor_non_t_nk_merged_id, mat_id=tumor_id, graph_id=tumor_id, lateral_gset_id = lateral_gset_id)
	
	
	# S5B
	mcsz = table(mc_t_nk@mc)
	mcsz2 = table(t_tcr$mc)
	f_tcr_mc = as.vector(mcsz2/mcsz)
	plt("CD8 - CD4", "%TCR cov", lfp_t, mc_t_nk@colors, x=lfp_t['CD8A', ] + lfp_t['CD8B', ] - lfp_t['CD4', ], y=f_tcr_mc, ofn=scfigs_fn(tumor_t_nk_id, "f_tcr_vs_CD8_CD4_diff"))
	
	# S5C -Cell cycle gene module gene cor heatmap
	cc_gset_id = "mel_cc_filt"
	f_cc_reg=1e-3
	cc_gset = scdb_gset(cc_gset_id)
	cc_genes = names(cc_gset@gene_set)
	.plot_start(scfigs_fn(tumor_t_nk_id, "cc_genes_cor_mc"), 150 + length(cc_genes)*12, 150+length(cc_genes)*12)
	cc_c = cor(t(lfp_t[cc_genes, ]))
	hc = hclust(dist(cor(cc_c)), method='ward.D2')
	diag(cc_c) = NA
	pheatmap(cc_c[hc$order, hc$order], cluster_rows=F, cluster_cols=F, breaks=seq(-1, 1, length=100), cellwidth=12, cellheight=12)
	dev.off()
	
	# S5B - distribution of CC score on all T cells
	f_cc = colSums(mat_t_nk@mat[cc_genes, ]) / colSums(mat_t_nk@mat)
	.plot_start(scfigs_fn(tumor_t_nk_id, "l_f_cc"), 300, 300)
	plot(density(log2(f_cc + f_cc_reg)), xlab="%prolif UMIs (log2)", lwd=2, col='navyblue', main="")
	abline(v=log2(f_cc_cutoff + f_cc_reg), lty=2, lwd=2)
	dev.off()

	
	# S5E - distribution of prolif score on top prolif bin vs. naive
	lf_cc_by_group = split(log2(f_cc + f_cc_reg), t_col2group[mc_t_nk@colors[mc_t_nk@mc]])
	
	.plot_start(scfigs_fn(tumor_t_nk_id, "lf_cc_ecdf_all_groups"), 400, 400)
	plot(ecdf(-lf_cc_by_group[['dysfunctional']]), do.points=F, lwd=2, col=t_group2col['dysfunctional'], xlab="- %prolif UMIs (log2)", ylab="ecdf", main="", xlim=c(3, 6), ylim=c(0, 0.15))
	for (g in names(lf_cc_by_group)) {
		lines(ecdf(-lf_cc_by_group[[g]]), do.points=F, lwd=2, col=t_group2col[g])
	}
	abline(v=-log2(f_cc_cutoff + f_cc_reg), lty=2, lwd=2)
	legend("topleft", legend=paste0(sprintf("%s (%.1f", names(lf_cc_by_group), 100*sapply(lf_cc_by_group, function(v) mean(v >= log2(f_cc_reg + f_cc_cutoff)))), "%)"), fill=t_group2col[names(lf_cc_by_group)], bty='n')
	dev.off()
	
	# compare cross-metacell cells membership between the metacells in the paper and the metacell object generated by build_metacells
	compare_mc_partitioning("mel_filt_Tumor_outClean_detailed_colors", "mel_filt_Tumor_new_outClean")
	
	# PBMC plots
	fig_s_pbmc_plots()
}

