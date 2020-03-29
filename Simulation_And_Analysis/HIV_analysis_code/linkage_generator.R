#Set start, end,p=0.05, gene_info, dpis
start =0
end =11
p=0.05
gene_info = list(c("Gag", "Pol", "Vif", "Vpr", "Tat", "Rev", "Vpu", "Env"), 
	c(262, 1551, 4507, 5025, 5296, 5435, 5529, 5692), 
	c(1758, 4562, 5085, 5315, 5510, 5510, 5774, 8274))
dpis = c(122, 562, 1084, 1254, 1282, 1393, 1861, 2303, 2578, 2639, 2922, 2996)
robj = rehh_maker(dir,start,end,reflen)
window_sizes = make_ehhs_graphs(robj$snp_pos, robj$hapobjs, dpis, egraph, dir)
make_window_graphs(window_sizes, gene_info,dir)
stat_tables = make_stat_tables(start, end, robj$hapobjs, p, xgraph)
cri_table = make_cri_table(robj$snp_pos, dpis, stat_tables, gene_info)
crx_table = make_crx_table(robj$snp_pos, dpis, stat_tables, gene_info)

make_all_haps = function(start, end, dir, reflen) {
	h = list()
	for (i in start:end) { # For each time-point, generate a haplotype object. 
		h[[i]] = data2haplohh(hap_file = paste(dir,i-1,".ms.txt", sep=""), position_scaling_factor = reflen) # Convert ms file to a haplohh object. 
	} 
	return(h)
}

make_ehhs_df = function(snp_pos, i, hapobjs, dpis) {
	df = data.frame(Position = snp_pos)
	index = match(snp_pos[i],hapobjs[[1]]@positions)	
	for (j in 1:length(hapobjs)) { # For each time-point...
		curr_ehhs = calc_ehhs(hapobjs[[j]], mrk = index, include_nhaplo = FALSE)
		df[toString(dpis[j])] = curr_ehhs$ehhs$EHHS
	}
	return(df)
}

make_ehhs_graphs = function(snp_pos, hapobjs, dpis, egraph, dir) {
	w = data.frame(DPI = dpis)
	for (i in 1:length(snp_pos)) { # For each SNP position... length(snp_pos)
		curr_ehhs = make_ehhs_df(snp_pos, i, hapobjs, dpis) # Make a master table of the EHHSes at each time-point.
		transform = melt(curr_ehhs, id.vars="Position")
		pos = round(snp_pos[i], digits = 0)
		w[toString(dpis[i])] = window_sizing(curr_ehhs,snp_pos)
		if (egraph) { # The goal is to make figures as well as the other thing.
			imgpath = paste(dir,"Position_",pos,"_EHHS",".png", sep="")
			g = ggplot(transform, aes(Position,value, col=variable, group=variable)) + geom_line() # + ylim(0,1000)
			g = g + scale_color_viridis(discrete = TRUE, name = "DPI")
			g = g + geom_vline(xintercept = pos, linetype="dashed", color = "black", size=1)
			g = g + ylab("EHHS") + xlab("Position (bp)") # + ggtitle(graph_name)
			ggsave(imgpath, plot = g, scale = 1, dpi = 1000, limitsize = TRUE)
		}
		print(paste(pos,"finished", sep=" "))
	}
	return(w)
}

# Fix legend label, remove title, change time-point to days
make_window_graphs = function(window_sizes, gene_info, dir) {
	transform = melt(window_sizes, id.vars="DPI")
	transform$variable = as.numeric(levels(transform$variable))[transform$variable] 
	for (i in 1:length(gene_info[[1]])) {
		imgpath = paste(dir,"/",gene_info[[1]][i],"_EHHS_WSize2",".pdf", sep="")
		pdf(file = pdf_path)
		gene_windows = subset(transform, variable >= gene_info[[2]][i] & variable <= gene_info[[3]][i])
		g = ggplot(gene_windows, aes(Time_Points,value, col=variable, group=variable)) + geom_line() # + ylim(0,1000)
		g = g + scale_color_viridis(name = "Focal SNP Position")
		g = g + ylab("EHHS Window Size (bp)") + xlab("DPI") # + ggtitle(graph_name)
		ggsave(imgpath, plot = g, scale = 1, dpi = 1000, limitsize = TRUE)
		print(g)
		dev.off()

	}
}

window_sizing = function(ehhs_table,snp_pos) {
# print(ehhs_table)
	w = integer()
	for (i in 2:(ncol(ehhs_table))) { # For each time point except the first...
		start = 1
		end = 1
		scan_to_end = FALSE
		scan_from_start = TRUE
		for (j in 1:(nrow(ehhs_table))) { # For each position take the first and last indices such that EHHS > 0.5. 
			if ((ehhs_table[j,i] >= 0.5) & (scan_from_start)) {
				start = round(snp_pos[j], digits = 0)
				scan_to_end = TRUE
				scan_from_start = FALSE
			}
			if ((ehhs_table[j,i] < 0.5) & (scan_to_end)) {
				end = round(snp_pos[j], digits = 0)
				scan_to_end = FALSE
				scan_from_start = TRUE
			}
		}
		if ((start != 1) & (end == 1)) {
	end = round(snp_pos[length(snp_pos)], digits = 0)
		}
		size = end - start + 1 # Subtract and save width.
		# print(paste(start,end,sep=" "))
		w = c(w, size)
	}
	return(w) # Add the position-average of a window across all time-points to a list
}

make_stat_tables = function(start, end, hapobjs, p, xgraph) { # window_sizes, 
	ihs_tables = list()
	# cri_tables = list() 
	cri_tables = list() # Originally called cr1
	xpehh_tables = list()
	crx_tables = list()
	thresh = abs(log10(as.numeric(p)))
	pdf(file = "XPEHH_Graphs.pdf")
	for (i in start:end) { # For each time-point, generate iHS. 
		ihh = scan_hh(hapobjs[[i]], threads=4, limehh = 0, discard_integration_at_border = FALSE) # Calculate the iHH and iHS for each time-point.
		ihs = ihh2ihs(ihh, freqbin = 0.1, min_maf = 0)
		# cri_tables[[i]] = calc_candidate_regions(ihs, threshold = thresh, pval = TRUE, window_size = window_sizes[i], min_n_extr_mrk = 1)
		cri_tables[[i]] = calc_candidate_regions(ihs, threshold = thresh, pval = TRUE, window_size = 1, min_n_extr_mrk = 1)
		if (i > 1) {
			xpehh = ies2xpehh(ihh,prev_ihh,popname1 = i, popname2 = (i - 1))
			xpehh_tables[[i]] = xpehh
			crx = calc_candidate_regions(xpehh, threshold = thresh, pval = TRUE, window_size = 1, min_n_extr_mrk = 1)
			crx_tables[[i]] = crx
			if (xgraph) {
			manhattanplot(xpehh, pval = TRUE, threshold = thresh, cr = crx)
			}
		}
		ihs_tables[[i]] = ihs
		prev_ihh = ihh
	}
	results = list(ihs = ihs_tables, cri = cri_tables, xpehh = xpehh_tables, crx = crx_tables) # cri = cri_tables, windows = window_sizes, 
	dev.off()
	return(results)	
}

make_cri_table = function(snp_pos, dpis, stat_tables, gene_info) {
	df = data.frame(matrix(ncol = 6, nrow = 0))
	names = c("Time Point", "SNP No.", "Genomic Pos.", "Gene", "iHS", "p-value") # "Window Size", "Tot. Markers", "XP-EHH", "p-value"
	colnames(df) = names
	for (i in 1:length(dpis)) { # For each time-point...
		time_point = dpis[i]
		crs = nrow(stat_tables$cri[[i]])
		if (crs > 0) {
			for (j in 1:nrow(stat_tables$cri[[i]])) { # For row of each candidate window...
				gpos = stat_tables$cri[i][[1]]$START[[j]]
				gene = identify_gene(gpos, gene_info[[1]], gene_info[[2]], gene_info[[3]])
				# wsize = stat_tables$cri[i][[1]]$END[[j]] - stat_tables$cri[i][[1]]$START[[j]] + 1
				# nmrk = stat_tables$cri[i][[1]]$N_MRK[[j]]
				
				index = match(gpos,stat_tables$ihs[[i]]$ihs$POSITION)	
				ihs = stat_tables$ihs[[i]]$ihs[index,3]
				p1 = 10^(-stat_tables$cri[i][[1]]$EXTR_MRK[[j]])
				# index = match(gpos,stat_tables$xpehh[[i]]$POSITION)	
				# xph = stat_tables$xpehh[[i]][index,3]
				# p2 = 10^(-stat_tables$xpehh[[i]][index,4])
				row = c(time_point, index, round(gpos, digits = 0), gene, ihs, p1) # wsize, nmrk, xph, p2
				df[i-1,] = row
			}
		}
	}
	return(df)
}

make_crx_table = function(snp_pos, dpis, stat_tables, gene_info) {
	df = data.frame(matrix(ncol = 6, nrow = 0))
	names = c("Time Point", "SNP No.", "Genomic Pos.", "Gene", "XP-EHH", "p-value") # "Window Size", "Tot. Markers", "iHS", "p-value", 
	colnames(df) = names
	for (i in 2:length(dpis)) { # For each time-point...
		time_point = dpis[i]
		crs = nrow(stat_tables$crx[[i]])
		if (crs > 0) {
			for (j in 1:nrow(stat_tables$crx[[i]])) { # For row of each candidate window...
				gpos = stat_tables$crx[i][[1]]$START[[j]]
				gene = identify_gene(gpos, gene_info[[1]], gene_info[[2]], gene_info[[3]])
				# wsize = stat_tables$cri[i][[1]]$END[[j]] - stat_tables$cri[i][[1]]$START[[j]] + 1
				# nmrk = stat_tables$crx[i][[1]]$N_MRK[[j]]
				
				# index = match(gpos,stat_tables$ihs[[i]]$ihs$POSITION)	
				# ihs = stat_tables$ihs[[i]]$ihs[index,3]
				# p1 = 10^(-stat_tables$cri[i][[1]]$MAX_MRK[[j]])
				index = match(gpos,stat_tables$xpehh[[i]]$POSITION)	
				xph = stat_tables$xpehh[[i]][index,3]
				p2 = 10^(-stat_tables$crx[i][[1]]$EXTR_MRK[[j]])
				row = c(time_point, index, round(gpos, digits = 0), gene, xph, p2) # wsize, nmrk, xph, p2
				df[i-1,] = row
			}
		}
	}
	return(df)
}

identify_gene = function(pos, names, starts, ends) {
	genes = ""
	first = TRUE
	for (i in 1:length(names)) {
		if ((pos >= starts[i]) && (pos <= ends[i])){
			if (first) {
				genes = names[i]
				first = FALSE
			} else {
				paste(genes,names[i], sep = ",", collapse = NULL)
			}
		}
	}
	if (genes == '') {
	genes = "NA"
	}
	return(genes)
}



rehh_maker = function(dir,start,end,reflen) {
	library(rehh)
	library(gap)
	library(ggplot2)
	library(reshape2)
	library(viridis)
	hapobjs = make_all_haps(start, end, dir, reflen) # A list of rehh haplotype objects for each time-point.
	#print (hapobjs)
	snp_pos = hapobjs[[1]]@positions # A vector of position names. 
	egraph = FALSE
	xgraph = TRUE

	results = list(hapobjs = hapobjs, snp_pos = snp_pos, dpis = dpis)
	return(results)
}
