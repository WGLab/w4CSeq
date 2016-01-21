# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>
#
#	This demo draw chromosome ideogram with padding between chromosomes, highlights, 
#	chromosome names, and tile plot. 
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Tile.Plot.Demo");
# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>



	#	Load RCircos library
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	library(RCircos);


	#	Load human cytoband data and gene expression data
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data(RCircos.Tile.Data);
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


	#	Setup RCircos core components:
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


	#	Open the graphic device (here a pdf file)
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	out.file <- "RCircos.Tile.Plot.Demo.pdf";
	pdf(file=out.file, height=8, width=8);

	RCircos.Set.Plot.Area();


	#	Draw chromosome ideogram
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Draw chromosome ideogram ...\n");

	RCircos.Chromosome.Ideogram.Plot();
	title("RCircos Tile Plot Demo");


	#	Tile plot. Note: tile plot data have chromosome locations only and
	#	each data file is for one track
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Add tile track ...\n");

	tile.data <- RCircos.Tile.Data;
	tile.colors <- rainbow(nrow(tile.data));
	tile.data["PlotColor"] <- tile.colors;

	track.num <- 9;
	RCircos.Tile.Plot(tile.data, track.num, "in");



	#	Close the graphic device and clear memory
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	dev.off();
	print("RCircos Tile Plot Demo Done!");

	rm(list=ls(all=T));
