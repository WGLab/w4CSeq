# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>
#
#	This demo draw chromosome ideogram with padding between chromosomes, highlights, 
#	chromosome names, and link lines. 
#
#	Usage:
#
#	library(RCircos);
#	demo("RCircos.Link.Plot.Demo");
# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>





	#	Load RCircos library
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	library(RCircos);


	#	Load human cytoband data and link data
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	data(RCircos.Link.Data);
	data(UCSC.HG19.Human.CytoBandIdeogram);
	cyto.info <- UCSC.HG19.Human.CytoBandIdeogram;


	#	Setup RCircos core components:
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	RCircos.Set.Core.Components(cyto.info, NULL, 10, 0);


	#	Open the graphic device (here a pdf file)
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	out.file <- "RCircos.Link.Plot.Demo.pdf";
	pdf(file=out.file, height=8, width=8);

	RCircos.Set.Plot.Area();



	#	Draw chromosome ideogram
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Draw chromosome ideogram ...\n");

	RCircos.Chromosome.Ideogram.Plot();
	title("RCircos Link Plot Demo");


	#	Link lines. Link data has only paired chromosome locations in
	#	each row and link lines are always drawn inside of chromosome 
	#	ideogram.
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Add link track ...\n");
	
	link.data <- RCircos.Link.Data;
	link.colors <- rep("blue", nrow(link.data));
	rows <- seq(1, nrow(link.data), by=5);
	link.colors[rows] <- "red";
	link.data["PlotColor"] <- link.colors;

	track.num <- 2;
	RCircos.Link.Plot(link.data, track.num, FALSE);


	#	Close the graphic device and clear memory
	#  	_________________________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	dev.off();	
	print("RCircos Link Plot Demo Done!");

	rm(list=ls(all=T));

