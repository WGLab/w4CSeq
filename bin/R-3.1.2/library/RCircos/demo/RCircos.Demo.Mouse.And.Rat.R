# ________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>
#
#	This demo draws a heatmap with both mouse and rat gene expresion data 
#	and links between same genes. 
#
#	We will put mouse data on right side and rat data on left side. Since we are 
#	going to draw the image with both mouse and rat chromosome ideogram, 
#	we deduct total points of each ideogram by half. This is done by modify the
#	base.per.unit from 3000 to 30000. Chromosome padding will be automatically
#	modified if we do not give it new value.
#
#	Last revised on May 8, 2013
#
#_________________________________________________________________________________________
# <><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><RCircos DEMO><>




	#	Load RCircos package and defined parameters. 
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
	
	library(RCircos);	


	#	Load chromosome cytoband data. We will generate 
	#	a pseudo karyotype by concatenating the second 
	#	ideogram data to the first one.
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Read in cytoband data ...\n");

	data(UCSC.Mouse.GRCm38.CytoBandIdeogram);
	data(UCSC.Baylor.3.4.Rat.cytoBandIdeogram);


	#	Load mouse and rat gene expression data
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Read in gene expression data ...\n");

	data(RCircos.Mouse.Expr.Data);
	data(RCircos.Rat.Expr.Data);


	#	Generate link data from the two datasets. These
	#	lines are regular R script code.
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	mouse.data <- RCircos.Mouse.Expr.Data;
	rat.data <- RCircos.Rat.Expr.Data;

	mouse.gene <- as.character(mouse.data$Gene);
	rat.gene   <- as.character(rat.data$Gene);

	mouse.set <- mouse.data[mouse.gene %in% rat.gene,];
	mouse.set <- mouse.set[order(mouse.set$Expr.Mean),];
	last <- nrow(mouse.set);
	first <- last - 49;
	mouse.set <- mouse.set[first:last,];

	mouse.set.gene <- as.character(mouse.set$Gene);
	rat.set <- rat.data[rat.gene %in% mouse.set.gene,];

	mouse.set <- mouse.set[order(mouse.set$Gene),];
	rat.set <- rat.set[order(rat.set$Gene),];

	link.data <- data.frame(mouse.set[,1:3], rat.set[,1:3]);


	#	Modify the datasets and ideogram so that they have
	#	same chromosome names. Please make sure the items
	#	in the two lists have same order.
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	species.list <- c("M", "R");

	cyto.list <- list(UCSC.Mouse.GRCm38.CytoBandIdeogram, 
			UCSC.Baylor.3.4.Rat.cytoBandIdeogram);
	RCircos.Multiple.Species.Core.Components(cyto.list, 
			species.list, NULL, 5, 0);

	data.list <- list(RCircos.Mouse.Expr.Data, 
			RCircos.Rat.Expr.Data);
	expr.data <- RCircos.Multiple.Species.Dataset(data.list, 
			species.list);

	link.data[,1] <- paste(species.list[1], link.data[,1], sep="");
	link.data[,4] <- paste(species.list[2], link.data[,4], sep="");


	#	Setup RCircos core components then change the
	#	base.per.unit to 30000 for fast drawing
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	params <- RCircos.Get.Plot.Parameters();
	params$base.per.unit <- 30000;
	RCircos.Reset.Plot.Parameters(params);



	#	Initialize graphic device (GUI or image file)
	#
	#	tiff(file="RCircos.Demo.Mouse.And.Rat.tif", height=8,
	#		width=9, unit="in", type="cairo", res=300);
	#
	# 	png(file="RCircos.Demo.Mouse.And.Rat.png", height=8,
	#		width=9, unit="in", type="cairo", res=300);
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Initialize graphic device ...\n");

	pdf(file="RCircos.Demo.Mouse.And.Rat.pdf", height=8, width=9);

	layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(8,1), 
					heights=8);
	RCircos.Set.Plot.Area();


	#	Plot chromosome ideogram, heatmap, and link lines
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Plot mouse and rat chromosome ideograms\n");
	RCircos.Chromosome.Ideogram.Plot();


	cat("Draw heatmap\n");
	data.col <- 5;
	track.num <- 1;
	RCircos.Heatmap.Plot(expr.data, data.col, track.num, "in");

	cat("Draw link lines\n");
	track.num <- 3;
	RCircos.Link.Plot(link.data, track.num, FALSE);


	#	Add a title and a legend
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Add title and legend\n");
	title("RCircos Demo: Mouse and Rat Gene Expression");
	legend(0.6, 2, legend=c("Right: GSE42081 Mouse", 
			"Left:  GSE42081 Rat"), cex=0.8);


	#	Add a color key
	#  	__________________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	cat("Add color key\n");  		

	ColorRamp <- RCircos.Get.Heatmap.ColorScales("BlueWhiteRed");
	ColorLevels <- seq(min(expr.data$Expr.Mean), 
				max(expr.data$Expr.Mean), 
				length=length(ColorRamp));

	par(mai=c(1.5, 0.1, 1.5, 0.5));
	image(1, ColorLevels, matrix(data=ColorLevels, 
			ncol=length(ColorLevels), nrow=1),
			col=ColorRamp,  
			xlab="",  ylab="",  xaxt="n");


	#	Done. Close the device if plot to image file...
	#  	_______________________________________________
	#	xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

	dev.off();
	cat("R Circos Demo with mouse and Rat data done ...\n\n");

	rm(list=ls(all=T));
