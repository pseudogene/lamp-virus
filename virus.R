#
# This file is part of lamp-virus. It provides a R function to map mutation on the
# provided alignment, calculate clusters/groups, provide the results of a PCA analysis,
# and draw a phylogenetic tree.
#
# Copyright (C) 2014-2016, Michaël Bekaert <michael.bekaert@stir.ac.uk>
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
library(adegenet);
library(ape);
library(parallel);

process_virus <- function( prefix, fasta, output, extra=FALSE ) {
    if (!dir.exists(output)) {
        dir.create(output);
    }
    virus <- fasta2genlight(fasta, chunk=10, saveNbAlleles=TRUE, quiet=TRUE, parallel=TRUE);

    cluster.strict <- find.clusters(virus, n.pca=200, choose=FALSE, n.iter=1e5, stat="BIC");
    cluster.loose <- find.clusters(virus, n.pca=200, n.clust=length(cluster.strict$size)*2, n.iter=1e5, stat="WSS");
    write.table(merge(cluster.strict$grp, cluster.loose$grp, by = "row.names", all = TRUE), paste(output, "/", prefix, ".groups.tsv", sep=""), sep="\t", col.names = FALSE, quote=FALSE, row.names = FALSE);

    png(paste(output,"/location_of_the_SNPs.png", sep=""));
    temp <- density(position(virus), bw=10);
    plot(temp, type="n", xlab="Position in the alignment", main="Location of the SNPs", xlim=c(0,10761));
    polygon(c(temp$x,rev(temp$x)), c(temp$y, rep(0,length(temp$x))), col=transp("blue",.3));
    points(position(virus), rep(0, nLoc(virus)), pch="|", col="blue");
    dev.off();

    pca1 <- glPca(virus, nf=3);
    
    if (extra==TRUE) {
        pdf(paste(output,"/PCA_colourscatterplot.pdf", sep=""));
        myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4);
        abline(h=0,v=0, col="grey");
        title("PCA of the virus data\n axes 1-2");
        dev.off();

        png(paste(output,"/PCA_colourscatterplot.png", sep=""));
        myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4);
        abline(h=0,v=0, col="grey");
        title("PCA of the virus data\n axes 1-2");
        dev.off();
    }
  
    png(paste(output,"/PCA_scatterplot.strict.png", sep=""));
    myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4);
    text(pca1$scores[,1], pca1$scores[,2], cluster.strict$grp);
    abline(h=0,v=0, col="grey");
    title("PCA of the virus data\n axes 1-2 (strict grouping)");
    dev.off();  

    png(paste(output,"/PCA_scatterplot.loose.png", sep=""));
    myCol <- colorplot(pca1$scores,pca1$scores, transp=TRUE, cex=4);
    text(pca1$scores[,1], pca1$scores[,2], cluster.loose$grp);
    abline(h=0,v=0, col="grey");
    title("PCA of the virus data\n axes 1-2 (loose grouping)");
    dev.off(); 

    write.table(pca1$scores,file=paste(output, "/", prefix, ".pca.tsv", sep=""),sep="\t");

    tre <- nj(dist(as.matrix(virus)));
    tre;
    
    if (extra==TRUE) {
        pdf(paste(output,"/NJ_tree.pdf", sep=""));
        plot(tre, typ="fan", show.tip=FALSE)
        tiplabels(pch=20, col=myCol, cex=4);
        title("NJ tree of the virus data");
        dev.off();
    
        pdf(paste(output,"/NJ_tree.named.pdf", sep=""));
        plot(tre, typ="fan", cex=0.7);
        title("NJ tree of the virus data");
        dev.off();

        png(paste(output,"/NJ_tree.png", sep=""));
        plot(tre, typ="fan", show.tip=FALSE)
        tiplabels(pch=20, col=myCol, cex=4);
        title("NJ tree of the virus data");
        dev.off();
    }

    png(paste(output,"/NJ_tree.strict.png", sep=""));
    plot(tre, typ="fan", show.tip=FALSE)
    tiplabels(pch=20, col=myCol, text=cluster.strict$grp, cex=4, bg=NULL, frame="none");
    title("NJ tree of the virus data\n(strict grouping)");
    dev.off();
    
    png(paste(output,"/NJ_tree.loose.png", sep=""));
    plot(tre, typ="fan", show.tip=FALSE)
    tiplabels(pch=20, col=myCol, text=cluster.loose$grp, cex=4, bg=NULL, frame="none");
    title("NJ tree of the virus data\n(loose grouping)");
    dev.off();

    write.tree(tre, file = paste(output, "/", prefix, ".tree.tre", sep=""));

    save(virus, tre, pca1, myCol, cluster.strict, cluster.loose, file=paste(output, "/", prefix, ".Rsave", sep=""));
}
