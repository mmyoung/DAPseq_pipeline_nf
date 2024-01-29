#!/usr/bin/env Rscript
### usage: Rscript input_depth.txt sample_name out_histogram.pdf
args = commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(dplyr)

depth_df <- 
      read.table(args[1],header=F) 

genome_cov <- 
      round(colSums(depth_df[which("V1">0),])[1]/colSums(depth_df)[1]*100,digits=2)

(depth_df %>%
filter(!V2==0) %>%
    ggplot(aes(x=V2,y=V1,group=1)) +
    geom_line()+
    labs(title=paste0(args[2],", coverage= ",genome_cov,"%"),
      x="reads coverage",y="loci number")+
    theme_classic()+
    theme(axis.text = element_text(size = 10,colour = "black"),
          plot.title = element_text(hjust = 0.5,vjust = 0.5)))%>%
    ggsave(filename=args[3],width=5,height=4)
