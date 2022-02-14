#library(Gviz)
library(stringr)

isEmpty <- function(x) { #This function checks if a data frame is empty or not
  return(length(x)==0)
}

system("sh ./bin/blaster.bash") ##With this script all multifatsa files al processed to get the blast against the custom database result and prokka resul files.

extract_table=data.frame(matrix(ncol=2))

files_totales=list.files("./Seqs/") 

input_acc=data.frame(matrix(ncol=1))
e=1

for(i in 1:length(files_totales)) ## Here we extract all the information for what sequence is include in part of each assembly project/genome.
{
  if(!isEmpty(grep(".fa$",files_totales[i])))
  {
    
    multifa_original=read.delim(paste("./Seqs",files_totales[i],sep = "/"),header=F)


    multifa_original=multifa_original[grepl("^>", multifa_original[,1]), ]

   #multifa_original=as.data.frame(multifa_original[ seq(1, nrow(multifa_original), by = 2),1])
    
    extract_table=rbind(extract_table,data.frame(cbind(files_totales[i],str_replace_all(gsub(" ", "", sapply(strsplit(sapply(strsplit(as.character(multifa_original), "\\ "), "[[", 1), "\\]"), "[[", 1), fixed = TRUE),">",""))))

    
  }
}

write.table(extract_table,"GCA_access.txt",sep="\t",row.names = F,col.names = F,quote = F) ##We save a table file with the information of what sequences are include in each assembly project


#########################################


multifa_original=read.delim("aux_file_prodigal_global_headers.txt",header=F) ##ORFs from Prodigal analysis are saved, taking care of the position information.

tab_aux_res=data.frame(matrix(ncol=3,nrow=length(multifa_original[,1])))

tab_aux_res[,1]=gsub(" ", "", sapply(strsplit(sapply(strsplit(as.character(multifa_original[,1]), "\\#"), "[[", 1), "\\]"), "[[", 1), fixed = TRUE)
tab_aux_res[,2]=as.numeric(gsub(" ", "", sapply(strsplit(sapply(strsplit(as.character(multifa_original[,1]), "\\#"), "[[", 2), "\\]"), "[[", 1), fixed = TRUE))
tab_aux_res[,3]=as.numeric(gsub(" ", "", sapply(strsplit(sapply(strsplit(as.character(multifa_original[,1]), "\\#"), "[[", 3), "\\]"), "[[", 1), fixed = TRUE))
tab_aux_res[,4]=tab_aux_res[,3]-tab_aux_res[,2]
tab_aux_res[,5]="+"

tab_aux_res[,1]=str_replace_all(tab_aux_res[,1],">","")



tab_GCA=read.delim("GCA_access.txt",header = F)

tab_blast=read.delim("blast_output_file",header = F) ###Blast result for each multifasta is loaded and the best hist from the most abundant hits are saved in each case

###########Aqui se puede meter código para filtar el blast, por ejemplo 90% minimo identidad

tab_blast$acc=sub('[_][^_]+$', '', as.character(tab_blast[,1]))

tab_blast$assigned=sapply(strsplit(sapply(strsplit(as.character(tab_blast[,2]), "\\?"), "[[", 3),"-"),"[[",1)


tab_blast$GCA=tab_GCA[match(tab_blast$acc,tab_GCA$V2),1]


o=1
u=1

tab_resultados=data.frame(matrix(ncol=dim(tab_blast)[2]))
tab_tipo_operon=data.frame(matrix(ncol=2))

filter_file=read.table("./Database/ETEC_CFs.txt",na.strings=c("","NA"),sep="\t",row.names=1)

for(e in 1:length(unique(tab_blast$GCA))) ##Within this loop information about what ETEC proteins are found for each assembly/project is saved in two files. "All_results_hit" will save the best hit for each ORF based on the most abundant one. "Composition_by_assembly" will save the information about the proteins included for each genome. 
{
  sub_blast=tab_blast[tab_blast$GCA==unique(tab_blast$GCA)[e],]
  
  for(i in 1:length(unique(sub_blast[,1])))
  {
    most_abundant=data.frame(table(sub_blast[sub_blast[,1]==unique(sub_blast[,1])[i],][,14]))[data.frame(table(sub_blast[sub_blast[,1]==unique(sub_blast[,1])[i],][,14]))[,2]==max(table(sub_blast[sub_blast[,1]==unique(sub_blast[,1])[i],][,14])),][1,]
    
    tab_resultados[o,]=sub_blast[sub_blast[,1]==unique(sub_blast[,1])[i],][sub_blast[sub_blast[,1]==unique(sub_blast[,1])[i],][,14]==most_abundant[1,1],][1,]
    
    o=o+1
  }
  
  tab_CS_aux=tab_resultados[tab_resultados$X15==unique(tab_blast$GCA)[e],] ###Here we only keep the major subunit result for the output! When taking into account the name for each protein we will name it after the information we get from the "major subunit"
  tab_CS_aux_not=NULL
  list_hist=unique(tab_CS_aux[,14])
  total_values_list=list()
  
  for(i in 1:length(list_hist))
  {
	values_look=as.character(na.omit(as.character(filter_file[rownames(filter_file)==list_hist[i],])))
	values_look[values_look=="character(0)"] <- NA
	values_look=na.omit(values_look)
	
	if(length(values_look)>0)
	{
	  total_values_list[list_hist[i]]=0
	  
	  for(u in 1:length(values_look))
	  {
	    if(!isEmpty(grep(values_look[u],tab_resultados[,2])))
	    {
	      total_values_list[list_hist[i]]=total_values_list[list_hist[i]][[1]]+1
	      #print(tab_resultados[grep(values_look[u],tab_resultados[,2]),14]) ##Si lo encuentra, guardar el resultado en una lista. Si el total de resultados es el mismo que el ttoal de valores que se buscan, se anota, si no no ACABAR
	    }
	  }
	  
	  if(total_values_list[list_hist[i]][[1]]!=length(values_look))
	  {
	    tab_CS_aux_not=rbind(tab_CS_aux_not,tab_CS_aux[tab_CS_aux[,14]==list_hist[i],])
	    tab_CS_aux=tab_CS_aux[tab_CS_aux[,14]!=list_hist[i],]
	  }
	}
	
  }
  
  ####WE ONLY KEEP THE VALUES OF CS WHEN MAJOR_SUBUNIT IS FOUND!!!
  
  #tab_CS_aux=rbind(tab_CS_aux[!grepl(pattern = "CS",x = tab_CS_aux[,2]),],tab_CS_aux[grepl(pattern = "major_subunit",x = tab_CS_aux[,2]),])
  
  
  tab_tipo_operon[e,]=c(unique(tab_blast$GCA)[e],paste(unique(tab_CS_aux[tab_CS_aux$X15==unique(tab_blast$GCA)[e],][,14]),collapse = ", "))
}



write.table(tab_tipo_operon,"Composition_by_assembly.txt",sep="\t",quote=F,row.names = F,col.names = F) ##Files are saved in this line
tab_resultados_print=tab_resultados
colnames(tab_resultados_print)=c("qseqid","sseqid","pident","length","mismatch","gapopen",
                                 "qstart","qend","sstart","send","evalue","bitscore","ACC","Protein","GCA")
write.table(tab_resultados_print,"All_results_hit.txt",sep="\t",quote=F,row.names = F)

colnames(tab_CS_aux_not)=c("qseqid","sseqid","pident","length","mismatch","gapopen",
                                 "qstart","qend","sstart","send","evalue","bitscore","ACC","Protein","GCA")
write.table(tab_CS_aux_not,"Incomplete_CFs.txt",sep="\t",quote=F,row.names = F)

system("perl ./bin/prokka_parser.pl") ##With this script we will parse Prokka results and get the information from the most abundant hit replacing the result from original prokka execution
system("mkdir ./Output/Aux_output_files/")
system("mv ./Output/*_prokka_out ./Output/Aux_output_files/")
system("rm -rf aux_file_* blast_output_file ./Seqs/aux_fil* ./Seqs/*txt GCA_access.txt") ##Temporary file are removed.

####GVIZ Analysis


#tab_aux_res[match(tab_resultados[,1],tab_aux_res[,1]),6]=tab_resultados[,2]
#tab_aux_res=tab_aux_res[!is.na(tab_aux_res[,6]),]

##tab_aux_res$acc=sapply(strsplit(as.character(tab_aux_res[,1]), "_"), "[[", 1)

##tab_aux_res[,7]=tab_GCA[match(sapply(strsplit(as.character(tab_aux_res[,1]), "_"), "[[", 1),tab_GCA[,2]),1]

#tab_aux_res[,7]=tab_GCA[match(sub('[_][^_]+$', '', as.character(tab_aux_res[,1])),tab_GCA[,2]),1]



#colnames(tab_aux_res)=c("chromosome","start","end","width","strand","symbol","genome")


##tab_aux_res[,1]=sapply(strsplit(as.character(tab_aux_res[,1]), "_"), "[[", 1)
#tab_aux_res[,1]=sub('[_][^_]+$', '', as.character(tab_aux_res[,1]))
##tab_aux_res[,8]=sapply(strsplit(as.character(tab_aux_res[,1]), "_"), "[[", 1)
#tab_aux_res[,8]=tab_aux_res[,1]

###tab_aux_res[,1]=paste("chr",match(tab_aux_res[,1],unique(sapply(strsplit(as.character(tab_aux_res[,1]), "_"), "[[", 1))),sep = "")



#for(e in 1:length(unique(tab_aux_res[,7])))
#{
 # pdf(paste(unique(tab_aux_res[,7])[e],"_genome_plot.pdf",sep = ""),width = 15,height = 5)
  
  #sub_blast=tab_aux_res[tab_aux_res[,7]==unique(tab_aux_res[,7])[e],]
  
  ##sub_blast[,1]=paste("chr",match(sub_blast[,1],unique(sapply(strsplit(as.character(sub_blast[,1]), "_"), "[[", 1))),sep = "")
  
  #for(i in 1:length(unique(sub_blast[,1])))
  #{
   # tab_aux_res2=sub_blast[sub_blast[,1]==unique(sub_blast[,1])[i],]
    
    #genome_GRange=makeGRangesFromDataFrame(tab_aux_res2)
    ##atrack <- AnnotationTrack(genome_GRange, name = chromosome)
    #gtrack <- GenomeAxisTrack()
	#options(ucscChromosomeNames=FALSE)
    #grtrack <- GeneRegionTrack(tab_aux_res2, name = tab_aux_res2[,8][1],transcriptAnnotation = "symbol",background.panel = "azure",
     #                          background.title = "darkblue",rstarts = "start",rends = "end")
    #plotTracks(list(gtrack, grtrack))
  #}
  
  #dev.off()
#}


#system(paste('/storage/parras/CMV/EggNOG/A_ver/eggnog-mapperpper.py -i /storage/parras/ProyectoAstrid/Prokka_prueba/Prueba_seq_unica/aux_file_prodigal_global_aa -o /storage/parras/ProyectoAstrid/Prokka_prueba/Prueba_seq_unica/JOJOJO/ --cpu 30 -m hmmer --data_dir /storage/parras/CMV/EggNOG/A_ver/eggnog-mapper-master/data/ -d Bacteria',sep=""))
