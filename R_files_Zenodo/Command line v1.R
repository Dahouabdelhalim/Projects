# Set working directory
setwd("~/Documents/Work/SDG_project/Phase2/Temp")

# Using the functions (all names and values are examples)
ReadFiles('output')
CoBibMatrix('output', 16)
Threshold_Eval('output', 16, 5, 30, 0.5)
Threshold_Cleaner('output', 'output_clean', 17)
Search_Function('output_clean', 'SDG_pubs', 'ramirez_thesaurus', 16)
SDG_Triads('SDG_pubs', 'output_clean', 'SDG_triads', 16)
SDG_Communities('SDG_pubs', 'output_clean', 0.4, 0.4, 'yes')

# Execute all lines manually
##########################################################################

################################
### STEP 1: DATA PREPARATION ###
################################

# Input values
output = 'output'

###
direct=getwd()
setwd(paste(direct,'Data',sep='/'))
lis=list.files(,pattern=".txt")

# Define the column that should be included
a=c('PT','UT','PY','TI','AB','ID','DE','SO','AU','C1','WC','SC','CR','Z9','DA','FU'); mm=c()

# Read data and paste into matrix 
for(k in lis){x=read.delim(k,skip=2, header = F,encoding="UTF-8", sep = "\\t",quote="")
  pos=c(which(substring(x[,1],1,2)=="PT"),dim(x)[1])
  m=matrix(,length(pos)-1,length(a));colnames(m)=a
  for(i in 1:dim(m)[1]){y=as.character(x[pos[i]:(pos[i+1]-1),1])
    b=substring(y,1,2)[substring(y,1,2)!='  ']
    m[i,]=trimws(unlist(lapply(a,function(j){if(length(y[substring(y,1,2)==j])!=0){if(j%in%c("CR",'AU')){trimws(substring(paste(trimws(y[which(substring(y,1,2)==j):(which(substring(y,1,2)%in%b)[which(substring(y,1,2)%in%b)>which(substring(y,1,2)==j)][1]-1)]),collapse='##'),3))}else{
    trimws(substring(paste(trimws(y[which(substring(y,1,2)==j):(which(substring(y,1,2)%in%b)[which(substring(y,1,2)%in%b)>which(substring(y,1,2)==j)][1]-1)]),collapse=' '),3))}}else{NA}})))}
    mm=rbind(mm,m)}

# Create folder 'Analysis' and write output
setwd(direct)
dir.create('Analysis'); setwd(paste(direct,'Analysis',sep='/'))
write.csv(unique(mm),paste(output,'.csv',sep=''),row.names=F)
setwd(direct)

#######################################
### STEP 2: CO-BIBLIOGRAPHY NETWORK ###
#######################################

# Input values
file='output'
cores=16

####
direct=getwd()
setwd(paste(direct,'Analysis',sep='/'))
data=read.csv(paste(file,'.csv',sep=''))

# Load libraries and extract citation data 
library(parallel);library(Matrix);library(igraph)
UT=unlist(mclapply(1:dim(data)[1],function(i){rep(as.character(data[i,"UT"]),length(trimws(toupper(unlist(strsplit(as.character(data[i,"CR"]),"##"))))))},mc.cores=cores))
CR=unlist(mclapply(1:dim(data)[1],function(i){trimws(toupper(unlist(strsplit(as.character(data[i,"CR"]),"##"))))},mc.cores=cores))
UT=substring(UT,1,30)

# Create matrix UTxUT with number of shared citations 
data=cbind(UT,CR)
data=data[data[,'CR']%in%rownames(as.matrix(table(data[,"CR"])[table(data[,"CR"])>1])),]
data=data[data[,'UT']%in%rownames(as.matrix(table(data[,"UT"])[table(data[,"UT"])>1])),]
M=Matrix(0,length(unique(data[,"UT"])),length(unique(data[,"UT"])))
colnames(M)=rownames(M)=sort(as.character(unique(data[,"UT"])))
a=c(round(seq(1,dim(M)[1],dim(M)[1]/100)),dim(M)[1])
for(i in 1:(length(a)-2)){
  x=mclapply((i+1):(length(a)-1),function(j){datas=data[data[,'UT']%in%c(colnames(M)[a[i]:a[i+1]],colnames(M)[a[j]:a[j+1]]),]
  n=Matrix(table(datas[,'UT'],datas[,'CR']))
  n=n%*%t(n)
  as.matrix(100*round(n/sqrt(diag(n)%*%t(diag(n))),3))},mc.cores=cores)
  for(j in 1:length(x)){
    b=as.data.frame(x[j])
    M[rownames(M)%in%rownames(b),rownames(M)%in%rownames(b)]=as.matrix(b)}}

# Create igraph network and write output
g=graph_from_adjacency_matrix(M, mode = "undirected",weighted = T, diag = F)
write.graph(g,paste(file,'network_full.graphml',sep='_'),format='graphml')
setwd(direct)


#################################################
### STEP 3: ANALYSING CO-BIBLIOGRAPHY NETWORK ###
#################################################

# Input values
file='output'
cores=10
min_val=5; max_val=30; step=0.5

####
direct=getwd()
setwd(paste(direct,'Analysis',sep='/'))

# Load packages and set up matrix
library(igraph);library(parallel)
g=read.graph(paste(file,'network_full.graphml',sep='_'),format='graphml')
a=seq(min_val,max_val,step)
m=matrix(0,length(a),14);rownames(m)=a

# Fill matrix as function of thresholds
x=mclapply(a,function(i){g2=delete.edges(g,E(g)[!E(g)$weight>i])
  g2=delete.vertices(g2,V(g2)[degree(g2)==0]);clu=cluster_louvain(g2)
  b1=paste(i,max(components(g2)$csize),sep="#")
  x=components(g2);g2=delete.vertices(g2,V(g2)[x$membership!=which(x$csize==max(x$csize))]);clu=cluster_louvain(g2)
  b2=paste(length(unique(clu$membership)),modularity(clu),length(table(clu$membership)[table(clu$membership)>10]),
         sum(table(clu$membership)[table(clu$membership)>10]),length(triangles(g2))/3,sep="#")
  paste(b1,b2,sep="#")},mc.cores=cores)
m=matrix(unlist(strsplit(unlist(x),"#")),byrow=T,ncol=7)
colnames(m)=c("Threshold","Nodes","#Com","Modularity","#Com.10","NodesCom.10","Triangles")

# Write output 
write.csv(m,paste(file,"T-Eval.csv",sep='_'),row.names=F)
setwd(direct)


#################################
### STEP 4: CLEANING THE DATA ###
#################################

# Input values
file='output'
output='output_cleaned'
value=17

####
direct=getwd()
setwd(paste(direct,'Analysis',sep='/'))
library(Matrix);library(igraph)
g=read.graph(paste(file,'network_full.graphml',sep='_'),format='graphml')
data=read.csv(paste(file,'.csv',sep=''))
data[,'UT']=substring(data[,'UT'],1,30)

# Remove the edges and nodes based on the given threshold
g=delete.edges(g,E(g)[!E(g)$weight>value]);x=components(g)
g=delete.vertices(g,V(g)[x$membership!=which(x$csize==max(x$csize))])

# Implement Louvain clustering algorithm
clu=cluster_louvain(g);modularity(clu)
clu$membership[clu$membership%in%rownames(as.matrix(table(clu$membership)[!table(clu$membership)>10]))]=0
g=delete.vertices(g,V(g)[clu$membership==0])
clu=cluster_louvain(g); modularity(clu) #mod = 0.9839205
g=set_vertex_attr(g,'COM',V(g),value=clu$membership)

# Add the communities to the data
datas=data[unlist(lapply(1:dim(data)[1],function(i){which(data[,"UT"]==names(V(g))[i])})),]
datas = datas[!duplicated(datas[,'UT']),]
datas=cbind(datas,paste('COM', clu$membership,sep=''))
colnames(datas)[dim(datas)[2]]="COM"
datas=datas[datas[,"COM"]!=0,]

# Write output
write.csv(datas,paste(output,'.csv',sep=''),row.names=F)
write.graph(g,paste(output,'.graphml',sep=''),format='graphml')
setwd(direct)


##############################
### STEP 5: KEYWORD SEARCH ###
##############################

# Input values
file = 'output'
output = 'data_searched_full'
thesaurus = 'ramirez_thesaurus'
cores = 16

###
direct=getwd()
setwd(paste(direct,'Analysis',sep='/'))
library(textstem);library(tm);library(parallel)
data=as.matrix(read.csv(paste(file,".csv",sep="")))
words=as.matrix(read.csv(paste(thesaurus,'.csv',sep='')))

# The search function
Search=function(j){
  if(words[j,'Mode']=='stemmed'){
    if(words[j,'type.1']==1){
      if(length(grep(paste("\\\\b", tolower(words[j,'stemmed']), "\\\\b", sep=""),datas_stemmed))!=0){paste(grep(paste("\\\\b", tolower(words[j,'stemmed']), "\\\\b", sep=""),datas_stemmed),j,
              unlist(lapply(datas_stemmed[grep(paste("\\\\b", tolower(words[j,'stemmed']), "\\\\b", sep=""),datas_stemmed)],function(k){
                              length(unlist(gregexpr(paste("\\\\b", tolower(words[j,'stemmed']), "\\\\b", sep=""),k)))})),sep='#')}}
    else{
      x=table(unlist(lapply(unlist(strsplit(words[j,'stemmed'],'#')),function(h){grep(paste("\\\\b",tolower(h),"\\\\b",sep=""),datas_stemmed)})))
      if(length(x[x==length(unlist(strsplit(words[j,'stemmed'],'#')))])!=0){paste(names(x[x==length(unlist(strsplit(words[j,'stemmed'],'#')))]),j,
              apply(matrix(unlist(lapply(unlist(strsplit(words[j,'stemmed'],'#')),function(h){unlist(lapply(datas_stemmed[as.numeric(names(x[x==length(unlist(strsplit(words[j,'stemmed'],'#')))]))],
              function(k){length(unlist(gregexpr(paste("\\\\b",h,"\\\\b",sep=""),k)))}))})),byrow=T,ncol=length(unlist(strsplit(words[j,'stemmed'],'#')))),1,min),sep='#')}}}
  else{
    if(words[j,'type.1']==1){
      if(length(grep(paste('\\\\b',tolower(words[j,'corrected']),sep=""),datas))!=0){paste(grep(paste('\\\\b',tolower(words[j,'corrected']),sep=""),datas),j,
              unlist(lapply(datas[grep(paste('\\\\b',tolower(words[j,'corrected']),sep=""),datas)],function(k){
                              length(unlist(gregexpr(paste('\\\\b',tolower(words[j,'corrected']),sep=""),k)))})),sep='#')}}
    else{
      x=table(unlist(lapply(unlist(strsplit(words[j,'corrected'],'#')),function(h){grep(paste('\\\\b',tolower(h),sep=""),datas)})))
      if(length(x[x==length(unlist(strsplit(words[j,'corrected'],'#')))])!=0){paste(names(x[x==length(unlist(strsplit(words[j,'corrected'],'#')))]),j,
              apply(matrix(unlist(lapply(unlist(strsplit(words[j,'corrected'],'#')),function(h){unlist(lapply(datas[as.numeric(names(x[x==length(unlist(strsplit(words[j,'corrected'],'#')))]))],
              function(k){length(unlist(gregexpr(paste('\\\\b',h,sep=""),k)))}))})),byrow=T,ncol=length(unlist(strsplit(words[j,'corrected'],'#')))),1,min),sep='#')}}}
}

# Prepare thesaurus and data
words[,'corrected']=unlist(lapply(words[,'corrected'], function(x) gsub('â€�', '\\\\\\\\b', x)))
words[,'corrected']=unlist(lapply(words[,'corrected'], function(x) gsub('\\\\*', '', x)))
datas=data[,c('TI','AB','DE')]
data[,'UT']=substring(data[,'UT'],1,30)
datas=gsub("  "," ",gsub('[[:digit:]]+', '', gsub('  ',' ',gsub('[[:punct:] ]+',' ',tolower(unlist(mclapply(1:dim(datas)[1],function(i){paste(datas[i,1],datas[i,2],datas[i,3])},mc.cores=cores)))))))
datas=gsub("  "," ",unlist(mclapply(datas,function(i){paste(unlist(strsplit(i,' '))[!unlist(strsplit(i,' '))%in%stopwords()],collapse=" ")},mc.cores=cores)))
datas_stemmed=gsub("  "," ",unlist(mclapply(datas,function(i){paste(stem_words(unlist(strsplit(i,' '))[!unlist(strsplit(i,' '))%in%stopwords()]),collapse=" ")},mc.cores=cores)))

# Search each publication for the keywords 
search=matrix(unlist(strsplit(unlist(mclapply(1:dim(words)[1],Search,mc.cores=cores)),"#")),byrow=T,ncol=3)
m=matrix(unlist(strsplit(unlist(lapply(1:dim(search)[1],function(i){paste(data[as.numeric(search[i,1]),'UT'], words[as.numeric(search[i,2]),2],words[as.numeric(search[i,2]),1],words[as.numeric(search[i,2]),3],search[i,3],sep='##')})),'##')),byrow=T,ncol=5)
colnames(m)=c('UT','KeyWord','SDG','Type','Freq')
write.csv(m,paste(output,"Words_Searched.csv",sep='_'),row.names=F)

# Determine the SDG(s) for each publication
datas=xtabs(as.numeric(m[,'Freq'])~m[,'UT']+m[,'SDG'])
datas=cbind(datas,rowSums(datas))
colnames(datas)=c(paste('SDG', 1:17, sep = '.'),'FreqT');class(datas)='numeric'
data=data[data[,'UT']%in%rownames(datas),]
datas=datas[unlist(mclapply(1:dim(data)[1],function(i){which(rownames(datas)==data[i,'UT'])},mc.cores=cores)),]
data=data.frame(data,datas)
SDG=unlist(mclapply(1:dim(data)[1],function(i){ifelse(sum(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T)[1])>.75,
                                                      paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1],collapse='#'),
                                                      ifelse(sum(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T)[1:2])>.60,
                                                             paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1:2],collapse='#'),
                                                             ifelse(sum(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T)[1:3])>.50,
                                                                    paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1:3],collapse='#'),
                                                                    paste(names(sort(data[i,paste('SDG.',1:17, sep='')]/data[i,'FreqT'],decreasing=T))[1:3],collapse='#'))))},mc.cores=cores))

data=cbind(data,SDG)

fields = c('UT', 'ACOUSTICS', 'ASTRONOMY & ASTROPHYSICS', 'CRYSTALLOGRAPHY', 'OPTICS')
delete=NULL
for (index in grep(paste(fields,collapse="|"),data$SC)){field=data[index,'SC']
if(length(unlist(strsplit(field,";")))==1){delete=append(delete,data[index,'UT'])}
else if(length(grep(paste(fields,collapse="|"),unlist(strsplit(field,";"))))==length(unlist(strsplit(field,";")))){
  delete=append(delete,data[index,'UT'])}}
data=subset(data,!UT%in%delete)

# Write output 
write.csv(data,paste(output,'.csv',sep=''),row.names=F)
setwd(direct)


#####################################
### STEP 6: TRIADS IN THE NETWORK ###
#####################################

# Input values
input_file='data_searched_full'
input_network='output_network_full'
output= 'SDG_triads'
cores=16

###
direct=getwd()
setwd(paste(direct,'Analysis',sep='/'))
SDG=c('TD','ST','ST','ST','TD','ST','ST','TD','ST','TD','ST','TD','TD','ST','TD','FC','FC')
library(igraph);library(parallel)
data=read.csv(paste(input_file,".csv",sep=""))
g=read.graph(paste(input_network,".graphml",sep=""),format="graphml")
'%notin%' <- Negate(`%in%`)
g=delete.vertices(g,names(V(g))[names(V(g)) %notin% data$UT])

# Prepare the data
a=b=c();for(i in 1:dim(data)[1]){a=c(a,paste(as.character(data[i,"UT"]),1:length(unlist(strsplit(as.character(data[i,"SDG"]),"#"))),sep="-"))
b=c(b,unlist(strsplit(as.character(data[i,"SDG"]),"#")))}
datas=cbind(a,b);colnames(datas)=c("UT","SDG")
g=as_edgelist(g);g[,1]=paste(g[,1],'-1',sep='');g[,2]=paste(g[,2],'-1',sep='')
for(i in 1:dim(data)[1]){for(j in 1:length(unlist(strsplit(as.character(data[i,"SDG"]),"#")))){
  if(j!=1){g1=g[g[,1]==paste(as.character(data[i,"UT"]),1,sep='-')|g[,2]==paste(as.character(data[i,"UT"]),1,sep='-'),]
  g1[g1==paste(as.character(data[i,"UT"]),1,sep='-')]=paste(as.character(data[i,"UT"]),j,sep="-")
  g=rbind(g,g1)}}}
g=graph_from_edgelist(g, directed = F)

# Triads for the whole network
tris=matrix(names(V(g))[triangles(g)],byrow=T,ncol=3)
Trigs=unlist(mclapply(1:dim(tris)[1],function(i){mar=sort(gsub("SDG.","",c(as.character(datas[datas[,"UT"]==as.character(tris[i,1]),"SDG"]),as.character(datas[datas[,"UT"]==as.character(tris[i,2]),"SDG"]),as.character(datas[datas[,"UT"]==as.character(tris[i,3]),"SDG"]))))
              mar2=sort(SDG[as.numeric(mar)]);paste(paste(mar[1],mar[2],mar[3],sep="_"),paste(mar2[1],mar2[2],mar2[3],sep="_"),sep="#")},mc.cores=cores))
M=matrix(unlist(strsplit(Trigs,"#")),byrow=T,ncol=2)
write.csv(M, paste(output,'full_network.csv',sep=''),row.names=F)  

# Triads per community
if ("COM" %in% colnames(data)){
  M=c()
  for(j in 1:length(unique(data[,"COM"]))){
    UTS=paste(sort(rep(as.character(data[data[,"COM"]==unique(data[,"COM"])[j],"UT"]),17)),1:17,sep="-")
    if(length(rownames(as.matrix(triangles(delete.vertices(g,rownames(as.matrix(V(g)))[!rownames(as.matrix(V(g)))%in%UTS])))))!=0){
      tris=matrix(rownames(as.matrix(triangles(delete.vertices(g,rownames(as.matrix(V(g)))[!rownames(as.matrix(V(g)))%in%UTS])))) ,byrow=T,ncol=3)
      Trigs=unlist(mclapply(1:dim(tris)[1],function(i){
        mar=sort(gsub("SDG.","",c(as.character(datas[datas[,"UT"]==as.character(tris[i,1]),"SDG"]),as.character(datas[datas[,"UT"]==as.character(tris[i,2]),"SDG"]),as.character(datas[datas[,"UT"]==as.character(tris[i,3]),"SDG"]))))
        mar2=sort(SDG[as.numeric(mar)])
        paste(paste(mar[1],mar[2],mar[3],sep="_"),paste(mar2[1],mar2[2],mar2[3],sep="_"),sep="#")},mc.cores=cores))
      m=matrix(unlist(strsplit(Trigs,"#")),byrow=T,ncol=2)
      M=rbind(M,cbind(m,rep(as.character(unique(data[,"COM"])[j]),dim(m)[1])))}}
  write.csv(M, paste(output,'communities.csv',sep=''),row.names=F)}
setwd(direct)
          

###############################
### STEP 7: SDG communities ###
###############################

# Input values 
file_sdg='SDG_pubs'
file_data='output_clean'
T4_Cut=0.4
Total_Cut=0.3
save_fig='yes'; figname='SDG_share'

###
direct=getwd(); dir.create('Images')
setwd(paste(direct,'Analysis',sep='/'))
library(dplyr);library(ggplot2)
sdg_pubs = read.csv(paste(file_sdg,'.csv',sep=""))
all_pubs <- read.csv(paste(file_data,'.csv',sep=""))

#SDG share per year 
pubs_yr = merge(data.frame(table(all_pubs$COM, all_pubs$PY)), data.frame(table(sdg_pubs$COM, sdg_pubs$PY)), by = c('Var1', 'Var2'), all.x = T)
pubs_yr[is.na(pubs_yr)]=0
pubs_yr$share.xy = pubs_yr$Freq.y/pubs_yr$Freq.x
pubs_yr = reshape(pubs_yr, idvar = "Var1", timevar = "Var2", direction = "wide", new.row.names=1:length(unique(all_pubs$COM))) #; names(pubs_years)= c('COM', paste('pubs', 2000:2020, sep = '_'))
names(pubs_yr)=gsub(x = names(pubs_yr), pattern = "Freq.x.", replacement = "pubs_"); names(pubs_yr)=gsub(x = names(pubs_yr), pattern = "Freq.y.", replacement = "sdgpubs_"); names(pubs_yr)=gsub(x = names(pubs_yr), pattern = "share.xy.", replacement = "share_"); names(pubs_yr)=gsub(x = names(pubs_yr), pattern = "Var1", replacement = "COM") 
pubs_yr$total_pubs = rowSums(pubs_yr[, c(paste('pubs', unique(all_pubs$PY), sep = '_'))]); pubs_yr$total_sdgpubs = rowSums(pubs_yr[, c(paste('sdgpubs', unique(all_pubs$PY), sep = '_'))])
pubs_yr$total_sdg_share=pubs_yr$total_sdgpubs/pubs_yr$total_pubs

sdg_share=pubs_yr[c('COM', paste('share', unique(all_pubs$PY), sep = '_'))]
sdg_share[is.na(sdg_share)]=0
fits <- lm.fit(cbind(1, seq_len(ncol(sdg_share[,-1]))), t(sdg_share[,-1]))

#SDG publications in communities over the years (increase/decrease), timeframes of 5 years 
pubs_tf = data.frame(table(all_pubs$COM, all_pubs$PY)) 

timeframes=cut(min(all_pubs$PY):max(all_pubs$PY),4)
intervals=unlist(lapply(levels(timeframes),function(x) regmatches(x, gregexpr("[[:digit:]]+", x))))

pubs_tf$Var2=gsub(paste(intervals[1]:intervals[2], collapse="|"), 't1', as.character(pubs_tf$Var2)); pubs_tf$Var2=gsub(paste((as.numeric(intervals[3])+1):intervals[4], collapse="|"), 't2', pubs_tf$Var2); pubs_tf$Var2=gsub(paste((as.numeric(intervals[5])+1):intervals[6], collapse="|"), 't3', pubs_tf$Var2); pubs_tf$Var2=gsub(paste((as.numeric(intervals[7])+1):intervals[8], collapse="|"), 't4', pubs_tf$Var2)
pubs_tf = pubs_tf %>% group_by(Var1, Var2) %>% summarise(Freq = sum(Freq))
sdg_pubs_tf = data.frame(table(sdg_pubs$COM, sdg_pubs$PY))
sdg_pubs_tf$Var2=gsub(paste(intervals[1]:intervals[2], collapse="|"), 't1', as.character(sdg_pubs_tf$Var2)); sdg_pubs_tf$Var2=gsub(paste((as.numeric(intervals[3])+1):intervals[4], collapse="|"), 't2', sdg_pubs_tf$Var2); sdg_pubs_tf$Var2=gsub(paste((as.numeric(intervals[5])+1):intervals[6], collapse="|"), 't3', sdg_pubs_tf$Var2); sdg_pubs_tf$Var2=gsub(paste((as.numeric(intervals[7])+1):intervals[8], collapse="|"), 't4', sdg_pubs_tf$Var2)
sdg_pubs_tf = sdg_pubs_tf %>% group_by(Var1, Var2) %>% summarise(Freq = sum(Freq))

pubs_tf = merge(pubs_tf, sdg_pubs_tf, by = c('Var1', 'Var2'), all.x = T)
pubs_tf$share.xy = pubs_tf$Freq.y/pubs_tf$Freq.x
pubs_tf = reshape(pubs_tf, idvar = "Var1", timevar = "Var2", direction = "wide", new.row.names=1:length(unique(all_pubs$COM))) #; names(pubs_tf)= c('COM', paste('pubs', 2000:2020, sep = '_'))
names(pubs_tf)=gsub(x = names(pubs_tf), pattern = "Freq.x.", replacement = "pubs_"); names(pubs_tf)=gsub(x = names(pubs_tf), pattern = "Freq.y.", replacement = "sdgpubs_"); names(pubs_tf)=gsub(x = names(pubs_tf), pattern = "share.xy.", replacement = "share_"); names(pubs_tf)=gsub(x = names(pubs_tf), pattern = "Var1", replacement = "COM") 
pubs_tf$total_pubs = rowSums(pubs_tf[, c(paste('pubs', paste('t', 1:4, sep=""), sep = '_'))]); pubs_tf$total_sdgpubs = rowSums(pubs_tf[, c(paste('sdgpubs', paste('t', 1:4, sep=""), sep = '_'))])
pubs_tf$total_sdg_share=pubs_tf$total_sdgpubs/pubs_tf$total_pubs

sdg_share_tf=pubs_tf[c('COM', paste('share', paste('t', 1:4, sep=""), sep = '_'), 'total_sdg_share')]
sdg_share_tf[is.na(sdg_share_tf)]=0
sdg_share_tf <- cbind(sdg_share_tf, t(coef(fits))[,2]); names(sdg_share_tf)[-(1:6)] <- c("Slope")

#Analyse cutoff points 
p_t4 = data.frame(sort(sdg_share_tf$share_t4)) %>% rename(x = sort.sdg_share_tf.share_t4.) %>%
  ggplot() + geom_bar(aes(x=seq_along(x),y=x), stat='identity', width=1, color='black') +
  ggtitle('T4 SDG share') + labs(x='', y='SDG pub share') +  theme_classic() + theme(plot.title=element_text(hjust=0.5, size=14), plot.subtitle=element_text(size=12,hjust=0.5), 
                                                                                     axis.text.x=element_text(angle=90), legend.position = 'none') + 
  geom_hline(yintercept=T4_Cut, color='red', linetype='dashed', size=0.5) + ylim(0, 1)

p_total = data.frame(sort(sdg_share_tf$total_sdg_share)) %>% rename(x = sort.sdg_share_tf.total_sdg_share.) %>%
  ggplot() +  geom_bar(aes(x=seq_along(x),y=x), stat='identity', width=1, color='black') +
  ggtitle('Total SDG share') + labs(x='', y='SDG pub share') + theme_classic() + theme(plot.title=element_text(hjust=0.5, size=14), plot.subtitle=element_text(size=12,hjust=0.5), 
                                                                                       axis.text.x=element_text(angle=90), legend.position = 'none') + 
  geom_hline(yintercept=Total_Cut, color='red', linetype='dashed', size=0.5) + ylim(0,1)

p=ggarrange(p_t4, p_total, ncol=2, nrow=1, align='h')
print(p)
if (save_fig == 'yes'){png(paste(direct,'/Images/',figname,'.png',sep=''), width=600, height=300);print(p);dev.off()}
    
#Extra communities based on different criteria 
sdg_communities = as.character(merge(sdg_share_tf[(sdg_share_tf$Slope > 0) & (sdg_share_tf$total_sdg_share > Total_Cut), ], 
                                     sdg_share_tf[(sdg_share_tf$share_t4 > sdg_share_tf$share_t3) & (sdg_share_tf$share_t4 > T4_Cut), ], all.x = T, all.y = T)[,'COM']); 
sdg_pubs$sdg_com = ifelse(sdg_pubs$COM %in% sdg_communities, 'yes', 'no')
write.csv(sdg_pubs, paste(file_sdg,'sdg_com.csv',sep="_"), row.names=F)

sdg_communities=unique(sdg_pubs[sdg_pubs$sdg_com == 'yes','COM'])
all_pubs$sdg_com = ifelse(all_pubs$COM %in% sdg_communities, 'yes', 'no')
all_pubs = plyr::join(all_pubs,sdg_pubs[,c('UT','SDG')])
write.csv(all_pubs, paste(file_data,'sdg_com.csv',sep="_"), row.names=F)

setwd(direct)


###############################
### STEP 8: HCPC clustering ###
###############################

file='output_clean_sdg_com'
ClustSDG=6
ClustNoSDG=6

direct=getwd()
setwd(paste(direct,'Analysis',sep='/'))
library(FactoMineR)
data=as.matrix(read.csv(paste(file,'.csv',sep='')))
datas=data[data[,'sdg_com']=='yes',]
m=matrix(unlist(strsplit(unlist(lapply(1:dim(datas)[1],function(i){paste(rep(datas[i,'COM'],length(unlist(strsplit(datas[i,'SDG'],'#')))),unlist(strsplit(datas[i,'SDG'],'#')),sep='#')})),'#')),byrow=T,ncol=2)
m=table(m[,1],m[,2])
m=m[,colnames(m)!='NA']%*%t(m[,colnames(m)!='NA'])
SC=HCPC(CA(m,graph=F),graph=F,nb.clust=ClustSDG)$data.clust
datas=data[data[,'sdg_com']=='no',]
m=matrix(unlist(strsplit(unlist(lapply(1:dim(datas)[1],function(i){paste(rep(datas[i,'COM'],length(unlist(strsplit(datas[i,'SC'],'; ')))),unlist(strsplit(datas[i,'SC'],'; ')),sep='#')})),'#')),byrow=T,ncol=2)
m=table(m[,1],m[,2])
m=m[,colnames(m)!='NA']%*%t(m[,colnames(m)!='NA'])
NC=HCPC(CA(m,graph=F),graph=F,nb.clust=ClustNoSDG)$data.clust
m=cbind(c(rownames(SC),rownames(NC)),c(paste('Y',SC[,'clust'],sep=''),paste('N',NC[,'clust'],sep='')))
Cluster=unlist(lapply(data[,'COM'],function(i){m[m[,1]==i,2]}))

data=cbind(data,Cluster)

l=list()
for (clust in unique(data[,'Cluster'])){
  dat=data[data[,'Cluster']==clust,]
  sc=sort(table(unlist(sapply(strsplit(dat[,'SC'], ";"), function(x) trimws(x)))),decreasing=T)
  l[[clust]]=data.frame(SC_1=names(sc[1]),SC_2=names(sc[2]),SC_3=names(sc[3]))}

m=cbind(names(l),as.data.frame(matrix(unlist(l), nrow=length(l), byrow=TRUE)))
names(m)=c('Cluster','SC1','SC2','SC3')
data=plyr::join(as.data.frame(data[,1:20]),m,by='Cluster')

write.csv(data,paste(file,'.csv',sep=''), row.names=F)


