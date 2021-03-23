library(ggplot2)

gene_plot<-function(mat,align_id='All',pos=16857,add_pos=NULL,mark='align_id',point_size=10,front=NULL,behind=NULL){
  sum_pos<-c(pos)
  if(nrow(mat)<1){
    p<-ggplot()+aes(x=1,y=1)+geom_hline(yintercept = 1,linetype='dashed')+
      theme_classic()+
      theme(axis.title = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())+
      scale_y_continuous(labels = 'Intron',breaks = 1,limits = c(1))
    return(p)
  }
  if(length(align_id)==1){
    if(align_id!='All'){
      mat<-mat[mat[,11]%in%align_id,]
    }
  }
  if(length(align_id)>1){
    mat<-mat[test[,11]%in%align_id,]
  }
  #make plot
  st<-max(as.numeric(mat[,3]))
  ed<-max(as.numeric(mat[,4]))
  d1<-merge(pos,c(1:nrow(mat)))
  p<-ggplot()+
    theme_classic()+
    geom_segment(aes(x=d1$x,xend=d1$x,y=d1$y,yend=d1$y+.5))+
    geom_point(aes(x=d1$x,y=d1$y+.5),color='red',size=point_size)+
    geom_hline(yintercept = c(1:nrow(mat)),linetype='dashed')+
    theme(axis.title = element_blank(),axis.text.x = element_blank(),
          axis.ticks = element_blank(),axis.line = element_blank())+
    geom_segment(aes(x=as.numeric(mat[,3]),xend=as.numeric(mat[,4]),
                     y=c(1:nrow(mat)),yend=c(1:nrow(mat))),size=1)
    # scale_y_continuous(labels = mat$hg19.knownGene.alignID,breaks = c(1:nrow(mat)))
  
  if(!is.null(add_pos)){
    add_pos1=strsplit(add_pos,',')[[1]]
    d2<-cbind.data.frame(x=as.numeric(rep(add_pos1,each=nrow(mat))),y=rep(c(1:nrow(mat)),length(add_pos1)))
    p<-p+geom_segment(aes(x=d2$x,xend=d2$x,y=d2$y,yend=d2$y+.5),size=1)+
      geom_point(aes(x=d2$x,y=d2$y+.5),color='red',size=point_size)
    sum_pos<-c(sum_pos,add_pos1)
  }
  
  if(front==0&behind==0){}else{
    rang1<-range(as.numeric(sum_pos))
    p<-p+coord_cartesian(xlim=c(rang1[1]-front,rang1[2]+behind))
  }
  
  if(mark=='align_id'){
    p<-p+scale_y_continuous(labels = mat[,11],breaks = c(1:nrow(mat)),limits = c(0,nrow(mat)+.5) )
  }
  if(mark=='symbol'){
    p<-p+scale_y_continuous(labels = mat[,12],breaks = c(1:nrow(mat)),limits = c(0,c(nrow(mat)+.5) ) )
  }
  if(mark=='sum_id'){
    p<-p+scale_y_continuous(labels = sprintf('%s (%s)',mat[,12],mat[,11]),breaks = c(1:nrow(mat)),limits = c(0,c(nrow(mat)+.5) ) )
  }

  for(i in 1:nrow(mat)){
    x1<-mat[,8][i]
    x1<-as.numeric(strsplit(x1,',')[[1]])
    x2<-mat[,9][i]
    x2<-as.numeric(strsplit(x2,',')[[1]])
    d=data.frame(x1=x1,x2=x2,y1=(i-.3),y2=(i+.3))
    p<-p+
      geom_rect(data = d,mapping = aes(xmin=x1,xmax=x2,ymin=y1,ymax=y2),color='black',alpha=1)
  }
  
  return(p)
}

# test<-read.table('./data1/2020Year/rs_finder/Plot/hg19_genelist',sep='\t',header=T,stringsAsFactors = F)
# test<-test[test$hg19.knownGene.chrom=='chr1',]
# pos=14361
# forw1=14361-1000
# reve1=14361+1000
# 
# test<-test[test$hg19.knownGene.txStart<pos&test$hg19.knownGene.txEnd>pos,]
# gene_plot(test,pos=pos,front = 1000,behind=1000)
# 
# mat=test;pos=pos;front = 1000;behind=1000;mark='align_id';point_size=10;add_pos='14361,14561,14061'

