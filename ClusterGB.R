library(R6)

ClusterGB <- R6Class("ClusterGB",
                public = list(
                  GbTable = NULL,
                  initialize = function(filetxt=NULL  ,fileName = NA) {
                    if(is.null(filetxt)){
                    private$fileText <- scan(fileName, what = "", sep = "\n", quiet = TRUE)
                    self$GbTable <- private$getGeneDataframe()
                    } else{
                      private$fileText <- filetxt
                      self$GbTable <- private$getGeneDataframe()
                    }
                    
                  },
                  SequenceList = function(string=T,toUpper=F){
                    CDS <- self$GbTable[["CDS"]]
                    seq <- private$getWholeSequence()
                    genes <- lapply(CDS,function(x,seq){
                      pos <-as.numeric(strsplit( 
                        gsub( '[^0-9\\.]', '', x ), '\\.\\.' )[[1]])
                      subseq <- substr(seq,pos[1],pos[2])
                      if(regexpr("complement",x)==-1){
                          return(subseq)
                      }else {
                          #reverse string
                          subseq <- sapply(lapply(strsplit(subseq, NULL), rev), 
                                           paste, collapse="")
                          return(chartr('atgc','tacg',subseq))
                        }
                    },seq=seq)
                    names(genes) <- self$GbTable[["db_xref"]]
                    if(toUpper) genes <- sapply(genes,toupper)
                    if(string) return(genes)
                    else return(sapply(genes,strsplit,split=NULL))
                    #else return(sapply())
                    
                  }
                  ),
                
                private = list(
                  fileText = NULL,
                  getGeneFeature = function(block){
                    
                    block <- unlist(block)
                    tag_start <- regexpr("/",block[grep("/",block)[1]])[1]
                    
                    
                    some_feature <- block[1]
                    value <- substr(some_feature,tag_start,nchar(some_feature))
                    FeatureValue <- list(value)
                    
                    geneFeature = c("/db_xref=\""," /function="," /product="," /gene=",
                                    " /protein_id="," /note=")
                    for(i in geneFeature){
                      
                      f_loc <- grep(i,block)[1]
                      
                      if(is.na(f_loc)){ 
                        FeatureValue <- append(FeatureValue,"no tag")
                      } else  {
                        value <- substr(block[f_loc],(tag_start+nchar(i)),
                                        nchar(block[f_loc])-1)
                        k <- 1
                        
                        if(f_loc<length(block)){
                          while((regexpr( "/",block[f_loc+k])==-1)&
                                (regexpr( "\\.\\.",block[f_loc+k])==-1 )){
                            value <- paste(value,(substr(block[f_loc+k],tag_start,
                                                         nchar(block[f_loc+k]))))
                            k=k+1
                            if((f_loc+k)>= length(block)) break
                          }
                        }
                        value <- gsub(pattern = '[;]',',',value)
                        value <- gsub(pattern = '[\n]',' ',value)
                        FeatureValue <- append(FeatureValue,value)
                      }
                    }
                    ifelse(regexpr("complement",FeatureValue[1])==-1,
                           FeatureValue <- append(FeatureValue,"+"),
                           FeatureValue <- append(FeatureValue,"-"))
                    
                    names(FeatureValue) = c("CDS","db_xref","Function","Product","Gene",
                                            "Protein_id","Note","Strand_of_DNA")
                    
                    return(FeatureValue)
                    
                  },
                  getGenBloks = function(){
                    edge_tags <- c("  CDS  ","  gene  ","  misc_feature  ","  regulatory  ")
                    edge_tags_pos <- unlist(sapply(edge_tags, grep, x=private$fileText))
                    seq_pos <- gsub( '[^0-9\\.<>]', '', private$fileText[edge_tags_pos])
                    seq_posU <- unique(seq_pos)
                    
                    bInex <- sapply(seq_posU,function(x,edge_tags_pos,seq_pos){
                      return( min(edge_tags_pos[x == seq_pos]))
                    },edge_tags_pos = edge_tags_pos, seq_pos = seq_pos )
                    
                    bInex <- c(bInex,grep("ORIGIN",private$fileText)[T])
                    
                    gBlocks<- lapply(1:(length(bInex)-1),function(x,bInex,clst){
                      return(clst[bInex[x]:(bInex[x+1]-1)]) },bInex=bInex,clst=private$fileText)
                    return(gBlocks)
                  },
                  getGeneDataframe = function(){
                    df <- sapply(private$getGenBloks(),private$getGeneFeature)
                    df1 <- data.frame(t(df))
                    FeatureValue <- c()
                    clusterFeatures <- c("ORGANISM","DEFINITION","LOCUS")
                    for(i in clusterFeatures){
                      
                      f_loc <- grep(i,private$fileText)[1]
                      
                      some_feature <- private$fileText[f_loc]
                      if(is.na(some_feature)){ 
                        FeatureValue <- append(FeatureValue,"no tag")
                      } else  {
                        value <- substr(some_feature,(regexpr(i,
                                                              some_feature)[1]+nchar(i)),
                                        nchar(some_feature))
                        FeatureValue <- append(FeatureValue,value)
                      }
                    }
                    
                    df1$DEFINITION <- FeatureValue[[2]]
                    df1$ORGANISM <- FeatureValue[[1]]
                    df1$LOCUS <- FeatureValue[[3]]
                    return(df1)
                  },
                  getWholeSequence =  function(){
                    if(!is.character(private$fileText)) {
                      stop("Error, file_string is not character!")
                      return(NULL)}
                    
                    firstN <- grep("^ {0,}ORIGIN", private$fileText) + 1
                    lastN <- which(private$fileText == "//") - 1
                    tmp <- gsub("[[:digit:] ]", "", private$fileText[firstN:lastN])
                    return(paste0(tmp,collapse = ''))
                  }
                )
)

ClusterSetsGB <- R6Class("ClusterSetsGB",
                     inherit = ClusterGB,
                     public = list(
                       ClusterGBList = list(),
                       GbTableList = list(),
                       fileList = list(),
                       initialize = function(fileList = NA){
                         self$fileList <- fileList
                         clusters <- private$getClusters()
                         self$ClusterGBList <- lapply(clusters,ClusterGB$new)
                         self$GbTable <- data.frame(t(
                           data.frame(sapply(self$ClusterGBList,function(x){
                             return(t(x$GbTable))
                           }))))
                       },
                       SequenceList = function(){
                         seq <- unlist(lapply(self$ClusterGBList,function(x){return(x$SequenceList())}))
                         names(seq) <- self$GbTable[,2]
                         return(seq)
                       }
                       ),
                     private = list(
                       
                       getClusters = function(){
                         clusters <- list()
                         for(fileName in self$fileList){
                           
                           fileText <- scan(fileName, what = "", sep = "\n", quiet = TRUE)
                           locusPos <- grep("LOCUS",fileText)
                           locusPos <- append(locusPos,length(fileText)+1)
                           
                           i <- 1
                           while(i<length(locusPos)){
                             clusters <- append(clusters,
                                                list(fileText[locusPos[i]:(locusPos[i+1]-1)]))
                             i <- i+1
                             
                           }
                         }
                         return(clusters)
                         
                       }
                       
                     )
)

