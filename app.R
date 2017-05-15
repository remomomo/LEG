library(shiny)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(stringr)

load("data/impute_ranks_SOR.Rdata")
load("data/impute_ranks_combined.Rdata")
load("data/160216_promoters_ensemblTSS_v83.Rdata")
load("data/160216_velements.Rdata")
load("data/top10000_SOR.Rdata")
load("data/top10000_combined.Rdata")

# allowing uploads of bedfiles up to 15MB
options(shiny.maxRequestSize=15*1024^2)

# function to find the overlapping tiles with highest scores for mouse-coordinate input
get_best_overlapping<-function(predictions, ranges, return.coordinates=FALSE ,...){
  olaps<-findOverlaps(ranges, predictions, ... )
  if(length(olaps)==0) stop("None of the ranges supplied are rankable.")
  scores<-data.table(data.frame(scores=mcols(predictions)[subjectHits(olaps),]),mm10.chr=as.character(seqnames(predictions)[subjectHits(olaps)]),mm10.start=start(predictions)[subjectHits(olaps)],qh=queryHits(olaps))
  setnames(scores, c("score", "mm10.chr", "mm10.start", "qh"))
  scores<-scores[,list("score"=max(score), "mm10.chr"=mm10.chr[which.max(score)],"mm10.start"=mm10.start[which.max(score)]),by=qh]
  if (return.coordinates ){
    scores$mm10.end<-scores$mm10.start+2000
  }
  return(scores)
}

# function to find the overlapping tiles with highest scores for human-coordinate input
get_best_overlapping_hg<-function(predictions, ranges, return.coordinates=FALSE ,...){
  olaps<-findOverlaps(ranges, predictions, ... )
  if(length(olaps)==0) stop("None of the ranges supplied are rankable.")
  olranges<-GenomicRanges::ranges(olaps,ranges(ranges),ranges(predictions))
  # this is fancy as fuck...
  scores<-data.table( qh=queryHits(olaps),mm10=predictions$name[ subjectHits(olaps) ], score=predictions$score[ subjectHits(olaps) ],size=width(olranges) )
  scores<-scores[ ,list("size"=sum(size)),by=c("qh","mm10","score") ]
  scores<-scores[ ,list(score=score[which.max(size)],longest=mm10[which.max(size)],longest.width=size[which.max(size)]),by="qh"]
  chr.start<-extr_ranges(scores$longest)
  scores$mm10.chr<-chr.start$chr
  scores$mm10.start<-chr.start$start
  rm(chr.start)
  scores[,longest:=NULL]
  setnames(scores, c("qh","score", "overlap", "mm10.chr", "mm10.start"))
  if (return.coordinates ){
    scores$mm10.end<-scores$mm10.start+2000
  }
  return(scores)
}

extr_ranges<-function(x){
  chr<-str_extract(string = x, pattern="chr[0-9]*[XYM]*")
  start<-as.numeric( sapply( str_split(string = str_extract(string = x, pattern="[0-9]+-[0-9]+"), pattern="-"), '[[', 1) )
  return(data.frame(chr=chr, start=start,stringsAsFactors=FALSE))
}

extract.ranges<-function(x){
  chr<-str_extract(string = x, pattern="chr[0-9]*[XY]*")
  start.end<-as.numeric( gsub(",","",unlist(strsplit(str_extract(pattern = '[0-9,]+-[0-9,]+', x),"-"))))
  return(GRanges(chr, IRanges(start.end[1],start.end[2])))
}

reduce_fragments<-function(f){
  s<-data.table(chr=as.factor(seqnames(f)),data.frame(ranges(f)),name=f$name)
  s<-s[,reduce_dt(.SD),by=c("name","chr"),]
  s<-GRanges(s$chr, IRanges(s$start, s$end), name=s$name)
  return(s)
}

reduce_dt<-function(x){
  x<-x[ order(x$start, x$end) ]
  gaps<- x$start[-1] - x$end[ -nrow(x) ]
  if ( all(gaps <= 100) | length(gaps) == 0 ){
    return( data.frame(start=min(x$start), end=max(x$end)) )
  }else{
    big.gap<-which.max(gaps)
    if( sum(x$width[1:big.gap]) > sum(x$width[(big.gap+1):nrow(x)]) ){
      return( reduce_dt( x[1:big.gap]))
    }else{ return( reduce_dt( x[(big.gap+1):nrow(x)] ) ) }
  }
}


# the server:
server <- function(input, output, session) { 
  
  overchain<-reactive({
    
    if(input$genome=="mm9"){
      
      return( import.chain("data/mm9ToMm10.over.chain") )
      
    }else if(input$genome=="hg19"){
      
      return( import.chain("data/hg19ToHg38.over.chain") )
      
    }else{return(NULL)}
    
  })
  
  observeEvent(input$genome,{
    
    if ( input$genome %in% c("hg19","hg38") ){
      
      updateRadioButtons(session,inputId = "analysis_type",choices = c("Score short region(s)"),selected = "Score short region(s)",inline = TRUE, label = "Analysis Type" )
      updateSelectInput(session, inputId ="range.chr", choices = paste0("chr",c(1:22,"X","Y")))
      
    }else{
      
      updateRadioButtons(session = session,inputId = "analysis_type",choices = c("Scan for top","Score short region(s)"),inline = TRUE, label = "Analysis Type", selected = input$analysis_type )
      updateSelectInput(session, inputId ="range.chr", choices = paste0("chr",c(1:18,"X","Y")))
      
    }
    
  })
  
  observeEvent(input$submit.button,{
    
    updateTabsetPanel(session, "tabs",selected = "Results")
    
  })
  
  outputs<-eventReactive(input$submit.button,{
    
    withProgress(message=NULL,expr={
      
      setProgress(0.3, message = "Importing coordinates...")
      
      # the input is evaluated and a GRanges is constructed for further analysis:
      validate( need( input$range.start!=""|!is.null(input$bedfile)|input$range.gb!="", "Please supply coordinates.") )
      
      # checking if a bedfile was selected, if yes -> import
      if ( !is.null(input$bedfile) ){
        bedfile<-input$bedfile
        ranges<-import.bed(bedfile$datapath)
        if( !all(end(ranges)>=start(ranges))) end(ranges)[ start(ranges) > end(ranges) ]<-start(ranges)[ start(ranges) > end(ranges) ]
      }else if(input$range.start!=""){
        # if no bedfile was given, take the text input...
        ranges<-GRanges(input$range.chr, IRanges(as.numeric(input$range.start), as.numeric(input$range.end)))
      }else{
        # if no text input was given, take the copy paste coordinates...
        ranges<-extract.ranges(input$range.gb)
      }
      # naming...
      if ( length(ranges)>1 ){
        if ( !all(is.na(ranges$name)) ){
          # if a mix of named and unnamed ranges is supplied, call those without name 'unnamed.X'
          ranges$name[is.na(ranges$name)]<-paste("unnamed",1:sum(is.na(ranges$name)),sep=".")
        }else{
          # if all ranges are unnamed call them 'unnamed.X'
          ranges$name<-paste("unnamed",1:length(ranges),sep=".")
        }
        if(any(duplicated(ranges$name))){
          dups<- duplicated(ranges$name, fromlast=T) | duplicated(ranges$name, fromlast=F)
          new_names<-data.table(dup=ranges$name[dups], i=1:sum(dups))
          new_names<-new_names[,list(name=paste0(dup,'.',1:length(dup)), i), by=new_names$dup]
          ranges$name[dups][new_names$i]<-new_names$name   
        }
      }else{
        if (is.null(ranges$name)){ ranges$name<-"unnamed.1"
        } else if (is.na(ranges$name)) ranges$name<-"unnamed.1"
      }
      
      orig.ranges<<-ranges
      if (input$analysis_type=="Score short region(s)"){
        
        # resizing regions if necessary...
        ranges[width(ranges)<1000]<-resize(ranges[width(ranges)<1000], width = 1000, fix="center")
        ranges[width(ranges)>10000]<-resize(ranges[width(ranges)>10000], width = 10000, fix="center")
        
      }
      
      if( input$genome %in% c("mm9","hg19") ){ 
        
        # if coordinates given are mm9 liftOver to mm10
        # if cooridnates given are hg19 liftOver to hg38
        
        ranges<-unlist(liftOver(ranges, chain=overchain()))
        
        if ( length(ranges)==0 ){ stop("Can't process regions you provided: liftOver failed for all regions.") }
        
        # sometimes liftOver fucks up the ranges pretty badly, this tries to fix it as good as possible:
        if ( !all( table(ranges$name)==1 ) ){
          
          fragments<-which( duplicated(ranges$name) | duplicated(ranges$name, fromLast = TRUE) )
          f<-ranges[fragments]
          ranges<-ranges[-fragments]
          
          f<-reduce_fragments(f)
          if (length(ranges$score)!=0) f$score<-0
          
          ranges<-c(ranges, f)
          ranges[width(ranges)<1000]<-resize(ranges[width(ranges)<1000], width = 1000, fix="center")
        }
        
      }
      
      names(ranges)<-ranges$name
      setProgress(0.5, message="Retrieving scores...")
      
      # depending on the analysis type, we submit the ranges to one of two methods, that calculate the outputs:
      
      ranges<-keepSeqlevels( ranges, value = paste0("chr",c(1:22, "X", "Y")) )
      
      if (input$analysis_type=="Scan for top"){return(find.highranking(ranges))
      }else return(rank.regions(ranges))
    })
    
  })
  
  
  find.highranking<-function(ranges){
    
    # this is the first method, which simply finds overlaps with the top 10'000 predictions
    
    if (input$Method == "Combined Model" ){ 
      top.preds<-top.10000.combined
      pred.type<-"predictions.combined"
    }else{ 
      top.preds<-top.10000.SOR
      pred.type<-"predictions.SOR"
    }
    
    message('ranges:')
    print(ranges)
    
    message('orig.ranges:')
    print(orig.ranges)
    
    olaps<-findOverlaps(ranges, top.preds)
    
    message('olaps:')
    print(olaps)
    
    validate( need( length(olaps)!=0 , "No high-ranking limb-enhancers overlap the regions you provided.") )
    
    # getting the overlapping top predictions and scores:
    hits<-top.preds[ subjectHits(olaps) ]
    
    if ( input$Method == "Combined Model" ){
      gw.pred<<-import.bedGraph("./data/predictions_Ridge.bdg.gz", which=hits)
    }else{
      gw.pred<<-import.bedGraph("./data/predictions_SOR.bdg.gz", which=hits)
    }
    
    message('hits:')
    print(hits)
    
    setProgress(0.8, message = "Ranking, constructing metadata...")
    
    s<-get_best_overlapping(predictions=gw.pred, hits, minoverlap=1999)
    message('s:')
    print(s)
    
    hits$score<-s$score
    hits$rank<-hits$reduced.rank
    hits$reduced.rank<-NULL
    
    qh<-queryHits(olaps)
    
    # small trick if it's mm9...
    if ( input$genome == "mm9" ){
      
      if (input$Method == "Combined Model" ){
        
        load("data/top10000_combined_mm9.Rdata")
        mm9.hits<-top.10000.combined.mm9[ match( hits$rank, top.10000.combined.mm9$reduced.rank) ]
        
      }else{
        
        load( "data/top10000_SOR_mm9.Rdata" )
        mm9.hits<-top.10000.SOR.mm9[ match( hits$rank, top.10000.SOR.mm9$reduced.rank) ]
        
      }
      
      message('mm9.hits:')
      print(mm9.hits)
      
      results_table<-data.frame(input.region = names(ranges)[qh], chr=seqnames(mm9.hits),start=start(mm9.hits),end=end(mm9.hits),data.frame(mcols(hits)))
      message('results table:')
      print(head(results_table,10))
      
    }else results_table<-data.frame(input.region = names(ranges)[qh], chr=seqnames(hits),start=start(hits),end=end(hits),data.frame(mcols(hits)))
    
    
    # making a summary
    # quantiles<-data.frame( "genome.wide.score.quantiles"=rev(quantile(gw.pred$score, c(0.5,0.9,0.95,0.99,0.999,1))) )
    best<-which.max(results_table$score)
    
    summary_table<-data.frame("n.regions.submitted"=length(orig.ranges), "n.overlapped"=length(hits), best.score=results_table$score[best], best.rank=results_table$rank[best], best.id=results_table$input.region[best] )
    row.names(summary_table)<-NULL
    
    results_table<-results_table[ order(results_table$score, na.last=TRUE, decreasing=TRUE), ]
    
    setProgress(1, message=NULL)
    # return(list(results_table,summary_table,quantiles))
    
    return(list(results_table,summary_table))
    
  }
  
  rank.regions<-function(ranges){
    
    # this is the second method, that ranks any set of ranges regardless if they overlap promoters, the training set, etc.
    # the genome wide ranks that are returned are 'imputed', based on the given top 10'000 ranks
    
    if (input$Method == "Combined Model" ){ 
      top.preds<-top.10000.combined
      pred.type<-"predictions.combined"
      imputefun<-impute.ranks.combined_pred
      if ( input$genome %in% c("mm9","mm10") ){ gw.pred<<-import.bedGraph("./data/predictions_Ridge.bdg.gz", which=ranges)
      }else gw.pred<<-import.bed("./data/hg38_predictions_Ridge.sorted.bed.gz", which=ranges)
      
    }else{ 
      top.preds<-top.10000.SOR
      pred.type<-"predictions.SOR"
      imputefun<-impute.ranks.SOR
      if ( input$genome %in% c("mm9","mm10") ){ gw.pred<<-import.bedGraph("./data/predictions_SOR.bdg.gz", which=ranges)
      }else gw.pred<<-import.bed("./data/hg38_predictions_SOR.sorted.bed.gz", which=ranges)
    }
    
    if ( input$genome %in% c("mm9","mm10") ){
      
      olaps.top<-findOverlaps(ranges, top.preds)
      olaps.promoter<-findOverlaps(ranges, promoters)
      olaps.training<-findOverlaps(ranges, velements)
      
      ranges$overlaps.top10000<-1:length(ranges) %in% queryHits(olaps.top)
      ranges$overlaps.promoter<-1:length(ranges) %in% queryHits(olaps.promoter)
      ranges$overlaps.training<-1:length(ranges) %in% queryHits(olaps.training)
      
      # getting the scores:
      hits<-ranges
      rm(ranges)
      
      scoring.index<-unique(queryHits(findOverlaps(hits, gw.pred, minoverlap=1000)))
      hits$score<-0
      
      best.overlaps<-get_best_overlapping(predictions=gw.pred, hits[scoring.index], minoverlap=1000, return.coordinates = TRUE)
      
      hits$mm10.chr<-NA
      hits$mm10.start<-NA
      hits$mm10.end<-NA
      
      hits$score[scoring.index]<-best.overlaps$score
      hits$mm10.chr[scoring.index]<-best.overlaps$mm10.chr
      hits$mm10.start[scoring.index]<-best.overlaps$mm10.start
      hits$mm10.end[scoring.index]<-best.overlaps$mm10.end
      
      setProgress(0.8, message="Ranking, constructing metadata...")
      hits$imputed.gw.rank<-sapply(hits$score, imputefun)
      hits$imputed.gw.rank[ hits$imputed.gw.rank == 10000 ]<-NA
      hits$rank.within.set<--1*(rank(hits$score)-length(hits)-1)
      
      # for mapping back to the original ranges... this should be included! ! ! ! !
      orig.ranges<-orig.ranges[match(names(hits), orig.ranges$name)]
      
      hits$name<-NULL
      
      print(head(hits))
      
      # constructing table:
      results_table<-data.frame(input.region = names(hits), chr=seqnames(orig.ranges), start=start(orig.ranges), end=end(orig.ranges), data.frame(mcols(hits))[,c("score","rank.within.set","imputed.gw.rank","mm10.chr","mm10.start","mm10.end","overlaps.top10000","overlaps.promoter","overlaps.training")])
      setnames(results_table, c("chr","start","end"), paste(input$genome, c("chr","start","end"), sep=".") )
      
      # making a summary
      best<-which.max(results_table$score)
      summary_table<-data.frame("n.regions.submitted"=length(hits), best.score=results_table$score[best], best.rank=results_table$imputed.gw.rank[best], best.id=results_table$input.region[best], average.score=mean(results_table$score) )
      row.names(summary_table)<-NULL
      
      results_table<-results_table[ order(results_table$score, na.last=TRUE, decreasing=TRUE), ]
      
    }else{
      
      hits<-ranges
      rm(ranges)
      
      scoring.index<-unique(queryHits(findOverlaps(hits, gw.pred)))
      hits$score<-0
      
      best.overlaps<-get_best_overlapping_hg(predictions=gw.pred, hits[scoring.index], return.coordinates = TRUE)
      
      hits$mm10.chr<-NA
      hits$mm10.start<-NA
      hits$mm10.end<-NA
      hits$lifted.bp<-NA
      
      hits$score[scoring.index]<-best.overlaps$score
      hits$mm10.chr[scoring.index]<-best.overlaps$mm10.chr
      hits$mm10.start[scoring.index]<-best.overlaps$mm10.start
      hits$mm10.end[scoring.index]<-best.overlaps$mm10.end
      hits$lifted.bp[scoring.index]<-best.overlaps$overlap
      
      hitsmm10<-GRanges(seqnames=as.character( best.overlaps$mm10.chr) ,ranges=IRanges(start=best.overlaps$mm10.start,end=best.overlaps$mm10.end))
      
      olaps.top<-findOverlaps(hitsmm10, top.preds)
      olaps.promoter<-findOverlaps(hitsmm10, promoters)
      olaps.training<-findOverlaps(hitsmm10, velements)
      
      hits$overlaps.top10000<-NA
      hits$overlaps.promoter<-NA
      hits$overlaps.training<-NA
      
      hits$overlaps.top10000[scoring.index]<-1:length(hitsmm10) %in% queryHits(olaps.top)
      hits$overlaps.promoter[scoring.index]<-1:length(hitsmm10) %in% queryHits(olaps.promoter)
      hits$overlaps.training[scoring.index]<-1:length(hitsmm10) %in% queryHits(olaps.training)
      
      setProgress(0.8, message="Ranking, constructing metadata...")
      hits$imputed.gw.rank<-sapply(hits$score, imputefun)
      hits$imputed.gw.rank[ hits$imputed.gw.rank == 10000 ]<-NA
      hits$rank.within.set<--1*(rank(hits$score)-length(hits)-1)
      
      # for mapping back to the original ranges:
      orig.ranges<-orig.ranges[match(names(hits), orig.ranges$name)]
      hits$name<-NULL
      
      # constructing table:
      results_table<-data.frame(input.region = names(hits), chr=seqnames(orig.ranges),start=start(orig.ranges),end=end(orig.ranges),data.frame(mcols(hits))[,c("lifted.bp","score","rank.within.set","imputed.gw.rank","mm10.chr","mm10.start","mm10.end","overlaps.top10000","overlaps.promoter","overlaps.training")])
      setnames(results_table, c("chr","start","end"), paste(input$genome, c("chr","start","end"), sep=".") )
      
      # making a summary
      best<-which.max(results_table$score)
      summary_table<-data.frame("n.regions.submitted"=length(hits), best.score=results_table$score[best], best.rank=results_table$imputed.gw.rank[best], best.id=results_table$input.region[best], average.score=mean(results_table$score) )
      row.names(summary_table)<-NULL
      
      results_table<-results_table[ order(results_table$score, na.last=TRUE, decreasing=TRUE), ]
      
    }
    
    setProgress(1, "The Genie has granted your wish!")
    return(list(results_table,summary_table))
    
  }
  
  output$table<-renderDataTable(outputs()[[1]], options = list(scrollX = TRUE))
  output$summary<-renderTable(outputs()[[2]],digits=8)
  
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
    filename = function() {
      paste("EnhancerGenieResults", input$filetype, sep = ".")
    },
    
    # This function writes data to a file given to it by the argument 'file'.
    content = function(file){
      if ( input$filetype=="bed" ){
        
        dat<-outputs()[[1]]
        
        # selecting the top 'n' predictions:
        if ( input$num == 0 ){ 
          i<-1:nrow(dat)
        }else i<-1:min(as.numeric(input$num),nrow(dat))
        
        dat<-dat[i,]
        
        # setting which coordinates should be used for output:
        if ( input$analysis_type == "Score short region(s)"){
          if ( input$coordinate.output == "Predictions (mm10)" ){
            dat<-dat[ !is.na(dat$mm10.start), ]
            out.bed<-GRanges(dat$mm10.chr, IRanges(dat$mm10.start, dat$mm10.end), score=dat$score)
          }else out.bed<-GRanges(dat[,paste0(input$genome,".chr")], IRanges(dat[,c(paste0(input$genome,".start"))], dat[,c(paste0(input$genome,".end"))]), score=dat$score)
        }else out.bed<-GRanges(dat$chr, IRanges(dat$start, dat$end), score=dat$score)
        out.bed$score<-round(1000*out.bed$score, 0)
        names(out.bed)<-dat$input.region
        export.bed( out.bed, file )
      }
      
      if ( input$filetype=="csv"){
        
        dat<-outputs()[[1]]
        
        if ( input$num == 0 ){
          i<-1:nrow(dat)
        }else i<-1:min(as.numeric(input$num),nrow(dat))
        write.csv(dat[i,], file, row.names=FALSE, quote=FALSE)
      }
    }
  )
}

######### the User Interface:
ui <- shinyUI(fluidPage(
  
  titlePanel("LEG - the Limb Enhancer Genie"),
  
  sidebarLayout(
    
    sidebarPanel(
      radioButtons(inputId = "analysis_type",choices = c("Scan for top","Score short region(s)"),selected = "Scan for top",inline = TRUE, label = "Analysis Type"),
      radioButtons(inputId = "genome",choices = c("mm10","mm9","hg38","hg19"),selected = "mm10",inline = TRUE, label = "Input Genome"),
      radioButtons(inputId = "Method",choices = c("Combined Model","Sum of Ranks"),selected = "Combined Model",inline = TRUE, label = "Method"),
      selectInput("range.chr","chromosome",choices = paste0("chr",c(1:18,"X","Y"))),
      fixedRow(column(width=6,textInput("range.start","start")),column(width=6,textInput("range.end","end"))),
      textInput("range.gb","... or copy-paste from Genome Browser:", value = "chr13:2,281,724-28,912,023"),
      fileInput(inputId = "bedfile",label = "... or upload BED"),
      actionButton(inputId = "submit.button", label = "Submit")
    ),
    
    mainPanel(
      tags$body(tags$script(src="iframeResizer.contentWindow.min.js")),
      tabsetPanel(id = "tabs",
                  tabPanel("Documentation",
                           includeMarkdown("title_page.Rmd")),
                  tabPanel("Results", dataTableOutput("table")),
                  # tabPanel("Plot", plotOutput("plot")),
                  # tabPanel("Summary", tableOutput("summary"), tableOutput("quantiles")),
                  tabPanel("Summary", tableOutput("summary")),
                  tabPanel("Download",
                           radioButtons("filetype", "File type:",
                                        choices = c("bed", "csv")),
                           conditionalPanel('input.filetype == "bed" && input.analysis_type == "Score short region(s)" ',
                                            radioButtons("coordinate.output", "BED-Coordinates:", c("Input Regions", "Predictions (mm10)"))),
                           numericInput("num", label = "Download best n (0 means all)", value = 0, min = 0, width = "120px", step = 25),
                           downloadButton('downloadData', 'Download'))
      )
    )
  )
))

shinyApp(ui = ui, server = server) # this launches your app