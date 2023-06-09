rm(list=ls(all=TRUE))
library(httr)
library(jsonlite)
library(uuid)

GRSFRequest="http://athena.ics.forth.gr:8890/sparql?default-graph-uri=&should-sponge=&format=application%2Fsparql-results%2Bjson&query=SELECT%20distinct%20%3Fuuid%20%3Furl%20%3Fsemantic_id%20%3Fshort_name%20%3Fgrsf_name%20%3Ftype%20%3Ftraceability_flag%20%3Fsdg_flag%20%3Fgeo_polygon%0AFROM%20%3Chttp%3A%2F%2Fgrsf%3E%0AWHERE%7B%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23type%3E%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2FStock%3E.%20%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2Fhas_uuid%3E%20%3Fuuid.%20%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2Fhas_public_catalog_url%3E%20%3Furl.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2FP1_is_identified_by%3E%20%3Fsemantic_id_iri.%0A%20%3Fsemantic_id_iri%20%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23type%3E%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2FE42_Identifier%3E.%0A%20%3Fsemantic_id_iri%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23label%3E%20%3Fsemantic_id.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2FP1_is_identified_by%3E%20%3Fgrsf_name_iri.%0A%20%3Fgrsf_name_iri%20%3Chttp%3A%2F%2Fwww.w3.org%2F1999%2F02%2F22-rdf-syntax-ns%23type%3E%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2FE41_Appellation%3E.%0A%20%3Fgrsf_name_iri%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23label%3E%20%3Fgrsf_name.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23label%3E%20%3Fshort_name.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2FP2_has_type%3E%20%3Ftype_iri.%0A%20%3Ftype_iri%20%3Chttp%3A%2F%2Fwww.w3.org%2F2000%2F01%2Frdf-schema%23label%3E%20%3Ftype.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2Fhas_traceability_flag%3E%20%3Ftraceability_flag.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2Fhas_sdg%3E%20%3Fsdg_flag.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2Fhas_polygon_from_dominant_record%3E%20%3Fgeo_polygon.%0A%20%3Fstock_iri%20%3Chttp%3A%2F%2Fwww.cidoc-crm.org%2Fcidoc-crm%2Fhas_status%3E%20%3Fstatus.%0A%20FILTER(%3Fstatus%3D%3Chttp%3A%2F%2Fwww.bluebridge-vres.eu%2Fgrsf%2Fstatus%2Fapproved%3E%20%7C%7C%20%3Fstatus%3D%3Chttp%3A%2F%2Fwww.bluebridge-vres.eu%2Fgrsf%2Fstatus%2Farchived%3E)%0A%20FILTER%20contains(%3Fsemantic_id%2C%22asfis%3A#CODE#%22)%0A%7D"

inputTable<-"ASFIS_sp_2019.txt"
#inputTable<-"ASFIS_sp_2019_redu.txt"
outputTable="ASFIS_sp_2019_geom.txt"
dataVessel<-read.csv(inputTable,header=T,sep=",")
polygons<-c()
counter = 0;
start_time <- Sys.time()



for (code in dataVessel$X3A_CODE){
  #code<-"PLE"
  status = (counter*100.0/length(dataVessel$X3A_CODE))
  if (counter%%10==0){
    end_time <- Sys.time()
    cat("status: ",status,"elapsed",(end_time-start_time),"\\n")
    start_time<-end_time
  }
  request = gsub("#CODE#", code, GRSFRequest)
  requestFile<-paste0("response_",UUIDgenerate(),".txt")
  got<-system(paste0("wget -q -O ",requestFile," ",request),ignore.stdout = T,ignore.stderr=T,intern = T, show.output.on.console=F)
  outJson<-read.csv(requestFile,header = F, sep = "#", stringsAsFactors=F, blank.lines.skip = TRUE, skipNul = T, as.is=T, quote="|")
  outJson<-paste(as.character(outJson$V1), collapse="",sep = "")
  answer<-fromJSON(txt=outJson)
  unlink(requestFile)
  nresults<-length(answer$results$bindings)
  if (nresults>0){
    cat("Received answer for",code,"\\n")
    polygonsWKT<-answer$results$bindings$geo_polygon$value
    multipol<-""
    for (pol in polygonsWKT){
      #POLYGON((1 1,5 1,5 5,1 5,1 1),(2 2, 3 2, 3 3, 2 3,2 2))
      polstring<-substr(pol,gregexpr(text=pol,pattern = "\\\\:")[[1]][2]+1,nchar(pol)-1)
      polstring<-gsub("\\\\],\\\\[", "#", polstring)
      polstring<-gsub(",", " ", polstring)
      polstring<-gsub("\\\\[", "", polstring)
      polstring<-gsub("\\\\]", "", polstring)
      polstring<-gsub("#", ",", polstring)
      if (nchar(multipol)>0)
        multipol<-paste0(multipol,",","(",polstring,")")
      else
        multipol<-paste0("(",polstring,")")
      
    }
    
    multipol<-paste0("POLYGON(",multipol,")")
    polygons<-c(polygons,multipol)
      
  }
  else{
    polygons<-c(polygons,"")
  }
  
  counter=counter+1
  
}

cat("Finished\\n")

dataVessel$geometry=polygons
write.csv(dataVessel,file=outputTable,row.names = F,na = "\\"\\"")
