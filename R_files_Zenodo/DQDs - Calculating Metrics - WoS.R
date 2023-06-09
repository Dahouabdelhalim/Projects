library(boot)        	        # for bootstrapping samples
library(countrycode)        	# for identifying countries' names
library(zipcode)              # for identifying the US's zip codes
library(ISOcodes)             # for identifying ISO codes
library(maps)                 # for identifying cities'names 
library(hunspell)             # for checking spelling in English


min_time <- 1980 
max_time <- 2017 


substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))}

substrLeft = function(text, num_char) {
  substr(text, 1, num_char)
}

trim_org = function(text) {
  left_part <- text
  last_char <- substrRight(text,1)
  while (last_char %in% c(","," "))
  {
    left_part <- substrLeft(left_part, nchar(left_part)-1)
    last_char <- substrRight(left_part,1)
  }
  return(left_part)
}

get_from_affiliation = function(affiliation) {
  
  temp_org <- NULL
  temp_country <- "" 
  
  if ((!anyNA(affiliation$name)) & (!is.null(affiliation$name)))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$name} 
    else {temp_org <- paste(temp_org, affiliation$name, sep = ", ")}
  }
  if (!anyNA(affiliation$city)&!is.null(affiliation$city))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$city} 
    else {temp_org <- paste(temp_org, affiliation$city, sep = ", ")}
  }
  if (!anyNA(affiliation$city_id)&!is.null(affiliation$city_id))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$city_id} 
    else {temp_org <- paste(temp_org, affiliation$city_id, sep = ", ")}
  }
  if (!anyNA(affiliation$country)&!is.null(affiliation$country))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$country} 
    else 
    {
      temp_org <- paste(temp_org, affiliation$country, sep = ", ")
      temp_country <- affiliation$country 
    }
  }
  if (!anyNA(affiliation$country_code)&!is.null(affiliation$country_code))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$country_code} 
    else {temp_org <- paste(temp_org, affiliation$country_code, sep = ", ")}
  }
  if (!anyNA(affiliation$state)&!is.null(affiliation$state))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$state} 
    else {temp_org <- paste(temp_org, affiliation$state, sep = ", ")}
  }
  if (!anyNA(affiliation$state_code)&!is.null(affiliation$state_code))
  {
    if (is.null(temp_org)) {temp_org <- affiliation$state_code} 
    else {temp_org <- paste(temp_org, affiliation$state_code, sep = ", ")}
  }
  
  if (length(temp_org)==0){temp_org <- "Empty value"}
  else {temp_org <- trim_org(temp_org)}
  
  my_list <- list("org" = temp_org, "country" = temp_country)
  return(my_list) 
}

istheUS <- function(name,string){
  data(zipcode)
  if (is.na(name)) {return(FALSE)} else
  {
    if ((length(state.abb[which(state.name == name)])==1)|(name %in% state.abb) | (name=='USA') | (name =='U.S.A') | (name =='U.S.A.') |(substrRight(string,13)=='United States'))  
    { 
      return(TRUE)} else
      {
        if (grepl(' ',name))
        {
          last_space_position <- gregexpr(' ',name)
          last_space_position_value <- last_space_position[[1]][length(last_space_position[[1]])]
          z_code <- trimws(substrRight(name, nchar(name)-last_space_position_value))
        }
        else
        {z_code <- name}
        
        zip_code <- clean.zipcodes(z_code)
  
        if (!is.na(match(zip_code,zipcode$zip)))
        {
          s_abb <- zipcode$state[match(zip_code,zipcode$zip)]
          city_name <- zipcode$city[match(zip_code,zipcode$zip)]
          state_name <- state.name[grep(s_abb, state.abb)]
          
          if (length(state_name)==0) 
          {(grepl(paste(s_abb,' ',z_code),name)) | grepl(city_name,string)}
          else
          {(grepl(paste(s_abb,' ',z_code),name)) | grepl(city_name,string) | grepl(state_name,string)}
        }
  
        else
        {
          if (gsub("[ ]","",gsub("[.]","",name)) %in% state.abb)
          {
            s_abb <- gsub("[ ]","",gsub("[.]","",name))
            city_name <- zipcode$city[match(s_abb,zipcode$city)]
            state_name <- state.name[grep(s_abb, state.abb)]
            if (!is.na(city_name))
            {
              grepl(city_name,string) | grepl(state_name,string)
            }
            else
            {
              grepl(state_name,string)
            }
          }
          else
          {
            grepl("Washington,DC",gsub("[.]","",gsub("[ ]","",string)))
          }
        }
      }
  }
}

istheUK <- function(name){
  UK_list <- c('UK', 'U.K', 'U.K.', 'England', 'ENGLAND', 'Scotland', 'SCOTLAND','Wales', 'WALES', 'Northern Ireland', 'NORTHERN IRELAND','Great Britain','GREAT BRITAIN')
  if (is.na(name)) {return(FALSE)}
  else {return((name %in% UK_list) | (substrRight(name,15)=='United Kingdoms'))}
}

isChina_Taiwan <- function(name, string){
  if (is.na(string))  { return('None')} 
  else {
    if (((substrRight(name,3)=='ROC') |(substrRight(name,5)=='R.O.C') | (substrRight(name,5)=='R.O.C.')) & (grepl('Taiwan',string))) {return ('TW')}         
    else {
      if ((substrRight(name,11)=='Taiwan, ROC')|(substrRight(name,13)=='Taiwan, R.O.C') | (substrRight(name,14)=='Taiwan, R.O.C.'))
      { return ('TW')}
      else { 
        if ((substrRight(name,5)==', PRC') | (substrRight(name,7)==', P.R.C') | (substrRight(name,9)==' PR China')|(substrRight(name,11)==' P.R. China')) {return ('CN')}
        else {return ('None')} 
      }         
    }      
  }  
}

affiliation_data_preparation <- function (my_string){
  
  country_name <- "" 
  country_name2 <- ""
  
  my_string <- gsub("Univ.", "University", my_string)
  my_string <- gsub("Sch.", "School", my_string)
  my_string <- gsub("Coll.", "College", my_string)
  my_string <- gsub("Dept.", "Department", my_string)
  my_string <- gsub("Lab.", "Laboratory,", my_string)
  my_string <- gsub("Acad.", "Academy", my_string)
  my_string <- gsub("Inst.", "Institute", my_string)
  
  if (grepl('vcard_country-name=',my_string))
  {
    position <- gregexpr('vcard_country-name=',my_string)
    position_value <- position[[length(position)]][1]
    right_part <- trimws(substrRight(my_string, nchar(my_string)-position_value-18))
    position <- gregexpr(',',right_part)
    position_value <- position[[1]][1]
    country_name <- trimws(substrLeft(right_part, position_value-1))
  }
  else
  {
    if ((grepl('email:',my_string,ignore.case=TRUE)) | (grepl('e-mail:',my_string,ignore.case=TRUE)) | (grepl('emails:',my_string,ignore.case=TRUE)) |(grepl('e-mails:',my_string,ignore.case=TRUE)))
    {
      position_value1 <- 0
      position_value2 <- 0
      position_value3 <- 0
      position_value4 <- 0
      if (grepl('email:',my_string,ignore.case=TRUE)) 
      {position1 <- gregexpr('email:',my_string,ignore.case=TRUE)
      position_value1 <- position1[[1]][length(position1[[1]])]}
      if (grepl('e-mail:',my_string,ignore.case=TRUE)) 
      {position2 <- gregexpr('e-mail:',my_string,ignore.case=TRUE)
      position_value2 <- position2[[1]][length(position2[[1]])]}
      if (grepl('emails:',my_string,ignore.case=TRUE)) 
      {position3 <- gregexpr('emails:',my_string,ignore.case=TRUE)
      position_value3 <- position3[[1]][length(position3[[1]])]}
      if (grepl('e-mails:',my_string,ignore.case=TRUE)) 
      {position4 <- gregexpr('e-mails:',my_string,ignore.case=TRUE)
      position_value4 <- position4[[1]][length(position4[[1]])]}
      position_value <- max(position_value1,position_value2,position_value3,position_value4)
      left_part <- trimws(substrLeft(my_string,position_value-1))
      last_char <- substrRight(left_part,1)
      
      while (last_char %in% special_string_endings)
      {
        left_part <- substrLeft(left_part, nchar(left_part)-1)
        last_char <- substrRight(left_part,1)
      }
      
      position <- gregexpr(',',left_part)
      position_value <- position[[1]][length(position[[1]])]
      country_name <- trimws(substrRight(left_part, nchar(left_part)-position_value))
    }
    else {   
      
      my_string <- gsub("#TAB#","",my_string) 
      my_string <- gsub(", EU","",my_string) 
      my_string <- gsub("#R#","",my_string) 
      my_string <- gsub("#N#","",my_string) 
      last_char <- substrRight(my_string,1)
      
      while (last_char %in% special_string_endings)
      {
        my_string <- substrLeft(my_string, nchar(my_string)-1)
        last_char <- substrRight(my_string,1)
      }
      
      position <-gregexpr(",",my_string)
      position_value <- position[[1]][length(position[[1]])]
      
      pcountry_value <- position_value
      
      if (pcountry_value[1]==-1)
      {
        pcountry <-gregexpr(";",my_string)
        pcountry_value <- pcountry[[1]][length(pcountry[[1]])]
      }
      if (pcountry_value[1]==-1)
      {
        pcountry <-gregexpr("\\\\|",my_string)
        pcountry_value <- pcountry[[1]][length(pcountry[[1]])]
      }
      
      country_name <- trimws(substrRight(my_string, nchar(my_string)-pcountry_value)) 
      
      if ((nchar(my_string)-pcountry_value)>0)
      {
        
        my_string_left <- trimws(substrLeft(my_string, pcountry_value-1))     
        pcountry2 <-gregexpr(",",my_string_left)
        pcountry2_value <- pcountry2[[1]][length(pcountry2[[1]])]
        
        if (pcountry2_value==-1)
        {
          pcountry2 <-gregexpr(";",my_string_left)
          pcountry2_value <- pcountry2[[1]][length(pcountry2[[1]])]
        }
        if (pcountry_value[1]==-1)
        {
          pcountry <-gregexpr("\\\\|",my_string)
          pcountry_value <- pcountry[[1]][length(pcountry[[1]])]
        }
        
        country_name2 <- trimws(substrRight(my_string_left, nchar(my_string_left)-pcountry2_value))
      }
    }
  }
  my_list <- list("string"=my_string, "name"=country_name, "name2"=country_name2)
  return(my_list)
}

string_matching <- function(country_name,my_string){
  
  country_code  <- NA 
  
  if (!is.na(countrycode(country_name, 'iso2c', 'country.name'))){
    country_code <- as.character(country_name)
  } else { country_code <- countrycode(country_name, 'country.name', 'iso2c')}
  
  if (is.na(country_code))      	
  {
    if (istheUS(country_name,my_string)) 
    {
      country_code <- 'US'
    } else
    {
      if (istheUK(country_name)) 
      {
        country_code <- 'GB'
      }
      else
      {
        China_Taiwan <- isChina_Taiwan(country_name,my_string) 
        if (!(China_Taiwan=='None') ){
          if (China_Taiwan=='CN') {country_code <- 'CN'}
          if (China_Taiwan=='TW') {country_code <- 'TW'}
        }
        else
        {
          country_code <- NA
        } 
      }
    }
  }
  return(country_code)
}


Wikidata_query <- function(country_name2, country_name,my_string){
  
  country_code  <- NA 
  
  my_query <- paste('https://query.wikidata.org/sparql?format=json&query=',RCurl::curlEscape(paste('PREFIX
                                                                                                   schema: <http://schema.org/> PREFIX wdt:
                                                                                                   <http://www.wikidata.org/prop/direct/> SELECT ?countryLabel WHERE
                                                                                                   {<https://en.wikipedia.org/wiki/',gsub(' ','_',gsub("\\\\|","",gsub("\\\\\\\\","",my_string))),'> schema:about
                                                                                                   ?datalink. ?datalink wdt:P17 ?country. SERVICE wikibase:label {
                                                                                                   bd:serviceParam wikibase:language "en" .}}',sep='')),sep='')
  
  
  mod2 <- try({
    country_label <- fromJSON(url(my_query))
  }, TRUE)
  
  if(isTRUE(class(mod2)=="try-error")) 
  {
    country_code  <- NA 
  } 
  else { 
    if(length(country_label$results$bindings$countryLabel$value) >0){
      country_name <- country_label$results$bindings$countryLabel$value[1]
      country_code <- countrycode(country_name, 'country.name', 'iso2c')
    }
    else
    {
      university_positions <- gregexpr('university',my_string,ignore.case=TRUE)
      comma_positions <- gregexpr(',',my_string,ignore.case=TRUE)
      university_list <-list()
      entity_list <-list()
      country_possibilities <-list()
      
      country_label_0 <- list () 
      
      k<-1
      if ((university_positions[[1]][1]>-1)&(comma_positions[[1]][1]>-1))
      {
        for (n in 1:length(university_positions[[1]]))
        {
          for (m in 1:length(comma_positions[[1]]))
          {
            if (university_positions[[1]][n]<comma_positions[[1]][m])
            {university_list[k] <- substrLeft(my_string, comma_positions[[1]][m]-university_positions[[1]][n])
            k<-k+1
            break
            }
          }
        }
      }
      entity_list <- university_list
      entity_list[k]<-country_name2
      entity_list[k+1]<-country_name
      
      for (l in 1: length(entity_list))
      {
        my_query_0 <- paste('https://query.wikidata.org/sparql?format=json&query=',RCurl::curlEscape(paste('PREFIX
                                                                                                           schema: <http://schema.org/> PREFIX wdt:
                                                                                                           <http://www.wikidata.org/prop/direct/> SELECT ?countryLabel WHERE
                                                                                                           {<https://en.wikipedia.org/wiki/',gsub(' ','_',gsub("\\\\|","",gsub("\\\\\\\\","",entity_list[l]))),'> schema:about
                                                                                                           ?datalink. ?datalink wdt:P17 ?country. SERVICE wikibase:label {
                                                                                                           bd:serviceParam wikibase:language "en" .}}',sep='')),sep='')
        mod2_0 <- try({
          country_label_0 <- fromJSON(url(my_query_0))
        }, TRUE)
        
        if ((isTRUE(class(mod2_0)=="try-error"))|(length(country_label_0$results$bindings$countryLabel$value)==0))
        {
          country_possibilities[l] <- NA   
        }
        else
        {
          country_possibilities[l] <- country_label_0$results$bindings$countryLabel$value[1]
        }
      }
      
      n_possibilities <-length(country_possibilities)
      
      if ((!is.na(country_possibilities[n_possibilities]))&(country_possibilities[n_possibilities] %in% country_possibilities[1:(n_possibilities-1)])){
        country_name <- as.character(country_possibilities[n_possibilities][[1]])
        country_code <- countrycode(country_name, 'country.name', 'iso2c')
      }
      else
      {
        if ((length(country_possibilities)>2)&(!is.na(country_possibilities[n_possibilities-1]))&(country_possibilities[n_possibilities-1] %in% country_possibilities[1:(n_possibilities-2)])){
          country_name <- as.character(country_possibilities[n_possibilities-1][[1]])
          country_code <- countrycode(country_name, 'country.name', 'iso2c')
        }
        else{
          country_code <- NA
        } 
      }
      
  }
}
  return (country_code)
  }


check_formed_datatype_literals <- function(my_time, my_string){ 
  
  time_check <- FALSE
  if (!is.na(my_time))
  {
    a_year <- gsub(' |,','',my_time)
    if (!is.na(as.numeric(a_year)))
    {
      if ((as.numeric(a_year)<max_time)&(as.numeric(a_year)>min_time))
      {
        time_check <- TRUE
      }
    }
  }
  
  spell_check <- FALSE
  word_list <- strsplit(my_string, " ")[[1]]
  word_list_length <- length(word_list)
  if (word_list_length<1)
  {
    print(paste("Error: ",my_string, "has the length of ", word_list_length))
  } else
  {
    for (i in 1:word_list_length)
    {
      spell_check <- spell_check | hunspell_check(word_list[i])
    }
  }
  
  formed_datatype_literals <- time_check & spell_check 
  return(formed_datatype_literals)
}


check_absent_inconsitent_uncondensed_form <- function(my_time, my_string, my_country){ 
  
  data_existed <- 0
  data_standard <- 0
  data_uniquely_used <- 0
  data_compact <- 0
  data_compact_ratio <- 0
  country_code <- ""
  compact_string <- ""
  
  country_code <- countrycode(my_country, 'country.name', 'iso2c')
  if (!is.na(country_code)) 
  {
    data_existed <- 1
    data_uniquely_used <- 1
    original_my_string <- my_country
    my_phrase <- my_country
    if (my_phrase==toupper(countrycode(my_phrase, 'country.name', 'iso.name.en'))) 
    {
      data_standard <- 1
    } 
    compact_string <- countrycode(country_code, 'iso2c','country.name')
    if (nchar(compact_string)==nchar(original_my_string))
    {
      data_compact <- 1 
      data_compact_ratio <- 1
    }
    else
    {
      data_compact <- 0 
      if (nchar(original_my_string)!=0)
      {
        data_compact_ratio <- nchar(compact_string)/nchar(original_my_string)
      }
    }
  }
  
  else
  {
    original_my_string <- my_string
    
    prior_character <- character()
    prior_character<-''
    
    last_char <- substrRight(my_string,1)
    while (last_char %in% special_string_endings)
    {
      my_string <- substrLeft(my_string, nchar(my_string)-1)
      last_char <- substrRight(my_string,1)
    }
    
    position_v1 <- 0
    position_v2 <- 0
    position_v3 <- 0
    position_v4 <- 0
    if (grepl(', ',my_string,ignore.case=TRUE)) 
    {position1 <- gregexpr(', ',my_string,ignore.case=TRUE)
    position_v1 <- position1[[1]][length(position1[[1]])]}
    if (grepl(';',my_string,ignore.case=TRUE)) 
    {position2 <- gregexpr(';',my_string,ignore.case=TRUE)
    position_v2 <- position2[[1]][length(position2[[1]])]}
    if (grepl('\\\\(',my_string,ignore.case=TRUE)) 
    {position3 <- gregexpr('\\\\(',my_string,ignore.case=TRUE)
    position_v3 <- position3[[1]][length(position3[[1]])]}
    if (grepl(',',my_string,ignore.case=TRUE)) 
    {position4 <- gregexpr(',',my_string,ignore.case=TRUE)
    position_v4 <- position4[[1]][length(position4[[1]])]}
    
    position_value <- max(position_v1,position_v2,position_v3,position_v4)
    
    if (position_v1==position_value) 
    {prior_character<-', '} else
    {
      if (position_v2==position_value) {prior_character<-';'}
      if (position_v3==position_value) {prior_character<-'\\\\('}
      if (position_v4==position_value) {prior_character<-','}  
    }
    
    my_phrase <- trimws(substrRight(my_string, nchar(my_string)-position_value))
    
    code_found <- countrycode(my_phrase, 'iso2c', 'iso2c')
    if (!is.na(code_found))
    {
      if(!((code_found %in% state.abb)&(grepl(state.name[match(code_found,state.abb)],my_string,ignore.case=TRUE))))
      {
        country_code <- as.character(my_phrase)
        data_existed <- 1
        data_standard <- 1 
        compact_string <- country_code
      if (isTRUE(prior_character==', ')){name_uniquely_used <- 1}
      }else 
      {
        #do nothing because this is the case of the US states'codes
      }
    } 
    else 
    { 
      country_code <- countrycode(my_phrase, 'country.name', 'iso2c')
      if (!is.na(country_code)) 
      {
        data_existed <- 1
        compact_string <- countrycode(country_code, 'iso2c','country.name')
        if (my_phrase==toupper(countrycode(my_phrase, 'country.name', 'iso.name.en'))) 
       
        {
          data_standard <- 1
         
        } 
        if (isTRUE(prior_character==', ')){data_uniquely_used <- 1} 
      }
    }
    
    if (data_existed==1) 
    {
      if (!grepl(compact_string,original_my_string, fixed=TRUE))
      {
        compact_string <- my_phrase
      }
      if (nchar(compact_string)==nchar(original_my_string))
      {
        data_compact <- 1 
        data_compact_ratio <- 1
      }
      else
      {
        data_compact <- 0 
        if (nchar(original_my_string)!=0)
        {
          data_compact_ratio <- nchar(compact_string)/nchar(original_my_string)
        }
      }
    }
  }

  if (is.na(as.numeric(my_time)))
  {
    data_existed <- 0
    data_standard <- 0
    data_uniquely_used <- 0
    data_compact <- 0
    data_compact_ratio <- 0
  }else
  {
    if ((as.numeric(my_time)>9999) | (as.numeric(my_time)<0))
    {
      data_standard <- 0
      data_uniquely_used <- 0
      data_compact <- 0
      data_compact_ratio <- 0
    }
  }
  
  my_list <- list("data_existed"=data_existed, "data_standard"=data_standard, "data_uniquely_used"=data_uniquely_used, "phrase_found"=compact_string,"country_code"=country_code, "data_compact"=data_compact, "data_compact_ratio"=data_compact_ratio) 
  return(my_list)
  
}


check_column_completeness <- function(my_time, my_string){ 
  
  if ((is.na(my_string))|(is.null(my_string))|(my_string=="Empty value")|(is.na(my_time))|(is.null(my_time)))
  {
    return(FALSE)
  }
  else
  {
    return(TRUE)
  }
}

check_population_completeness <- function(my_time, country_code, Population_presence){ 
  if (!is.na(country_code))
  {
    for (i in 1: country_count)
    {
      if (country_code ==ISO_3166_1[i,1])
      {
        Population_presence$Country_presence[i] <- 1
      }
    }
  }
  if (!is.na(my_time))
  {
    a_year <- gsub(' |,','',my_time)
    if (!is.na(as.numeric(a_year)))
    {
      if ((as.numeric(a_year)<max_time)&(as.numeric(a_year)>min_time))
      {
        for (j in min_time:max_time)
        {
          if (j==as.numeric(a_year))
          {
            Population_presence$Year_presence[j-min_time+1] <- 1
          }
        }
        
      }
    }
  }
  return(Population_presence)
}

check_appropriate_amount_data <- function(my_time, my_string){ 
  
  country_name <- "" 
  country_name2 <- "" 
  
  affiliation_p <- affiliation_data_preparation (my_string)
  country_code <- string_matching(affiliation_p$name,affiliation_p$string)
  
  if (is.na(country_code))
  {
    country_code <- Wikidata_query(affiliation_p$name2, affiliation_p$name, affiliation_p$string)
  }
  
  a_year <- gsub(' |,','',my_time)
  if (!is.na(as.numeric(a_year)))
  {
    if ((as.numeric(a_year)<max_time)&(as.numeric(a_year)>min_time))
    {
      #do nothing
    }else
    {
      a_year <- NA
    }
  }else
  {
    a_year <- NA
  }
  my_list <- list("a_year"=a_year, "country_code"=country_code)
  return(my_list)
}

check_affiliation_vocabularies <- function(my_string){
  
  presence_affiliation_vocabularies <- FALSE
  
  for (i in 1:length(affiliation_vocabularies))
  {
    presence_affiliation_vocabularies <- grepl(affiliation_vocabularies[i], my_string, ignore.case=TRUE)
    if (presence_affiliation_vocabularies==TRUE) {break}
  }
  
  return(presence_affiliation_vocabularies)
}

check_English_spelling <- function(my_string){ 
  
  list_of_unknown <- hunspell(my_string)
  
  if (length(list_of_unknown[[1]])>0)
  {
    for (i in 1: length(list_of_unknown[[1]]))
    {
      if (is.element(list_of_unknown[[1]][i],list_of_cities)) 
      {
        list_of_unknown[[1]] <- list_of_unknown[[1]][list_of_unknown[[1]]!=list_of_unknown[[1]][i]]
      }
    }
  }
  
  if (length(list_of_unknown[[1]])!=0)
  {
    return(FALSE)
  }
  else
  {
    return(TRUE)
  }
}

My_Calculation <- function(MyData, starting_position){
  
  print(Sys.time())
  
  n_row <- nrow(data.frame(MyData)) 
  
  Year_presence <- rep(c(0),times=(max_time-min_time+1))
  Country_presence <- rep(c(0),times=country_count)
  Population_presence <- list("Year_presence"=Year_presence, "Country_presence"=Country_presence)

  Missing_values_counter <- 0
  Absent_data_counter <- 0
  Inconsistent_data_counter <- 0
  Misused_data_counter <- 0
  Inappropriate_amount_data_counter <- 0
  Absent_affiliation_vocabularies_counter <- 0
  Incorrect_English_spelling_counter <- 0
  Malformed_datatype_literals_counter <- 0 
  Uncondensed_Form_counter <- 0 
  Compact_Form_ratio <- 0 
  
  Accuracy_Malformed_Datatype_Free <- 0 
  Concise_Compact_Form <- 0 
  Concise_Compact_Form_ratio <- 0 
  
  Completeness_Property <- 0
  Completeness_Population <- 0
  Appropriate_Data_Explicitly <- 0
  Appropriate_Data_Implicitly <- 0
  Consistency_Standard <- 0
  Consistency_Syntax <- 0
  EoU_Presence_Relevant_Vocabularies <- 0
  EoU_Correct_Spelling <- 0
  
  for (j in 1:n_row)
  {

    my_time <-as.character(MyData[j,1])
    my_string <-as.character(MyData[j,2])
    my_country <- as.character(MyData[j,3])

    if (check_column_completeness(my_time, my_string)==FALSE) 
    {
      Missing_values_counter <- Missing_values_counter +1 
    } 
    else
    {

      data_report <- check_absent_inconsitent_uncondensed_form(my_time, my_string, my_country)
      if (data_report$data_existed==0) 
      {
        Absent_data_counter <- Absent_data_counter +1
        
        AA_data <- check_appropriate_amount_data(my_time, my_string) 
        if ((is.na(AA_data$a_year)) | (is.na(AA_data$country_code))) 
        {            
          Inappropriate_amount_data_counter <- Inappropriate_amount_data_counter +1
        }
        Population_presence <- check_population_completeness(AA_data$a_year, AA_data$country_code, Population_presence) 
        if (check_formed_datatype_literals(my_time, my_string)==FALSE) 
        {
          Malformed_datatype_literals_counter <- Malformed_datatype_literals_counter + 1
        } 
      }
      else 
      {
        Population_presence <- check_population_completeness(my_time, data_report$country_code, Population_presence) 
        
        if ((data_report$data_standard==0))
        {
          Inconsistent_data_counter <- Inconsistent_data_counter +1
        }
        if (data_report$data_uniquely_used==0)
        {
          Misused_data_counter <- Misused_data_counter +1
        }
        if (data_report$data_compact==0)
        {
          Uncondensed_Form_counter <- Uncondensed_Form_counter + 1
        }
        Compact_Form_ratio <- Compact_Form_ratio + data_report$data_compact_ratio
      }  
      
      if (check_affiliation_vocabularies(my_string)==FALSE)
      {
        Absent_affiliation_vocabularies_counter <- Absent_affiliation_vocabularies_counter + 1
      }
      if (check_English_spelling(my_string)==FALSE)
      {
        Incorrect_English_spelling_counter <- Incorrect_English_spelling_counter + 1
      }  
  }
  }
  
  print(Population_presence$Year_presence)
  print(Population_presence$Country_presence)
  
  #Accuracy DQD
  Accuracy_Malformed_Datatype_Free <- (1- Malformed_datatype_literals_counter/(n_row-Missing_values_counter)) 
  
  #Completeness DQD
  Completeness_Property <- (1- Missing_values_counter/n_row)
  Completeness_Population <- sum(Population_presence$Year_presence)/(max_time-min_time+1) * sum(Population_presence$Country_presence)/country_count  
  
  #Appropriate amount of data
  Appropriate_Data_Explicitly <- (1- Absent_data_counter/(n_row-Missing_values_counter))
  Appropriate_Data_Implicitly <- (1- Inappropriate_amount_data_counter/Absent_data_counter)
  
  #Representational Consistency
  Consistency_Standard <- (1- Inconsistent_data_counter/(n_row - Missing_values_counter - Absent_data_counter)) 
  Consistency_Syntax <- (1- Misused_data_counter/(n_row - Missing_values_counter - Absent_data_counter)) 
  
  #Concise Representation
  Concise_Compact_Form <- (1- Uncondensed_Form_counter/(n_row - Missing_values_counter - Absent_data_counter)) 
  Concise_Compact_Form_ratio <- Compact_Form_ratio/(n_row - Missing_values_counter - Absent_data_counter) 
  
  #Ease of Understanding
  EoU_Presence_Relevant_Vocabularies <- (1- Absent_affiliation_vocabularies_counter/(n_row-Missing_values_counter))
  EoU_Correct_Spelling <- (1- Incorrect_English_spelling_counter/(n_row-Missing_values_counter))
  
  write.table(data.frame(
    Missing_values_counter,
    Absent_data_counter,
    sum(Population_presence$Year_presence),
    sum(Population_presence$Country_presence),
    Inappropriate_amount_data_counter,
    Malformed_datatype_literals_counter,
    Absent_affiliation_vocabularies_counter,
    Incorrect_English_spelling_counter,
    Uncondensed_Form_counter,
    Compact_Form_ratio,
    Inconsistent_data_counter,
    Misused_data_counter,
    stringsAsFactors = F),
    file = paste("DQDs - Calculated Metrics - Wos - Counts.csv"), sep=",",append = T,row.names = F, col.names = F)

  write.table(data.frame(
    Completeness_Property, #M1
    Completeness_Population, #M2
    Appropriate_Data_Explicitly, #M3
    Appropriate_Data_Implicitly, #M4
    Accuracy_Malformed_Datatype_Free, #M5
    EoU_Presence_Relevant_Vocabularies, #M6
    EoU_Correct_Spelling, #M7
    Concise_Compact_Form, #M8
    Concise_Compact_Form_ratio,
    Consistency_Standard, #M9
    Consistency_Syntax, #M10
    stringsAsFactors = F),
    file = paste("DQDs - Calculated Metrics - Wos - Metrics.csv"), sep=",",append = T,row.names = F, col.names = F)

}

special_string_endings <- c(',','.',' ','\\\\(','\\\\)',';','-') 
affiliation_vocabularies <- c('university','univ.','school','sch.','college','coll.','deparment','dep.','laboratory','lab.','labs','faculty','academy','acad.','institute','inst.','center','centre','company','comp.','corporation','corp.','group','team','research','science','service','technology','Ltd.','co.','Inc.')

list_of_cities <- maps::world.cities[[1]] 

country_count <- length(ISO_3166_1[,1])
my_block <- 10000
my_sample_iteration <- 40

input_file <- "DQDs - Listed relevant attributes - WoS.csv"
file_length <- length(count.fields(input_file, sep = "\\n"))

WoS <- read.csv(input_file, sep=",")
set.seed(0) 
sample_index <- sample(file_length, my_block*my_sample_iteration, TRUE)

for (i in 1:my_sample_iteration)
{
  print(paste("Calculating metrics of Sample number ", i))
  
  AllData <-NULL
  for (j in 1:my_block)
  {
    AllData <- rbind(AllData, WoS[(i-1)*my_block+j,])  
  }
  
  MyData <- as.data.frame(AllData[,c(2,3,4)]) 
  
  My_Calculation(MyData, 0) 
  
}


min_time <- 1980 
max_time <- 2017 

for (i in 1:my_sample_iteration)
{
  print(paste("Exporting the Year sum of Sample number ", i))
  
  AllData <-NULL
  for (j in 1:my_block)
  {
    AllData <- rbind(AllData, WoS[(i-1)*my_block+j,])  
  }
  
  MyData <- as.data.frame(AllData[,c(2,3,4)]) 
  
  total_count <- 0
  Year_sum <- rep(c(0),times=(max_time-min_time+1))
  
  print(nrow(MyData))
  for (j in 1:nrow(MyData))
  {
    my_time <- MyData[j,1]
    if (!is.na(my_time))
    {
      is_year <- as.numeric(my_time)
      Year_sum[is_year-min_time+1] <- Year_sum[is_year-min_time+1] + 1
      
      if (total_count>my_block){break}
      total_count <- total_count +1
    }else
    {
      print(MyData[j,])
    }
  }
  print(total_count)
  for (k in min_time:max_time)
  {
    write.table(data.frame(k,Year_sum[k-min_time+1], stringsAsFactors = F), file="Year_Sum_WoS.csv", sep=",",append = T,row.names = F, col.names = F)
  }
  
}