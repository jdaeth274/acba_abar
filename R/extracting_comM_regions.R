###############################################################################
## checking out the comM results ##############################################
###############################################################################

require(dplyr)

get_input <- function(){
  args <- commandArgs(trailingOnly = TRUE)
  
  if(length(args) != 3){
    print("Not correct number of inputs, need 3: <comM_hits_csv> <contig_loc> <out_name>")
    quit(save = "no", status = 1, runLast = FALSE)
  }
  
  left_csv <- args[1]
  right_csv <- args[2]
  contig_dir <- args[3]
  
  args <- c(left_csv, right_csv, contig_dir)
  return(args)
}


## lets loop through the assemblies and try to merge them
contig_location_adder <- function(isolate, contig_dir, hit_table, iter){
  ## Function to take an input isolate and add in the contig 
  ## numbers for the different hit locations 
  #if(iter == 71) browser()
  contig_dir <- gsub("/$","",contig_dir)
  contig_file <- paste(contig_dir , "/" , isolate , "#contig_bounds.csv", sep = "")
  
  contig_table <- read.csv(contig_file, stringsAsFactors = FALSE)
  hit_table$contig <- NA
  for(row in 1:nrow(hit_table)){
    current_hit <- hit_table[row,]
    current_pos <- c(current_hit$sstart, current_hit$send)
    contig <- which((contig_table[,1] <= (min(current_pos) + 15)) & ((contig_table[,2] >= (max(current_pos) - 15))))
    if(length(contig) == 1)
      hit_table$contig[row] <- contig
  }
  
  return(hit_table)  
  
}



test_merger <- function(comm_broken, contig_dir){
  
  assembly_num <- dplyr::count(comm_broken, assembly)
  cat("\n","\n")
  combined_hits <- NULL
  tot_isos <- nrow(assembly_num)
  tot_chars <- nchar(tot_isos)
  
  
  for(k in 1:nrow(assembly_num)){
    
    current_assembly <- assembly_num[k,1]
    if(assembly_num[k,2] < 2){
      next
    }
    
    
    
    
    cat_num <- k
    cat_char <- nchar(k)
    nchar_dif <- tot_chars - cat_char
    if(nchar_dif == 0){
      cat("\r",as.character(k), " of ", as.character(tot_isos), " done", sep = "")
    }else{
      cat("\r",paste(rep("0",nchar_dif), collapse = ""),k, " of ", as.character(tot_isos), " done", sep = "")
    }
    
    ## Get the current splits 
    current_comm <- comm_broken[comm_broken$assembly == current_assembly,] %>%
      arrange(qstart) %>% mutate(index = row_number())
    
    current_comm <- contig_location_adder(isolate = current_comm[1,"subject"],
                                          contig_dir = contig_dir,
                                          hit_table = current_comm, 
                                          iter = k)
    
    ## While loop through the fragments till we get a whole one
    align_length <- current_comm[1,"align"]
    current_hit <- 1
    comp_hit <- FALSE
    concat_hits <- current_comm[1,]
    n_hits <- nrow(current_comm)
    next_gap <- 1
    current_gap <- 1
    while(align_length < 1450){
      next_hit <- current_hit + next_gap
      curr_qs <- current_comm[current_hit, "qstart"]
      curr_qe <- current_comm[current_hit, "qend"]
      next_qs <- current_comm[next_hit, "qstart"]
      next_qe <- current_comm[next_hit, "qend"]
      ## give 20 bp overlap potential, really need to know some biochem to back this up
      if((next_qs >= (curr_qe - 20)) & (next_qs <= (curr_qe + 20))){
        align_length <- align_length + current_comm[next_hit, "align"]
        concat_hits <- bind_rows(concat_hits, current_comm[next_hit,])
        next_gap <- 1
        current_gap <- 1
      }else{
        ## try again with the next hit if longer in length 
        if(current_comm[next_hit, "align"] > current_comm[current_hit, "align"]){
          concat_hits <- current_comm[next_hit,]
          align_length <- current_comm[next_hit, "align"]
          next_gap <- 1
          current_gap <- 1
        }else{
          next_gap <- current_comm[next_hit, "index"] - current_comm[current_hit, "index"] + 1
          current_gap <- 0
        }
      }
      current_hit <- current_hit + current_gap
      if((current_hit + next_gap) > n_hits){
        break
      }
    }
    
    if(nrow(concat_hits) > 1){
      ## We do have a potential combiner so lets do some initial qc
      ## first check if all the oris the same, if not we'll lose this one I'm afraid.
      if(length(unique(concat_hits$ori)) != 1){
        ## drop these hits 
        concat_hits <- NULL
      }
      ## Check on same contig 
      if(length(unique(concat_hits$contig)) != 1){
        ## drop these hits 
        concat_hits <- NULL
      }
      
    }else{
      concat_hits <- NULL
    }
    if(align_length < 1450){
      concat_hits <- NULL
    }
    
    if(!is.null(concat_hits)){
      ## Merge into one hit with start and end 
      ## Have id hit start hit end ori contig 
      concat_hits <- concat_hits %>% arrange(qstart)
      out_row <- data.frame(matrix(nrow = 1, ncol = 6))
      colnames(out_row) <- c("id","hit_start","hit_end","ori","contig", "inversion")
      out_row$id <- concat_hits$subject[1]
      out_row$ori <- concat_hits$ori[1]
      out_row$contig <- concat_hits$contig[1]
      if(concat_hits$ori[1] == "forward"){
        out_row$hit_start <- concat_hits$send[1]
        out_row$hit_end <- concat_hits$sstart[nrow(concat_hits)]
        out_row$inversion <- "no"
        if(concat_hits$send[1] > concat_hits$sstart[nrow(concat_hits)]){
          out_row$inversion <- "yes"
          out_row$hit_start <- concat_hits$sstart[nrow(concat_hits)]
          out_row$hit_end <- concat_hits$send[1]
        } 
        
      }else{
        out_row$hit_start <- concat_hits$send[1]
        out_row$hit_end <- concat_hits$sstart[nrow(concat_hits)]
        out_row$inversion <- "no"
        if(concat_hits$send[1] < concat_hits$sstart[nrow(concat_hits)]){
          out_row$inversion <- "yes"
          out_row$hit_start <- concat_hits$send[nrow(concat_hits)]
          out_row$hit_end <- concat_hits$sstart[1]
        } 
        
      }
      ## check if an inversion in the seqs
      
      combined_hits <- bind_rows(combined_hits, out_row)  
    }
    
    
    
    
  }
  
  return(combined_hits)
}


main <- function(){
  
  input_args <- get_input()
  
  comM_hits <- read.csv(input_args[1],
                        stringsAsFactors = FALSE, header = FALSE)
  colnames(comM_hits) <- c("query","subject","pid","align","gap","mismatch",
                           "qstart","qend","sstart","send","eval","bit")
  comM_hits <- comM_hits %>%  mutate(assembly = sub("\\..*$","",subject,perl = TRUE))
  
  ## first remove the complete hits any align >= 1450
  
  comm_broken <- comM_hits %>% filter(align < 1450) %>%
    mutate(ori = ifelse(sstart > send, "reverse","forward"))
  
  combo_hits <- test_merger(comm_broken, input_args[2])
  
  if(!grepl("\\.csv$",input_args[3])){
    out_csv_name <- paste(input_args[3], ".csv",sep = "")
  }else{
    out_csv_name <- input_args[3]
  }
  
  write.csv(combo_hits, out_csv_name, row.names = FALSE, quote = FALSE)
  
}

main()
cat("\n","\n", "Done")
