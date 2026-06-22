
##################################################################################
###############################
# Main taxon detection function
# handles data and auxiliary functions
taxon_detect <- function(mpse = NULL, groupVar = NULL, abdVar = NULL, dbPath = NULL, speciesSearch = NULL){
    if(require(dplyr) != T){
        install.packages(dplyr)
    }
    library(dplyr)
    if(require(tidyr) != T){
        install.packages(tidyr)
    }
    library(tidyr)
    # set dplyr.summarise.inform to avoid annoying messages from dplyr
    dplyr.summarise.opt_orig <- getOption("dplyr.summarise.inform")
    on.exit(options(dplyr.summarise.inform = dplyr.summarise.opt_orig), add = TRUE)
    options(dplyr.summarise.inform = FALSE)
    
    #calculate relative abundance
    rel_abd_tab <- calc_rel_abd(mpse = mpse, groupVar = groupVar, sumAbdVar = {{ abdVar }})
    
    #search db for species of interest
    ## the speciesSearch must contain a column named "species_name"
    speciesSearch.matched <- species_db_match(dbPath = dbPath, speciesSearch = speciesSearch)
    
    #process the rel abd tab for species detections
    detections_tab = species_detections(speciesDbMatch = speciesSearch.matched, relAbdTab = rel_abd_tab, groupVar = groupVar)
    return(list(db_matches = speciesSearch.matched, detections_tab = detections_tab))
    
}
##################################################################################
#taxon_detect_list.Site = taxon_detect(mpse=mpse2, groupVar = c("Site"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)
#taxon_detect_list.Site_Lure = taxon_detect(mpse=mpse2, groupVar = c("Site", "Lure"), abdVar = Abundance, dbPath = dbPath, speciesSearch = spp_search_terms)
##################################################################################


##################################################################################
##################################################################################
calc_rel_abd <- function(mpse = mpse, groupVar = groupVar, sumAbdVar = sumAbdVar){
    
    data.frame(mpse) %>% dplyr::group_by(across(all_of(groupVar))) %>% dplyr::summarize(Abundance = sum( {{ sumAbdVar }}) ) -> groupTots
    taxaLvls = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    if(sum(colnames(data.frame(mpse)) %in% taxaLvls) != length(taxaLvls)){
        missingTaxa = paste(taxaLvls[!taxaLvls %in% colnames(data.frame(mpse))])
        stop("Taxonomic levels are missing: ", missingTaxa) 
    }
    data.frame(mpse) %>% dplyr::group_by(across(all_of(c(taxaLvls, groupVar)))) %>% dplyr::summarize(Abundance = sum( {{ sumAbdVar }}) ) -> sppTots
    #print(groupTots)
    
    sppTots %>% tidyr::pivot_wider(names_from = all_of(groupVar), names_sep = "__nameSep__", values_from = Abundance) -> sppTots.wide
    cbind(sppTots.wide[,1:7], #grabs the taxonomy
        apply(
            sppTots.wide[,8:ncol(sppTots.wide)], 
            1, 
            FUN = function(x) x/groupTots$Abundance*100 #note the vector in the denominator is recycled, so when traversing columns with apply the group totals are recycled for each species. Output gives the original cols in rows so we t() to original format
        ) %>% t
    ) -> sppTots.wide.RA
    if(length(groupVar)>1){
        sppTots.wide.RA %>% pivot_longer(cols = where(is.numeric), names_sep = "__nameSep__", names_to = groupVar, values_to = "relAbd") -> sppTots.long.RA
    }else if(length(groupVar)==1){
        sppTots.wide.RA %>% pivot_longer(cols = where(is.numeric), names_to = groupVar, values_to = "relAbd")  -> sppTots.long.RA
    }
    return(sppTots.long.RA)
}
##################################################################################
#relAbdFoo = calc_rel_abd(mpse = mpse2, groupVar = c("Site"), sumAbdVar = Abundance)
#relAbdFoo = calc_rel_abd(mpse = mpse2, groupVar = c("Site", "Lure"), sumAbdVar = Abundance)
##################################################################################


##################################################################################
##################################################################################
species_db_match <- function(dbPath = dbPath, speciesSearch = speciesSearch){
    #note this function depends on external software including perl and core utilities like grep 
    # the speciesSearch must contain a column names "species_name"
    
    # extract headers from db file (grep) and extract the genus and species info
    genComm <- paste("grep \'>\'",dbPath,"| perl -lne 'print $1 if /g__(.*);/' | sort | uniq | grep -v 'gen_Incertae_sedis' ") 
    dbGenera <- system(genComm, intern = TRUE) 
    #print(dbGenera)
    sppComm <- paste("grep \'>\'",dbPath,"| perl -lne 'print $1 if /s__(.*)/' | grep -v '.*_sp$' | sort | uniq ") 
    dbSpecies <- system(sppComm, intern = TRUE) 
    #print(dbSpecies)

    speciesSearch$speciesInDb <- speciesSearch$species_name %in% dbSpecies
    speciesSearch$genusInDb <- unlist(lapply(strsplit(speciesSearch$species_name, "_"), `[[`, 1)) %in% dbGenera

    return(speciesSearch)
}
##################################################################################
#sppFoo = species_db_match(dbPath = dbPath, speciesSearch = spp_search_terms)
##################################################################################


##################################################################################
species_detections <- function(speciesDbMatch = speciesDbMatch, relAbdTab = relAbdTab, groupVar = groupVar){

    # set up the search cols
    # lapply(strsplit(sub("s__", "",foo), "_g__"), '[', 2)
    #relAbdTab$Genus = sub("g__", "", relAbdTab$Genus)
    relAbdTab$gen_sp <- ifelse(
        # MPSE does subs of NA values to it's own text, so if we have a spp NA we wil get "s__un_" followed by the genus preceded with "g__"
        # Note that in this case the family level is also pasted onto the genus so we should not rely on that unless we really explore the format
        # (and the sp level actually gets the full genus+family pasted on so we can get something like "s__un_g__Geosmithia_g__Geosmithia_f__Dothioraceae"
        # foo = "s__un_g__Geosmithia_g__Geosmithia_f__Dothioraceae"
        # foo = "s__un_g__Geosmithia"
        # lapply(strsplit(sub("s__", "",foo), "_g__"), '[', 2)
        # paste0(lapply(strsplit(sub("s__", "",foo), "_g__"), '[', 2), "_NA")
        grepl("s__un_", relAbdTab$Species) == T, 
        paste0(lapply(strsplit(sub("s__", "",relAbdTab$Species), "_g__"), '[', 2), "_NA"), 
        paste(sub("g__", "", relAbdTab$Genus), sub("s__", "", relAbdTab$Species), sep = "_")
    )
    
    # to hold search results
    sppDetectList = list()
    
    #we set up a new df to populate for spp that have no detects
    relAbdTab$gen_sp %>% unique -> detected_sp
    noDetectsTable = relAbdTab %>% ungroup %>% filter(gen_sp == detected_sp[1]) %>% select(all_of(c(groupVar, "relAbd")))
    noDetectsTable$relAbd = "not detected"
    noDetectsTableGen = noDetectsTable
    noDetectsTableGen$relAbd = "sp. not in db; no unkn. spp. match gen."
    
    for(i in 1:nrow(speciesDbMatch)){
        searchGen = unlist(strsplit(speciesDbMatch$species_name[i], "_"))
        if(speciesDbMatch$species_name[i] %in% relAbdTab$gen_sp){
            #print("foo")
        #}
            relAbdTab %>% 
                filter(gen_sp == speciesDbMatch$species_name[i]) %>%
                ungroup %>%
                select(all_of(c(groupVar, "relAbd"))) -> tempTab
            
            for(j in 1:nrow(tempTab)){
                if(tempTab$relAbd[j] == 0){tempTab$relAbd[j] = "not detected"}
                else if(as.numeric(tempTab$relAbd[j]) > 0 & as.numeric(tempTab$relAbd[j]) < 0.01){tempTab$relAbd[j] = "present <0.01%"}
                else if(as.numeric(tempTab$relAbd[j]) >= 0.01 & as.numeric(tempTab$relAbd[j]) < 0.1){tempTab$relAbd[j] = "present 0.01-0.1%"}
                else if(as.numeric(tempTab$relAbd[j]) >= 0.1 & as.numeric(tempTab$relAbd[j]) < 1){tempTab$relAbd[j] = "present 0.1-1%"}
                else if(as.numeric(tempTab$relAbd[j]) >= 1 & as.numeric(tempTab$relAbd[j]) < 10){tempTab$relAbd[j] = "present 1-10%"}
                else if(as.numeric(tempTab$relAbd[j]) >= 10 & as.numeric(tempTab$relAbd[j]) <= 100){tempTab$relAbd[j] = "present >10%"}
            }
            
            sppDetectList[[ speciesDbMatch$species_name[i] ]] <-  tempTab
            
        }else if(speciesDbMatch$speciesInDb[i] == F & paste0(searchGen[1], "_NA") %in% relAbdTab$gen_sp){
            relAbdTab %>% 
                filter(gen_sp == paste0(searchGen[1], "_NA")) %>%
                group_by(across(all_of(c("gen_sp", groupVar)))) %>%
                summarize(relAbd = sum(relAbd)) %>%
                ungroup %>% select(-gen_sp) -> sppDetectList[[ speciesDbMatch$species_name[i] ]]
                
                sppDetectList[[ speciesDbMatch$species_name[i] ]]$relAbd[sppDetectList[[ speciesDbMatch$species_name[i] ]]$relAbd > 0] = "unknown spp. in gen. detected; sp. not in db"
                sppDetectList[[ speciesDbMatch$species_name[i] ]]$relAbd[sppDetectList[[ speciesDbMatch$species_name[i] ]]$relAbd == 0] = "sp. not in db; no unkn. spp. match gen."
        }else if(speciesDbMatch$speciesInDb[i] == F & !paste0(searchGen[1], "_NA") %in% relAbdTab$gen_sp){
                sppDetectList[[ speciesDbMatch$species_name[i] ]] <- noDetectsTableGen
        }else if(speciesDbMatch$speciesInDb[i] == T){
                sppDetectList[[ speciesDbMatch$species_name[i] ]] <- noDetectsTable
        }
    }
    
    full_tab = bind_rows(sppDetectList, .id = "Species")
    colnames(full_tab)[ncol(full_tab)] = "Detections"
    return(full_tab)
}
##################################################################################
#species_detections(speciesDbMatch = sppFoo, relAbdTab = relAbdFoo, groupVar = c("Site"))
##################################################################################
