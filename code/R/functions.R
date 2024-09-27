################# Data preprocessing
# Describe_directory----
# Purpose: to describe files with sizes in a folder/directory
# Note: Will sort in ascending order according to size 
describe_directory <- function(dir) {
  x = tibble(filenames = list.files(dir, full.names = F)) %>% 
    mutate(size_Mb = file.size(list.files(dir, full.names = T)) / 1024^2) %>% 
    arrange(size_Mb)
  return(x)
}

# Variables in directory----
# Purpose: to list variables in all files in a directory with duplicates for variables in several files (i.e., keys)
# Note: Throws a warning due to silly workadound (setwd) - this is ok

variables_in_directory <- function(dir, delim = "|") {
  files = describe_directory(dir)
  setwd(dir) 
  var_names = data.frame()
  for (i in 1:nrow(files)) {
    temp_data = vroom(files$filenames[i], delim = delim, .name_repair = "minimal",
                      locale = locale(encoding = "latin1"), n_max = 1000, guess_max = Inf, show_col_types = FALSE)
    var_names = rbind(
      var_names, 
      tibble(name = colnames(temp_data)) %>% 
        mutate(source = files$filenames[i]))
  }
  return(var_names)
}


# Convert_to_parquet----
# Purpose: Remove redundant variables and convert all files in a directory to parquet
# Note: Time consuming if large dataset (several hours)

convert_to_parquet <- function(file_list, var_list) {
  for (i in 1:nrow(file_list)) {
    # determining keep variables for each file
    temp_var_list <- var_list %>% filter(source == file_list$filenames[i] & keep)
    # reading only those variables
    temp_file <- vroom(file_list$filenames[i], delim = "|", col_select = c(temp_var_list$name),
                       .name_repair = "unique", guess_max = Inf, locale = locale(encoding = "latin1"), show_col_types = FALSE)
    # This is needed due to a misspecification in some files
    temp_file <- temp_file %>%
      filter(if_any(matches("Alias"), ~(!is.na(.x) & nchar(.x) %in% c(7,8)))) 
    
    # remove the .csv suffix and save as .parquet
    write_parquet(temp_file, sink = sub('\\.csv$', '\\.parquet', paste("data/converted_parquet/", 
                                                                       raw_files$filenames[i], 
                                                                       sep = "")))
  }
}





################## Microbiology data cleaning
# Clean LIMS file ----------------------
# purpose: to clean a file originating from LIMS system, with the purpose to merge with a file from wwLab
# Note: Heavily dependent on local variable names and formatting, everything needs to be double checked if new data

clean_LIMS_file <- function(raw_LIMS_file){
  cleaned_LIMS_file <- raw_LIMS_file %>% 
    janitor::clean_names() %>%
    dplyr::mutate(
      finding = ifelse(
        fynd == "Ingen växt", "Negative", ifelse(mikroorganism == "Terminerad",  "Negative", mikroorganism)),
      ttp = round(as.numeric(period_to_seconds(period(str_replace_all(omslagstid, c("m" = "M", "h" = "H")))) / 3600),1), # workaround to get TTP in hours
      bottle = paste("bottle", substr(flaska, 19, 20), sep = ""), # this retrieves the bottle number
      sex = ifelse(kon == "MAN", "Male", "Female"), 
      incubation_start = as_datetime(NA), # this can be found in wwLab file, adds NA here which will be filled later
      sample_date = as.Date(provtagningsdatum, format = "%Y-%m-%d"),
      sample_time = as_hms(format(provtagningsdatum, format = "%H:%M:%S")),
      sample_datetime = provtagningsdatum,
      .keep = "unused") %>% 
    dplyr::select(
      -c(ar, 
         ankomstdatum, 
         svarsdatum, 
         undersokning, 
         material, 
         kund,
         flasknummer, 
         test_andrad_av, 
         test_resultatdatum, 
         first_prel_forsta_preliminarsvar), 
      -where(is.logical)) %>% # need to remove these first as the ab columns will be sorted then 
    dplyr::select(
      sample_id = mikrobiologi_prov_alias, 
      patient_id = rs_pat_alias, 
      labnr, 
      sample_date, 
      sample_time,
      sample_datetime,
      incubation_start, 
      age = alder, 
      sex, 
      bottle, 
      ttp, 
      finding, 
      sample_site = kundkod, 
      hospital = ort,
      sort(colnames(.))
    )
  return(cleaned_LIMS_file)
}

# Clean wwLab file -------
# purpose: to clean wwLab file and prepare for merging
# Note: wwLab file has one row per set of bottles, why pivoting is needed, this is sensitive to error
# In addition, highly dependent on local circumstances

clean_wwLab_file <- function(raw_wwLab_file){
  cleaned_wwLab_file <- raw_wwLab_file %>% 
    janitor::clean_names() %>% 
    dplyr::rename(
      bottle1_ttp = ttd, 
      bottle2_ttp = ttd_1) %>%  # these need to be renamed before pivoting
    dplyr::mutate(
      finding = ifelse(fynd == "Ingen växt", "Negative", fynd),
      labnr = paste(str_sub(ar, 3,4), avd, avdnr, sep = ""), # making one labnr as this was several variables in raw data
      bottle1_find = case_when(
        str_detect(result, "Negativ") ~"Negative", # harmonise with LIMS
        str_detect(finding, "Negative") ~ "Negative",
        TRUE ~ finding),
      bottle2_find = case_when(
        str_detect(result_1, "Negativ") ~"Negative", # harmonise with LIMS
        str_detect(finding, "Negative") ~ "Negative",
        TRUE ~ finding),
      sample_date = as.Date(ifelse(!is.na(provdatum), provdatum, ankomstdatum)),
      sample_time = provtid,
      sample_datetime = as.POSIXct(paste(sample_date, sample_time, sep = " "), format = "%Y-%m-%d %H:%M:%S"),
      sex = ifelse(kon == "Man", "Male", "Female"),
      .keep = "unused") %>%
    tidyr::pivot_longer(cols = c(bottle1_ttp, bottle2_ttp, bottle1_find, bottle2_find), # detta är viktigt eftersom LIMS har en rad per flaska men wwLab har en per set
                        names_to = c("bottle", ".value"), 
                        names_sep = "_") %>% 
    dplyr::mutate(finding = if_else(find == "Negative", "Negative", finding)) %>% # samma variabler som ovan
    dplyr::select(
      sample_id = mikrobiologi_prov_alias, 
      patient_id = rs_pat_alias, 
      labnr, 
      sample_date, 
      sample_time,
      sample_datetime,
      incubation_start = start_dtm_1, 
      age = alder, 
      sex, 
      bottle, 
      ttp, 
      finding, 
      sample_site = mgkod, 
      hospital = city, 
      contains(c("sir"))
    )
  return(cleaned_wwLab_file)
}

# Fix SIR data -----
# purpose: To fix antibiotic susceptibility (SIR) data, to work for both wwLab or LIMS
# col_num_first_ab is hardcoded and must be double checked
# suffices must be double checked using names (df)
# Ablist 2 is adaptation to different naming schemes in LIMS and wwLab - otherwise redundant

fix_sir_data <- function(df, first_ab_col = col_num_first_ab){ 
  abcol_numbers = first_ab_col:ncol(df) # this is the colnumbers of antibiotics cols
  abcols <- df %>% select(all_of(abcol_numbers))
  ablist <- unique(str_remove_all(colnames(abcols), "_mic|_buljong|_2|sir_|mic_sir_|13210|13010|13510|13220|30ec|_b|_res_monster")) #These should be removed from names - double check if needed
  for (i in 1:length(ablist)) {
    temp_ab = ablist[i]
    temp_df <- data.frame(df %>% select(matches(temp_ab)))
    temp_df[,temp_ab] <- temp_df[cbind(1:nrow(temp_df), max.col(!is.na(temp_df), ties.method = "last"))] # ties last to prioritize MIC (which is last alphabetically) change to first if zones are better
    df[,temp_ab] <- temp_df[,temp_ab]
  }  
  
  # after loop, renaming old abs to same - for consistency across files
  ablist2 <- case_when(str_detect(ablist, "klavulans") ~ "amoxicillin_klavulansyra",
                       str_detect(ablist, "tericin") ~ "amphotericin_b",
                       str_detect(ablist, "ceftazidimavibak") ~ "ceftazidim_avibaktam",
                       str_detect(ablist, "ceftolozan") ~ "ceftolozan_tazobaktam",
                       str_detect(ablist, "piperacillin") ~ "piperacillin_tazobaktam",
                       str_detect(ablist, "trimetoprimsulfa") ~ "trimetoprim_sulfametoxazol",
                       str_detect(ablist, "penicillinv") ~ "penicillin_v",
                       TRUE ~ ablist)
  
  df <- df %>% 
    select(-colnames(abcols), 1:first_ab_col - 1, all_of(ablist)) %>% 
    setnames(old = ablist, new = ablist2) %>% 
    select(1:first_ab_col - 1, sort(colnames(.)))
  
  # then we fill antibiotics as there is only susceptibility data on one bottle / row
  df <- df %>% 
    group_by(patient_id, sample_date, finding) %>% 
    fill(c(col_num_first_ab:ncol(.)), .direction = "downup")
  
  return(df)
}


# Removing sir data  ----
# Purpose: to remove SIR data from merged file
# Note: need to check which variables.

remove_sir_data <- function(df, first_ab_col = col_num_first_ab) {
  sir_data <- df %>% 
    filter(finding != "Negative") %>% 
    select(patient_id, sample_date, finding, first_ab_col:ncol(df)) %>% 
    pivot_longer(cols = 4:ncol(.), names_to = "antibiotic", values_to ="sir")  %>% 
    filter(!is.na(sir))
  
  df <- df %>% select(-(first_ab_col:ncol(df)))
  write_parquet(sir_data, here("data", "processed_data", "sir_data.parquet"))
  return(df)
}
# Clean merged files -----
# purpose: to further clean merged file

clean_merged_file <- function(df) {
  df <- df %>% 
    mutate(hospital = case_when(hospital == "MALMÖ" ~ "MALMÖ",
                                hospital == "LUND" ~ "LUND",
                                hospital == "HELSINGBORG" ~ "HELSINGBORG",
                                hospital == "KRISTIANSTAD" ~ "KRISTIANSTAD",
                                hospital == "YSTAD" ~ "YSTAD",
                                hospital == "BROBY" ~ "KRISTIANSTAD",
                                hospital == "BUNKEFLOSTRAND" ~ "MALMÖ",
                                hospital == "HÖGANÄS" ~ "HELSINGBORG",
                                hospital == "HÖLLVIKEN" ~ "MALMÖ",
                                hospital == "KLIPPAN" ~ "HELSINGBORG",
                                hospital == "LIMHAMN" ~ "MALMÖ",
                                hospital == "LOMMA" ~ "MALMÖ",
                                hospital == "STAFFANSTORP" ~ "LUND",
                                hospital == "KÄVLINGE" ~ "LUND",
                                hospital == "SVALÖV" ~ "LUND",
                                hospital == "ÅHUS" ~ "KRISTIANSTAD",
                                hospital == "EKEBY" ~ "HELSINGBORG",
                                TRUE ~ hospital),
           has_cabinet = hospital %in% c("HELSINGBORG", "KRISTIANSTAD", "LUND", "MALMÖ", "YSTAD")) %>% 
    relocate(has_cabinet, .after = hospital) %>% 
    filter(sample_date %within% interval(start = start_date, end = end_date)) %>% # some samples are just outside the period
    arrange(patient_id) # for simplicity
}


# Fix TTP ----------------
# Purpose: to set ttp to NA if finding negative, to stratify TTP
fix_ttp <- function(df, ttp_cutoff = chosen_ttp_cutoff){
  df <- df %>% 
    mutate(
      ttp = ifelse(is.na(ttp) | finding == "Negative", NA, ttp),
      ttp_strata = case_when(ttp <= ttp_cutoff ~"<= 10 hours",
                             ttp > ttp_cutoff ~ "> 10 hours",
                             is.na(ttp) & finding != "Negative" ~ NA,
                             finding == "Negative" ~ "Negative"))
}
# Rename bacteria -----
# Purpose: merge names from LIMS and wwLab, fix language issues
rename_bacteria <- function(df) {
  ## renaming and translating species - this takes height for all misspellings etc we have encountered in Bactius, wwLab, and LIMS
  df <- df %>%
    mutate(finding = tolower(finding), # this is just for harmonising, in the older systems, names were fully CAPITALIZED
           species = case_when(str_detect(finding, "actinomycetemcomitans") ~ "aggregatibacter actinomycetemcomitans",
                               str_detect(finding, "acinetobacter calcoaceticus-komplexet") ~ "acinetobacter calcoaceticus",
                               str_detect(finding, "dubliniensis") ~ "candida dubliniensis",
                               str_detect(finding, "acinetobacter baumannii") ~ "acinetobacter baumannii",
                               str_detect(finding, "acinetobacter-art") ~ "acinetobacter species",
                               str_detect(finding, "aeromonas hydrophila") ~ "aeromonas hydrophila",
                               str_detect(finding, "aeromonas-art") ~ "aeromonas species",
                               str_detect(finding, "anaerob blandflora") ~ "anaerobic mixed growth",
                               str_detect(finding, "anaeroba gramnegativa kocker") ~ "anaerobic gramnegative cocci",
                               str_detect(finding, "anaeroba bakterier") ~ "anaerobic bacteria",
                               str_detect(finding, "anaerob gram-negativ och gram-positiv blandflora") ~ "anaerobic gramnegative and grampositive mixed flora",
                               str_detect(finding, "anaeroba gram-negativa kocker") ~ "anaerobic gramnegative cocci",
                               str_detect(finding, "anaerob gram negativ blandflora") ~ "anaerobic gramnegative mixed growth",
                               str_detect(finding, "anaeroba gramnegativa stavar") ~ "anaerobic gramnegative rods",
                               str_detect(finding, "anaeroba gram-negativa stavar") ~ "anaerobic gramnegative rods",
                               str_detect(finding, "anaeroba grampositiva kocker") ~ "anaerobic grampositive cocci",
                               str_detect(finding, "anaeroba gram-positiva kocker") ~ "anaerobic grampositive cocci",
                               str_detect(finding, "anaerob gram positiv blandflora") ~ "anaerobic grampositive mixed growth",
                               str_detect(finding, "anaeroba grampositiva stavar") ~ "anaerobic grampositive rods",
                               str_detect(finding, "anaeroba gram-positiva stavar") ~ "anaerobic grampositive rods",
                               str_detect(finding, "anaerob växt") ~ "anaerobic growth",
                               str_detect(finding, "blandflora av anaeroba bakt.") ~ "anaerobic mixed growth",
                               str_detect(finding, "anaerob blandflora") ~ "anaerobic mixed growth",
                               str_detect(finding, "anaeroba streptokocker") ~ "anaerobic streptococci",
                               str_detect(finding, "bacillus cereus-komplexet") ~ "bacillus cereus",
                               str_detect(finding, "bacillus-art") ~ "bacillus species",
                               str_detect(finding, "bacillus") ~ "bacillus species",
                               str_detect(finding, "bacteroides tillh. b.fragilis gruppen") ~ "bacteroides fragilis",
                               str_detect(finding, "bacteroides-art") ~ "bacteroides species",
                               str_detect(finding, "brevibacterium-art") ~ "brevibacterium species",
                               str_detect(finding, "campylobacter jejuni") ~ "campylobacter jejuni",
                               str_detect(finding, "chryseobacterium-art") ~ "chryseobacterium species",
                               str_detect(finding, "citrobacter-art") ~ "citrobacter species",
                               str_detect(finding, "clostridium-art") ~ "clostridium species",
                               str_detect(finding, "corynebacterium-art") ~ "corynebacterium species",
                               str_detect(finding, "delftia-art") ~ "delftia species",
                               str_detect(finding, "desulfovibrio-art") ~ "desulfovibrio species",
                               str_detect(finding, "eggerthella-art") ~ "aegerthella species",
                               str_detect(finding, "enterobacter cloacae-komplexet") ~ "enterobacter cloacae",
                               str_detect(finding, "enterobacter-art") ~ "enterobacter species",
                               str_detect(finding, "enterococcus gallinarum-gruppen") ~ "enterococcus gallinarum",
                               str_detect(finding, "enterokock-art") ~ "enterococcus species",
                               str_detect(finding, "e.coli") ~ "escherichia coli",
                               str_detect(finding, "e.coli") ~ "escherichia coli",
                               str_detect(finding, "eubacterium-art") ~ "eubacterium species",
                               str_detect(finding, "exiguobacterium-art") ~ "exiguobacterium species",
                               str_detect(finding, "finegoldia magna") ~ "finegoldia magna",
                               str_detect(finding, "fusobacterium-art") ~ "fusobacterium species",
                               str_detect(finding, "gramnegativa kocker") ~ "gramnegative cocci",
                               str_detect(finding, "gram negativa kocker, ej artbestämda") ~ "gramnegative cocci",
                               str_detect(finding, "gramnegativ miljöbakterie") ~ "gramnegative environmental bacteria",
                               str_detect(finding, "gram negativ miljöbakterie") ~ "gramnegative environmental bacteria",
                               str_detect(finding, "gram negativ blandflora, bl.a. proteus") ~ "gramnegative mixed flora",
                               str_detect(finding, "gramnegativ blandflora") ~ "gramnegative mixed flora",
                               str_detect(finding, "gramnegativa stavar") ~ "gramnegative rods",
                               str_detect(finding, "gram-negativa stavar") ~ "gramnegative rods",
                               str_detect(finding, "gram neg stavar, ej artbestämda") ~ "gramnegative rods",
                               str_detect(finding, "gramnegativ tarmbakterie") ~ "gramnegative rods (enterobacterales)",
                               str_detect(finding, "gram-negativ stav tillhö. enterobacteriaceae") ~ "gramnegative rods (enterobacterales)",
                               str_detect(finding, "grampositiva kocker") ~ "grampositive cocci",
                               str_detect(finding, "grampositiva kocker, troligen strepto- eller enterokocker") ~ "grampositive cocci",
                               str_detect(finding, "gram-positiva kocker") ~ "grampositive cocci",
                               str_detect(finding, "gram positiv kock,ej artbestämd") ~ "grampositive cocci",
                               str_detect(finding, "grampositiva kocker, troligen stafylokocker") ~ "grampositive cocci",
                               str_detect(finding, "grampositiv blandflora") ~ "grampositive mixed flora",
                               str_detect(finding, "grampositiva stavar") ~ "grampositive rods",
                               str_detect(finding, "gram-positiva stavar") ~ "grampositive rods",
                               str_detect(finding, "gramlabila stavar") ~ "grampositive rods",
                               str_detect(finding, "haemophilus influenzae serotyp/biotyp") ~ "haemophilus influenzae",
                               str_detect(finding, "enterobacter aerogenes") ~ "klebsiella aerogenes  ",
                               str_detect(finding, "klebsiella oxytoca-gruppen") ~ "klebsiella oxytoca",
                               str_detect(finding, "klebsiella pneumoniae") ~ "klebsiella pneumoniae",
                               str_detect(finding, "klebsiella-art") ~ "klebsiella species",
                               str_detect(finding, "lactobacillus-art") ~ "lactobacillus species",
                               str_detect(finding, "leptotrichia-art") ~ "leptotrichia species  ",
                               str_detect(finding, "leuconostoc-art") ~ "leuconostoc species",
                               str_detect(finding, "listeria monocytogenes serotyp") ~ "listeria monocytogenes",
                               str_detect(finding, "mikroaerofila gram positiva kocker") ~ "microaerophilic grampositive cocci",
                               str_detect(finding, "mikroaerofila gram-positiva stavar") ~ "microaerophilic grampositive rods",
                               str_detect(finding, "mikrokock-art") ~ "micrococcus species",
                               str_detect(finding, "mikroaerofila streptokocker") ~ "microaerophilic streptococci",
                               str_detect(finding, "blandflora") ~ "mixed flora",
                               str_detect(finding, "blandflora") ~ "mixed flora",
                               str_detect(finding, "blandflora av tarmbakterier") ~ "mixed growth of enterobacterales",
                               str_detect(finding, "moraxella-art") ~ "moraxella species",
                               str_detect(finding, "mycobacterium fortuitum-komplexet") ~ "mycobacterium fortuitum",
                               str_detect(finding, "myroides species f.d flavobacterium") ~ "myroides species",
                               str_detect(finding, "neisseria meningitidis") ~ "neisseria meningitidis",
                               str_detect(finding, "neisseria-art") ~ "neisseria species",
                               str_detect(finding, "troligen neisseria") ~ "neisseria species",
                               str_detect(finding, "paenibacillus-art") ~ "pantoea species",
                               str_detect(finding, "pantoea-art") ~ "pantoea species",
                               str_detect(finding, "peptostreptococcus-art") ~ "peptostreptococcus species",
                               str_detect(finding, "porphyromonas-art") ~ "porphyromonas species",
                               str_detect(finding, "bivia") ~ "prevotella bivia",
                               str_detect(finding, "buccalis") ~ "prevotella buccalis",
                               str_detect(finding, "melaninogenica") ~ "prevotella melaninogenica",
                               str_detect(finding, "prevotella-art") ~ "prevotella species",
                               str_detect(finding, "propionibacterium-art") ~ "proprionebacterium species",
                               str_detect(finding, "proteus mirabilis") ~ "proteus mirabilis",
                               str_detect(finding, "proteus vulgaris-gruppen") ~ "proteus vulgaris",
                               str_detect(finding, "providencia-art") ~ "providencia species",
                               str_detect(finding, "pseudomonas-art") ~ "pseudomonas species",
                               str_detect(finding, "psychrobacter-art") ~ "psychrobacter species",
                               str_detect(finding, "ralstonia-art") ~ "ralstonia species",
                               str_detect(finding, "ornithinolytica") ~ "raoultella ornithinolytica",
                               str_detect(finding, "planticola") ~ "raoultella planticola",
                               str_detect(finding, "rothia-art") ~ "rothia species  ",
                               str_detect(finding, "salmonella-art") ~ "salmonella species",
                               str_detect(finding, "selenomonas-art") ~ "selenomonas species",
                               str_detect(finding, "serratia-art") ~ "serratia species",
                               str_detect(finding, "hudflora") ~ "skin flora",
                               str_detect(finding, "meticillinresistent s.aureus") ~ "staphylococcus aureus",
                               str_detect(finding, "troligen staphylococcus aureus") ~ "staphylococcus aureus",
                               str_detect(finding, "staphylococcus species, ej s.aureus") ~ "staphylococcus species",
                               str_detect(finding, "stafylokocker ?, ej s.aureus") ~ "staphylococcus species",
                               str_detect(finding, "koagulas-negativa stafylokocker") ~ "staphylococcus species",
                               str_detect(finding, "staphylococcus species") ~ "staphylococcus species",
                               str_detect(finding, "streptokocker grupp b") ~ "streptococcus agalactiae",
                               str_detect(finding, "streptokocker grupp b serotyp ii") ~ "streptococcus agalactiae",
                               str_detect(finding, "streptokocker grupp b serotyp iii") ~ "streptococcus agalactiae",
                               str_detect(finding, "streptococcus agalactiae") ~ "streptococcus agalactiae",
                               str_detect(finding, "streptococcus anginosus-komplexet") ~ "streptococcus anginosus",
                               str_detect(finding, "streptococcus anginosus-gruppen") ~ "streptococcus anginosus",
                               str_detect(finding, "streptococcus bovis-komplexet") ~ "streptococcus bovis",
                               str_detect(finding, "streptococcus bovis-gruppen") ~ "streptococcus bovis",
                               str_detect(finding, "streptococcus species tillhörande s.bovis gruppen.") ~ "streptococcus bovis",
                               str_detect(finding, "betastreptokocker grupp g") ~ "streptococcus dysgalactiae",
                               str_detect(finding, "streptococcus dysgalactiae") ~ "streptococcus dysgalactiae",
                               str_detect(finding, "betastreptokocker grupp c") ~ "streptococcus dysgalactiae",
                               str_detect(finding, "streptococcus dysgalactiae") ~ "streptococcus dysgalactiae",
                               str_detect(finding, "beta-streptokocker grupp g") ~ "streptococcus dysgalactiae",
                               str_detect(finding, "beta-streptokocker grupp c") ~ "streptococcus dysgalactiae",
                               str_detect(finding, "beta-streptokocker") ~ "streptococcus dysgalactiae",
                               str_detect(finding, "streptococcus mitis-gruppen") ~ "streptococcus mitis",
                               str_detect(finding, "streptococcus mitis-komplexet") ~ "streptococcus mitis",
                               str_detect(finding, "streptococcus mitis-oralis") ~ "streptococcus mitis",
                               str_detect(finding, "streptococcus pneumoniae") ~ "streptococcus pneumoniae",
                               str_detect(finding, "streptococcus pneumoniae typ") ~ "streptococcus pneumoniae",
                               str_detect(finding, "pneumokocker") ~ "streptococcus pneumoniae",
                               str_detect(finding, "gram-positiva kocker, prel pneumokocker") ~ "streptococcus pneumoniae",
                               str_detect(finding, "beta-streptokocker grupp a") ~ "streptococcus pyogenes",
                               str_detect(finding, "beta-streptokocker grupp a") ~ "streptococcus pyogenes",
                               str_detect(finding, "streptococcus pyogenes") ~ "streptococcus pyogenes",
                               str_detect(finding, "streptococcus salivarius-komplexet") ~ "streptococcus salivarius",
                               str_detect(finding, "streptococcus salivarius-gruppen") ~ "streptococcus salivarius",
                               str_detect(finding, "streptococcus sanguinis-komplexet") ~ "streptococcus sanguinis",
                               str_detect(finding, "alfa-streptokocker") ~ "streptococcus species",
                               str_detect(finding, "streptokock-art") ~ "streptococcus species",
                               str_detect(finding, "streptokocker av alfa-streptokock typ") ~ "streptococcus species",
                               str_detect(finding, "veillonella-art") ~ "veillonella species",
                               str_detect(finding, "jästsvamp") ~ "yeast",
                               str_detect(finding, "jästsvamp") ~ "yeast",
                               str_detect(finding, "jästsvamp, ej c.albicans") ~ "yeast",
                               str_detect(finding, "jästsvamp?") ~ "yeast",
                               str_detect(finding, "streptococcus mutans-komplexet") ~ "streptococcus mutans",
                               str_detect(finding, "maltophilia") ~ "stenotrophomonas maltophilia",
                               str_detect(finding, "sporosarcina-art") ~ "sporosarcina species",
                               str_detect(finding, "serratia liquefaciens-komplexet") ~ "serratia liquefaciens",
                               str_detect(finding, "mucilaginosus") ~ "rothia mucilaginosa",
                               str_detect(finding, "radiobacter") ~ "rhizobium radiobacter",
                               str_detect(finding, "asaccharolytica") ~ "porphyromonas asaccharolytica",
                               str_detect(finding, "agglomerans") ~ "pantoea agglomerans",
                               str_detect(finding, "morganii") ~ "morganella morganii",
                               str_detect(finding, "laktobaciller ,") ~ "lactobacillus species",
                               str_detect(finding, "laktobaciller") ~ "lactobacillus species",
                               str_detect(finding, "fusarium solani.") ~ "fusarium solani",
                               str_detect(finding, "pneumosintes") ~ "dialister pneumosintes",
                               str_detect(finding, "vesicularis") ~ "brevundimonas vesicularis",
                               str_detect(finding, "blandflora av anaeroba gram.neg. stavar") ~ "anaerobic mixed growth",
                               str_detect(finding, "acinetobacter johnsonii") ~ "acinetobacter johnsonii",
                               str_detect(finding, "alfastreptokocker") ~ "streptococcus species",
                               str_detect(finding, "bifidobacterium-art") ~ "bifidobacterium species",
                               str_detect(finding, "brachybacterium-art") ~ "brachybacterium species  ",
                               str_detect(finding, "burkholderia cepacia-komplexet") ~ "burkholderia cepacia",
                               str_detect(finding, "comamonas-art") ~ "comamonas species",
                               str_detect(finding, "difteroida stavar") ~ "difteroid rods",
                               str_detect(finding, "trådsvamp") ~ "fungus",
                               str_detect(finding, "munflora") ~ "oral flora",
                               str_detect(finding, "stafylokocker") ~ "staphylococcus species",
                               str_detect(finding, "strept. gallolyticus ssp gallolyticus") ~ "streptococcus gallolyticus",
                               str_detect(finding, "streptococcus equi subsp zooepidemicus") ~ "streptococcus equi",
                               str_detect(finding, "svamp species") ~ "fungus",
                               str_detect(finding, "troligen pneumokocker") ~ "streptococcus pneumoniae",
                               str_detect(finding, "wickerhamomyces anomalus") ~ "candida pelliculosa",
                               str_detect(finding, "clavispora lusitaniae") ~ "candida lusitaniae",
                               str_detect(finding, "kluyveromyces marxianus") ~ "candida kefyr",
                               str_detect(finding, "ingen växt") ~ "negative",
                               TRUE ~ finding),
           # this creates a classification whith is more clinically relevant but still quite granular
           genus = case_when(str_detect(species, "escherichia coli") ~ "escherichia coli",
                             str_detect(species, "staphylococcus aureus") ~ "staphylococcus aureus",
                             str_detect(species, "lugdunensis") ~ "staphylococcus lugdunensis",
                             str_detect(species, "staphylococcus") ~ "staphylococcus species",
                             str_detect(species, "bacteroides") ~ "bacteroides species",
                             str_detect(species, "klebsiella") ~ "klebsiella species",
                             str_detect(species, "streptococcus pneumoniae") ~ "streptococcus pneumoniae",
                             str_detect(species, "dysgalactiae|Betastreptokocker grupp g|agalactiae|pyogenes") ~ "beta-hemolytic streptococci",
                             str_detect(species, "enterococcus") ~ "enterococcus species",
                             str_detect(species, "pseudomonas") ~ "pseudomonas species",
                             str_detect(species, "proteus") ~ "proteus species",
                             str_detect(species, "streptococcus") ~ "alpha-hemolytic streptococci",
                             str_detect(species, "candida") ~ "candida species",
                             str_detect(species, "clostridium") ~ "clostridium species",
                             str_detect(species, "enterobacter") ~ "enterobacter species",
                             str_detect(species, "citrobacter") ~ "citrobacter species",
                             str_detect(species, "listeria") ~ "listeria species",
                             str_detect(species, "neisseria") ~ "neisseria species",
                             str_detect(species, "fusobacterium") ~ "fusobacterium species",
                             str_detect(species, "salmonella") ~ "salmonella species",
                             str_detect(species, "aerococcus") ~ "aerococcus species",
                             TRUE ~ species),
  species = StrCap(species, "first"),
  genus = StrCap(genus, "first"),
  finding = StrCap(finding, "first"))
  return(df)
}




# Flag contaminants-----
# Purpose: to flag contaminants before deduplication

flag_contaminants <- function(df){
    df <- df %>% 
    group_by(patient_id, sample_date, finding) %>% 
    mutate(contaminant = ifelse(is.na(potential_contaminant), NA, 
                                ifelse(!potential_contaminant, FALSE, 
                                       ifelse(n_distinct(labnr) == 1, TRUE, FALSE)))) # NOTE, contmainant if < 3 bottles pos this works since we grouped by sample_date
  df <- df %>% relocate(c(potential_contaminant, contaminant), .after = "bacterial_class")
  return(df)
}

# Determine episodes ----
# Purpose: to determine how many culturing episodes each individual has had during the study period
# Note 1: the first culturing date is index date, regardless of results, so a negative followed by positive will be considered negative
# Note 2: Depends on hardcoded duplicate time indays
# Note 3: This functions sorts df on id and date - to make the loop work

determine_episodes <- function(df, dedup_time = duplicate_time) {
  temp <- df %>% 
    arrange(patient_id, sample_date) %>% 
    group_by(patient_id) %>% 
    mutate(days_from_first_sample = as.numeric(difftime(sample_date, min(sample_date), units = "days")),
           episode_no = NA) %>% 
    ungroup()
  
  # loop determining whether which episode each sample belongs to
  for (i in 1:nrow(temp)) {
    if (i == 1) { 
      temp$episode_no[i] <- 1
    } 
    else if ((temp$patient_id[i]) != (temp$patient_id[i - 1])) { # om det är en ny person sätts episod = 1
      temp$episode_no[i] <- 1
    } 
    else if (temp$days_from_first_sample[i] - temp$days_from_first_sample[i - 1] > dedup_time) { # om det inte är samma vårdtillfälle men har gått mer än 30 dagar = ny episod
      temp$episode_no[i] <- temp$episode_no[i - 1] + 1
    } 
    else {
      temp$episode_no[i] <- temp$episode_no[i - 1]
    }
  }
  
# make unique episode_id
  temp <- temp %>% 
    mutate(episode_id = paste(temp$patient_id, episode_no, sep = "")) %>% 
    select(-days_from_first_sample, -episode_no) %>% 
    relocate(c(episode_id), .after = patient_id)
  
  return(temp)
}
  
# Flag polymicrobial ----
# Purpose: to identify those with multiple different findings on the first day of each culturing episode
# Note: creates several interim variables that are not needed in the final product
flag_polymicrobial <- function(df) {
  df <- df %>% 
    group_by(episode_id) %>% 
    mutate(sample_taken_on_index_day = sample_date == min(sample_date),
           relevant_finding = ifelse((contaminant | finding == "Negative"), NA, finding),
           contaminant_finding = ifelse((!contaminant | finding == "Negative"), NA, finding))
  
  df <- df %>%     
    group_by(episode_id, sample_date) %>% 
    mutate(no_of_unique_relevant_findings = n_distinct(relevant_finding, na.rm = T), # number of different relevant findings for each individual and day
         no_of_unique_contaminants = n_distinct(contaminant_finding, na.rm = T), # number of different contaminants for each individual and day
         finding2 = ifelse(finding == "Negative", finding, 
                           ifelse(no_of_unique_relevant_findings > 1, "polymicrobial", # two different relevant findings = polymicrobial
                           ifelse(no_of_unique_contaminants > 1, "contamination", # two different contaminants = contamination
                                  finding))),
         no_of_unique_findings = n_distinct(finding))
  return(df)
}


# Deduplicate data -----
# Purpose: to deduplicate data 
# Note: in case of several similar findings on the index date, the one with shortest TTP is chosen
deduplicate_data <- function(df) {
  df <- df %>%
    filter(sample_taken_on_index_day) %>% 
    group_by(episode_id) %>% 
    filter(!(no_of_unique_findings > 1 & finding == "Negative")) %>%   # Filter away negatives if there are both positive and negatives
    filter(!(no_of_unique_relevant_findings > 0 & contaminant)) %>% # Filter away contaminants if there are both relevant and contaminants
    group_by(episode_id, finding2) %>% # group also polymicrobial etc
    slice_min(ttp, n = 1, with_ties = FALSE) %>% # Use the shortest TTP upon deduplication
  return(df)
}

# Clean episode data----
# Purpose: to clean deduplicated data

clean_episode_data <- function(episode_df) {
  episode_df <- episode_df %>% 
    mutate(finding = ifelse(finding2 == "polymicrobial", finding2, finding),
         genus = ifelse(finding2 == "polymicrobial", finding2, genus),
         bacteria_topten = ifelse(finding2 == "polymicrobial", finding2, bacteria_topten),
         bacterial_class = ifelse(finding2 == "polymicrobial", finding2, 
                                  ifelse(contaminant, "contaminant", bacterial_class))) %>% 
    select(-c(sample_taken_on_index_day, relevant_finding, contaminant_finding, 
              no_of_unique_relevant_findings, no_of_unique_contaminants, finding2, no_of_unique_findings))

return(episode_df)
}

#################### Other data

##################### Merging with other datasets
# Clean hospitalisation data -----
# Purpose: to clean hospitalisation data
# Note: takes two files as input - may need to change in other setting
clean_hosp_data <- function(hf1, hf2) {
  temp <- bind_rows(
    read_parquet(hf1) %>% 
      janitor::clean_names() %>% 
      select(patient_id = rs_pat_alias, 
             hosp_id = vardtillfalle_alias, 
             hosp_start = vardtillfalle_start_datum, 
             hosp_stop = vardtillfalle_slut_datum,
             hosp_site = vardtillfalle_avdelning,
             hosp_type = vardtillfalle_vardform_text),
    read_parquet(hf2) %>% 
      janitor::clean_names() %>%
      select(patient_id = rs_pat_alias, 
             hosp_id = vardtillfalle_alias, 
             hosp_start = vardtillfalle_start_datum, 
             hosp_stop = vardtillfalle_slut_datum,
             hosp_site =vardtillfalle_avdelning,
             hosp_type = vardtillfalle_vardform_text)) %>% 
    distinct()  %>%  
    filter(!is.na(hosp_stop),
           hosp_type == "Slutenvård") %>% 
    select(-hosp_type) %>% 
    arrange(patient_id, hosp_start) # for easier inspection only 
  write_parquet(temp, here("data", "processed_data", "hosp_data.parquet"))
}


# Add hospitalisation data ----
add_hosp_data <- function(df, cleaned_hosp_file, max_win = 3){
  temp <- left_join(x = df, 
                    y = read_parquet(cleaned_hosp_file), 
                    relationship = "many-to-many") %>% 
    filter(sample_date %within% interval(hosp_start - days(max_win),
                                         hosp_stop)) %>% 
    group_by(episode_id) %>% 
    mutate(hosp_start = min(hosp_start, na.rm = T),
           hosp_stop = max(hosp_stop, na.rm = T)) %>% # this is so that sequential hospitalisation episodes will be treated as one - THIS HAS TO BE CHANGED IF ER IS OF INTEREST
    slice_max(hosp_stop, with_ties = FALSE) %>%  # prioritising last if sequential
    left_join(x = df, y = .)
  return(temp)
}


# Add recent hospitalisation ----
add_recent_hospitalisation <- function(original_df, hosp_path, period = 90) {
  df <- left_join(original_df %>% 
                    select(patient_id, episode_id, sample_date),
                  read_parquet(hosp_path),
                  by ="patient_id") %>% 
    mutate(days_before = as.numeric(difftime(sample_date, hosp_stop, units = "days")),
           hosp_before = days_before > 1 & days_before <= period & !str_detect(hosp_site, "Aku|aku")) %>% # to avoid ER visits
    filter(hosp_before) %>% 
    group_by(episode_id) %>% 
    slice_head(n = 1) %>% 
    select(episode_id, hosp_before) %>% 
    left_join(original_df, y = .) %>% 
    mutate(hosp_before = ifelse(is.na(hosp_before), FALSE, hosp_before))
  return(df)
}


# Clean ceiling data ----
# this function cleans ceiling data in a directory with parquet files
# NOTE: filenames may have to be changed in the function
# NOTE2: calls read_ceiling below.

clean_ceiling_data <- function(dir) {
  ceiling_df <- list.files(dir,
                         pattern = "Melior_OppenVardtillfalle_Fritext.parquet|Melior_OppenVardtillfalle_Flerval.parquet|Melior_SlutenVardtillfalle_Fritext.parquet|Melior_SlutenVardtillfalle_Flerval.parquet", full.names = TRUE) %>% 
    lapply(read_ceiling_data) %>% 
    dplyr::bind_rows() %>% 
    mutate(ceiling_decision = case_when(str_detect(ceiling_decision, "ntensivvård") ~ 2,
                                        str_detect(ceiling_decision, "HLR") ~ 1,
                                        str_detect(ceiling_decision, "alliativ") ~ 3),
           ceiling_decision = factor(
             ceiling_decision, 
             levels = c(1,2,3), 
             labels = c("No CPR", "No CPR or ICU", "Palliative care"))
    ) %>% 
    filter(!is.na(ceiling_decision))

  VTF <- list.files(dir, 
                  pattern = "Melior_OppenVardtillfalle.parquet|Melior_SlutenVardtillfalle.parquet", 
                  full.names = TRUE) %>% 
  lapply(read_parquet) %>% 
    dplyr::bind_rows() %>% 
    janitor::clean_names() %>% 
    select(patient_id = rs_pat_alias, vardtillfalle_alias) %>% 
    distinct()
  
  ceiling_df <- left_join(ceiling_df, VTF) %>% 
    select(patient_id, hosp_id = vardtillfalle_alias, ceiling_date, ceiling_decision) %>% 
    distinct() %>% 
    mutate(ceiling_date = as.POSIXct(ceiling_date))
  
  write_parquet(ceiling_df, here("data", "processed_data", "ceiling_data.parquet"))
  
}

# read_ceiling 
# function that is called from clean ceiling data
# Note: this is just because of many files
read_ceiling_data <- function(path) {  
  df <- read_parquet(path) %>% janitor::clean_names() 
  
  if("fritext_modifierad_datum" %in% colnames(df)) {
    df$ceiling_date <-  df$fritext_modifierad_datum
    df$ceiling_decision <- df$fritext_varde 
  }
  else if ( "flerval_modifierad_datum" %in% colnames(df)) {
    df$ceiling_date <- df$flerval_modifierad_datum
    df$ceiling_decision <- df$flerval_varde
  }
  
  return(df)
}


# Add ceiling of care ------
  # NOTE: THis adds ceiling of care data
  # NOTE: prioritizes the highest level of ceiling and the first date
  add_ceiling <- function(df, path) {
    df <- left_join(df, read_parquet(path)) %>% 
      group_by(episode_id) %>% 
      slice_min(ceiling_date, n = 1) %>% 
      slice_max(ceiling_decision, n = 1) %>% 
      mutate(ceiling_decision = ifelse(is.na(ceiling_decision), 0, ceiling_decision),
             ceiling_decision = factor(ceiling_decision, levels = c(0,1,2,3),
                                       labels = c("No ceiling of care", "No CPR", "No CPR or ICU", "Palliative care"))) %>% 
      select(-ceiling_date)
    return(df)
  }




# Clean dialysis data ----
# Purpose: to clean data on dialysis
clean_dialysis_data <- function(dialysis_file_raw) {
  dialysis_file_clean <- read_parquet(dialysis_file_raw) %>% janitor::clean_names() %>% 
    rename(dialysis_date = aktivitet_registrerad_datum,
           patient_id = rs_pat_alias) %>% 
    select(patient_id, dialysis_date) %>% 
    arrange(patient_id, dialysis_date) %>% 
    distinct()
  write_parquet(dialysis_file_clean, here("data", "processed_data", "dialysis_data.parquet"))
}

# Add dialysis data----
add_dialysis_data <- function(df, dialysis_file_clean) {
  temp <- left_join(df, read_parquet(dialysis_file_clean),
                    by = "patient_id",
                    relationship = "many-to-many") %>%
    group_by(episode_id) %>% 
    slice_min(dialysis_date, n = 1) %>% 
    mutate(dialysis = ifelse(!is.na(dialysis_date) & dialysis_date <= sample_date, 1, 0)) %>% 
    select(-dialysis_date)
  return(temp)
}

# Clean diagnosis data ----
clean_diagnosis_data <- function() {
  VTF <- list.files(here("data", "converted_parquet"), 
                    pattern = "Melior_OppenVardtillfalle.parquet|Melior_SlutenVardtillfalle.parquet", 
                    full.names = TRUE) %>% 
    lapply(read_parquet) %>% 
    dplyr::bind_rows() %>% 
    janitor::clean_names() %>% 
    select(patient_id = rs_pat_alias, hosp_id = vardtillfalle_alias) %>% 
    distinct()
  
  diagnosis_df1 <- list.files(here("data", "converted_parquet"),
                              pattern = "Melior_OppenVardtillfalle_PatientDiagnos|Melior_SlutenVardtillfalle_PatientDiagnos", full.names = TRUE) %>% 
    lapply(read_parquet) %>% 
    dplyr::bind_rows() %>% 
    janitor::clean_names() %>% 
    select(hosp_id = vardtillfalle_alias,
           activity = aktivitet,
           diagnosis_type = diagnostyp,
           diagnosis_text = patient_diagnos_beskrivning,
           diagnosis_code = patient_diagnos_kod,
           diagnosis_date = patient_diagnos_modifierad_datum) %>% 
    distinct()
  
  diagnosis_df2 <- read_parquet(here("data", "converted_parquet", "Melior_PatientDiagnos_365InnanProv.parquet")) %>% 
    janitor::clean_names() %>% 
    select(patient_id = rs_pat_alias,
           activity = aktivitet,
           diagnosis_type = diagnostyp,
           diagnosis_text = patient_diagnos_beskrivning,
           diagnosis_code = patient_diagnos_kod,
           diagnosis_date = patient_diagnos_modifierad_datum) %>% 
    distinct()
  
  diagnosis_df3 <- read_parquet(here("data", "converted_parquet", "PMO_PatientDiagnos_365InnanProv.parquet")) %>% 
    janitor::clean_names() %>% 
    mutate(activity = as.character(NA),
           diagnosis_type = as.character(NA)) %>% 
    select(patient_id = rs_pat_alias,
           activity,
           diagnosis_type,
           diagnosis_text = patient_diagnos_namn,
           diagnosis_code = patient_diagnos_kod,
           diagnosis_date = datarad_version_data_datum) %>% 
    distinct()
  
  test <- left_join(diagnosis_df1, VTF) %>% 
    select(patient_id, hosp_id, diagnosis_date, 
           diagnosis_code, diagnosis_text, diagnosis_type, activity)
  
  test2 <- bind_rows(diagnosis_df2, diagnosis_df3) %>% mutate(hosp_id = as.numeric(NA)) %>% 
    select(patient_id, hosp_id, diagnosis_date, 
           diagnosis_code, diagnosis_text, diagnosis_type, activity)
  
  test3 <- bind_rows(test, test2) %>% 
    arrange(patient_id)
  
  write_parquet(test3, here("data", "processed_data", "diagnosis_data.parquet"))
}

# Add main diagnosis ----

add_main_diagnosis <- function(original_df, diagnosis_file){
  df1 <- read_parquet(diagnosis_file) %>% 
    filter(str_detect(activity, "Epikris|Läk|Op-berättelse|Överflyttningsanteckning"),
           str_detect(diagnosis_type, "Huvuddiagnos|huvuddiagnos"),
           !is.na(hosp_id)) %>% 
    group_by(hosp_id) %>% 
    slice_max(diagnosis_date, n = 1, with_ties = FALSE) %>% 
    select(hosp_id, diagnosis_code, diagnosis_text)
  
  temp <- left_join(original_df, df1, by = "hosp_id") 
  return(temp)
}




# Add infection focus -----

add_infection_focus <- function(original_df, id_cat_file) {
  temp <- partial_join(
    x = original_df %>% select(episode_id, diagnosis_code), 
    y  = read_excel(id_cat_file), 
    by_x = "diagnosis_code", 
    pattern_y = "icd_code") %>% 
    arrange(episode_id) %>% 
    select(-icd_code) %>% 
    left_join(x = original_df, y = .) %>% 
    mutate(id_category = ifelse(is.na(id_category), "No Infection", id_category))
  
  return(temp)
}


# Add Charlson score----

add_charlson <- function(episode_file, diagnosis_file) {
  df1 <- read_parquet(diagnosis_file)
  df2 <- left_join(
    (episode_file %>% select(episode_id, patient_id, sample_date, hosp_id)), df1, 
    by = "patient_id") %>% 
    mutate(days_between_sample_and_diagnosis = 
             as.integer(difftime(sample_date, as.Date(diagnosis_date), units = "days"))) %>% 
    filter((is.na(hosp_id.y) & days_between_sample_and_diagnosis %in% c(0:365)) | hosp_id.x == hosp_id.y) %>% 
    select(episode_id, diagnosis_code) %>% 
    mutate(ami_diagnosis_before = str_starts(diagnosis_code, "I21|I22|I252"),
           heart_failure_before = str_starts(diagnosis_code, "I110|I130|I132|I255|I420|I426|I427|I428|I429|I43|I50"),
           peripheral_vascular_before = str_starts(diagnosis_code, "I70|I71|I731|I738|I739|I771|I790|I792|K55"),
           cerebrovascular_before = str_starts(diagnosis_code, "G45|I60|I61|I62|I63|I64|I67|I69"),
           copd_before = str_starts(diagnosis_code, "J43|J44"),
           other_pulmonary_before = str_starts(diagnosis_code, "J41|J42|J45|J46|J47|J60|J61|J62|J63|J64|J65|J66|J67|J68|J69|J70"),
           rheumatic_disease_before = str_starts(diagnosis_code, "M05|M06|M123|M070|M071|M072|M073|M08|M13|M30|M313|M314|M315|M316|M32|M33|M34|M350|M351|M353|M45|M46"),
           dementia_before = str_starts(diagnosis_code, "F00|F01|F02|F03|F051|G30|G311|G318"),
           hemiplegia_before = str_starts(diagnosis_code, "G114|G80|G81|G82|G830|G831|G832|G833|G838"),
           diabetes_before = str_starts(diagnosis_code, "E100|E101|E110|E111|E120|E121|E130|E131|E140|E141"),
           diabetes_end_organ_before = str_starts(diagnosis_code, "E102|E103|E104|E105|E107|E112|E113|E114|E115|E116|E117|E122|E123|E124|E125|E126|E127|
                                                E132|E133|E134|E135|E136|E137|E142|E143|E144|E145|E146|E147"),
           moderate_severe_kidney_before = str_starts(diagnosis_code, "N032|N033|N034|N035|N036|N037|N052|N053|N054|N055|N056|N057|N11|N18|N19|N250|
                                             I120|I131|Q611|Q612|Q613|Q614|Z49|Z940|Z992"),
           mild_liver_before = str_starts(diagnosis_code, "B15|B16|B17|B18|B19|K703|K73|K746|K703|K754"),
           moderate_severe_liver_before = str_starts(diagnosis_code, "R18|I850|I859|I982|I983"),
           ulcer_before = str_starts(diagnosis_code, "K25|K26|K27|K28"),
           malignancy_before = str_starts(diagnosis_code, "C") & str_starts(diagnosis_code, "C77|C78|C79|C80", negate = TRUE),
           metastatic_before = str_starts(diagnosis_code, "C77|C78|C79|C80"),
           aids_before = str_starts(diagnosis_code, "B20|B21|B22|B23|B24|F024|O987|R75|Z114|Z219Z711"))
  
  df3 <- df2 %>% 
    group_by(episode_id) %>% 
    mutate(ami = ifelse(sum(ami_diagnosis_before, na.rm = T) > 0, 1, 0),
           heart_failure = ifelse(sum(heart_failure_before, na.rm = T) > 0, 1, 0),
           peripheral_vascular= ifelse(sum(peripheral_vascular_before, na.rm = T) > 0, 1, 0),
           cerebrovascular = ifelse(sum(cerebrovascular_before, na.rm = T) > 0, 1, 0),
           copd = ifelse(sum(copd_before, na.rm = T) > 0, 1, 0),
           other_pulmonary = ifelse(sum(other_pulmonary_before, na.rm = T) > 0, 1, 0),
           rheumatic_disease = ifelse(sum(rheumatic_disease_before, na.rm = T) > 0, 1, 0),
           dementia = ifelse(sum(dementia_before, na.rm = T) > 0, 1, 0),
           hemiplegia = ifelse(sum(hemiplegia_before, na.rm = T) > 0, 1, 0),
           diabetes = ifelse(sum(diabetes_before, na.rm = T) > 0, 1, 0),
           diabetes_complications = ifelse(sum(diabetes_end_organ_before, na.rm = T) > 0, 1, 0),
           moderate_severe_kidney = ifelse(sum(moderate_severe_kidney_before, na.rm = T) > 0, 1, 0),
           mild_liver = ifelse(sum(mild_liver_before, na.rm = T) > 0, 1, 0),
           moderate_severe_liver = ifelse(sum(moderate_severe_liver_before, na.rm = T) > 0, 1, 0),
           ulcer = ifelse(sum(ulcer_before, na.rm = T) > 0, 1, 0),
           malignancy = ifelse(sum(malignancy_before, na.rm = T) > 0 & sum(metastatic_before, na.rm = T) == 0, 1, 0),
           metastatic = ifelse(sum(metastatic_before, na.rm = T) > 0, 1, 0),
           aids = ifelse(sum(aids_before, na.rm = T) > 0, 1, 0),
           cci_points = ami+ heart_failure + peripheral_vascular + 
             cerebrovascular + copd + other_pulmonary + rheumatic_disease +
             dementia + 2 * hemiplegia + diabetes + 2 * diabetes_complications + 
             2 * moderate_severe_kidney + mild_liver + 3 * moderate_severe_liver + 
             ulcer + 2 * malignancy + 6 * metastatic + 6 * aids) %>% # to not double count metastatic malignancies 
    distinct(across(episode_id), .keep_all = TRUE) %>% 
    select(-c(diagnosis_code, ends_with("before")))
  
  output <- left_join(episode_file , df3 %>% select(episode_id, cci_points))
  write_parquet(df3 %>% select(-cci_points), here("data", "processed_data", "cci_data.parquet"))
  return(output)
}







# Add specific comorbidities ----

add_comorbidities <- function(original_df, diagnosis_file) {
  temp <- read_parquet(diagnosis_file) %>% 
    mutate(diagnosis_group = case_when(
      str_starts(diagnosis_code, "C|Z511") ~ "Malignancy", # including treatment with cytotoxic drugs
      str_starts(diagnosis_code, "D5|D60|D61|D62|D63|D64") ~ "Anemia",
      str_starts(diagnosis_code, "B20|B21|B22|B23|B24|D70|D80|D81|D82|D83|D84|Z94") ~ "Immunodeficiency",
      str_starts(diagnosis_code, "E10|E11|E12|E13|E14") ~ "Diabetes",
      str_starts(diagnosis_code, "F") ~ "Psychiatric disorder",
      str_starts(diagnosis_code, "G|I6") & str_starts(diagnosis_code, "G0|G47", negate = T) ~ "Neurologic disease",
      str_starts(diagnosis_code, "I10|I11|I12|I13|I14|I15") ~ "Hypertension",
      str_starts(diagnosis_code, "I05|I06|I07|I08|I09|I20|I21|I22|I23|I24|I25|I34|I35|I36|I37|I42|I43|I44|I45|I47|I48|I49|I50|Z95") ~ "Cardiac disease", #arrythmias, heart failure, ischemic heart disease, valve disease
      str_starts(diagnosis_code, "I7|I8") ~ "Peripheral vascular disease",
      str_starts(diagnosis_code, "I26|I27|I28|J4|J84") ~ "Pulmonary disease",
      str_starts(diagnosis_code, "B1|K7") ~ "Hepatic disease",
      str_starts(diagnosis_code, "L") & str_starts(diagnosis_code, "L0", negate = T) ~ "Skin disease",
      str_starts(diagnosis_code, "M") & str_starts(diagnosis_code, "M01|M02|M03|M549|M545|M791|M796|M81", negate = T) ~ "Musculoskeletal disease",
      str_starts(diagnosis_code, "N0|N1|N2|N3|N4|Z49") & str_starts(diagnosis_code, "N10|N17|N30|N39|N136", negate = T) ~ "Genitourinary disease"
      ))
  
  temp2 <- left_join(
    (original_df %>% select(episode_id, patient_id, sample_date, hosp_id)), temp, by = "patient_id") %>% 
    mutate(days_between_sample_and_diagnosis = as.integer(difftime(sample_date, as.Date(diagnosis_date), units = "days"))) %>% 
    filter((is.na(hosp_id.y) & days_between_sample_and_diagnosis %in% c(0:365)) | hosp_id.x == hosp_id.y) %>% 
    filter(!is.na(diagnosis_group)) %>% 
    select(episode_id, diagnosis_group) %>% 
    distinct() %>% 
    mutate(value = 1) %>% # just an indicator for positivity
    pivot_wider(id_cols = episode_id, names_from = diagnosis_group, values_from = value) %>% 
    left_join(x = original_df %>% select(episode_id), y = .) %>% 
    mutate(across(everything(), ~replace_na(.,0))) %>% 
    janitor::clean_names()
  
  output_df <- left_join(original_df, temp2, by = "episode_id")
  return(output_df)
}

# Clean medications data------

clean_meds_before <- function() {
  df1 <- read_parquet(here("data", "converted_parquet", "Melior_Ordination_Utdelning_2MI2MEB_Ankomstar_2021.parquet")) %>% 
    janitor::clean_names() %>% 
    distinct(across(-mikrobiologi_prov_alias)) %>% 
    select(patient_id = rs_pat_alias,
           drug_atc = ordination_atc_kod,
           drug_name = ordination_preparat_namn,
           drug_date_in = ordination_insatt_datum)
  
  df2 <- read_parquet(here("data", "converted_parquet", "Melior_Ordination_Utdelning_2MI2MEB_Ankomstar_2022.parquet")) %>% 
    janitor::clean_names() %>% 
    distinct(across(-mikrobiologi_prov_alias)) %>% 
    select(patient_id = rs_pat_alias,
           drug_atc = ordination_atc_kod,
           drug_name = ordination_preparat_namn,
           drug_date_in = ordination_insatt_datum)
  
  df3 <- read_parquet(here("data", "converted_parquet", "Melior_Ordination_Utdelning_2MI2MEB_Ankomstar_2023.parquet")) %>% 
    janitor::clean_names() %>% 
    distinct(across(-mikrobiologi_prov_alias)) %>% 
    select(patient_id = rs_pat_alias,
           drug_atc = ordination_atc_kod,
           drug_name = ordination_preparat_namn,
           drug_date_in = ordination_insatt_datum)
  
  df4 <- read_parquet(here("data", "converted_parquet", "vw_LakemforskrAllaPerioder_2MI2MEB.parquet")) %>% 
    janitor::clean_names() %>% 
    distinct(across(-mikrobiologi_prov_alias)) %>% 
    select(patient_id = rs_pat_alias,
           drug_atc = atc_kod,
           drug_name = atckodstext,
           drug_date_in = kopdatum)
  
  df5 <- bind_rows(df1, df2, df3, df4) %>% 
    distinct() %>% 
    arrange(patient_id, drug_date_in)
  
  rm(df1,df2,df3,df4)
  
  write_parquet(df5, here("data", "processed_data", "meds_before_data.parquet"))
}

# add meds before data -----

add_ab_immuno_before <- function(original_df, path){
  df1 <- left_join(original_df %>% select(episode_id, patient_id, sample_date),
                 read_parquet(path), by = "patient_id", relationship = "many-to-many") %>% 
  mutate(days_between = as.numeric(difftime(sample_date, drug_date_in, units = "days")),
         ab_before = str_starts(drug_atc, "J01") & 
           days_between > 0 & days_between < 14,
         immunosupp_before = str_starts(drug_atc, "L01|L04") &
           days_between > 0 & days_between < 60) %>% 
  filter(days_between > 0)


  df2 <- df1 %>% select(episode_id, ab_before) %>% filter(ab_before) %>% 
  slice_head(n = 1)
  
  df3 <- df1 %>% select(episode_id, immunosupp_before) %>% filter(immunosupp_before) %>% 
  slice_head(n = 1)
  
  output_df <- left_join(original_df, df2) %>% 
  mutate(ab_before = ifelse(is.na(ab_before), FALSE, ab_before))
  
  output_df <- left_join(output_df, df3) %>% 
  mutate(immunosupp_before = ifelse(is.na(immunosupp_before), FALSE, immunosupp_before))
  
  return(output_df)
}



# Fix culture time (needs other data) ----------------------------
# Purpose: to approximate the time of culturing, in case this is not provided in raw data
# Note: priority order as follows
# 1. culture time was provided 
# 2: blood tests taken same day
# 3: hospitalisation started same day

fix_culture_time <- function(original_df, lab_file_parquet) {
  temp_lab <- read_parquet(lab_file_parquet) %>% 
    select(patient_id, lab_time)
  
  temp <- left_join(original_df, temp_lab, relationship = "many-to-many") %>% 
    mutate(culture_source_priority = ifelse(!is.na(sample_time) & !strftime(sample_time, format="%H:%M:%S") == "00:00:00",0,
                                            ifelse(!is.na(lab_time) & sample_date == as.Date(lab_time),1,
                                                   ifelse(!is.na(hosp_start) & sample_date == as.Date(hosp_start),2, 
                                                          3))))  %>% 
    group_by(episode_id) %>% 
    slice_min(culture_source_priority) %>% 
    distinct(across(-lab_time), .keep_all = TRUE)
  
  temp <- temp %>% 
    mutate(sample_datetime = as.POSIXct(ifelse(
      culture_source_priority == 0, sample_datetime,
      ifelse(culture_source_priority == 1, lab_time,
             ifelse(culture_source_priority == 2, hosp_start, NA))), tz = "UTC")) %>% 
    select(-lab_time) %>% 
    mutate(culture_source_priority = factor(culture_source_priority, levels = c(0,1,2,3),
                                            labels = c("Culture time provided",
                                                       "From venepuncture time",
                                                       "From hospitalisation start",
                                                       "None of the above")))
  
  return(temp)
}






# Clean outcome data ----
# Purpose: to combine data on icu and mortality outcomes

clean_outcome_data <- function(deceased_file, icu_file) {
  df1 <- read_parquet(deceased_file) %>% 
    janitor::clean_names() %>% 
    select(patient_id = rs_pat_alias,
         deceased = avliden,
         deceased_date = avliden_datum) %>% 
    arrange(patient_id)
  
  df2 <- read_parquet(icu_file) %>% 
    janitor::clean_names() %>% 
    select(patient_id = rs_pat_alias,
           hosp_id = vardtillfalle_alias,
           icu = inskriven_iva,
           icu_date = inskrivning_iva) %>%
    distinct() %>% 
    arrange(patient_id)
  
  df3 <- left_join(df1, df2, by = "patient_id")
  write_parquet(df3, here("data", "processed_data", "outcome_data.parquet"))
}

# Add vital status--------

add_vital_status <- function(original_df, outcome_file_path) {
  df <- left_join(
    original_df, 
    read_parquet(outcome_file_path) %>%
      filter(patient_id %in% episode_data$patient_id) %>% 
      select(patient_id, deceased, deceased_date) %>% 
      distinct()) %>% 
    mutate(
      in_hosp_mortality = if_else(!deceased, FALSE, deceased_date %within% interval(hosp_stop - days(1), hosp_stop + days(1))), # this is used to obtain a date and time for death in case of in-hospital mortality (otherwise only date is registered, which also may be incorrect by one day if death is near midnight)
      deceased_date = if_else(!is.na(in_hosp_mortality) & in_hosp_mortality, hosp_stop, deceased_date),
      fup_date_deceased = if_else(!deceased, as.Date("2024-03-12"), deceased_date),
      fup_time_deceased = as.numeric(difftime(fup_date_deceased, sample_date, units = "days")),
      fup_time_deceased = if_else(fup_time_deceased <= -1, NA, fup_time_deceased), #not interested in those in ICU before
      deceased_30 = if_else(deceased & fup_time_deceased <= 30,1,0),
      deceased_90 = if_else(deceased & fup_time_deceased <= 90,1,0),
      deceased_365 = if_else(deceased & fup_time_deceased <= 365,1,0)) %>%
    select(-fup_date_deceased) %>% 
    group_by(patient_id) %>% 
    mutate(deceased_30 = ifelse(is.na(deceased_30) & episode_id != max(as.numeric(episode_id)), # if hospitalised again after 30 days,30 day mortality is negative
                                0, deceased_30))
  return(df)
}

# Add icu ----------

add_icu <- function(original_df, outcome_file) {
  
  temp <- original_df %>% select(episode_id, patient_id, hosp_id)

  temp2 <- read_parquet(outcome_file) %>% 
    filter(patient_id %in% episode_data$patient_id,
           hosp_id %in% episode_data$hosp_id) %>%
    select(patient_id, hosp_id, icu, icu_date) %>% 
    left_join(original_df, y = .)  %>% 
    mutate(
      icu = if_else(is.na(icu) & is.na(hosp_id), 0 , icu), #if ICU is missing but not hospitalised -> no ICU
      ref_date = if_else(!is.na(icu) & icu == 1, icu_date, hosp_stop), # this is the reference date for ICU
      fup_time_icu = as.numeric(difftime(ref_date, sample_date, units = "days")),
      icu = if_else(!is.na(icu) & icu == 1 & !is.na(fup_time_icu) & fup_time_icu <= -1, 0, icu), # not ICU before sampling
      fup_time_icu = if_else(fup_time_icu <= -1, NA, fup_time_icu), #not interested in those in ICU before
      icu_30 = case_when(is.na(icu) | (icu == 1 & is.na(fup_time_icu)) ~ NA,
                         fup_time_icu <= 30 & icu == 1 ~ 1,
                         TRUE ~ 0),
      dead_or_icu = case_when(is.na(deceased) & is.na(icu) ~ NA,
                              deceased ~ 1,
                              icu == 1 ~ 1,
                              TRUE ~ 0),
      dead_or_icu2 = case_when(is.na(deceased) | is.na(icu) ~ NA,
                              deceased ~ 1,
                              icu == 1 ~ 1,
                              TRUE ~ 0),
      fup_time_both = case_when(is.na(fup_time_icu) & is.na(fup_time_deceased) ~ NA,
                                !is.na(fup_time_icu) & icu == 1 ~ fup_time_icu, # ICU admisison is always before death so if ICU admitted use this
                                !is.na(fup_time_deceased) ~ fup_time_deceased, # unless fup_time deceased is missing use this
                                !is.na(fup_time_icu) ~ fup_time_icu),
      dead_or_icu_30 = case_when(is.na(dead_or_icu) | is.na(fup_time_both) ~ NA,
                                   dead_or_icu == 1 & fup_time_both <= 30 ~ 1,
                                   TRUE ~ 0),
      dead_or_icu_30_2 = case_when(is.na(dead_or_icu2) | is.na(fup_time_both) ~ NA,
                                   dead_or_icu2 == 1 & fup_time_both <= 30 ~ 1,
                                   TRUE ~ 0))
      
      
  return(temp2)
}



# Clean lab data -----
# Purpose: to clean laboratory data
# Note: requires a list of laboratory values, this should be double-checked vs raw data to not miss any misspellings etc.

clean_lab_data <- function(original_df, raw_lab_data, lab_list) {
  clean_lab_data <- read_parquet(raw_lab_data) %>% 
    janitor::clean_names() %>%
    filter(str_detect(labanalys_beskrivning, lab_list) &
             str_detect(labanalys_beskrivning, "Ext|X", negate = T)) %>%
    mutate(
      lab_value = as.numeric(ifelse(str_detect(analyssvar_varde, ">|<"), 
                                    str_sub(analyssvar_varde,2), 
                                    analyssvar_varde)),
      lab_name = case_when(
        str_detect(labanalys_beskrivning, "P-CRP|B-CRP") ~ "CRP",
        str_detect(labanalys_beskrivning, "P-Kreatinin") ~ "Creatinine",
        str_detect(labanalys_beskrivning, "P-Laktat|B-Laktat") ~ "Lactate",
        str_detect(labanalys_beskrivning, "Prokalcitonin") ~ "Procalcitonin",
        str_detect(labanalys_beskrivning, "B-Trombocyter") ~ "Platelets",
        str_detect(labanalys_beskrivning, "P-Bilirubin") ~ "Bilirubin")) %>% 
    select(patient_id = rs_pat_alias,
           lab_time = analyssvar_provtagning_datum,
           lab_name,
           lab_value) %>% 
    distinct() %>% 
    filter(!is.na(lab_value),
           !is.na(lab_name),
           lab_value > 0) %>% 
    arrange(patient_id)
  
  clean_lab_data <- left_join(
    x = original_df %>% select(
      patient_id, episode_id, sample_datetime, ttp, ttp_strata, has_cabinet), y = clean_lab_data) %>% 
    mutate(hours_from_culture = interval(sample_datetime, lab_time) %/% hours(1))
  
  write_parquet(clean_lab_data, here("data", "processed_data", "lab_data.parquet"))
}

# Add baseline lab data ------
# addition of baseline lab to wide format file

add_baseline_lab <- function(original_df, lab_file, time_frame = 24) {
  temp <- read_parquet(lab_file) %>% 
    filter(between(hours_from_culture, -time_frame, time_frame)) %>% 
    group_by(episode_id, lab_name) %>% 
    slice_min(abs(hours_from_culture), with_ties = FALSE) %>% 
    select(episode_id, lab_name, lab_value) %>% 
    pivot_wider(id_cols = episode_id, names_from = lab_name, values_from = lab_value)
  
  output_df <- left_join(original_df, temp, by = "episode_id")
  return(output_df)
}



# Clean NEWS data ----
# Purpose: to clean raw NEWS data
# Note: distinct across removes 75% - which are duplicates - file is large
# At least 3 digits is needed for temperature
# NEWS values are written with comma as decimal point
# There are quite a few implausible values that are corrected already here in the last part of the function
# There is also a misspecification of the Obstetric news. all blood pressures with capital O are systolic, all with
# small o are diastolic

clean_news_data <- function(original_df, raw_news_data) {
  clean_news <- read_parquet(raw_news_data) %>% 
    janitor::clean_names() %>% 
    distinct(across(-mikrobiologi_prov_alias), .keep_all = TRUE) %>% 
    rename(patient_id = rs_pat_alias,
           NEWS_time = tal_modifierad_datum,
           NEWS_name = term_namn) %>% 
    mutate(NEWS_value = as.numeric(scan(text = tal_varde, dec = ",", sep = ".")), # changes comma to dot separator
           NEWS_name = case_when(str_detect(NEWS_name, "totalp") ~ "NEWS total score",
                                 str_detect(NEWS_name, "temper") ~ "Temperature",
                                 str_detect(NEWS_name, "blodtryck2") ~ "Diastolic blood pressure",
                                 str_detect(NEWS_name, "blodtryck") ~ "Systolic blood pressure",
                                 str_detect(NEWS_name, "puls") ~ "Heart rate",
                                 str_detect(NEWS_name, "medvetandeg") ~ "Mental alteration",
                                 str_detect(NEWS_name, "blodtryck") ~ "Systolic blood pressure",
                                 str_detect(NEWS_name, "syremättn") ~ "Oxygen saturation",
                                 str_detect(NEWS_name, "andningsfre") ~ "Respiratory rate"),
           NEWS_value = case_when(str_detect(NEWS_name, "NEWS total score") & NEWS_value %in% c(0:25) ~ NEWS_value, # these are plausible values
                                  str_detect(NEWS_name, "Temperature") & NEWS_value <=45 & NEWS_value >= 25 ~ NEWS_value,
                                  str_detect(NEWS_name, "Diastolic blood pressure") & NEWS_value %in% c(20:200) ~ NEWS_value,
                                  str_detect(NEWS_name, "Systolic blood pressure") & NEWS_value %in% c(30:300) ~ NEWS_value,
                                  str_detect(NEWS_name, "Heart rate") & NEWS_value %in% c(10:300) ~ NEWS_value,
                                  str_detect(NEWS_name, "Mental alteration") & NEWS_value %in% c(0:8) ~ NEWS_value,
                                  str_detect(NEWS_name, "Oxygen saturation") & NEWS_value %in% c(40:100) ~ NEWS_value,
                                  str_detect(NEWS_name, "Respiratory rate") & NEWS_value %in% c(5:80) ~ NEWS_value,
                                  TRUE ~ NA), #setting to NA if not plausible values
           NEWS_value = ifelse(NEWS_name == "Mental alteration" & NEWS_value != 0, 1, NEWS_value)) %>% # dichotomising this as different were used
    filter(!NEWS_name %in% c("Diastolic blood pressure"), !is.na(NEWS_value)) %>%
    select(patient_id, NEWS_time, NEWS_name, NEWS_value) %>% 
    left_join(x = original_df %>% select(patient_id, episode_id, sample_datetime, ttp, ttp_strata, has_cabinet), y =.) %>% 
    mutate(hours_from_culture = interval(sample_datetime, NEWS_time) %/% hours(1))
   
    
  write_parquet(clean_news, here("data", "processed_data", "news_data.parquet"))
}

# Add baseline NEWS data -----
add_baseline_news <- function(original_df, news_file, time_frame = 24) {
  temp <- read_parquet(news_file) %>% 
    filter(between(hours_from_culture, -time_frame, time_frame)) %>% 
    group_by(episode_id, NEWS_name) %>% 
    slice_min(abs(hours_from_culture), with_ties = FALSE) %>% 
    select(episode_id, NEWS_name, NEWS_value) %>% 
    pivot_wider(id_cols = episode_id, names_from = NEWS_name, values_from = NEWS_value)
  
  output_df <- left_join(original_df, temp, by = "episode_id") %>% 
    janitor::clean_names()
  return(output_df)
}


## tabl fix --------

tblfix <- function(data, indicator = "Yes", decimal){ 
  if(is.numeric(data)){
    paste(round(median(data, na.rm = T),decimal), " (",
          round(quantile(data, 0.25, na.rm = T),decimal), "-",
          round(quantile(data, 0.75, na.rm = T),decimal), ")", sep = "")
  }
  else {
    paste(sum(data == indicator, na.rm = T), " (",
          round((sum(data == indicator, na.rm = T) / sum(!is.na(data))*100),0), "%)", sep = "")
    
  }
}

## this is the same function but it provides denominators if n != N
tblfix2 <- function(data, indicator = "Yes", decimal=1){ 
  if(is.numeric(data)){
    paste(round(median(data, na.rm = T),decimal), " (",
          round(quantile(data, 0.25, na.rm = T),decimal), "-",
          round(quantile(data, 0.75, na.rm = T),decimal), ")", sep = "")
  }
  else {
    paste(sum(data == indicator, na.rm = T), 
          paste(ifelse(sum(is.na(data)) == 0, "",  
                       paste("/",sum(!is.na(data)), sep = ""))),
          " (",
          round((sum(data == indicator, na.rm = T) / sum(!is.na(data))*100),0), "%)", sep = "")
    
  }
}


# ggplot theme -----
theme_GT <- function(base_size=8, base_family="sans") {
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5, margin = margin(0,0,20,0)),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =2),
           axis.title.x = element_text(vjust = -2),
           axis.text = element_text(size = rel(0.8)), 
           axis.line = element_line(colour = "black", size = 0.2),
           axis.ticks = element_line(),
           panel.grid.major = element_line(colour="#f0f0f0"),
           panel.grid.minor = element_blank(),
           legend.key = element_rect(colour = NA),
           legend.position = "right",
           #legend.direction = "horizontal",
           legend.box = "vetical",
           legend.key.size= unit(0.5, "cm"),
           #legend.margin = unit(0, "cm"),
           legend.title = element_text(),
           plot.margin=unit(c(10,5,4,5),"mm"),
           strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
           strip.text = element_text(size = rel(0.9))
   ))
  
}


# Make map function ----
make_map <- function(){
  library(sf)
  library(pxweb)
  library(ggspatial)
  # reading in borders for communities (kommun) and regions (lan)
  lan <- read_sf(here("data", "shape_svenska_210505", "LanSweref99TM", "Lan_Sweref99TM_region.shp"))
  kommun <- read_sf(here("data","shape_svenska_210505","KommunSweref99TM","Kommun_Sweref99TM_region.shp"))
  
  # Download population data by kommun from api
  px_data_kn <- pxweb_get(url = "https://api.scb.se/OV0104/v1/doris/sv/ssd/BE/BE0101/BE0101C/BefArealTathetKon",
                          query = list("Region"=kommun$KnKod,"Kon" = c("1+2"),
                                       "ContentsCode"=c("BE0101U1"), "Tid"=c("2019")))
  
  # Convert to data.frame, namning to KnNamn for matching
  kn_df <- as.data.frame(px_data_kn, column.name.type = "text", variable.value.type = "text") %>%
    rename(KnNamn = region,
           `Population density` = `Invånare per kvadratkilometer`)
  
  # joining geometry with population density
  kn_geo_df <- left_join(kn_df, kommun)
  
  # this is needed to position labels for lan
  lan_points <- cbind(lan, st_coordinates(st_centroid(lan$geometry)))
  
  # cities with hospitals
  cities <- data.frame(city = c("Malmö", "Lund", "Helsingborg", "Ängelholm", "Kristianstad", "Hässleholm", 
                                "Simrishamn", "Ystad", "Trelleborg", "Landskrona")) %>%
    # hard coded coordinates to hospitals
    mutate(X = c(374153.53, 386736.96, 358000, 367321, 448483, 424078, 458407, 425088, 383460, 364972),
           Y = c(6162214.35, 6175417.96, 6213609, 6234800, 6209756, 6224538, 6156524, 6143690, 6138736, 6194624),
           `On-site blood culture cabinet` = c(TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE))
  
  cities
  
  figS1 <- ggplot(lan)+
    geom_sf(data = kn_geo_df, aes(geometry = geometry,fill = `Population density`), colour = "black", lwd = 0.05)+
    geom_sf(colour = "dark blue", size = 0.5, alpha = 0.1)+
    geom_point(data = cities, aes(x = X, y = Y, col = `On-site blood culture cabinet`), size = 3)+
    geom_text(data = cities, aes(x = X, y = Y+6000, label = city), size = 3, col = "dark blue")+
    geom_text(data = lan_points, aes(x = X, y = Y, label = LnNamn),
              col = "dark blue", fontface = "bold", size = 4)+
    coord_sf(xlim = c(300000, 490000), ylim = c(6100000, 6270000), expand = FALSE)+
    annotation_scale(location = "bl", width_hint = 0.5,
                     bar_cols = c("dark blue", "white"),
                     height = unit(0.15, "cm"),
                     text_col = "dark blue",
                     line_col = "dark blue")+
    annotation_north_arrow(pad_y = unit(0.22, "in"), 
                           style = north_arrow_fancy_orienteering(line_col = "dark blue",
                                                                  fill = c("white", "dark blue"),
                                                                  text_col = "dark blue",
                                                                  text_size = 9))+
    annotate(geom = "text", x =450000 , y = 6120000, label = "the Baltic sea",
             fontface = "italic", size = 4, col = "dark blue")+
    scale_fill_fermenter(breaks = c(0,10,100,1000),
                         limits = c(0,10000),
                         direction = 1
    )+
    labs(fill = "Population density\n(inh. per square km)")+
    scale_color_manual(values = c("dark blue", "red"))+
    theme(panel.grid.major = element_line(color = gray(0.6), linetype = "dashed",
                                          size = 0.1),
          panel.background = element_rect(fill = "aliceblue"),
          axis.title = element_blank(),
          axis.line = element_line(colour = "dark blue", size = 0.25),
          axis.text = element_text(colour = "dark blue"),
          axis.ticks = element_line(colour = "dark blue", size = 0.25),
          legend.background = element_rect(color = "dark blue", size = 0.25),
          legend.text = element_text(color = "dark blue", size = 9),
          legend.title = element_text(colour = "dark blue", size = 9))
  
  return(figS1)
}

