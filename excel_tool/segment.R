if(!"pacman" %in% rownames(installed.packages())) {
  install.packages("pacman", repos = "https://cran.rstudio.com/")
}
if(!"Biobase" %in% rownames(installed.packages())) {
  install.packages("BiocManager", repos = "https://cran.rstudio.com/")
  BiocManager::install("Biobase")
}
pacman::p_load(
  dplyr, tibble, stringr, haven, assertthat, purrr, furrr, tidyr, rlang, readr,
  NMF, archetypes, update = TRUE
)
#' Function for creating segment assignments in batch form, used to interface 
#' with Excel tool. 
#' @param data_path String indicating the windows file path to the 
#' data file (using forward slashes)
#' @param input_vars String with space delimited variable names.  
#' @param id_var String indicating the variable name of the id variable
#' @param algorithms Character vector indicating what methods by which to 
#' segment.
#' @param save_path Path to excel file name to where to save the results
segment <- function(data_path, input_vars, id_var, kmin, kmax, algorithms, 
                    save_path) {
  assert_that(file.exists(data_path))
  raw_data <- read_sav(data_path) %>% 
    rowid_to_column(".rowid")
  
  input_vars_vec <- str_split(c(id_var, input_vars), " ") %>% 
    purrr::list_c() %>% 
    trimws()
  assert_that(
    all(input_vars_vec %in% names(raw_data)), 
    msg = paste(
      "Variables:",
      paste(!input_vars %in% names(raw_data), collapse = ","), 
      "do not exist in the data provided"
    )
  )
  na_check <- map_lgl(c(id_var, input_vars_vec), ~anyNA(raw_data[[.x]]))
  assert_that(
    all(!na_check), 
    msg = paste(
      "Variables:", 
      paste(input_vars_vec[na_check], collapse = ","), 
      "contain missing values, and cannot serve as inputs"
    )
  )
  assert_that(
    !anyDuplicated(raw_data[[id_var]]), 
    msg = "The id variable specified contains duplicates and cannot be used"
  )
  id <- select(raw_data, all_of(id_var))
  inputs <- select(raw_data, all_of(input_vars_vec))
  splan <- expand_grid(
    algorithm = keep(algorithms, ~.x != "twostep"), k = kmin:kmax
  )
  has_twostep <- "twostep" %in% algorithms
  
  plan(multisession)
  
  if("wardsquared" %in% splan$algorithm) {
    dist <- dist(select(inputs, -1))
    hclust <- hclust(dist)
  } else {
    hclust <- NULL
  }
  results <- future_pmap_dfc(
    .l = as.list(splan), 
    .f = function(algorithm, k, d = select(inputs, -1), hc = hclust) {
      strip_newlines <- function(string) {
        str_replace_all(string, "[\\r\\n]", "") %>% 
          str_replace_all("  ", "")
      }
      archetypes <- function(data, k) {
        at <- archetypes::archetypes(data, k)
        labels <- 1:k %>% set_names(paste("Segment", 1:k))
        label <- str_glue("Segment archetypes: {k} cluster solution")
        tibble(max.col(at$alphas)) %>% 
          set_names(str_glue("archetypes_{k}")) %>% 
          mutate_all(~labelled_spss(., labels = labels, label = label))
      }
      kmeans <- function(data, k) {
        kms <- stats::kmeans(
          x = data, 
          centers = k, 
          iter.max = 1000L, 
          nstart = 50, 
          algorithm = "Lloyd"
        )
        labels <- 1:k %>% set_names(paste("Segment", 1:k))
        label <- str_glue("Segment kmeans: {k} cluster solution")
        tibble(kms$cluster) %>% 
          set_names(str_glue("kmeans_{k}")) %>% 
          mutate_all(~labelled_spss(., labels = labels, label = label))
      }
      
      wardsquared <- function(data, k) {
        labels <- 1:k %>% set_names(paste("Segment", 1:k))
        label <- str_glue("Segment ward squared: {k} cluster solution")
        cutree(hc, k = k) %>% 
          tibble() %>% 
          set_names(str_glue("wardsquared_{k}")) %>% 
          mutate_all(~labelled_spss(., labels = labels, label = label))
      }
      extract_nmf <- function(nmf, method) {
        labels <- 1:k %>% set_names(paste("Segment", 1:k))
        label <- str_glue("Segment NMF - {method}: {k} cluster solution")
        NMF::basis(nmf) %>% 
          max.col() %>% 
          tibble() %>% 
          set_names(str_glue("nmf{method}_{k}")) %>% 
          mutate_all(~labelled_spss(., labels = labels, label = label))
      }
      nmflee <- function(data, k) {
        NMF::nmf(data, rank = k, method = "lee", nrun = 30) %>% 
          extract_nmf("lee")
      }
      nmfbrunet <- function(data, k) {
        NMF::nmf(data, rank = k, method = "brunet", nrun = 30) %>% 
          extract_nmf("brunet")
      }
      nmfns <- function(data, k) {
        NMF::nmf(data, rank = k, method = "nsNMF", nrun = 30) %>% 
          extract_nmf("nonsmooth")
      }
      nmfls <- function(data, k) {
        NMF::nmf(
          data, rank = k, method = "ls-nmf", nrun = 30, weight = rep(1, times = nrow(data))
        )  %>%
          extract_nmf("leastsquares")
      }
      nmfoffset <- function(data, k) {
        NMF::nmf(data, rank = k, method = "offset", nrun = 30) %>% 
          extract_nmf("offset")
      }
      
      result <- eval(call2(algorithm, data = d, k = k))
      result
    }, 
    .options = furrr_options(
      packages = c("dplyr", "tibble", "stringr", "stats", "archetypes", "NMF")
    )
  ) %>% 
    bind_cols(id, .) %>% 
    (function (x) {
      raw_data %>% 
        left_join(x, by = id_var)
    }) %>% 
    select(-.rowid)
  write_sav(results, save_path)
  if(has_twostep) {
    sps_file <- suppressWarnings(normalizePath(file.path("~", "segtool.sps")))
    spj_file <- suppressWarnings(normalizePath(file.path("~", "prod.spj")))
    sps <- map_chr(kmin:kmax, function(k) {
      first <- str_glue("
TWOSTEP CLUSTER 
/CONTINUOUS VARIABLES={input_vars}
/NUMCLUSTERS FIXED={k}
/HANDLENOISE 0
/MEMALLOCATE 64
/CRITERIA INITHRESHOLD(0) MXBRANCH(8) MXLEVEL(3)
/VIEWMODEL DISPLAY=NO
/SAVE VARIABLE=twostep_{k}.
VARIABLE LABELS twostep_{k} 'Segment twostep: {k} cluster solution'. \n
    ")
      second <- map_chr(1:k, ~str_glue("VALUE LABELS twostep_{k} {.x} 'Segment {.x}'. ")) %>% 
        paste(collapse = "\n ") %>% 
        str_glue()
      str_glue(paste(first, second))
    }) %>% 
      paste(collapse = "\n ") %>% 
      str_glue()
    sps <- str_glue(
      paste0("GET FILE=@datafile. \n ", sps, "\n SAVE OUTFILE=@datafile")
    )
    spj <- str_glue(
      '<?xml version="1.0" encoding="UTF-8"?><job codepageSyntaxFiles="false" print="false" syntaxErrorHandling="continue" syntaxFormat="interactive" unicode="true" xmlns="http://www.ibm.com/software/analytics/spss/xml/production" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.ibm.com/software/analytics/spss/xml/production http://www.ibm.com/software/analytics/spss/xml/production/production-1.4.xsd"><locale charset="UTF-8" country="US" language="en"/><output outputFormat="viewer" outputPath="{normalizePath("~")}\\out.spv"/><syntax syntaxPath="{normalizePath("~")}\\segtool.sps"/><symbol name="datafile" quote="true"/></job>'
    )
    write_lines(sps, sps_file)
    write_lines(spj, spj_file)
    spss <- "C:\\PROGRA~1\\IBM\\SPSS\\Statistics\\"
    if(!dir.exists(spss)) {
      stop("Could not find SPSS")
    }
    version <- list.files(spss) %>% 
      parse_number() %>% 
      max()
    command <- str_glue('{spss}\\{version}\\stats.exe {spj_file} -production silent -symbol @datafile {save_path}')
    shell(command)
    file.remove(sps_file)
    file.remove(spj_file)
  }
}
