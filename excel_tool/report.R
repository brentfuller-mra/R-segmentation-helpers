if(!"pacman" %in% rownames(installed.packages())) {
  install.packages("pacman", repos = "https://cran.rstudio.com/")
}
pacman::p_load(
  dplyr, tibble, stringr, haven, assertthat, purrr, furrr, tidyr, rlang, readr,
  update = TRUE
)
options(dplyr.summarise.inform = FALSE)
#' Function for creating segmentation report from Excel tool
#' @param data_path String with path to SPSS data.  Should contain segment 
#' variables
#' @param report_save Path to save Excel report
#' @param seg_autofine Boolean; should segments be automatically detected in the
#' data? Finds segments with spss label starting with "Segment "
#' @param other_seg_vars String with space delimited variable names indicating
#' segment assignments, group delimited by ';'
#' @param weight_var If there is a weight, String indicating weight variable
#' @param profile Boolean; Whether to include profile sheet in report
#' @param profile_vars Space delimited string with variable names to use for 
#' profiling, group delimited with ';'.
#' @param migration_tables Boolean; Whether to include migration tables in 
#' report
#' @param inputs Boolean; Whether to include input data set in report
#' @param input_vars Space delimited string with input variables
#' @param id_var String indicating id variable, needed for input vars and 
#' assignments
#' @param assignments Boolean; Should segment assignments be included in the 
#' report

report <- function(data_path, report_save, seg_autofind = TRUE, 
                   other_seg_vars = NULL, weight_var = NULL, profile = FALSE, 
                   profile_vars = NULL, migration_tables = FALSE, inputs = TRUE, 
                   input_vars = NULL, id_var = NULL, assignments = FALSE) {
  
  wb <- createWorkbook()
  sav_data <- read_sav(data_path) %>% rowid_to_column(".rowid")
  data_labels <- map_dfc(sav_data, function(.x) {
    ifelse(
      test = is.null(attributes(.x)$label), 
      yes = "", 
      no = as.character(attributes(.x)$label[1])
    )
  }) %>% 
    pivot_longer(everything(), names_to = "vars", values_to = "label")
  other_seg_vars <- str_split(other_seg_vars, ";") %>% 
    purrr::list_c() %>% 
    trimws() %>% 
    map(~purrr::list_c(str_split(.x, " ")))
  if(seg_autofind) {
    seg_var_list <- data_labels %>% 
      filter(str_detect(label, "^Segment ")) %>% 
      mutate(group = str_extract(label, "(?<=^Segment ).*(?=\\:)")) %>% 
      group_by(group) %>% 
      group_split() %>% 
      map(~pluck(.x, "vars")) %>% 
      append(other_seg_vars)
  } else {
    seg_vars_list <- other_seg_vars
  }
  seg_vars_flat <- purrr::list_c(seg_var_list)
  if(is.null(weight_var)){
    weights <- rep(1, times = nrow(sav_data))
  } else{
    weights <- sav_data[[weight_var]]
  }
  assert_that(length(seg_vars_flat) > 0)
  mborder_right <- createStyle(
    border = "right", 
    borderStyle = "medium"
  )
  mborder_left <- createStyle(
    border = "left", 
    borderStyle = "medium"
  )
  mborder_top <- createStyle(
    border = "top", 
    borderStyle = "medium"
  )
  mborder_bottom <- createStyle(
    border = "bottom", 
    borderStyle = "medium"
  )
  tborder_right <- createStyle(
    border = "right", 
    borderStyle = "thin"
  )
  
  mborder_around <- function(wb, sheet, first_row, last_row, first_col, 
                             last_col) {
    addStyle(wb, sheet, mborder_left, first_row:last_row, first_col, 
             stack = TRUE)
    addStyle(wb, sheet, mborder_right, first_row:last_row, last_col, 
             stack = TRUE)
    addStyle(wb, sheet, mborder_top, first_row, first_col:last_col, 
             stack = TRUE)
    addStyle(wb, sheet, mborder_bottom, last_row, first_col:last_col, 
             stack = TRUE)
  }
  vh_center <- createStyle(halign = "center", valign = "center")
  h_center <- createStyle(halign = "center")
  twodec <- createStyle(numFmt = "#,#0.00")
  integer <- createStyle(numFmt = "#,###")
  prop <- createStyle(numFmt = "0.0%")
  yellow_fill <- createStyle(fgFill = "yellow")
  black_fill <- createStyle(fgFill = "black")
  bold <- createStyle(textDecoration = "bold")
  full_positive <- createStyle(
    bgFill = rgb(146, 208, 80, maxColorValue = 255)
  )
  half_positive <- createStyle(
    bgFill = rgb(214, 237, 189, maxColorValue = 255)
  )
  half_negative <- createStyle(
    bgFill = rgb(230, 185, 184, maxColorValue = 255)
  )
  full_negative <- createStyle(
    bgFill = rgb(194, 80, 77, maxColorValue = 255)
  )
  if(inputs) {
    input_vars_vec <- str_split(c(id_var, input_vars), " ") %>% 
      purrr::list_c() %>% 
      trimws()
    addWorksheet(wb, "Inputs")
    input_data <- select(sav_data, all_of(input_vars_vec))
    writeData(wb, "Inputs", input_data, headerStyle = bold)
  }
  if(assignments) {
    addWorksheet(wb, "Segment Assignments")
    seg_data <- select(sav_data, all_of(c(id_var, seg_vars_flat)))
    writeData(wb, "Segment Assignments", seg_data, headerStyle = bold)
  }
  if(profile) {
    assert_that(!is.null(profile_vars))
    profile_list <- str_split(profile_vars, ";") %>% 
      purrr::list_c() %>% 
      trimws() %>% 
      map(~purrr::list_c(str_split(.x, " ")))
    addWorksheet(wb, "Profile")
    writeData(wb, "Profile", 0.3, colNames = FALSE)
    addStyle(wb, "Profile", yellow_fill, 1, 1, stack = TRUE)
    mborder_around(wb, "Profile", 1,1,1,1)
    writeData(wb, "Profile", "<-- Formatting threshold", 2, 1, colNames = F)
    mborder_around(wb, "Profile", 1,1,2,2)
    writeData(wb, "Profile", "Variable Info", 1, 3, colNames = F)
    addStyle(wb, "Profile", bold, 3, 1, stack = TRUE)
    addStyle(wb, "Profile", vh_center, 3, 1, stack = TRUE)
    mergeCells(wb, "Profile", 1:4, 3:5)
    mborder_around(wb, "Profile", 3, 5, 1, 4)
    writeData(wb, "Profile", "Variable", 1, 6)
    writeData(wb, "Profile", "Label", 2, 6)
    writeData(wb, "Profile", "N-Size", 3, 6)
    writeData(wb, "Profile", "Total", 4, 6)
    mborder_around(wb, "Profile", 6, 6, 1, 4)
    addStyle(wb, "Profile", bold, 6, 1:4, stack = TRUE)
    addStyle(wb, "Profile", h_center, 6, 1:4, stack = TRUE)
    setColWidths(wb, "Profile", 2, 90)
    letters_xl <- expand_grid(
      a = c("", LETTERS), b = c("", LETTERS), c = c("", LETTERS)
    ) %>% 
      unite(cmb, c(a,b,c), sep = "") %>% 
      filter(cmb != "") %>% 
      filter(!duplicated(.)) %>% 
      pluck(1)
    
    g <- 0
    chunks <- map_dfr(profile_list, function(.x) {
      g <<- g + 1
      select(sav_data, all_of(.x), all_of(seg_vars_flat)) %>% 
        bind_cols(tibble(.weights = weights)) %>% 
        mutate_all(as.numeric) %>% 
        pivot_longer(cols = all_of(.x), names_to = "profvar") %>% 
        add_column(group = g)
    }) %>% 
      pivot_longer(
        cols = c(-.weights, -profvar, -value, -group), names_to = "seg", 
        values_to = "assignment"
      ) %>% 
      group_by(profvar, group, seg, assignment) %>% 
      summarize(cell = weighted.mean(value, .weights, na.rm = T)) %>% 
      mutate(cell = ifelse(is.nan(cell), NA, cell)) %>% 
      unite(seg, c(seg, assignment)) %>% 
      pivot_wider(names_from = seg, values_from = cell) %>% 
      mutate(rank = which(profvar == purrr::list_c(profile_list))) %>% 
      arrange(group, rank) %>% 
      group_by(group) %>% 
      group_split() %>% 
      map(function(rowchunk) {
        rowchunk <- rowchunk %>% select(-group) 
        map(c(seg_vars_flat), function(.x) {
          select(rowchunk, profvar, matches(str_glue("^{.x}_")))
        })
      })
    use_prop <- purrr::list_c(profile_list) %>% 
      tibble(vars = .) %>% 
      mutate(
        prop = map_lgl(vars, function(.x) {
          all(sav_data[[.x]] %in% c(0, 1, NA), na.rm = TRUE)
        })
      ) %>% 
      rowid_to_column(".rowid")
    info <- map(profile_list, function(.x) {
      filter(data_labels, vars %in% .x) %>% 
        mutate(
          N = map_dbl(vars, function(vars) {
            sum(weights[!is.na(sav_data[[vars]])])
          }), 
          Total = map_dbl(vars, function(vars){
            weighted.mean(sav_data[[vars]], weights, na.rm = TRUE)
          }), 
          order = map_dbl(vars, function(.y) which(.x %in% .y))
        ) %>% 
        arrange(order) %>% 
        select(-order)
    })
    
    headers <- map(seg_vars_flat, function(s) {
      seg <- sav_data[[s]]
      map(1:max(seg), function(.x) {
        N <- sum(weights[seg == .x])
        perc <- N/sum(weights)
        rbind(N, perc)
      }) %>% 
        reduce(cbind) %>% 
        `colnames<-`(1:ncol(.))
    }) %>% 
      set_names(seg_vars_flat)
    
    .row <- 7
    pwalk(list(chunks, info), function(chunks, info) {
      writeData(
        wb, "Profile", info, startRow = .row, colNames = FALSE, rowNames = FALSE
      )
      prop_rows <- use_prop %>% 
        filter(vars %in% info$vars, prop) %>% 
        pluck(".rowid") %>% 
        `+`(6)
      num_rows <- use_prop %>% 
        filter(vars %in% info$vars, !prop) %>% 
        pluck(".rowid") %>% 
        `+`(6)
      if(length(prop_rows) > 0) {
        addStyle(wb, "Profile", prop, prop_rows, 4, stack = TRUE)
      }
      if(length(num_rows) > 0) {
        addStyle(wb, "Profile", twodec, num_rows, 4, stack = TRUE)
      }
      
      mborder_around(wb, "Profile", .row, nrow(info) + .row - 1, 1, 4)
      addStyle(
        wb, "Profile", h_center, .row:(nrow(info) + .row - 1), 3:4, 
        stack = TRUE, gridExpand = TRUE
      )
      addStyle(
        wb, "Profile", integer, .row:(nrow(info) + .row - 1), 3, stack = TRUE
      )
      .col <- 6
      i <<- 1
      walk(chunks, function(chunk) {
        chunk <- select(chunk, -1)
        if(.row == 7) {
          writeData(wb, "Profile", names(headers)[i], .col, 3)
          writeData(
            wb, "Profile", headers[[i]], .col, 4, colNames = FALSE, 
            rowNames = FALSE
          )
          writeData(
            wb, "Profile", t(matrix(colnames(headers[[i]]))), .col, 6, 
            colNames = FALSE, rowNames = FALSE
          )
          mergeCells(wb, "Profile", .col:(.col + ncol(headers[[i]]) - 1), 3)
          mborder_around(
            wb, "Profile", 3, 3, .col, (.col + ncol(headers[[i]]) - 1)
          )
          mborder_around(
            wb, "Profile", 3, 6, .col, (.col + ncol(headers[[i]]) - 1)
          )
          mborder_around(
            wb, "Profile", 6, 6, .col, (.col + ncol(headers[[i]]) - 1)
          )
          addStyle(
            wb, "Profile", bold, c(3, 6), .col:(.col + ncol(headers[[i]]) - 1), 
            stack = TRUE, gridExpand = TRUE
          )
          addStyle(
            wb, "Profile", integer, 4, .col:(.col + ncol(headers[[i]]) - 1), 
            stack = TRUE, gridExpand = TRUE
          )
          addStyle(
            wb, "Profile", integer, 4, .col:(.col + ncol(headers[[i]]) - 1), 
            stack = TRUE, gridExpand = TRUE
          )
          addStyle(
            wb, "Profile", prop, 5, .col:(.col + ncol(headers[[i]]) - 1), 
            stack = TRUE, gridExpand = TRUE
          )
          addStyle(
            wb, "Profile", h_center, 5, .col:(.col + ncol(headers[[i]]) - 1), 
            stack = TRUE, gridExpand = TRUE
          )
          addStyle(
            wb, "Profile", tborder_right, 4:6, 
            .col:(.col + ncol(headers[[i]]) - 2), stack = TRUE, gridExpand = TRUE
          )
          addStyle(
            wb, "Profile", h_center, 3:6, 
            .col:(.col + ncol(headers[[i]]) - 1), stack = TRUE, gridExpand = TRUE
          )
          setColWidths(wb, "Profile", .col - 1, 2)
          i <<- i + 1
        }
        writeData(wb, "Profile", chunk, .col, .row, colNames = FALSE)
        if(length(prop_rows) > 0) {
          addStyle(
            wb, "Profile", prop, prop_rows, 
            .col:(.col + ncol(chunk) - 1), stack = TRUE, gridExpand = TRUE
          )
        }
        if(length(num_rows) > 0) {
          addStyle(
            wb, "Profile", twodec, num_rows, 
            .col:(.col + ncol(chunk) - 1), stack = TRUE, gridExpand = TRUE
          )
        }
        addStyle(
          wb, "Profile", tborder_right, .row:(nrow(chunk) + .row - 1), 
          .col:(ncol(chunk) - 1 + .col), gridExpand = TRUE, stack = TRUE
        )
        conditionalFormatting(
          wb, "Profile", .col:(ncol(chunk) - 1 + .col), 
          .row:(nrow(chunk) + .row - 1), 
          str_glue("({letters_xl[.col]}{.row}/$D{.row} - 1) <= $A$1*-1/2"), 
          half_negative
        )
        conditionalFormatting(
          wb, "Profile", .col:(ncol(chunk) - 1 + .col), 
          .row:(nrow(chunk) + .row - 1), 
          str_glue("({letters_xl[.col]}{.row}/$D{.row} - 1) >= $A$1*1/2"), 
          half_positive
        )
        conditionalFormatting(
          wb, "Profile", .col:(ncol(chunk) - 1 + .col), 
          .row:(nrow(chunk) + .row - 1), 
          str_glue("({letters_xl[.col]}{.row}/$D{.row} - 1) <= $A$1*-1"), 
          full_negative
        )
        conditionalFormatting(
          wb, "Profile", .col:(ncol(chunk) - 1 + .col), 
          .row:(nrow(chunk) + .row - 1), 
          str_glue("({letters_xl[.col]}{.row}/$D{.row} - 1) >= $A$1"), 
          full_positive
        )
        addStyle(
          wb, "Profile", h_center, .row:(nrow(chunk) + .row - 1), 
          .col:(ncol(chunk) - 1 + .col), gridExpand = TRUE, stack = TRUE
        )
        mborder_around(
          wb, "Profile", .row, .row + nrow(chunk) - 1, .col, .col + ncol(chunk) - 1
        )
        .col <<- .col + ncol(chunk) + 1
      })
      .row <<- .row + nrow(info)
    })
  }
  if(migration_tables){
    addWorksheet(wb, "Migration Tables")
    .row <- 1
    walk(seg_var_list, function(sg) {
      walk2(sg[-length(sg)], sg[-1], function(.x, .y) {
        tbl <- wtd.table(as.numeric(sav_data[[.x]]), as.numeric(sav_data[[.y]]), weights) %>% 
          as.matrix()
        writeData(
          wb, "Migration Tables", str_glue("{.x} vs {.y}"), 1, .row, 
          colNames = FALSE, rowNames = FALSE
        )
        mergeCells(wb, "Migration Tables", 1:(ncol(tbl) + 1), .row)
        addStyle(
          wb, "Migration Tables", h_center, .row, 1:(ncol(tbl) + 1), 
          stack = TRUE
        )
        addStyle(
          wb, "Migration Tables", bold, .row:(.row + 1), 1:(ncol(tbl) + 1), 
          stack = TRUE, gridExpand = TRUE
        )
        mborder_around(wb, "Migration Tables", .row, .row, 1, ncol(tbl) + 1)
        .row <<- .row + 1 
        mborder_around(wb, "Migration Tables", .row, .row, 1, ncol(tbl) + 1)
        addStyle(
          wb, "Migration Tables", black_fill, .row, 1, gridExpand = TRUE, 
          stack = TRUE
        )
        writeData(
          wb, "Migration Tables", tbl, 1, .row, colNames = TRUE, rowNames = TRUE
        )
        addStyle(
          wb, "Migration Tables", bold, .row:(.row + nrow(tbl) + 1), 1, 
          gridExpand = TRUE, stack = TRUE
        )
        .row <<- .row + 1
        mborder_around(
          wb, "Migration Tables", .row, .row + nrow(tbl) - 1, 1, ncol(tbl) + 1
        )
        walk(.row:(nrow(tbl) + .row), function(.r) {
          conditionalFormatting(
            wb, "Migration Tables", 2:(ncol(tbl) + 1), 
            .r, style = c("white", "#63BE7B"), 
            type = "colourScale"
          )
          .row <<- .row + 1
        })
      })
    })
  }
  saveWorkbook(wb, report_save, overwrite = TRUE)
  openXL(report_save)
}
