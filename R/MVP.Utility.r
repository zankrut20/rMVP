# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#' Print MVP Banner
#'
#' Build date: Aug 30, 2017
#' Last update: Dec 12, 2018
#' 
#' @author Lilin Yin, Haohao Zhang, and Xiaolei Liu
#' 
#' @param width the width of the message
#' @param verbose whether to print detail.
#' 
#' @return version number.
#'
#' @examples
#' MVP.Version()
MVP.Version <- function(width=65, verbose=TRUE) {
    welcome <- "Welcome to MVP"
    title   <- "an R package for Memory-efficient, Visualization-enhanced and Parallel-accelerated genome-wide association study"
    authors <- c("Design & Maintain: Lilin Yin, Haohao Zhang, and Xiaolei Liu", 
                 "Contributors: Zhenshuang Tang, Jingya Xu, Dong Yin, Zhiwu Zhang, Xiaohui Yuan, Mengjin Zhu, Shuhong Zhao, Xinyun Li")
    contact <- "Mailto: xiaoleiliu@mail.hzau.edu.cn, ylilin@mail.hzau.edu.cn"
    logo_s  <- c(" __  __  __   __  ___",
                 "|  \\/  | \\ \\ / / | _ \\",
                 "| |\\/| |  \\ V /  |  _/",
                 "|_|  |_|   \\_/   |_|")

    version <- print_info(welcome = welcome, title = title, logo = logo_s, authors = authors, contact = contact, linechar = '=', width = width, verbose = verbose)
    return(invisible(version))
}


#' Print progress bar
#'
#' @param i the current loop number
#' @param n max loop number
#' @param type type1 for "for" function
#' @param symbol the symbol for the rate of progress
#' @param tmp.file the opened file of "fifo" function
#' @param symbol.head the head for the bar
#' @param symbol.tail the tail for the bar
#' @param fixed.points whether use the setted points which will be printed
#' @param points the setted points which will be printed
#' @param symbol.len the total length of progress bar
#'
#' @keywords internal
print_bar <- function(i,
                    n,
                    type = c("type1", "type3"),
                    symbol = "-",
                    tmp.file = NULL,
                    symbol.head = ">>>",
                    symbol.tail = ">" ,
                    fixed.points = TRUE,
                    points = seq(0, 100, 1),
                    symbol.len = 48,
                    verbose = TRUE
) {
    switch(
        match.arg(type), 
        "type1"={
            if(fixed.points){
                point.index <- points
                point.index <- point.index[point.index > floor(100*(i-1)/n)]
                if(floor(100*i/n) %in% point.index){
                    if(floor(100*i/n) != max(point.index)){
                        print.len <- floor(symbol.len*i/n)
                        logging.log(
                            paste("\r", 
                                  paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                                  paste(rep(" ", symbol.len-print.len), collapse=""),
                                  sprintf("%.2f%%", 100*i/n)
                                  , sep=""),
                            verbose = verbose
                        )
                    }else{
                        print.len <- floor(symbol.len*i/n)
                        logging.log(
                            paste("\r", 
                                  paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                                  sprintf("%.2f%%", 100*i/n), "\n"
                                  , sep=""),
                            verbose = verbose
                        )
                    }
                }
            }else{
                if(i < n){
                    print.len <- floor(symbol.len*i/n)
                    logging.log(
                        paste("\r", 
                              paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                              paste(rep(" ", symbol.len-print.len), collapse=""),
                              sprintf("%.2f%%", 100*i/n)
                              , sep=""),
                        verbose = verbose
                    )
                }else{
                    print.len <- floor(symbol.len*i/n)
                    logging.log(
                        paste("\r", 
                              paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                              sprintf("%.2f%%", 100*i/n), "\n"
                              , sep=""),
                        verbose = verbose
                    )
                }
            }
        },
        "type3"={
            progress <- readBin(tmp.file, "double") + 1
            writeBin(progress, tmp.file)
            print.len <- round(symbol.len * progress / n)
            if(fixed.points){
                if(progress %in% round(points * n / 100)){
                    logging.log(
                        paste("\r", 
                              paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                              paste(rep(" ", symbol.len-print.len), collapse=""),
                              sprintf("%.2f%%", progress * 100 / n)
                              , sep=""),
                        verbose = verbose
                    )
                }
            }else{
                logging.log(
                    paste("\r", 
                          paste(c(symbol.head, rep("-", print.len), symbol.tail), collapse=""), 
                          paste(rep(" ", symbol.len-print.len), collapse=""),
                          sprintf("%.2f%%", progress * 100 / n)
                          , sep=""),
                    verbose = verbose
                )
            }
        }
    )
}


print_accomplished <- function(width = 60, verbose = TRUE) {
    logging.log(make_line("MVP ACCOMPLISHED", width = width, linechar = '='), "\n", verbose = verbose)
}

#' Print R Package information, include title, short_title, logo, version, authors, contact
#'
#' Build date: Oct 22, 2018
#' Last update: Oct 22, 2018
#' 
#' @keywords internal
#' @author Haohao Zhang
#' 
#' @param welcome welcome text, for example: "Welcom to <Packagename>"
#' @param title long text to introduct package
#' @param short_title short label, top-left of logo
#' @param logo logo
#' @param version short label, bottom-right of logo
#' @param authors list of author names to display
#' @param contact email or website
#' @param width banner width
print_info <- function(welcome=NULL, title=NULL, short_title=NULL, logo=NULL, version=NULL, authors=NULL, contact=NULL, linechar = '=', width=NULL, verbose=TRUE) {
    msg <- c()
    # width
    if (is.null(width)) { width <- getOption('width') }
    # version
    if (is.null(version)) {
        if (getPackageName() == ".GlobalEnv") {
            version <- "devel"
        } else {
            version <- as.character(packageVersion(getPackageName()))
        }
    }
    # welcome
    if (is.null(welcome)) { 
        if (getPackageName() == ".GlobalEnv") {
            welcome <- ""
        } else {
            welcome <- paste0("Welcome to ", getPackageName())
        }
    }
    msg <- c(msg, make_line(welcome, linechar = linechar, width = width))
    # title
    if (!is.null(title)) {
        msg <- c(msg, rule_wrap(string = title, width = width, align = "center"))
    }
    
    # align logo
    logo_width <- max(sapply(logo, nchar))
    for (i in 1:length(logo)) {
        l <- paste0(logo[i], paste(rep(" ", logo_width - nchar(logo[i])), collapse = ""))
        l <- make_line(l, width)
        msg <- c(msg, l)
    }
    
    # paste short_title label to logo top-left
    if (!is.null(short_title)) {
        i <- length(msg) - length(logo) + 1
        msg[i] <- paste_label(msg[i], paste0(short_title), side = "left")
    }
    
    # paste version label to logo bottom-right
    msg[length(msg)] <- paste_label(msg[length(msg)], paste0("Version: ", version), side = "right")
    
    # authors
    if (!is.null(authors)) {
        msg <- c(msg, rule_wrap(string = authors, align = "left", linechar = " ", width = width))
    }
    # contact
    if (!is.null(contact)) {
        msg <- c(msg, rule_wrap(string = contact, align = "left", linechar = " ", width = width))
    }
    # bottom line
    msg <- c(msg, paste0(rep(linechar, width), collapse = ''))
    
    logging.log(msg, sep = "\n", verbose = verbose)
    
    return(version)
}

#' make line
#' 
#' Build date: Dec 12, 2018
#' Last update: Dec 12, 2018
#' 
#' @keywords internal
#' @author Haohao Zhang
make_line <- function(string, width, linechar = " ", align = "center", margin = 1) {
    string <- paste0(paste0(rep(" ", margin), collapse = ""),
                     string,
                     paste0(rep(" ", margin), collapse = ""))
    
    if (align == "center") {
        if (width > nchar(string)) {
            left_width <- (width - nchar(string)) %/% 2
            right_width <- width - nchar(string) - left_width
            string <-
                paste0(paste0(rep(linechar, left_width), collapse = ""),
                       string,
                       paste0(rep(linechar, right_width), collapse = ""))
        }
    } else if (align == "left") {
        if (width > nchar(string)) {
            string <-
                paste0(linechar,
                       string,
                       paste0(rep(linechar, width - nchar(string) - 1), collapse = ""))
        }
    }
    return(string)
}

#' wrap text to multiple line, align left, right or center.
#' 
#' Build date: Oct 22, 2018
#' Last update: Dec 12, 2018
#' by using base::strwrap.
#' 
#' @keywords internal
#' @author Haohao Zhang
rule_wrap <- function(string, width, align = "center", linechar = " ") {
    # define
    msg <- c()
    lines <- strwrap(string, width = width - 4)

    # wrap
    for (i in 1:length(lines)) {
        l <- make_line(lines[i], width = width, linechar = linechar, align = align)
        msg <- c(msg, l)
    }
    return(msg)
}

#' Paste label to a line
#' 
#' Build date: Oct 22, 2018
#' Last update: Oct 22, 2018
#' 
#' @param line long text
#' @param label short label
#' @param side "right" or "left"
#' @param margin default 2
#' 
#' @keywords internal
#' @author Haohao Zhang
paste_label <- function(line, label, side = "right", margin = 2) {
    if (side == "right") {
        end   <- nchar(line) - margin
        start <- end - (nchar(label) - 1)
    } else {
        start <- 1 + margin
        end   <- start + (nchar(label) - 1)
    }
    substr(line, start, end) <- label
    return(line)
}

#' format time
#' 
#' @param x seconds
#' 
#' @keywords internal
#' @author Lilin Yin
format_time <- function(x) {
    h <- x %/% 3600
    m <- (x %% 3600) %/% 60
    s <- ((x %% 3600) %% 60)
    index <- which(c(h, m, s) != 0)
    num <- c(h, m, s)[index]
    num <- round(num, 0)
    char <- c("h", "m", "s")[index]
    return(paste0(num, char, collapse = ""))
}


load_if_installed <- function(package) {
    if (!identical(system.file(package = package), "")) {
        do.call('library', list(package))
        return(TRUE)
    } else {
        return(FALSE) 
    }
}


.safe_solve <- function(X) {
    result <- tryCatch(solve(X), error = function(e) MASS::ginv(X))
    result
}


mkl_env <- function(exprs, threads = 1) {
    if (load_if_installed("RevoUtilsMath")) {
        math.cores <- eval(parse(text = "getMKLthreads()"))
        eval(parse(text = "setMKLthreads(threads)"))
    }else{
        math.cores <- blas_get_num_procs()
        blas_set_num_threads(threads)
    }
    result <- exprs
    if (load_if_installed("RevoUtilsMath")) {
        eval(parse(text = "setMKLthreads(math.cores)"))
    }else{
        blas_set_num_threads(math.cores)
    }
    return(result)
}


remove_bigmatrix <- function(x, desc_suffix=".geno.desc", bin_suffix=".geno.bin") {
    name <- basename(x)
    path <- dirname(x)
    
    descfile <- paste0(x, desc_suffix)
    binfile  <- paste0(x, bin_suffix)
    
    remove_var <- function(binfile, envir) {
        for (v in ls(envir = envir)) {
            if (is(get(v, envir = envir), "big.matrix")) {
                desc <- describe(get(v, envir = envir))@description
                if (desc$filename == binfile) {
                    rm(list = v, envir = envir)
                    gc()
                }
            }
        }
    }
    
    # remove_var(binfile, globalenv())
    remove_var(binfile, as.environment(-1L))
    
    if (file.exists(descfile)) {
        file.remove(descfile)
    }
    if (file.exists(binfile)) {
        file.remove(binfile)
    }
}


# ---------------------------------------------------------------------------
# Phenotype utility functions for logistic regression support
# Added in Phase 1 of the logistic-regression feature branch.
# ---------------------------------------------------------------------------

#' Detect phenotype family (binary vs. continuous)
#'
#' Inspects the unique non-missing values of a numeric phenotype vector and
#' returns \code{"binomial"} when exactly two values are found and they belong
#' to one of the recognised binary codings (\{0,1\}, \{1,2\}, or \{-1,1\}).
#' Otherwise \code{"gaussian"} is returned.
#'
#' Build date: 2026-06-17
#'
#' @param y Numeric vector of phenotype values (NAs allowed).
#' @param verbose Logical. Print detection result when \code{TRUE} (default).
#'
#' @return Character scalar: \code{"binomial"} or \code{"gaussian"}.
#'
#' @examples
#' detect_family(c(0, 1, 1, 0, NA))   # "binomial"
#' detect_family(c(1, 2, 2, 1))        # "binomial"
#' detect_family(c(-1, 1, -1, 1))      # "binomial"
#' detect_family(c(1.2, 3.5, 2.8))     # "gaussian"
#'
#' @export
detect_family <- function(y, verbose = TRUE) {
    uvals <- sort(unique(y[!is.na(y)]))
    is_binary <- (length(uvals) == 2L) &&
        (identical(uvals, c(0, 1))  ||
         identical(uvals, c(1, 2))  ||
         identical(uvals, c(-1, 1)))
    family <- if (is_binary) "binomial" else "gaussian"
    if (verbose) {
        logging.log(
            sprintf("Phenotype detected as '%s' (unique values: %s)\n",
                    family, paste(uvals, collapse = ", "))
        )
    }
    return(family)
}


#' Recode phenotype to standard \{0, 1\} binary coding
#'
#' Converts common binary phenotype codings to the \{0,1\} standard required by
#' logistic regression:
#' \itemize{
#'   \item \{0,1\} - returned unchanged.
#'   \item \{1,2\} - recoded as \code{y - 1}.
#'   \item \{-1,1\} - recoded as \code{(y + 1) / 2}.
#' }
#' The original coding is stored in the \code{"original_coding"} attribute of
#' the returned vector.  NAs are preserved in their original positions.
#'
#' Build date: 2026-06-17
#'
#' @param y Numeric vector of phenotype values.
#' @param from_coding Optional character string specifying the input coding
#'   (\code{"01"}, \code{"12"}, or \code{"-11"}).  When \code{NULL} (default)
#'   the coding is detected automatically.
#' @param verbose Logical. Print recoding action when \code{TRUE} (default).
#'
#' @return Numeric vector coded as \{0,1\} with attribute
#'   \code{"original_coding"} recording the detected input scheme.
#'
#' @examples
#' recode_phenotype(c(1, 2, 1, 2, NA))   # -> c(0, 1, 0, 1, NA)
#' recode_phenotype(c(-1, 1, -1, 1))      # -> c(0, 1, 0, 1)
#' recode_phenotype(c(0, 1, 1, 0))        # unchanged
#'
#' @export
recode_phenotype <- function(y, from_coding = NULL, verbose = TRUE) {
    uvals <- sort(unique(y[!is.na(y)]))

    # Determine coding scheme
    if (is.null(from_coding)) {
        if (identical(uvals, c(0, 1)))  from_coding <- "01"
        else if (identical(uvals, c(1, 2)))  from_coding <- "12"
        else if (identical(uvals, c(-1, 1))) from_coding <- "-11"
        else stop("recode_phenotype: unrecognised binary coding. ",
                  "Unique values found: ", paste(uvals, collapse = ", "))
    }

    y_new <- switch(from_coding,
        "01"  = y,
        "12"  = y - 1,
        "-11" = (y + 1) / 2,
        stop("recode_phenotype: 'from_coding' must be one of '01', '12', '-11'.")
    )

    attr(y_new, "original_coding") <- from_coding

    if (verbose) {
        action <- switch(from_coding,
            "01"  = "No recoding needed ({0,1} already standard)",
            "12"  = "Recoded {1,2} -> {0,1} (subtracted 1)",
            "-11" = "Recoded {-1,1} -> {0,1} ((y+1)/2)"
        )
        logging.log(action, "\n")
    }

    return(y_new)
}


#' Validate binary phenotype for logistic regression
#'
#' Checks that the phenotype column of an \eqn{n \times 2}{n x 2} matrix
#' (ID, Y) is suitable for binary analysis:
#' \enumerate{
#'   \item No NA values.
#'   \item Exactly two unique values.
#'   \item Minimum number of cases and controls.
#' }
#' A warning is emitted (but \code{TRUE} is still returned) when the case or
#' control count falls below the respective minimum.  \code{FALSE} is returned
#' when a hard requirement fails.
#'
#' Build date: 2026-06-17
#'
#' @param phe An \eqn{n \times 2}{n x 2} matrix or data frame with columns
#'   \code{[ID, Y]}.  Column 2 must be numeric.
#' @param min_cases Integer. Minimum number of cases (Y == 1 after recoding)
#'   required.  Default \code{10}.
#' @param min_controls Integer. Minimum number of controls (Y == 0 after
#'   recoding) required.  Default \code{10}.
#' @param verbose Logical. Print validation messages when \code{TRUE} (default).
#'
#' @return Logical scalar: \code{TRUE} if the phenotype passes all hard
#'   requirements, \code{FALSE} otherwise.
#'
#' @examples
#' phe <- data.frame(ID = 1:20, Y = c(rep(0, 10), rep(1, 10)))
#' validate_binary_phenotype(phe)  # TRUE
#'
#' @export
validate_binary_phenotype <- function(phe, min_cases = 10, min_controls = 10,
                                      verbose = TRUE) {
    y <- as.numeric(phe[, 2])

    # Hard check 1: no NAs
    if (anyNA(y)) {
        logging.log("validate_binary_phenotype: NA values found in phenotype.\n")
        return(FALSE)
    }

    uvals <- sort(unique(y))

    # Hard check 2: exactly 2 unique values
    if (length(uvals) != 2L) {
        logging.log(
            sprintf("validate_binary_phenotype: expected 2 unique values, found %d (%s).\n",
                    length(uvals), paste(uvals, collapse = ", "))
        )
        return(FALSE)
    }

    # Recode to {0,1} for counting
    y_01 <- tryCatch(
        recode_phenotype(y, verbose = FALSE),
        error = function(e) {
            logging.log("validate_binary_phenotype: unrecognised coding - ", conditionMessage(e), "\n")
            return(NULL)
        }
    )
    if (is.null(y_01)) return(FALSE)

    n_cases    <- sum(y_01 == 1L, na.rm = TRUE)
    n_controls <- sum(y_01 == 0L, na.rm = TRUE)

    if (verbose) {
        logging.log(
            sprintf("Binary phenotype: %d cases, %d controls.\n",
                    n_cases, n_controls)
        )
    }

    # Soft checks: warn but don't fail
    if (n_cases < min_cases) {
        warning(sprintf(
            "validate_binary_phenotype: only %d cases (minimum recommended: %d).",
            n_cases, min_cases
        ))
    }
    if (n_controls < min_controls) {
        warning(sprintf(
            "validate_binary_phenotype: only %d controls (minimum recommended: %d).",
            n_controls, min_controls
        ))
    }

    return(TRUE)
}
