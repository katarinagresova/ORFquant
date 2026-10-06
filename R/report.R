# ORFquant, a software to perform splice-aware 
# quantification of ORF translation using Ribo-seq data
#
# Authors: 
# Lorenzo Calviello (calviello.l.bio@gmail.com)
# Uwe Ohler (Uwe.Ohler@mdc-berlin.de)
#
# This software is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software. If not, see
# <http://www.gnu.org/licenses/>.

# The HTML report of the results.


#' Create an html report summarizing ORFquant results
#'
#' This function creates an html report showing summary statistics for ORFquant-detected ORFs.
#' 
#' @keywords ORFquant
#' @author Lorenzo Calviello, \email{calviello.bio@@gmail.com}
#' 
#' @param input_files Character vector with full paths to plot files (*ORFquant_plots_RData) 
#' generated with \code{plot_ORFquant_results}. 
#' Must be of same length as \code{input_sample_names}.
#' 
#' @param input_sample_names Character vector containing input names.
#' Must be of same length as \code{input_files}.
#' 
#' @param output_file String; full path to html report file.
#' @return The function saves the html report file with the file path \code{output_file}.
#' @details This function creates the html report visualizing final ORFquant results. \cr \cr
#' Input are two lists of the same length: \cr \cr
#' a) \code{input_files}: list of full paths to one or multiple input files 
#' (*ORFquant_plots_RData files generated with \code{plot_ORFquant_results}) and \cr \cr
#' b) \code{input_sample_names}: list of corresponding names describing the file content (these are used as names in the report). \cr \cr
#' For the report, a RMarkdown file is rendered as html document, saved as \code{output_file}. \cr \cr
#' @seealso \code{\link{plot_ORFquant_results}}, \code{\link{run_ORFquant}}
#' @export

create_ORFquant_html_report <- function(input_files, input_sample_names, output_file){
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Package \"rmarkdown\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  # get input and output file paths
  input_files <- paste(normalizePath(dirname(input_files)),basename(input_files),sep="/")
  output_file <- paste(normalizePath(dirname(output_file)),basename(output_file),sep="/")
  
  # get path to RMarkdown file (to be rendered)
  rmd_path <- paste(system.file(package="ORFquant"),"/rmd/ORFquant_template.Rmd",sep="")
  
  # render a copy in a new temporary directory, so that nothing is written into the installed package
  rmd_dir <- tempfile("ORFquant_report_")
  dir.create(rmd_dir)
  on.exit(unlink(rmd_dir, recursive = TRUE), add = TRUE)
  file.copy(rmd_path, rmd_dir)
  rmd_path <- paste(rmd_dir,"ORFquant_template.Rmd",sep="/")
  
  sink(file = paste(output_file,"_ORFquant_report_output.txt",sep = ""))
  on.exit(sink(), add = TRUE)
  # render RMarkdown file > html report
  suppressWarnings(rmarkdown::render(rmd_path, 
                          params = list(input_files = input_files,
                                        input_sample_names = input_sample_names),
                          output_file = output_file))
  invisible()
}
