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

# The checks of the input, which run before the long steps.


#checks of the input, done before the long steps start, so that a wrong input stops the run with a clear message

check_files_exist<-function(files){
  files<-files[!is.na(files)]
  missing_files<-files[file.access(files,0)==-1]
  if(length(missing_files)>0){
    stop("The following files don't exist:\n",paste(missing_files,collapse="\n"),"\n",call. = FALSE)
  }
}

check_n_cores<-function(n_cores){
  if(missing(n_cores) || !is.numeric(n_cores) || length(n_cores)!=1 || is.na(n_cores) || n_cores<1){
    stop("n_cores must be a number, 1 or more!",call. = FALSE)
  }
}

#the output files are written at the end of a run: their directory must exist and be writable
check_output_dir<-function(prefix){
  out_dir<-dirname(prefix)
  if(!dir.exists(out_dir) || file.access(out_dir,2)==-1){
    stop("Cannot write the output files, ",out_dir," is not a directory that exists and can be written to (",prefix,")",call. = FALSE)
  }
}
