library(optparse)
library(tidyverse)
set.seed(1234)
option_list = list(
  make_option(c("-i", "--input"),
              type="character",
              default=NULL,
              help="path to data folder (e.g. xenium output folder: /diskmnt/primary/Xenium/data/20240626__183233__20240626_5k_PDAC-bone-PKD/output-XETG00122__0037704__HT910P1-S1H1Fp1U1__20240626__183338/)",
              metavar="character")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
if (is.null(opt$input)){
  print_help(opt_parser)
  stop("Path to data is required (--input).n", call.=FALSE)
}
input = opt$input
PATH = paste0(input,"/", "transcripts.parquet")
OUTPUT <- gsub('\\.parquet$', '.csv', PATH)
if (!(file.exists(paste0(OUTPUT,".gz")))) {
  library(arrow)
  CHUNK_SIZE <- 1e6
  parquet_file <- arrow::read_parquet(PATH, as_data_frame = FALSE)
  start <- 0
  while(start < parquet_file$num_rows) {
    end <- min(start + CHUNK_SIZE, parquet_file$num_rows)
    chunk <- as.data.frame(parquet_file$Slice(start, end - start))
    data.table::fwrite(chunk, OUTPUT, append = start != 0)
    start <- end
  }
  if(require('R.utils', quietly = TRUE)) {
    R.utils::gzip(OUTPUT)
  }
}
