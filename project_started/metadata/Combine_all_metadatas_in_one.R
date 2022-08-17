# combine Sennet several metadata files into one
suppressMessages(library(dplyr))
suppressMessages(library(tibble))
suppressMessages(library(data.table))
suppressMessages(library(stringr))
suppressMessages(library(optparse))
suppressMessages(library(googlesheets4))


option_list = list(
  make_option(c("-i", "--input.folder"),
              type="character",
              default="/diskmnt/Projects/SenNet_analysis/Main.analysis/cell_typing/all_metadata", 
              help="output folder path",
              metavar="character"),
  make_option(c("-v", "--version"),
              type="character",
              default="v1.0", 
              help="version of the metadata",
              metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

input.folder <- opt$input.folder
ver <- opt$version

meta.files <- list.files(path = input.folder, pattern = '*data', full.names = T)
meta.files <- meta.files[!grepl('All_', meta.files)]

# this script will take metadata tables from regular ATAC and combo objects and create a single metadats file
gs4_deauth()
samples <- read_sheet("https://docs.google.com/spreadsheets/d/1VeWme__vvVHAhHaQB3wCvAGuq-w3WrhZ5cT-7Mh7Sr0/edit#gid=0", 
                      sheet = "Patient single cells data", trim_ws = T) %>%
  dplyr::filter(`Include in downstream` == 'Yes') %>%
  dplyr::select(`Patient ID`,`Sample name`, `Data type`, Age, Tissue) %>%
  as.data.frame() %>% rename(Sample = `Sample name`,
                             Patient_ID = `Patient ID`,
                             Data_type = `Data type`) %>%
  mutate(Data_type = str_replace(string = Data_type, pattern = ' RNA', replacement = ''),
         Data_type = str_replace(string = Data_type, pattern = ' ATAC', replacement = '')) %>% 
  distinct() %>%
  dplyr::filter(Data_type != "5' TCR-seq")
head(samples)
print(colnames(samples))

samples.id <- samples$Sample %>% as.character()
length(samples.id)

total.metadata <- lapply(meta.files, function (x) {
  fread(x, header = TRUE)})

total.colnames <- lapply(meta.files, function (x) {
  colnames(fread(x, header = TRUE))})

common.columns <- Reduce('intersect', total.colnames)
print(common.columns)

total.metadata <- lapply(total.metadata, FUN = function(x) {
  x[,common.columns, with = F]
})

total.metadata <- rbindlist(total.metadata)
total.metadata$Sample <- total.metadata$orig.ident

print(total.metadata$orig.ident %>% unique)
total.metadata <- merge(total.metadata, samples, by = 'Sample', all.x = TRUE) %>% 
  data.frame(row.names = paste(.$Sample, .$V1, sep = '_')) %>% dplyr::select(-V1)

fwrite(total.metadata, paste0('All_',length(samples.id),'_samples_metadata_data_freeze_', ver,'.tsv'), sep = '\t', row.names = TRUE)
