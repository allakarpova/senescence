# written by Austin Southard-Smith
# some modification/feedback provided by Evan Chien-Wei Peng in v6
# 2025-02-26
# based on https://satijalab.org/seurat/articles/spatial_vignette_2#mouse-brain-10x-genomics-xenium-in-situ
# addresses error mentioned https://github.com/satijalab/seurat/issues/7491
# update 2024-07-01: added in v5 of the script on 2024-07-01 to generate the `transcripts.csv.gz`file that is missing from the output of the Xenium onboard Analysis v3.0 software. I assume that Seurat will be made compatible with this in the future but right now they have not added any sort of translation functionality.
# update 2024-12-31: added in v6 of the script:
        # Output plots are now rasterized to (300 dpi for dotplots and 600 dpi for spatial image plots). This decreases the file size for plots of larger datasets
        # Added the Xenium.with.snvs.microbiome and Xenium.with.microbiome assays and their respective fovs. With v1 of the Xenium 5k panel some probes were generated to target transcripts expressed by the microbiome. Like the competing Variant/Reference probes targeting snvs and indels these are now removed from the Default Xenium assay present in the object and added as part of a separate assay so they are not normalized with the rest of the expression data during processing and analysis.
# update 2025-02-26: added in v8 of the script:
        # prior version the output objects were a mix of seurat v3 and seurat v5 assays. By default the script now only creates seurat v5 assays. I also added an option that when specified results in only v3 assays. To support legacy code/projects the flag also supports a "legacy" option which will keep the object assays as they were originally created. In this case the first assay that is created is a seurat v5 asssay (this varies depending on the commands used in v6) while all other assays were seurat v3 assays.
        # added option to provide a whitelist of barcodes when creating the object. This option accepts a file that is just a list of barcodes. If this option is specified and a file is provided, then script will remove all cells not present in that list. This is being added because in some cases different tissue sections are placed too close together on the Xenium slide. The result is that when fovs are selected for different tissue pieces one of the two output datasets for those two samples will have a small amount of tissue from a different sample in it. In order to maintain project and sample separation and not contaminate analysis with tissue from different cases it is necessary to remove the cells that belong to one sample from the other sample's dataset.
        # added option to crop the FOV cell segmentations and molecules using a custom convex hull polygon. This option is intended to be used alongside the barcode subsetting option. If this option is specified and a file is provided, then the script will remove all transcript spots from all of the fovs in the object that fall outside of the polygon. This option takes a table of X, Y coordinates. First column must be the X coordinate. Second column must be the Y coordinate. It is compatible with the X,Y selection coordinates file that is exported from XeniumExplorer.
        # S Q U A R E S   Y O U R   P L O T: by order of the S^2 all DimPlots and FeaturePlots shall now fit in the square hole. https://www.youtube.com/watch?v=cUbIkNUFs-4
        # Have you ever opened the plots after running this script and thought to yourself: "Oh god! Stop it! That's disgusting!" because your ImageDimPlots and ImageFeaturePlots have unsightly lines crossing them? Well no more!. Through countless seconds of arduous coding we have removed these unsightly blemishes creating a seemless black finish to your plots and in doing so added a bespoke custom black border as well. This is a feature nobody asked me for and I did it anyways because fuck it.
            # If you don't like this feature then replace each instance of the code in the plot with print(p*) where p* is the plot number.
                # grid:::grid.newpage()
                # grid:::grid.rect(gp=grid:::gpar(fill="#000000", col="#000000"))
                # ggplot2:::print.ggplot(p5, newpage=FALSE)
            # If you don't like this feature then you will also need to remove the parts from the theme() function call: panel.background = element_rect(fill = "#000000", color = NA), plot.background = element_rect(fill = "#000000", color = NA), panel.border = element_rect(fill = NA, color = "#000000"),  text = element_text(color = "white")
# important consideration for running SCTransform in Seurat5 on an object with multiple layers: "if you run SCTransform on your object as a multi-layer object (i.e. not with layers all joined, and not as a list of separate objects), then by default SCTransform will identify variable features that are generally shared across the different samples. Please let me know if you continue to run into any issues!" - https://github.com/satijalab/seurat/issues/9341#issuecomment-2479410705 
        # **DECIDED AGAINST. DO NOT DO THIS.** I'm debating adding an option to create a single object from multiple input datasets. In some cases due to tissue detachment a single section is split into two different datasets when Xenium data is generated. This is because we do not want to include the detached regions in the selected fovs of the dataset if the majority of the region has detached (it will result in lower quality data and problems with image alignment). Because Xenium Instrument software version 3.0.2.0 requires that a region's dataset be contiguous when it is generated on the instrument we must merge the datasets either after or during object generation. When merged this will result in a single counts matrix and additional fovs per additional dataset. These fovs will be named after the additional sample that is being added.
        # if this were done then it would be done by creating two separate objects first and then merging them together at which point the multiple SCTransform layers may come into play depending on how the objects are handled.
        # currently I am aware of only one case where this would be applicable across all of the datasets that we have generated.        
        # i've decided this is a bad idea precisely because it would result in objects that are the result of multiple samples. This would solve my current problem but it would introduce multiple downstream problems. For instance I would now have a initial object that is formatted differently that would then have to be accounted for every single time I write a script. If I had a way to register the spatial data then this could be avoided but I . When it comes to data upload/sharing at time of publication we would have to split the datasets anyways because the two input datasets are different. Keeping track of it would mean I would have to have duplicate values for input paths at some level in any project tracking sheets I used. 
# important consideration for running SCTransfrom in Seurat5 when trying to regress out multiple variables:
        # version 5.0.0, 5.0.1, and 5.0.2 have a bug in the SCTransfrom function that results in regression of multiple variables to be no different from if you did not regress out variables at all (see https://github.com/satijalab/seurat/issues/8148). The work-around was to add `ncells = ncol(seurat_object)` in the SCTransform() function call. This was supposedly fixed in Seurat version 5.0.3. 
        # The SCTransfrom function in this script does not rely on this solution because we do not need to regress out mitochondrial, ribosomal, or cell cycle related stuff in the Xenium spatial data. 
library(Seurat)
library(tidyverse)
library(future)
library(optparse)
library(magrittr)
library(ggrastr)
library(data.table)
options(future.globals.maxSize= +Inf)
set.seed(1234)
option_list = list(
    make_option(c("-i", "--input"),
                type="character",
                default=NULL,
                help="path to data folder (e.g. xenium output folder: /diskmnt/primary/Xenium/data/20230803__195918__20230803_firstsamples/output-XETG00122__0010660__HT270P1-S1H1A1Fp1Us1_1__20230803__200008/)",
                metavar="character"),
    make_option(c('-s', "--sampleID"),
                type="character",
                default=NULL,
                help="SampleID that output file names will be based on",
                metavar="character"),
    make_option(c("-o", "--output"),
                type="character",
                default="./",
                help="output folder path",
                metavar="character"),
    make_option(c("--mols.qv.threshold"),
                type="integer",
                default=20,
                help="threshold for the qv score cutoff threshold. Transcripts above this will be included.",
                metavar="integer"),
    make_option(c("--pc_num"),
                type="integer",
                default=30,
                help="number of principal components to use",
                metavar="integer"),
    make_option(c("--counts_cutoff"),
                type="integer",
                default=0,
                help="retain only cells with nCounts greater than the value specified by the counts_cutoff",
                metavar="integer"),
    make_option(c("--plot_features"),
                type='character',
                default=NA,
                help = "comma separated list of features to plot when the object is generated (for pancreas this looks like 'KRT7,AMY2A,TFF2,CFTR,EPCAM,MKI67,TOP2A,PPARG,PECAM1,VWF,MS4A1,CD3D,CD3E,CD4,MYC,INS,GCG,C7,PDGFRA,SRPX,MAMDC2,VCAN,MEDAG,MFAP5,GPRC5A,CAV1,PTPRC,IL7R,TMC5,MET,EHF,CAPN8,CD163,MS4A6A,PDGFRB,MYH11,KCNMA1'",
                metavar="character"),
    make_option(c("--with_snvs"),
                action="store_true",
                default=FALSE,
                help="Use this flag when the input xenium data also contains SNVs present. By default if this flag is used then the assay will only remove those Xenium probes with the pattern '_WT' or '_ALT' in the probe name. To specify other snv probes to remove use the `--snv_probes` flag. These snv probes need to be filtered out prior to any normalization or subsequent marker analysis. When used this flag results in the resulting object having two raw assays one is the 'Xenium' assay and the other is the 'Xenium.with.snvs' assay. The 'Xenium' assay has no 'snv' target probes in it and is used for subsequent SCT normalization and can be used for future non-mutation related analysis. There are then two whole slide FOVs generated. One is the normal 'fov' which is intended for use with the 'Xenium' and 'SCT' assay. The other is the 'fov.with.snvs' image which is intended for use with the 'Xenium.with.snvs' assay. In order to visualize the localization of the SNV probes on the section you need to use the 'fov.with.snvs' image. If this option is not specified then only the 'Xenium' assay and 'fov' image will be present in the resulting seurat object. If this flag is used with the --with_microbiome flag then an additional raw assay 'Xenium.with.snvs.microbiome' and full section FOV 'fov.with.snvs.microbiome' will be generated."),
    make_option(c("--with_microbiome"),
                action="store_true",
                default=FALSE,
                help="Use if the gene panel includes both SNVS and microbiome probes. By default if this flag is used then the assay will only remove those Xenium probes with the pattern 'Bacteroides', 'Escherichia', 'Fusobacterium', and 'papillomavirus' in the probe name. To specify other snv probes to remove use the `--microbiome_probes` flag. These snv probes need to be filtered out prior to any normalization or subsequent marker analysis. When used this flag results in the resulting object having two raw assays one is the 'Xenium' assay and the other is the 'Xenium.with.microbiome' assay. The 'Xenium' assay has no 'microbiome' target probes in it and is used for subsequent SCT normalization and can be used for future non-mutation related analysis. There are then two whole slide FOVs generated. One is the normal 'fov' which is intended for use with the 'Xenium' and 'SCT' assay. The other is the 'fov.with.microbiome' image which is intended for use with the 'Xenium.with.microbiome' assay. In order to visualize the localization of the SNV probes on the section you need to use the 'fov.with.microbiome' image. If this option is not specified then only the 'Xenium' assay and 'fov' image will be present in the resulting seurat object. If this flag is used with the --with_snvs flag then an additional raw assay 'Xenium.with.snvs.microbiome' and full section FOV 'fov.with.snvs.microbiome' will be generated."),
    make_option(c("--snv_probes"),
                type='character',
                default=NA,
                help="Only use this flag in combination with the `--with_snvs`. This flag does not do anything if it is used and `--with_snvs` is not also used. Use this flag followed by the path to a file containing a list of SNV probes that need to be removed (all probes to remove need to be present in this file) from the seurat object when the object is being generated. (e.g. '--with_snvs --snv_probes /diskmnt/Projects/Users/austins2/tools/xenium_snvs/Xenium_human_custom_panel_v2_snv_probe_names.tsv'). Any probe/gene specified in this file will be removed from the 'Xenium` assay and 'fov' image. These probes will still be retained and accessible from the 'Xenium.with.snvs' assay and 'fov.with.snvs' image."),
    make_option(c("--microbiome_probes"),
                type='character',
                default=NA,
                help="Only use this flag in combination with the `--with_snvs_microbiome`. This flag does not do anything if it is used and `--with_microbiome` is not also used. Use this flag followed by the path to a file containing a list of microbiome probes that need to be removed (all probes to remove need to be present in this file) from the seurat object when the object is being generated. (e.g. '--with_microbiome --microbiome_probes /diskmnt/Projects/Users/austins2/tools/Xenium_human_panel_5k_v1_DHYHMW_microbiome_probe_names.tsv'). Any probe/gene specified in this file will be removed from the 'Xenium` assay and 'fov' image. These probes will still be retained and accessible from the 'Xenium.with.microbiome' assay and 'fov.with.microbiome' image."),
    make_option(c("-w","--barcode_whitelist"),
                type='character',
                default=NULL,
                help="Cells that are subset out when this option is provided will be retained in the <Sample_ID>_raw.rds intermediate file. This option is used to specify the path to a file containing a list of barcodes where the first column (file can have any number of additional columns) is all cell barcodes that should be retained in the object. The file should have a column header. (e.g. barcode subsetting list: /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/individual_object/seurat_5.0.1_5k_v7/tests/HT227P1-S1H1L1U2_half_section_Xenium_Explorer_region_selection.csv)"),
    make_option(c("-c","--coordinate_whitelist"),
                type='character',
                default=NULL,
                help="Molecules that are cropped out when this option is provided will be retained in the <Sample_ID>_raw.rds intermediate file. This option is used to specify the path to a file containing a list of X and Y coordinates where the first column is the X coordinate and the second column is the Y coordinate (file can have any number of additional columns). These coordinates will be used to generate a convex hull polygon where only molecules present inside the polygon will be retained in the fov of the object. The file should have a column header. (e.g. coordinate cropping list: /diskmnt/Projects/HTAN_analysis_2/PDAC/xenium/individual_object/seurat_5.0.1_5k_v7/tests/HT227P1-S1H1L1U2_half_section_Selection_1_coordinates.csv)"),
    make_option(c("-a","--assay_version"),
                type='character',
                default="v5",
                help="Specify the version of all of the assays when creating the Seurat object. Options are: 'v3', 'v5', or 'legacy' (used like: -a v5). 'v3' results in all of the assays in the output seurat object being v3 assays. These are normal sparse matrices like those prior to the release of Seurat version 5. 'v5' results in all of the assays in the output seurat object being v5 assays. These are the updated layers that were released in Seurat v5. 'legacy' results in object with a combination of v3 and v5 assays mimicing the behavior of the script prior to the v7 of this script. When the 'legacy' option is used only the assay that is made by calling CreateSeuratObject() is a v5 assay and all others are v3 assays. The problem with using the 'legacy' option is that the output assays in the output objects may not be the same version if the script is not called in the same manner every single time (say one time you use the --with_snvs while another time you use both `--with_snvs` and `--with_microbiome`). This can complicate downstream analysis and for that reason it is recommended to use either `v5` or `v3`")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# complain if there's no data
if (is.null(opt$input)){
    print_help(opt_parser)
    stop("Path to data is required (--input).n", call.=FALSE)
}
if (is.null(opt$sampleID)){
    print_help(opt_parser)
    stop("sampleID is required (--sampleID).n", call.=FALSE)
}

input = opt$input
sampleID = opt$sampleID
out_path = paste0(opt$output,"/",sampleID,"/")
dir.create(out_path)
num_pcs = opt$pc_num
print(sampleID)
# added in v5 of the script on 2024-07-01 to generate the `transcripts.csv.gz`file that is missing from the output of the Xenium onboard Analysis v3.0 software.
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

probe_hyphens <- function(probe_list) {
    probe_list <- gsub(":", "-", probe_list)
    probe_list <- gsub("\\.", "-", probe_list)
    probe_list <- gsub("_", "-", probe_list)
    probe_list <- gsub("\\+", "-plus-", probe_list)
    return(probe_list)
}

######
# capitalize this you filthy casual
opt$assay_version <- tolower(opt$assay_version)
# defining functions used for the generation of Seurat objects in this script and setting options based on the specified version of the output object
if ((opt$assay_version == "v5" )) {
    options(Seurat.object.assay.version = "v5")
    CreateAssay <- SeuratObject::CreateAssay5Object
    object_first_assay_version = "v5"
    subsequent_assay_version = "v5"
} else if (opt$assay_version == "v3") {
    options(Seurat.object.assay.version = "v3")
    CreateAssay <- SeuratObject::CreateAssayObject
    object_first_assay_version = "v3"
    subsequent_assay_version = "v3"
} else if (opt$assay_version == "legacy") {
    options(Seurat.object.assay.version = "v5")
    CreateAssay <- SeuratObject::CreateAssayObject
    object_first_assay_version = "v5"
    subsequent_assay_version = "v3"
} else { 
    # default is v5. warn the user their input does not match one of the required options so we are defaulting to using v5 objects
    print("The option passed using '-a' or '--assay_version' does not match one of the required 3 options 'v5', 'v3', or 'legacy'. Using the default version option of 'v5'.")
    options(Seurat.object.assay.version = "v5")
    CreateAssay <- SeuratObject::CreateAssay5Object
    object_first_assay_version = "v5"
    subsequent_assay_version = "v5"
}
print(paste0("The first assay in the output seurat object will be ", object_first_assay_version))
print(paste0("All other assays in the seurat object will be ", subsequent_assay_version))
######

create_assay_and_fov <- function(counts, pixels, segmentations, list.to.remove, assay_name) {
    counts_working <- counts
    counts.no.remove <- counts_working[(!rownames(counts_working) %in% c(list.to.remove)),]
    rm(counts_working)
    assay.no.remove <- CreateAssay(counts = counts.no.remove)
    pixels_working <- pixels
    # print(dim(pixels))
    # print(head(pixels))
    pixels.indices.remove <- which(pixels_working[,"gene"] %in% list.to.remove, arr.ind = T)
    # print(head(pixels.indices.remove))
    # print(dim(pixels.indices.remove))
    pixels.no.remove <- pixels_working[(!rownames(pixels_working) %in% c(pixels.indices.remove)),]
    # print(dim(pixels.no.remove))
    # print(head(pixels.no.remove))
    rm(pixels_working)
    fov.no.remove <- CreateFOV(
        coords = segmentations,
        type = c("segmentation", "centroids"),
        molecules = pixels.no.remove,
        assay = assay_name
    )
    rm(counts.no.remove, pixels.no.remove)
    gc()
    return(list(assay.no.remove, fov.no.remove))
}

plan("multicore", workers = 30)
# The below errors out because it is reading in the Cell IDs incorrectly per https://github.com/satijalab/seurat/issues/7491
# xenium.obj <- LoadXenium(data.dir=input, fov = "fov", assay = "Xenium") #this might be fixed in Seurat5. # this appears to be fixed along with a few other things in v5.2.0
xenium.list <- ReadXenium(data.dir=input, outs = c("matrix", "microns"), type = c("centroids", "segmentations"), mols.qv.threshold = opt$mols.qv.threshold)
plan("sequential")

# fixing input matrix names
rownames(xenium.list$matrix$`Gene Expression`) <- gsub(":","-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub(":","-",xenium.list$microns$gene)
rownames(xenium.list$matrix$`Gene Expression`) <- gsub("\\.","-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub("\\.","-",xenium.list$microns$gene)
rownames(xenium.list$matrix$`Gene Expression`) <- gsub("_","-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub("_","-",xenium.list$microns$gene)
rownames(xenium.list$matrix$`Gene Expression`) <- gsub("\\+","-plus-",rownames(xenium.list$matrix$`Gene Expression`))
xenium.list$microns$gene <- gsub("\\+","-plus-",xenium.list$microns$gene)
# extracting the list of features that variant/reference and microbiome probes will be filtered from
targets <- rownames(xenium.list$matrix[["Gene Expression"]])

if (opt$with_snvs && opt$with_microbiome) { 
    # compiling the list of snvs features that will be removed
    if (is.na(opt$snv_probes)) {
        wt_targets <- targets[grepl("-WT",targets)]
        alt_targets <- targets[grepl("-ALT",targets)]
        snv.row.names.remove <- c(wt_targets, alt_targets)
    } else {
        snv.row.names.remove <- readLines(opt$snv_probes)
    }
    snv.row.names.remove <- probe_hyphens(snv.row.names.remove)
    print("Removing the following SNV probes from the default Xenium matrix")
    print(snv.row.names.remove)
    
    if (is.na(opt$microbiome_probes)) {
        bacteroides_targets  <- targets[grepl("Bacteroides", targets)]
        escherichia_targets  <- targets[grepl("Escherichia", targets)]
        fusobacterium_targets  <- targets[grepl("Fusobacterium", targets)]
        hpv_targets  <- targets[grepl("papillomavirus", targets)]
        microbiome.row.names.remove <- c(bacteroides_targets, escherichia_targets, fusobacterium_targets, hpv_targets)
    } else {
        microbiome.row.names.remove <- readLines(opt$microbiome_probes)
    }
    microbiome.row.names.remove <- probe_hyphens(microbiome.row.names.remove)
    print("Removing the following Microbiome probes from the default Xenium matrix")
    print(microbiome.row.names.remove)
    
    all.row.names.remove <- c(snv.row.names.remove, microbiome.row.names.remove)
    # each assay and FOV will be associated with the same segmentation data since we aren't subsetting cells yet.
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    # generating the object with the initial assay with everything
    assay <- "Xenium.with.snvs.microbiome"
    fov <- "fov.with.snvs.microbiome" #If you are using Seurat v4.3.0 do not rename this as it will result in you being unable to plot molecules.
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay) # This will follow the version option specified above
    xenium.obj[["BlankCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- all.coords
    # generating the assay and fov with genes and variant/reference probes
    assay <- "Xenium.with.snvs"
    fov <- "fov.with.snvs"
    with.snvs.list <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                           xenium.list$microns,
                                           segmentations.data,
                                           microbiome.row.names.remove,
                                           assay)
    xenium.obj[[assay]] <- with.snvs.list[[1]]
    xenium.obj[[fov]] <- with.snvs.list[[2]]
    # generating the assay and fov with genes and microbiome probes
    assay <- "Xenium.with.microbiome"
    fov <- "fov.with.microbiome"
    with.microbiome.list <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                           xenium.list$microns,
                                           segmentations.data,
                                           snv.row.names.remove,
                                           assay)
    xenium.obj[[assay]] <- with.microbiome.list[[1]]
    xenium.obj[[fov]] <- with.microbiome.list[[2]]
    # generating the assay and fov with only genes
    assay <- "Xenium"
    fov <- "fov"
    xenium.list.sub <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                                 xenium.list$microns,
                                                 segmentations.data,
                                                 all.row.names.remove,
                                                 assay)
    xenium.obj[[assay]] <- xenium.list.sub[[1]]
    xenium.obj[[fov]] <- xenium.list.sub[[2]]
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
} else if (opt$with_microbiome) {
    # compiling the list of snvs features that will be removed
    if (is.na(opt$microbiome_probes)) {
        bacteroides_targets  <- targets[grepl("Bacteroides", targets)]
        escherichia_targets  <- targets[grepl("Escherichia", targets)]
        fusobacterium_targets  <- targets[grepl("Fusobacterium", targets)]
        hpv_targets  <- targets[grepl("papillomavirus", targets)]
        microbiome.row.names.remove <- c(bacteroides_targets, escherichia_targets, fusobacterium_targets, hpv_targets)
    } else {
        microbiome.row.names.remove <- readLines(opt$microbiome_probes)
    }
    microbiome.row.names.remove <- probe_hyphens(microbiome.row.names.remove)
    print("Removing the following Microbiome probes from the default Xenium matrix")
    print(microbiome.row.names.remove)
    # each assay and FOV will be associated with the same segmentation data since we aren't subsetting cells yet.
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    # generating the object with the initial assay with everything
    assay <- "Xenium.with.microbiome"
    fov <- "fov.with.microbiome" #If you are using Seurat v4.3.0 do not rename this as it will result in you being unable to plot molecules.
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay) # This will follow the version option specified above
    xenium.obj[["BlankCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- all.coords
    
    # generating the assay and fov with only genes
    assay <- "Xenium"
    fov <- "fov"
    xenium.list.sub <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                                 xenium.list$microns,
                                                 segmentations.data,
                                                 microbiome.row.names.remove,
                                                 assay)
    xenium.obj[[assay]] <- xenium.list.sub[[1]]
    xenium.obj[[fov]] <- xenium.list.sub[[2]]
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
} else if (opt$with_snvs) {
    # compiling the list of snvs features that will be removed
    if (is.na(opt$snv_probes)) {
        wt_targets <- targets[grepl("-WT",targets)]
        alt_targets <- targets[grepl("-ALT",targets)]
        snv.row.names.remove <- c(wt_targets, alt_targets)
    } else {
        snv.row.names.remove <- readLines(opt$snv_probes)
    }
    snv.row.names.remove <- probe_hyphens(snv.row.names.remove)
    print("Removing the following SNV probes from the default Xenium matrix")
    print(snv.row.names.remove)
    # each assay and FOV will be associated with the same segmentation data since we aren't subsetting cells yet.
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    # generating the object with the initial assay with everything
    assay <- "Xenium.with.snvs"
    fov <- "fov.with.snvs"
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay) # This will follow the version option specified above
    xenium.obj[["BlankCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Probe"]])
    xenium.obj[[fov]] <- all.coords
    # generating the assay and fov with only genes
    assay <- "Xenium"
    fov <- "fov"
    xenium.list.sub <- create_assay_and_fov(xenium.list$matrix[["Gene Expression"]], 
                                                 xenium.list$microns,
                                                 segmentations.data,
                                                 snv.row.names.remove,
                                                 assay)
    xenium.obj[[assay]] <- xenium.list.sub[[1]]
    xenium.obj[[fov]] <- xenium.list.sub[[2]]
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
} else {
    assay <- "Xenium"
    segmentations.data <- list(
        "centroids" = CreateCentroids(xenium.list$centroids),
        "segmentation" = CreateSegmentation(xenium.list$segmentations)
    )
    all.coords <- CreateFOV(
        coords = segmentations.data,
        type = c("segmentation", "centroids"),
        molecules = xenium.list$microns,
        assay = assay
    )
    xenium.obj <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay) # This will follow the version option specified above
    xenium.obj[["BlankCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Unassigned Codeword"]])
    xenium.obj[["ControlCodeword"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Codeword"]])
    xenium.obj[["ControlProbe"]] <- CreateAssay(counts = xenium.list$matrix[["Negative Control Probe"]])
    fov <- "fov" #If you are using Seurat v4.3.0 do not rename this as it will result in you being unable to plot molecules.
    xenium.obj[[fov]] <- all.coords
    DefaultFOV(xenium.obj, assay='Xenium') <- 'fov'
}

#check to make sure that plotting with transcript/cell coordinates works
pdf(paste(out_path,"VlnPlot_QC_",sampleID,".pdf", sep=''),useDingbat=FALSE, width=8, height=6)
VlnPlot(xenium.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
dev.off()
# pdf(paste("ImageDimPlot_QC_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #This shows all of the dots on one plot
# ImageDimPlot(xenium.obj, fov = "fov", molecules = c("TFF2", "KRT7", "VCAN", "CFTR"), nmols = 20000)
# dev.off()
# pdf(paste("ImageDimPlot_QC_features_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# ImageFeaturePlot(xenium.obj, features = c("TFF2", "KRT7", "VCAN", "CFTR"), max.cutoff = c(20, 20, 15, 5), size = 0.25, cols = c("lightgrey", "blue", "yellow"), border.size= NA)
# dev.off()
# pdf(paste("ImageDimPlot_QC_features_2color_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# ImageFeaturePlot(xenium.obj, features = c("TFF2", "KRT7", "VCAN", "CFTR"), max.cutoff = c(20, 20, 15, 5), size = 0.25, cols = c("blue", "yellow"), border.size= NA)
# dev.off()
# save the raw
print("Saving the raw object")
saveRDS(xenium.obj, paste(out_path,sampleID,"_raw.rds",sep=''))

# check if there is a whitelist of barcodes that was provided and if so remove them.
source("/diskmnt/Projects/SenNet_analysis/Main.analysis/scripts/project_started/spatial/subset_obj_seurat.R")
# modified function originally from GATK to check if there are header lines in the cell barcode file that start with a '#'.
# function originally checked for and removed '@': https://github.com/broadinstitute/gatk/blob/master/src/main/resources/org/broadinstitute/hellbender/tools/copynumber/plotting/CNVPlottingLibrary.R
# changed function name to ReadCSV from ReadTSV
# function will mostly be used for reading in csv files, since the table output by Xenium Explorer is a csv file if you are using that to manually select cells belonging to a specific sample and that is the current expected use case.
# also allow for handling of CSV, TSV, and space separated values files. current implementation will break if the table relies on a combination of separateors in different fields but the input files are not anticipated to be complex and in 99% of cases will only use a single field separator.
ReadCSV = function(tsv_file) {
    # We need to filter out header lines beginning with '#';
    # however, the standard 'fread("grep ...")' causes issues with the default [gatk] Docker container, so we use a temporary file.
    # See https://github.com/broadinstitute/gatk/issues/4140.
    temp_file = tempfile() # if you are running this inside a docker container with the tmp folder mapped to a different system volume where to docker user does not have write access or with permission of the container so that the user does not have write access to the tmp folder then this will fail. You can sometimes avoid this by doing the following with your docker run commands: docker run -it -u $(id -u):$(id -g) <docker_repo/image:tag_goes_here>
    system(sprintf('grep -v ^# "%s" > %s', tsv_file, temp_file))
    return_table <- suppressWarnings(fread(temp_file, sep=",", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE)) 
    # this section of if statement allows me to check if the file is a tsv, csv, or spaces separated file. In the event that the file only has one column and no field separators then it will still return a 1 column dataframe.
    if (ncol(return_table) == 1) {
        # this is here in case someone provides a tsv instead of a csv.
        return_table <- suppressWarnings(fread(temp_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE))
    }
    if (ncol(return_table) == 1) {
        # on the off chance that someone is using a file with values separated by spaces then include this here. I don't expect this to every actually be needed.
        return_table <- suppressWarnings(fread(temp_file, sep=" ", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE, data.table=FALSE, showProgress=FALSE, verbose=FALSE))
    }
    file.remove(temp_file)
    return(return_table)
}

xenium.obj$original_Xenium_barcode <- rownames(xenium.obj@meta.data)
if (!(is.null(opt$barcode_whitelist))){
    writeLines(paste0("Filtering the Seurat object cell matrices using the provided barcode whitelist:\n",opt$barcode_whitelist))
    whitelist_table <- ReadCSV(opt$barcode_whitelist)
    whitelist_barcodes <- whitelist_table[,1]
    xenium.sub <- subset_opt(xenium.obj, subset = original_Xenium_barcode %in% c(whitelist_barcodes))
} else {
    xenium.sub <- xenium.obj
}
if (!(is.null(opt$coordinate_whitelist))){
    # alikhuseynov uploaded a custom croping script here: 
    # https://github.com/alikhuseynov/add-on_R/blob/develop/R/crop_seurat_v1.R
    # I found this after I already had built the above barcode whitelist functions so I'm just going to keep that as is and focus on subsetting the molecules. The Centroids and Segmentations are already handled by subsetting the cell ids.
    # # method 1
    # decide to use method 2 since it relies on native Seurat functions which are more likely to be supported moving forwards.
    # this method relies on partial solutions shared on github and stackoverflow and the documentation for spatial polygons:
    # https://github.com/satijalab/seurat/issues/8457 # whoever this person is "alikhuseynov" "they are a seurat object wizard.
    # https://stackoverflow.com/questions/37538230/how-to-change-class-from-data-frame-to-spatial-polygon
    # https://mhallwor.github.io/_pages/basics_SpatialPolygons 
    # print(paste0("Filtering the Seurat object fovs using the provided barcode whitelist:\n",opt$coordinate_whitelist))
    # whitelist_table <- ReadCSV(opt$coordinate_whitelist)
    # whitelist_table <- whitelist_table[,c(1,2)]
    # colnames(whitelist_table) <- c("x","y")
    # coordinate_mask <- sp::SpatialPolygons(list(sp::Polygons(list(sp::Polygon(whitelist_table)), ID="segmentation")))
    # Image_list <- Images(xenium.sub)
    # for (image_name in Image_list) {
    #     xenium.sub[[paste0(image_name,".test.1")]]  <- Overlay(xenium.sub[[image_name]], coordinate_mask, invert = FALSE)
    # }
    # method 2
    # this method is pulled from https://github.com/satijalab/seurat/issues/9257 
    print(paste0("Filtering the Seurat object fovs using the provided barcode whitelist:\n",opt$coordinate_whitelist))
    whitelist_table <- ReadCSV(opt$coordinate_whitelist)
    whitelist_table <- whitelist_table[,c(1,2)] # x and y coordinate. first column must be x and second column must be y.
    colnames(whitelist_table) <- c("x","y")
    whitelist_table$cell <- rep("segmentation", nrow(whitelist_table))
    whitelist_table <- whitelist_table[,c('x','y','cell')] #making sure everything is ordered properly.
    segmentation_mask <- CreateSegmentation(whitelist_table)
    Image_list <- Images(xenium.sub)
    for (image_name in Image_list) {
        xenium.sub[[image_name]] <- Overlay(xenium.sub[[image_name]], segmentation_mask, invert = FALSE)
    }
}

# some cells will be empty so you need to remove them
# works fine Seruat v4.3.0.1 but not in v5.0.1
# counts <- xenium.sub@assays$Xenium@counts
# empty_cells <- colnames(counts)[(colSums(counts) == 0)]
print("Filtering the Seurat object based on specified quality metrics")
if (sum(xenium.sub$nCount_Xenium <= opt$counts_cutoff) > 0) {
    xenium.empty <- subset_opt(xenium.sub, subset = nCount_Xenium <= opt$counts_cutoff)
    empty_cells = Cells(xenium.empty)
} else {
    empty_cells = c()
}
write.table(empty_cells,paste(out_path,"Empty_cells_",sampleID,'.tsv',sep=''),sep="\t",quote=FALSE)
rm(xenium.empty)
# cells_to_keep <- colnames(counts)[(colSums(counts) > opt$counts_cutoff)]
# xenium.sub <- subset(xenium.sub, cells = cells_to_keep)
#xenium.sub <- CreateSeuratObject(counts = xenium.list$matrix[["Gene Expression"]], assay = assay)
xenium.sub <- subset_opt(xenium.sub, subset = nCount_Xenium > opt$counts_cutoff)


print("Normalizing")
# normalize,umap,cluster, SCTv2 handles low number of genes very well compared to LogNormalize (see Ilya's results)
xenium.sub <- SCTransform(xenium.sub, assay = "Xenium", return.only.var.genes = F, vst.flavor="v2")
xenium.sub <- RunPCA(xenium.sub, npcs = num_pcs, features = rownames(xenium.sub))
xenium.sub <- RunUMAP(xenium.sub, dims = 1:num_pcs, reduction.name = paste("umap.",num_pcs,"PC", sep = ""), reduction.key = paste("UMAP",num_pcs,"PC_",sep=""))
xenium.sub <- FindNeighbors(xenium.sub, reduction = "pca", dims = 1:num_pcs)
xenium.sub <- FindClusters(xenium.sub, resolution = 0.3)
DefaultFOV(xenium.sub, assay='SCT') <- 'fov'
# the above line results in the following warning and may cause problems down the line. to do this properly you would copy the current 'fov' slot and save it to the FOV slot for the SCT assay I think. Not 100% sure how to do this yet.
# Warning: FOV 'fov' currently associated with assay 'Xenium', changing to 'SCT'

# save the processed object.
Idents(xenium.sub) <- "seurat_clusters"
# saveRDS(xenium.sub, file='HT270P1-S1H1A1Fp1Us1_1_processed.rds')
# adding in cell labels (here they are based on the xenium browser clusters)
# can't use the traditional method for seurat cluster adding and I think it has to do with the xenium browser clusters not being a factor.
xenium_browser_clusters <- read.table(paste(input,"/analysis/clustering/gene_expression_graphclust/clusters.csv",sep=''), sep = ",", header=T)
xenium_browser_clusters <- xenium_browser_clusters[(xenium_browser_clusters$Barcode %in% Cells(xenium.sub)),]
rownames(xenium_browser_clusters) <- xenium_browser_clusters$Barcode
xenium_browser_clusters$Barcode <- NULL
xenium.sub <- AddMetaData(object = xenium.sub, metadata = xenium_browser_clusters, col.name = "xenium_browser_graph_based_clusters")
# adding in cell labels (here they are based on the xenium browser clusters)
# can't use the traditional method for seurat cluster adding and I think it has to do with the xenium browser clusters not being a factor.
# could also be due to the xenium object not containing the right cell type labels.
#cell_types <- c("your","clusters","go","here")
#xenium.sub$cell_type_individual <- NA
#xenium_clusters = sort(unique(xenium.sub$xenium_browser_graph_based_clusters))
#for (i in 1:length(xenium_clusters)) {
#    value = xenium_clusters[i]
#    cell_type = cell_types[i]
#    xenium.sub$cell_type_individual[xenium.sub$xenium_browser_graph_based_clusters == value] <- cell_type
#}
# Setting the default values to be included in plots
DefaultAssay(xenium.sub) <- "SCT"
Idents(xenium.sub) <- "seurat_clusters"
DefaultFOV(xenium.sub, assay='SCT') <- 'fov'
DefaultBoundary(xenium.sub[["fov"]]) <- "segmentation"
# save the processed object.
saveRDS(xenium.sub, paste(out_path,sampleID,"_processed.rds",sep=''))
tmp_df <- xenium.sub@meta.data[,c("original_Xenium_barcode","seurat_clusters")]
colnames(tmp_df) <- c("cell_id","seurat_clusters")
write.table(tmp_df, paste0(out_path,sampleID,"_seurat_clusters.tsv"),row.names=FALSE, sep="\t", quote=FALSE)
colnames(tmp_df) <- c("cell_id","group")
write.table(tmp_df, paste0(out_path,sampleID,"_seurat_clusters.csv"),row.names=FALSE, sep=",", quote=FALSE)
print("generating_QC_plots")

# S Q U A R E S   Y O U R   P L O T 
coordinates <- Embeddings(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""))
dim_names <- colnames(coordinates)
dim_xlims <- c(floor(x = min(coordinates[, dim_names[1]])), ceiling(x = max(coordinates[, dim_names[1]])))
dim_ylims <- c(floor(min(coordinates[, dim_names[2]])), ceiling(x = max(coordinates[, dim_names[2]])))
sorted_dim_ylims <- sort(dim_ylims)
sorted_dim_xlims <- sort(dim_xlims)
x_length <- sorted_dim_xlims[2] - sorted_dim_xlims[1]
y_length <- sorted_dim_ylims[2] - sorted_dim_ylims[1]
xy_ratio <- x_length/y_length


pdf(paste(out_path,"DimPlots_QC_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6)
p1 <- rasterize(DimPlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), group.by = "seurat_clusters", label=TRUE, label.size=4, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers='Point', dpi=300)
p2 <- rasterize(DimPlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), group.by = "xenium_browser_graph_based_clusters", label=TRUE, label.size=4, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers='Point', dpi=300)
p3 <- rasterize(FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = "nCount_Xenium", label=TRUE, label.size=4, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers='Point', dpi=300)
p4 <- rasterize(FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = "nFeature_Xenium", label=TRUE, label.size=4, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers='Point', dpi=300)
print(p1)
print(p2)
print(p3)
print(p4)
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.color=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "segmentation"`
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.size=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "centroids"`
p5 <- rasterize(ImageFeaturePlot(xenium.sub, features = "nCount_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("blue", "yellow"), border.color = NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#000000", color = NA), plot.background = element_rect(fill = "#000000", color = NA), panel.border = element_rect(fill = NA, color = "#000000"),  text = element_text(color = "white")) + coord_fixed(ratio=1), layers='Polygon', dpi=600) #genes with high expression can have low misleading background so we set all zero value to lightgrey and max.cutoff is 90th percentile
grid:::grid.newpage()
grid:::grid.rect(gp=grid:::gpar(fill="#000000", col="#000000"))
ggplot2:::print.ggplot(p5, newpage=FALSE)
# print(p5)
p6 <- rasterize(ImageFeaturePlot(xenium.sub, features = "nFeature_Xenium", max.cutoff = 'q95', size = 0.15, cols = c("blue", "yellow"), border.color = NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#000000", color = NA), plot.background = element_rect(fill = "#000000", color = NA), panel.border = element_rect(fill = NA, color = "#000000"),  text = element_text(color = "white")) + coord_fixed(ratio=1), layers='Polygon', dpi=600) #max.cutoff is 90th percentile
grid:::grid.newpage()
grid:::grid.rect(gp=grid:::gpar(fill="#000000", col="#000000"))
ggplot2:::print.ggplot(p6, newpage=FALSE)
# print(p5)
p7 <- rasterize(ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", border.color=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#000000", color = NA), plot.background = element_rect(fill = "#000000", color = NA), panel.border = element_rect(fill = NA, color = "#000000"),  text = element_text(color = "white")) + coord_fixed(ratio=1), layers='Polygon', dpi=600)
grid:::grid.newpage()
grid:::grid.rect(gp=grid:::gpar(fill="#000000", col="#000000"))
ggplot2:::print.ggplot(p7, newpage=FALSE)
# print(p5)
#ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", border.color = "white", border.size=0.001)
p8 <- rasterize(ImageDimPlot(xenium.sub, fov = "fov", group.by="xenium_browser_graph_based_clusters", border.color=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#000000", color = NA), plot.background = element_rect(fill = "#000000", color = NA), panel.border = element_rect(fill = NA, color = "#000000"),  text = element_text(color = "white")) + coord_fixed(ratio=1), layers='Polygon', dpi=600)
grid:::grid.newpage()
grid:::grid.rect(gp=grid:::gpar(fill="#000000", col="#000000"))
ggplot2:::print.ggplot(p8, newpage=FALSE)
# print(p5)
#ImageDimPlot(xenium.sub, fov = "fov", group.by="xenium_browser_graph_based_clusters", border.color = "white", border.size=0.001)
dev.off()
if (is.na(opt$plot_features)){
    print("skipping plotting of features as no features were provided")
} else {
    genes_filtered <-c()
    plot_features=str_split(opt$plot_features, pattern=',')[[1]]
    plot_features=unique(plot_features)
    sample_genes <- rownames(x = xenium.sub@assays$SCT@counts)
    
    for (gene in plot_features) {
        if (gene %in% sample_genes) {
            genes_filtered <- c(genes_filtered,gene)
        }
    }
    print("the following genes are unique and are present in the SCT assay for inclusion in plots")
    print(genes_filtered)
    pdf(paste(out_path,"FeaturePlots_markers_of_interest_SCT_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6)
    for (marker in genes_filtered) {
        print(rasterize(FeaturePlot(xenium.sub, reduction = paste("umap.",num_pcs,"PC", sep = ""), features = marker, label=TRUE, label.size=4, raster=FALSE) + scale_x_continuous(limits = dim_xlims) + scale_y_continuous(limits = dim_ylims) + coord_fixed(ratio=xy_ratio), layers='Point', dpi=300))
    }
    dev.off()
    pdf(paste(out_path,"ImageDimPlot_QC_features_",sampleID,".pdf",sep=''), useDingbat=FALSE, width=8, height=6) #features splits out the different genes to separate plots.
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.color=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "segmentation"`
# If you are plotting an ImageDimPlot() then use ImageDimPlot(border.size=NA) if you are trying to generate a full section image and the cell boundaries are displayed as `DefaultBoundary(obj[["fov"]]) <- "centroids"`
    for (marker in genes_filtered){
        # print(rasterize(ImageFeaturePlot(xenium.sub, features = marker, max.cutoff = 'q95', size = 0.15, cols = c("lightgrey", "darkblue", "yellow"), border.color= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + coord_fixed(ratio=1), layers='Polygon', dpi=600)) #genes with high expression can have low misleading background so we set all zero value to lightgrey and max.cutoff is 90th percentile
        p9 <- rasterize(ImageFeaturePlot(xenium.sub, features = marker, max.cutoff = 'q90', size = 0.15, cols = c("darkblue", "yellow"), border.color= NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#000000", color = NA), plot.background = element_rect(fill = "#000000", color = NA), panel.border = element_rect(fill = NA, color = "#000000"),  text = element_text(color = "white")) + coord_fixed(ratio=1), layers='Polygon', dpi=600) #max.cutoff is 90th percentile
        grid:::grid.newpage()
        grid:::grid.rect(gp=grid:::gpar(fill="#000000", col="#000000"))
        ggplot2:::print.ggplot(p9, newpage=FALSE)
        # print(p9)
    }
    # uncomment this section if you want to check that Overlay is properly cropping the fov molecules of each image just like the cell barcode whitelist is doing.
    # new_colors = {set.seed(123); sample(grDevices::colors()[!(grepl("grey",grDevices::colors())) & !(grepl("gray",grDevices::colors()))])}[1:length(genes_filtered)]
    # all_images = Images(xenium.sub)
    # for (image in all_images) {
    #     print(image)
    #     DefaultBoundary(xenium.sub[[image]]) <- "segmentation"
    #     p10 <- rasterize(ImageDimPlot(xenium.sub, fov = image, #axes = TRUE,
    #                                  border.color = "white", border.size = 0.02, cols = rep("#000000", times = length(unique(xenium.sub$seurat_clusters))), #can only have 1 border color or need to provide a vector of border colors for each cell individually
    #                                  mols.size = 0.1,
    #                                  molecules = genes_filtered, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
    #                                  axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "#000000", color = NA), plot.background = element_rect(fill = "#000000", color = NA), panel.border = element_rect(fill = NA, color = "#000000"),  text = element_text(color = "white")) + labs(title=image) + coord_fixed(ratio=1), layers='Polygon', dpi=600)
    #     grid:::grid.newpage()
    #     grid:::grid.rect(gp=grid:::gpar(fill="#000000", col="#000000"))
    #     ggplot2:::print.ggplot(p10, newpage=FALSE)
    # }
    dev.off()
}

# Further examples of how to crop to new FOVs and 
# 
# options(future.globals.maxSize= +Inf)
# new_markers = c("CD68","CD4","FGG","KRT7","KRT19","ACTA2","PECAM1","CD8A")
# new_colors = c("magenta","cyan","yellow","pink","orange","red","green","white")
# Idents(nano.sub) <- "orig.ident"
# DefaultBoundary(xenium.sub[["fov"]]) <- "segmentation"
# pdf('test7_6.pdf',useDingbats = F, width = 8, height = 6)
# #ImageDimPlot(nano.obj, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", #size = 0.5, 
#              border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", #size = 0.5, 
#              #border.size = NA,
#              border.color = "white", border.size=0.001,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "fov", #axes = TRUE, 
#              #border.size = NA, 
#              border.color = "white", border.size = 0.001, cols = NA,
#              mols.size = 0.1,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes = F, stroke=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()

# #Zoom1
# # c(0, 3000), y = c(0, 3000)
# cropped.coords <- Crop(xenium.sub[["fov"]], x = c(0, 3000), y = c(0, 3000), coords = "plot")
# xenium.sub[["zoom1"]] <- cropped.coords
# DefaultBoundary(xenium.sub[["zoom1"]]) <- "segmentation"
# source("/diskmnt/Projects/Users/austins2/tools/subset_obj_seurat.R")
# pdf('test7_zoom1.pdf',useDingbats = F, width = 8, height = 6)
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "zoom1", group.by="seurat_clusters", #size = 0.5,
#              border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "zoom1", group.by="seurat_clusters", #size = 0.5,
#              #border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "zoom1", #axes = TRUE,
#              #border.size = NA,
#              #border.color = "white", border.size = 0.2, cols = "black",
#              mols.size = 0.5,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# xenium.sub = subset_opt(xenium.sub, subset = seurat_clusters %in% c(0,1,2))
# ImageDimPlot(xenium.sub, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "zoom1", #axes = TRUE,
#              #border.size = NA,
#              #border.color = "white", border.size = 0.2, cols = "black",
#              mols.size = 0.5,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes = F) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()

# # Zoom2 
# # x = c(1000, 4000), y = c(1000, 4000)
# cropped.coords <- Crop(xenium.sub[["fov"]], x = c(1000, 4000), y = c(1000, 4000), coords = "plot")
# xenium.sub[["zoom2"]] <- cropped.coords
# DefaultBoundary(xenium.sub[["zoom2"]]) <- "segmentation"
# pdf('test7_zoom2.pdf',useDingbats = F, width = 8, height = 6)
# #ImageDimPlot(nano.obj, fov = "fov", group.by="seurat_clusters", size = 0.15, border.size=NA)
# ImageDimPlot(xenium.sub, fov = "zoom2", group.by="seurat_clusters", #size = 0.5, 
#              border.size = NA,
#              #border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot(xenium.sub, fov = "zoom2", group.by="seurat_clusters", #size = 0.5, 
#              #border.size = NA,
#              border.color = "white", border.size=0.001,
#              axes = FALSE) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# ImageDimPlot.FIX(xenium.sub, fov = "zoom2", axes = TRUE,
#              #border.size = NA, 
#              border.color = "white", border.size = 0.2, cols = "black",
#              mols.size = 0.3,
#              molecules = new_markers, coord.fixed = FALSE, mols.cols = new_colors, nmols = 10000000,
#              axes=F, stroke=NA) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
# dev.off()
