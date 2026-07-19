################################################################################
#  Annotation & Pathway Enrichment Pipeline
#  Run in Google Colab R kernel AFTER the WGCNA pipeline
#
#  This script:
#    1. Annotates WGCNA hub metabolites (KEGG + HMDB)
#    2. Maps RNA hub genes to human orthologs
#    3. Annotates iRF stable biomarkers
#    4. Runs pathway enrichment via gprofiler2
#    5. Generates summary tables for presentation
################################################################################

# ============================================================================
# PART 0: Install & Load
# ============================================================================
install.packages(c("gprofiler2", "ggplot2", "dplyr", "tidyr"),
                 repos = "https://cloud.r-project.org")

library(gprofiler2)
library(ggplot2)
library(dplyr)
library(tidyr)

# ============================================================================
# PART 1: Set File Paths
# ============================================================================
# --- Annotation reference files (from instructor) ---
neg_kegg_ann_path <- "polar_neg_pkl_to_kegg_annotations.tsv"
neg_hmdb_ann_path <- "polar_neg_pkl_to_hmdb_annotations.tsv"
pos_kegg_ann_path <- "polar_pos_pkl_to_kegg_annotations.tsv"
pos_hmdb_ann_path <- "polar_pos_pkl_to_hmdb_annotations.tsv"
rna_hsa_map_path  <- "rna_dma_to_hsa_mappings.tsv"
rna_hsa_gn_path   <- "rna_dma_to_hsa_gn_mappings.tsv"
rna_dme_map_path  <- "rna_dma_to_dme_mappings.tsv"

# --- WGCNA output files ---
wgcna_dir <- "wgcna_results"
neg_hub_path <- file.path(wgcna_dir, "neg_metab_05_hub_features.csv")
pos_hub_path <- file.path(wgcna_dir, "pos_metab_05_hub_features.csv")
rna_hub_path <- file.path(wgcna_dir, "rna_05_hub_features.csv")
neg_mod_path <- file.path(wgcna_dir, "neg_metab_module_assignments.csv")
pos_mod_path <- file.path(wgcna_dir, "pos_metab_module_assignments.csv")
rna_mod_path <- file.path(wgcna_dir, "rna_module_assignments.csv")

# --- iRF stable biomarkers ---
neg_irf_path <- "neg_metab_stable_biomarkers_for_annotation.csv"
pos_irf_path <- "pos_metab_stable_biomarkers_for_annotation.csv"
rna_irf_path <- "rna_stable_biomarkers_for_annotation.csv"

# Output
ann_dir <- "annotation_results"
dir.create(ann_dir, showWarnings = FALSE)

# ============================================================================
# PART 2: Load Annotation Reference Files
# ============================================================================
cat("\n>>> Loading annotation reference files...\n")

# Metabolite annotations (format: feature_id \t mz \t KEGG/HMDB_ids)
load_metab_ann <- function(filepath) {
  df <- read.delim(filepath, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("FeatureID", "mz", "IDs"))
  df$IDs[is.na(df$IDs)] <- ""
  return(df)
}

neg_kegg <- load_metab_ann(neg_kegg_ann_path)
neg_hmdb <- load_metab_ann(neg_hmdb_ann_path)
pos_kegg <- load_metab_ann(pos_kegg_ann_path)
pos_hmdb <- load_metab_ann(pos_hmdb_ann_path)

# RNA mappings (format: daphnia_gene \t human_ids)
load_rna_map <- function(filepath) {
  df <- read.delim(filepath, header = FALSE, stringsAsFactors = FALSE,
                   col.names = c("DaphniaID", "MappedIDs"))
  df$MappedIDs[is.na(df$MappedIDs)] <- ""
  return(df)
}

rna_hsa    <- load_rna_map(rna_hsa_map_path)
rna_hsa_gn <- load_rna_map(rna_hsa_gn_path)
rna_dme    <- load_rna_map(rna_dme_map_path)

cat("  Annotation files loaded successfully\n")

# ============================================================================
# PART 3: Helper Functions for Annotation
# ============================================================================

#' Annotate metabolite features with KEGG and HMDB IDs
annotate_metabolites <- function(feature_ids, kegg_ref, hmdb_ref) {
  results <- data.frame(FeatureID = feature_ids, stringsAsFactors = FALSE)
  
  # Match KEGG
  kegg_match <- kegg_ref[match(feature_ids, kegg_ref$FeatureID), ]
  results$KEGG_IDs <- ifelse(is.na(kegg_match$IDs), "", kegg_match$IDs)
  results$mz <- ifelse(is.na(kegg_match$mz), NA, kegg_match$mz)
  
  # Match HMDB
  hmdb_match <- hmdb_ref[match(feature_ids, hmdb_ref$FeatureID), ]
  results$HMDB_IDs <- ifelse(is.na(hmdb_match$IDs), "", hmdb_match$IDs)
  
  # Has annotation?
  results$HasAnnotation <- (results$KEGG_IDs != "" | results$HMDB_IDs != "")
  
  return(results)
}

#' Map Daphnia gene IDs to human orthologs
map_rna_to_human <- function(gene_ids, hsa_ref, hsa_gn_ref) {
  results <- data.frame(DaphniaID = gene_ids, stringsAsFactors = FALSE)
  
  # Ensembl IDs
  hsa_match <- hsa_ref[match(gene_ids, hsa_ref$DaphniaID), ]
  results$HSA_Ensembl <- ifelse(is.na(hsa_match$MappedIDs), "", hsa_match$MappedIDs)
  
  # Gene names
  gn_match <- hsa_gn_ref[match(gene_ids, hsa_gn_ref$DaphniaID), ]
  results$HSA_GeneName <- ifelse(is.na(gn_match$MappedIDs), "", gn_match$MappedIDs)
  
  results$HasMapping <- (results$HSA_Ensembl != "")
  
  return(results)
}

# ============================================================================
# PART 4: Annotate WGCNA Hub Features
# ============================================================================
cat("\n>>> Annotating WGCNA hub features...\n")

# --- 4.1 Negative metabolomics hubs ---
neg_hubs <- read.csv(neg_hub_path, stringsAsFactors = FALSE)
neg_hubs_ann <- annotate_metabolites(neg_hubs$Feature, neg_kegg, neg_hmdb)
neg_hubs_ann <- cbind(neg_hubs, neg_hubs_ann[, c("mz", "KEGG_IDs", "HMDB_IDs", "HasAnnotation")])
write.csv(neg_hubs_ann, file.path(ann_dir, "neg_metab_hubs_annotated.csv"), row.names = FALSE)
cat("  Neg metab hubs:", sum(neg_hubs_ann$HasAnnotation), "/", nrow(neg_hubs_ann), "annotated\n")

# --- 4.2 Positive metabolomics hubs ---
pos_hubs <- read.csv(pos_hub_path, stringsAsFactors = FALSE)
pos_hubs_ann <- annotate_metabolites(pos_hubs$Feature, pos_kegg, pos_hmdb)
pos_hubs_ann <- cbind(pos_hubs, pos_hubs_ann[, c("mz", "KEGG_IDs", "HMDB_IDs", "HasAnnotation")])
write.csv(pos_hubs_ann, file.path(ann_dir, "pos_metab_hubs_annotated.csv"), row.names = FALSE)
cat("  Pos metab hubs:", sum(pos_hubs_ann$HasAnnotation), "/", nrow(pos_hubs_ann), "annotated\n")

# --- 4.3 RNA hubs ---
rna_hubs <- read.csv(rna_hub_path, stringsAsFactors = FALSE)
rna_hubs_ann <- map_rna_to_human(rna_hubs$Feature, rna_hsa, rna_hsa_gn)
rna_hubs_ann <- cbind(rna_hubs, rna_hubs_ann[, c("HSA_Ensembl", "HSA_GeneName", "HasMapping")])
write.csv(rna_hubs_ann, file.path(ann_dir, "rna_hubs_annotated.csv"), row.names = FALSE)
cat("  RNA hubs:", sum(rna_hubs_ann$HasMapping), "/", nrow(rna_hubs_ann), "mapped to human\n")

# ============================================================================
# PART 5: Annotate iRF Stable Biomarkers
# ============================================================================
cat("\n>>> Annotating iRF stable biomarkers...\n")

# --- 5.1 Negative metabolomics biomarkers ---
neg_irf <- read.csv(neg_irf_path, stringsAsFactors = FALSE)
neg_irf_ann <- annotate_metabolites(neg_irf$FeatureID, neg_kegg, neg_hmdb)
neg_irf_ann <- cbind(neg_irf, neg_irf_ann[, c("mz", "KEGG_IDs", "HMDB_IDs", "HasAnnotation")])
write.csv(neg_irf_ann, file.path(ann_dir, "neg_irf_biomarkers_annotated.csv"), row.names = FALSE)
cat("  Neg iRF biomarkers:", sum(neg_irf_ann$HasAnnotation), "/", nrow(neg_irf_ann), "annotated\n")

# --- 5.2 Positive metabolomics biomarkers ---
pos_irf <- read.csv(pos_irf_path, stringsAsFactors = FALSE)
pos_irf_ann <- annotate_metabolites(pos_irf$FeatureID, pos_kegg, pos_hmdb)
pos_irf_ann <- cbind(pos_irf, pos_irf_ann[, c("mz", "KEGG_IDs", "HMDB_IDs", "HasAnnotation")])
write.csv(pos_irf_ann, file.path(ann_dir, "pos_irf_biomarkers_annotated.csv"), row.names = FALSE)
cat("  Pos iRF biomarkers:", sum(pos_irf_ann$HasAnnotation), "/", nrow(pos_irf_ann), "annotated\n")

# --- 5.3 RNA biomarkers ---
rna_irf <- read.csv(rna_irf_path, stringsAsFactors = FALSE)
rna_irf_ann <- map_rna_to_human(rna_irf$FeatureID, rna_hsa, rna_hsa_gn)
rna_irf_ann <- cbind(rna_irf, rna_irf_ann[, c("HSA_Ensembl", "HSA_GeneName", "HasMapping")])
write.csv(rna_irf_ann, file.path(ann_dir, "rna_irf_biomarkers_annotated.csv"), row.names = FALSE)
cat("  RNA iRF biomarkers:", sum(rna_irf_ann$HasMapping), "/", nrow(rna_irf_ann), "mapped\n")

# ============================================================================
# PART 6: Pathway Enrichment for KEY Modules (RNA)
# ============================================================================
cat("\n>>> Running pathway enrichment for key RNA modules...\n")

# Load RNA module assignments
rna_mods <- read.csv(rna_mod_path, stringsAsFactors = FALSE)

# KEY MODULES from heatmap analysis:
# ME4: SiteNum cor=0.38 (up downstream), IsDownstream cor=0.30
# ME5: Is10x cor=-0.36 (down at 10x)
# ME7: Is1x cor=-0.33 (down at 1x)
# ME6: Is10x cor=0.25 (up at 10x)

key_rna_modules <- c("4", "5", "7", "6")
# Note: module assignments use color names, not numbers.
# We need to check what colors correspond to these MEs.
# The module number in ME4 = the 4th module label = look at module_assignments

# Actually, WGCNA labels2colors maps numeric labels to colors.
# Let's work with the color assignments directly.
# We need to figure out which color = which ME number.
# The safest way: load the RData and check, OR just enumerate modules.

# Let's get all unique modules
all_modules <- unique(rna_mods$Module)
cat("  RNA modules found:", paste(all_modules, collapse = ", "), "\n")

#' Run gprofiler2 enrichment for a set of Daphnia genes
#' Maps to human first, then queries
run_enrichment <- function(daphnia_ids, module_name, hsa_ref, output_dir) {
  # Map to human Ensembl IDs
  mapped <- hsa_ref[hsa_ref$DaphniaID %in% daphnia_ids, ]
  mapped <- mapped[mapped$MappedIDs != "", ]
  
  if (nrow(mapped) == 0) {
    cat("    No human orthologs found for module", module_name, "\n")
    return(NULL)
  }
  
  # Expand semicolon-separated IDs
  all_ensembl <- unique(unlist(strsplit(mapped$MappedIDs, ";")))
  all_ensembl <- all_ensembl[all_ensembl != ""]
  
  cat("    Module", module_name, ":", length(daphnia_ids), "Daphnia genes →",
      length(all_ensembl), "human Ensembl IDs\n")
  
  if (length(all_ensembl) < 5) {
    cat("    Too few genes for enrichment, skipping\n")
    return(NULL)
  }
  
  # Run gprofiler2
  tryCatch({
    gost_res <- gost(
      query = all_ensembl,
      organism = "hsapiens",
      ordered_query = FALSE,
      significant = TRUE,
      user_threshold = 0.05,
      correction_method = "fdr",
      sources = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "WP")
    )
    
    if (is.null(gost_res$result) || nrow(gost_res$result) == 0) {
      cat("    No significant enrichment found\n")
      return(NULL)
    }
    
    # Save results
    res <- gost_res$result[, c("source", "term_id", "term_name", 
                                "p_value", "term_size", "intersection_size",
                                "query_size")]
    res <- res[order(res$p_value), ]
    
    write.csv(res, file.path(output_dir, 
              paste0("rna_module_", module_name, "_enrichment.csv")),
              row.names = FALSE)
    
    cat("    Found", nrow(res), "significant terms\n")
    
    # Plot top 15 terms
    top_terms <- head(res, 15)
    top_terms$term_short <- substr(top_terms$term_name, 1, 50)
    top_terms$neg_log_p <- -log10(top_terms$p_value)
    top_terms$term_short <- factor(top_terms$term_short, 
                                    levels = rev(top_terms$term_short))
    
    p <- ggplot(top_terms, aes(x = neg_log_p, y = term_short, fill = source)) +
      geom_col(width = 0.7) +
      scale_fill_brewer(palette = "Set2") +
      theme_bw(base_size = 11) +
      labs(title = paste("Module", module_name, "- Pathway Enrichment"),
           x = "-log10(adjusted p-value)", y = "", fill = "Source") +
      theme(axis.text.y = element_text(size = 9))
    
    ggsave(file.path(output_dir, 
           paste0("rna_module_", module_name, "_enrichment_plot.pdf")),
           p, width = 10, height = 6)
    
    return(res)
    
  }, error = function(e) {
    cat("    gprofiler2 error:", conditionMessage(e), "\n")
    return(NULL)
  })
}

# Run enrichment for each key module color
# First, we need to identify which colors correspond to the key modules.
# Since we don't have the color-to-number mapping loaded, let's run
# enrichment for ALL non-grey modules (safe approach)

enrichment_results <- list()
for (mod_color in all_modules) {
  if (mod_color == "grey") next  # Skip unassigned
  
  mod_genes <- rna_mods$Feature[rna_mods$Module == mod_color]
  cat("\n  Processing RNA module:", mod_color, "(", length(mod_genes), "genes)\n")
  
  res <- run_enrichment(mod_genes, mod_color, rna_hsa, ann_dir)
  if (!is.null(res)) {
    enrichment_results[[mod_color]] <- res
  }
}

# ============================================================================
# PART 7: Metabolite Pathway Lookup via KEGG IDs
# ============================================================================
cat("\n>>> Extracting KEGG compound IDs for metabolite pathway analysis...\n")

# For metabolites, gprofiler2 doesn't work directly.
# Instead, we extract KEGG compound IDs for use in:
#   - KEGG Mapper: https://www.genome.jp/kegg/tool/map_pathway1.html
#   - IMPaLA: http://impala.molgen.mpg.de/
#   - MetaboAnalyst: https://www.metaboanalyst.ca/

# Function to extract all KEGG IDs from a module
extract_kegg_for_module <- function(module_features, kegg_ref, module_name) {
  matched <- kegg_ref[kegg_ref$FeatureID %in% module_features, ]
  matched <- matched[matched$IDs != "", ]
  
  all_kegg <- unique(unlist(strsplit(matched$IDs, ";")))
  all_kegg <- all_kegg[all_kegg != ""]
  
  cat("  Module", module_name, ":", length(module_features), "features →",
      nrow(matched), "annotated →", length(all_kegg), "unique KEGG IDs\n")
  
  return(all_kegg)
}

# Get module assignments for metabolites
neg_mods <- read.csv(neg_mod_path, stringsAsFactors = FALSE)
pos_mods <- read.csv(pos_mod_path, stringsAsFactors = FALSE)

# Extract KEGG IDs for all non-grey modules
cat("\n  --- Negative metabolomics ---\n")
neg_kegg_by_module <- list()
for (mod_color in unique(neg_mods$Module)) {
  if (mod_color == "grey") next
  mod_features <- neg_mods$Feature[neg_mods$Module == mod_color]
  kegg_ids <- extract_kegg_for_module(mod_features, neg_kegg, 
                                       paste0("neg_", mod_color))
  if (length(kegg_ids) > 0) {
    neg_kegg_by_module[[mod_color]] <- kegg_ids
    writeLines(kegg_ids, file.path(ann_dir, 
               paste0("neg_module_", mod_color, "_kegg_ids.txt")))
  }
}

cat("\n  --- Positive metabolomics ---\n")
pos_kegg_by_module <- list()
for (mod_color in unique(pos_mods$Module)) {
  if (mod_color == "grey") next
  mod_features <- pos_mods$Feature[pos_mods$Module == mod_color]
  kegg_ids <- extract_kegg_for_module(mod_features, pos_kegg, 
                                       paste0("pos_", mod_color))
  if (length(kegg_ids) > 0) {
    pos_kegg_by_module[[mod_color]] <- kegg_ids
    writeLines(kegg_ids, file.path(ann_dir, 
               paste0("pos_module_", mod_color, "_kegg_ids.txt")))
  }
}

# Also extract KEGG IDs for iRF biomarkers
cat("\n  --- iRF biomarker KEGG IDs ---\n")
neg_irf_kegg <- extract_kegg_for_module(neg_irf$FeatureID, neg_kegg, "neg_iRF")
pos_irf_kegg <- extract_kegg_for_module(pos_irf$FeatureID, pos_kegg, "pos_iRF")

if (length(neg_irf_kegg) > 0) 
  writeLines(neg_irf_kegg, file.path(ann_dir, "neg_irf_kegg_ids.txt"))
if (length(pos_irf_kegg) > 0) 
  writeLines(pos_irf_kegg, file.path(ann_dir, "pos_irf_kegg_ids.txt"))

# ============================================================================
# PART 8: Generate Summary Report
# ============================================================================
cat("\n>>> Generating summary report...\n")

# --- 8.1 Key module summary with top annotated hub features ---
generate_module_summary <- function(hubs_annotated, dataset_name) {
  # Only annotated features
  ann <- hubs_annotated[hubs_annotated$HasAnnotation == TRUE | 
                        hubs_annotated$HasMapping == TRUE, ]
  
  if (nrow(ann) == 0) return(NULL)
  
  # For metabolites
  if ("KEGG_IDs" %in% colnames(ann)) {
    summary_df <- ann %>%
      group_by(Module) %>%
      arrange(desc(abs(kME))) %>%
      slice_head(n = 5) %>%
      select(Module, Feature, kME, mz, KEGG_IDs, HMDB_IDs) %>%
      ungroup()
  } else {
    # For RNA
    summary_df <- ann %>%
      group_by(Module) %>%
      arrange(desc(abs(kME))) %>%
      slice_head(n = 5) %>%
      select(Module, Feature, kME, HSA_GeneName) %>%
      ungroup()
  }
  
  write.csv(summary_df, file.path(ann_dir, 
            paste0(dataset_name, "_key_hubs_summary.csv")), row.names = FALSE)
  
  return(summary_df)
}

neg_summary <- generate_module_summary(neg_hubs_ann, "neg_metab")
pos_summary <- generate_module_summary(pos_hubs_ann, "pos_metab")
rna_summary <- generate_module_summary(rna_hubs_ann, "rna")

# --- 8.2 Overlap between iRF biomarkers and WGCNA hub features ---
cat("\n>>> Checking overlap: iRF biomarkers vs WGCNA hub features...\n")

check_overlap <- function(irf_ids, hub_ids, dataset_name) {
  overlap <- intersect(irf_ids, hub_ids)
  cat("  ", dataset_name, ": iRF has", length(irf_ids), 
      "biomarkers, WGCNA has", length(hub_ids), "hub features,",
      length(overlap), "overlap\n")
  return(overlap)
}

neg_overlap <- check_overlap(neg_irf$FeatureID, neg_hubs$Feature, "Neg metab")
pos_overlap <- check_overlap(pos_irf$FeatureID, pos_hubs$Feature, "Pos metab")
rna_overlap <- check_overlap(rna_irf$FeatureID, rna_hubs$Feature, "RNA")

# Save overlaps
overlap_df <- data.frame(
  Dataset = c(rep("neg_metab", max(1, length(neg_overlap))),
              rep("pos_metab", max(1, length(pos_overlap))),
              rep("rna", max(1, length(rna_overlap)))),
  Feature = c(if(length(neg_overlap)>0) neg_overlap else NA,
              if(length(pos_overlap)>0) pos_overlap else NA,
              if(length(rna_overlap)>0) rna_overlap else NA),
  stringsAsFactors = FALSE
)
overlap_df <- overlap_df[!is.na(overlap_df$Feature), ]
write.csv(overlap_df, file.path(ann_dir, "irf_wgcna_overlap.csv"), row.names = FALSE)


# ============================================================================
# PART 9: Instructions for Manual Pathway Tools
# ============================================================================
cat("\n========================================================\n")
cat("  ANNOTATION & ENRICHMENT COMPLETE\n")
cat("========================================================\n")
cat("\nOutput files in 'annotation_results/':\n")
cat("  *_hubs_annotated.csv        — WGCNA hub features with IDs\n")
cat("  *_irf_biomarkers_annotated.csv — iRF biomarkers with IDs\n")
cat("  rna_module_*_enrichment.csv  — GO/KEGG/Reactome enrichment\n")
cat("  rna_module_*_enrichment_plot.pdf — Enrichment bar plots\n")
cat("  *_module_*_kegg_ids.txt     — KEGG IDs for manual lookup\n")
cat("  *_key_hubs_summary.csv      — Top annotated hubs per module\n")
cat("  irf_wgcna_overlap.csv       — Features found by BOTH methods\n")

cat("\n========================================================\n")
cat("  NEXT STEPS (Manual)\n")
cat("========================================================\n")
cat("\n1. KEGG Pathway Mapper (for metabolite modules):\n")
cat("   https://www.genome.jp/kegg/tool/map_pathway1.html\n")
cat("   → Paste KEGG compound IDs from *_kegg_ids.txt files\n")
cat("   → Select 'Compound' and submit\n")
cat("   → Screenshot highlighted pathways for your slides\n")
cat("\n2. IMPaLA (multi-omics pathway, if you have time):\n")
cat("   http://impala.molgen.mpg.de/\n")
cat("   → Upload both KEGG metabolite IDs + human gene IDs\n")
cat("   → This shows pathways enriched across BOTH omics levels\n")
cat("\n3. For presentation slides:\n")
cat("   → Use module-trait heatmaps as the main result figure\n")
cat("   → Show eigengene patterns for the top 2-3 modules\n")
cat("   → Present enrichment bar plots for RNA modules\n")
cat("   → Summarise key metabolites from hub annotation tables\n")
cat("   → Discuss the non-linear dose response (hormesis)\n")
