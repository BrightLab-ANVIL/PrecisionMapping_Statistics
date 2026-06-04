## 1. Libraries ----
library(RNifti)
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(rstatix) 

## 2. Safe loader + Dice coefficient helpers ----
read_nifti_3d <- function(path, label) {
  if (!file.exists(path)) {
    stop(label, " file not found: ", path)
  }
  img <- RNifti::readNifti(path)
  d <- dim(img)
  if (length(d) != 3L) {
    stop(label, " is not 3D (dim = ", paste(d, collapse = "x"), "): ", path)
  }
  img
}

dice_coeff <- function(a, b) {
  if (!is.array(a) || !is.array(b) || length(dim(a)) != 3L || length(dim(b)) != 3L) {
    stop("dice_coeff(): inputs must both be 3D arrays")
  }
  if (!all(dim(a) == dim(b))) {
    stop("dice_coeff(): arrays must have matching dimensions")
  }
  
  a <- as.logical(a)
  b <- as.logical(b)
  
  a_sum <- sum(a)
  b_sum <- sum(b)
  
  if (a_sum == 0L || b_sum == 0L) return(NA_real_)
  
  (2 * sum(a & b)) / (a_sum + b_sum)
}

## 3. Subjects, sessions, base paths ----
## ---- Subject IDs and paths ----
subs <- c("sub-01")  # add/remove as needed

# Root directory for all subjects
base_root <- ""

# Session mapping: specify the session for each subject
session_map <- tibble::tribble(
  ~sub_id, ~session,
  "sub-01", "ses-02"
)

get_session <- function(sub) {
  sess <- session_map %>% filter(sub_id == !!sub) %>% pull(session)
  if (length(sess) != 1L) stop("Session not found or ambiguous for subject: ", sub)
  sess
}

audio_dir <- function(sub) {
  sess <- get_session(sub)
  file.path(base_root, sub, sess, "func", "AUDIO")
}

## ---- Folds ----
# fold_letter is what changes in the folder names (a, b, c, d)
folds <- tibble::tribble(
  ~fold_id, ~fold_letter,
  "foldA",  "a",
  "foldB",  "b",
  "foldC",  "c",
  "foldD",  "d"
)
fold_ids <- folds$fold_id

## 4. Methods: folder patterns and file suffixes ----
# folder_pattern: a string with %s where the fold_letter goes
# suffix: the part of the filename AFTER the subject ID

method_specs <- tibble::tribble(
  ~method,      ~folder_pattern,                                      ~suffix,
  # FDR and Top5 share the "output.GLM ..." folders
  "FDR_p05",        "output.GLM_ICA_ConsOrth_concat3%sruns_sm_1ses",      "AUDIO_TaskReg_tstat_fdrp05_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p01",        "output.GLM_ICA_ConsOrth_concat3%sruns_sm_1ses",      "AUDIO_TaskReg_tstat_fdrp01_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p001",        "output.GLM_ICA_ConsOrth_concat3%sruns_sm_1ses",      "AUDIO_TaskReg_tstat_fdrp001_wholeBrain_func2anat_bin.nii.gz",
  "Top5",       "output.GLM_ICA_ConsOrth_concat3%sruns_sm_1ses",      "AUDIO_TaskReg_tstat_%s_thrP95_func2anat_bin.nii.gz",
  "Top10",       "output.GLM_ICA_ConsOrth_concat3%sruns_sm_1ses",      "AUDIO_TaskReg_tstat_%s_thrP90_func2anat_bin.nii.gz",
  "ClusterFWE_p005", "Cluster_Audio-0_ICA_concat3%sruns_sm_1ses",          "Cluster_Audio-0_clusters_bcoef_p005_a05_func2anat_bin.nii.gz",
  "ClusterFWE_p001", "Cluster_Audio-0_ICA_concat3%sruns_sm_1ses",          "Cluster_Audio-0_clusters_bcoef_p001_a05_func2anat_bin.nii.gz"
)

methods <- method_specs$method

# Define each method's "gold" file, which is their 40-minute concatenated data cluster
method_region_gold <- tibble::tribble(
  ~method,          ~region,      ~folder_pattern_gold,                                       ~gold_filename,
  "FDR_p05",        "thalamus",   "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp05_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p05",        "cerebellum", "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp05_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p05",        "brainstem",  "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp05_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p01",        "thalamus",   "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp01_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p01",        "cerebellum", "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp01_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p01",        "brainstem",  "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp01_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p001",        "thalamus",   "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp001_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p001",        "cerebellum", "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp001_wholeBrain_func2anat_bin.nii.gz",
  "FDR_p001",        "brainstem",  "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",              "AUDIO_TaskReg_tstat_fdrp001_wholeBrain_func2anat_bin.nii.gz",
  "Top5",          "thalamus",   "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",                "AUDIO_TaskReg_tstat_thalamus_thr50_thrP95_func2anat_bin.nii.gz",
  "Top5",          "cerebellum", "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",                "AUDIO_TaskReg_tstat_cerebellum_thrP95_func2anat_bin.nii.gz",
  "Top5",          "brainstem",  "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",                "AUDIO_TaskReg_tstat_brainstem_thrP95_func2anat_bin.nii.gz",
  "Top10",          "thalamus",   "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",                "AUDIO_TaskReg_tstat_thalamus_thr50_thrP90_func2anat_bin.nii.gz",
  "Top10",          "cerebellum", "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",                "AUDIO_TaskReg_tstat_cerebellum_thrP90_func2anat_bin.nii.gz",
  "Top10",          "brainstem",  "output.GLM_ICA_ConsOrth_concat4runs_sm_1ses",                "AUDIO_TaskReg_tstat_brainstem_thrP90_func2anat_bin.nii.gz",
  "ClusterFWE_p005",          "thalamus",   "Cluster_Audio-0_ICA_ConsOrth_concat_sm",                "Cluster_Audio-0_clusters_bcoef_p005_a05_func2anat_bin.nii.gz",
  "ClusterFWE_p005",          "cerebellum", "Cluster_Audio-0_ICA_ConsOrth_concat_sm",                "Cluster_Audio-0_clusters_bcoef_p005_a05_func2anat_bin.nii.gz",
  "ClusterFWE_p005",          "brainstem",  "Cluster_Audio-0_ICA_ConsOrth_concat_sm",                "Cluster_Audio-0_clusters_bcoef_p005_a05_func2anat_bin.nii.gz",
  "ClusterFWE_p001",          "thalamus",   "Cluster_Audio-0_ICA_ConsOrth_concat_sm",                "Cluster_Audio-0_clusters_bcoef_p001_a05_func2anat_bin.nii.gz",
  "ClusterFWE_p001",          "cerebellum", "Cluster_Audio-0_ICA_ConsOrth_concat_sm",                "Cluster_Audio-0_clusters_bcoef_p001_a05_func2anat_bin.nii.gz",
  "ClusterFWE_p001",          "brainstem",  "Cluster_Audio-0_ICA_ConsOrth_concat_sm",                "Cluster_Audio-0_clusters_bcoef_p001_a05_func2anat_bin.nii.gz"
)

# Define the 3 regions
regions <- tibble::tribble(
  ~region,      ~suffix_replace,        ~gold_filename,                                    ~mask_filename,
  "thalamus",   "thalamus_thr50",       "Sitek-sub-invivo_MNI_MGN_combined_2mm_bin_stand2func_bin_func2anat_bin.nii.gz", "harvardoxford-subcortical_prob_CombinedThalamus_thr50_CORRECT_dilateSitekMGN_bin_stand2func_bin_func2anat_bin.nii.gz",
  "cerebellum", "cerebellum",     "cerebellum_mnifnirt_prob_Combined_VIIbVIIIa_thr50_2mm_bin_stand2func_bin_func2anat_bin.nii.gz", "cerebellum_mnifnirt_prob_Combined_VIIbVIIIa_thr50_2mm_bin_stand2func_bin_func2anat_bin_dilate2vox.nii.gz",
  "brainstem",  "brainstem",      "Sitek-sub-invivo_MNI_IC_combined_2mm_bin_stand2func_bin_func2anat_bin.nii.gz", "harvardoxford-subcortical_Brainstem_thr25_midbrain-IC-CORRECT_bin_stand2func_bin_func2anat_bin.nii.gz"
)

## 5. Region-specific paths (parameterized) ----
get_region_gold_path <- function(sub, method, region) {
  row <- method_region_gold %>%
    dplyr::filter(method == !!method, region == !!region)
  if (nrow(row) != 1L) {
    stop("No unique gold path for ", method, " / ", region)
  }
  
  folder    <- row$folder_pattern_gold
  file_name <- paste0(sub, "_", row$gold_filename)
  
  file.path(audio_dir(sub), folder, file_name)
}

load_region_gold <- function(sub, method, region) {
  path <- get_region_gold_path(sub, method, region)
  read_nifti_3d(path, paste0(region, " gold (", sub, ", ", method, ")"))
}

load_region_mask <- function(sub, region) {
  r <- regions %>% dplyr::filter(region == !!region)
  if (nrow(r) != 1L) stop("Unknown region: ", region)
  
  path <- file.path(audio_dir(sub), "output.reg_2mm", r$mask_filename)
  img  <- read_nifti_3d(path, paste0(region, " mask (", sub, ")"))
  as.array(img) > 0
}

# Updated map_path to handle region-specific suffixes
map_path <- function(sub, method, fold_id, region) {
  frow <- folds %>% filter(fold_id == !!fold_id)
  mrow <- method_specs %>% filter(method == !!method)
  r <- regions %>% filter(region == !!region)
  
  if (nrow(frow) != 1L || nrow(mrow) != 1L || nrow(r) != 1L) {
    stop("Bad fold/method/region: ", fold_id, "/", method, "/", region)
  }
  
  suffix <- mrow$suffix
  if (grepl("%s", suffix)) {
    suffix <- sprintf(suffix, r$suffix_replace)
  }
  
  folder <- sprintf(mrow$folder_pattern, frow$fold_letter)
  file_name <- paste0(sub, "_", suffix)
  file.path(audio_dir(sub), folder, file_name)
}

load_map_img <- function(sub, method, fold_id, region) {
  read_nifti_3d(map_path(sub, method, fold_id, region),
                paste0("Map (", sub, ", ", method, ", ", region, ", ", fold_id, ")"))
}

## 6. Dice computation loop (MULTI-REGION, LEFT/RIGHT HEMISPHERES) ----
dice_results <- purrr::map_dfr(regions$region, function(region) {
  purrr::map_dfr(subs, function(sub) {
    
    roi <- load_region_mask(sub, region)
    
    purrr::map_dfr(methods, function(meth) {
      
      gold_img_obj <- load_region_gold(sub, meth, region)
      gold_arr     <- as.array(gold_img_obj)
      
      if (!all(dim(roi) == dim(gold_arr))) {
        warning("Dim mismatch for ", sub, "/", meth, "/", region)
        return(tibble())
      }
      
      dims <- dim(gold_arr)
      
      mid_x <- floor(dims[1] / 2)
      hemi_list <- list(
        left  = 1:mid_x,
        right = (mid_x + 1):dims[1]
      )
      
      purrr::map_dfr(names(hemi_list), function(hemi_name) {
        x_idx <- hemi_list[[hemi_name]]
        
        roi_hemi <- roi
        roi_hemi[-x_idx, , ] <- FALSE
        
        gold_bin  <- gold_arr > 0
        gold_hemi <- gold_bin
        gold_hemi[!roi_hemi] <- FALSE
        
        purrr::map_dfr(fold_ids, function(fid) {
          
          map_img_obj <- load_map_img(sub, meth, fid, region)
          map_arr     <- as.array(map_img_obj)
          
          if (!all(dim(map_arr) == dim(gold_arr))) {
            warning("Map/gold dim mismatch: ", sub, "/", region, "/", meth, "/", fid)
            return(tibble())
          }
          
          map_bin  <- map_arr > 0
          map_hemi <- map_bin
          map_hemi[!roi_hemi] <- FALSE
          
          dice <- dice_coeff(gold_hemi, map_hemi)
          
          tibble(
            sub_id      = sub,
            region      = region,
            hemi        = hemi_name,
            method      = meth,
            fold        = fid,
            dice        = dice
          )
        })
      })
    })
  })
})

## 7. Summary data for bars + points ----
dice_results <- dice_results %>% filter(!is.na(dice))

dice_bar_plot <- dice_results %>%
  group_by(region, method, hemi) %>%
  summarise(
    dice_mean = mean(dice, na.rm = TRUE),
    dice_sd   = sd(dice, na.rm = TRUE),
    n         = n(),
    dice_se   = dice_sd / sqrt(n),
    .groups   = "drop"
  ) %>%
  mutate(
    method = factor(method,
                    levels = c("FDR_p05", "FDR_p01", "FDR_p001", "ClusterFWE_p005", "ClusterFWE_p001", "Top5", "Top10"),
                    labels = c("VoxFDR_p05", "VoxFDR_p01", "VoxFDR_p001", "ClustFWE_p005", "ClustFWE_p001", "Top5%_tstats", "Top10%_tstats")
    ),
    region = factor(region,
                    levels = c("thalamus", "cerebellum", "brainstem"),
                    labels = c("Thalamus", "Cerebellum", "Midbrain")
    ),
    hemi = factor(hemi,
                  levels = c("left", "right"),
                  labels = c("Left", "Right")
    )
  )

dice_points <- dice_results %>%
  mutate(
    method = factor(method,
                    levels = c("FDR_p05", "FDR_p01", "FDR_p001", "ClusterFWE_p005", "ClusterFWE_p001", "Top5", "Top10"),
                    labels = c("VoxFDR_p05", "VoxFDR_p01", "VoxFDR_p001", "ClustFWE_p005", "ClustFWE_p001", "Top5%_tstats", "Top10%_tstats")
    ),
    region = factor(region,
                    levels = c("thalamus", "cerebellum", "brainstem"),
                    labels = c("Thalamus", "Cerebellum", "Midbrain")
    ),
    hemi = factor(hemi,
                  levels = c("left", "right"),
                  labels = c("Left", "Right")
    ),
    fold = factor(fold, levels = c("foldA", "foldB", "foldC", "foldD"))
  )

## 8. Dice Plot: bars + error bars + filled subject/fold points ----
ggplot() +
  geom_col(
    data = dice_bar_plot,
    aes(x = method, y = dice_mean, fill = hemi),
    position = position_dodge(width = 0.7),
    width = 0.6,
    alpha = 0.4
  ) +
  geom_point(
    data = dice_points,
    aes(x = method, y = dice, group = hemi, color = sub_id, shape = fold),
    position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.7),
    size = 1.5,
    stroke = 0.6
  ) +
  geom_errorbar(
    data = dice_bar_plot,
    aes(x = method, ymin = dice_mean - dice_se, ymax = dice_mean + dice_se, group = hemi),
    position = position_dodge(width = 0.7),
    width = 0.2
  ) +
  facet_wrap(~ region, nrow = 1, scales = "free_y") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, 1)) +
  scale_shape_manual(
    name = "Fold",
    values = c("foldA" = 16, "foldB" = 17, "foldC" = 15, "foldD" = 18)
  ) +
  scale_fill_manual(
    name = "Hemisphere",
    values = c("Left" = "grey80", "Right" = "grey20")
  ) +
  labs(
    x = "Method",
    y = "Dice coefficient",
    fill = "Hemisphere",
    color = "Subject",
    shape = "Fold",
    title = "Left and Right Dice by Subregion"
  ) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 15)
  ) +
  guides(color = guide_legend(ncol = 2), shape = guide_legend(ncol = 2))
