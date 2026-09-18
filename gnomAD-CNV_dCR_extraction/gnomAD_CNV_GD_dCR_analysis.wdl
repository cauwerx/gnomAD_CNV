version 1.0

## gnomAD_CNV_GD_dCR_analysis.wdl
##
## Terra/Cromwell translation of the original R script that:
##   1) builds a BED of genomic-disorder (GD) regions from the GD catalog xlsx
##   2) filters the gnomAD-SV manifest to high-confidence, releasable samples
##   3) for each gCNV batch: pulls the batch's dCR matrix, tabix-extracts the
##      GD regions, subsets to HC samples in that batch, floors/ceils dCR
##      values, and computes per-sample/per-GD summary stats (mean/MAD/SD/N)
##   4) merges all per-batch summaries into one final table
##
## Design notes vs. the original script:
##   - The per-batch loop becomes a `scatter`, so batches run in parallel
##     instead of serially.
##   - Batch IDs are NOT a workflow input. They're derived at runtime from
##     PrepareManifest's output (the unique gcnv_batch values among HC,
##     releasable samples) via read_lines() on its batches.txt output,
##     then scattered over directly.
##   - dCR file/index paths are NOT workflow inputs either. Since all ~1000
##     files live at a fixed, predictable location
##     (gs://.../dcr_reheader/BATCH.dcr.bed.gz[.tbi]), each scatter shard
##     builds its own File path by string-interpolating the batch ID into
##     `dcr_bucket_prefix`. WDL/Cromwell coerces a computed String into a
##     File and localizes it like any other input — no manual enumeration
##     of all 1000 files needed.
##   - The final merged table is copied to a fixed GCS destination via a
##     dedicated CopyToFinalDestination task (uses google/cloud-sdk:slim,
##     so the main R/tabix image doesn't need gsutil bundled in).
##   - Everything else (column selection, gsub logic, dCR ceiling/flooring,
##     the findOverlaps + by() summary logic) is preserved as-is in R,
##     just split across tasks instead of one script.

workflow gnomAD_CNV_GD_dCR_analysis {
  input {
    File   gd_catalog_xlsx      # GenomicDisorderRegions_hg38_CAuwerx-....xlsx
    File   manifest             # GNOMAD_V4.6_merged_manifest.txt.gz

    # Fixed location pattern for per-batch dCR files: <prefix>BATCH.dcr.bed.gz[.tbi]
    String dcr_bucket_prefix = "gs://fc-712bf694-df47-4018-8788-bfdc120cdd67/dcr_reheader/"

    # Where the final merged table gets copied to
    String final_output_path = "gs://fc-712bf694-df47-4018-8788-bfdc120cdd67/CAuwerx/GenomicDisorders/extract_dCR/output_data/GD_dCR_summary_gnomAD_CNV_v4.txt.gz"

    String docker         = "us.gcr.io/YOUR_PROJECT/gnomad-gd-dcr:latest"  # see Dockerfile
    String gsutil_docker  = "google/cloud-sdk:slim"
    Int    process_batch_disk_gb = 30
  }

  call PrepareGDCatalog {
    input:
      gd_catalog_xlsx = gd_catalog_xlsx,
      docker = docker
  }

  call PrepareManifest {
    input:
      manifest = manifest,
      docker = docker
  }

  # Batch IDs derived from step 2's output — the unique gcnv_batch values
  # among release/HC/PASS samples. No batch list is supplied by the user.
  Array[String] batch_ids = read_lines(PrepareManifest.batches_list)

  scatter (b in batch_ids) {
    call ProcessBatch {
      input:
        batch_id        = b,
        dcr_file        = dcr_bucket_prefix + b + ".dcr.bed.gz",
        dcr_index       = dcr_bucket_prefix + b + ".dcr.bed.gz.tbi",
        gd_catalog_bed  = PrepareGDCatalog.gd_catalog_bed,
        gd_catalog_full = PrepareGDCatalog.gd_catalog_full,
        hc_samples_tsv  = PrepareManifest.hc_samples_tsv,
        disk_gb         = process_batch_disk_gb,
        docker          = docker
    }
  }

  call MergeBatches {
    input:
      summary_files = ProcessBatch.summary_file,
      docker = docker
  }

  call CopyToFinalDestination {
    input:
      local_file       = MergeBatches.merged_summary,
      destination_path = final_output_path,
      docker           = gsutil_docker
  }

  output {
  }
}

# ------------------------------------------------------------------------
# STEP 1 (original): build the GD region BED + cleaned catalog table
# ------------------------------------------------------------------------
task PrepareGDCatalog {
  input {
    File   gd_catalog_xlsx
    String docker
  }

  command <<<
    set -euo pipefail
    Rscript --vanilla - <<'EOF'
    library(openxlsx)
    library(data.table)

    gd_catalog <- read.xlsx("~{gd_catalog_xlsx}")
    gd_catalog <- gd_catalog[, c(1:4, 9)]
    colnames(gd_catalog) <- c("chr", "start", "end", "gd_id", "SD-flanked")
    gd_catalog$gd_id <- gsub("DUP_", "", gsub("DEL_", "", gd_catalog$gd_id))
    gd_catalog$gd_label <- gsub("GD_", "", gsub("_chr.*", "", gd_catalog$gd_id))
    gd_catalog <- unique(gd_catalog)
    gd_catalog <- gd_catalog[which(gd_catalog$chr != "chrX"), ]

    fwrite(gd_catalog[, c(1:3)], "gd_catalog.bed",
           col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    fwrite(gd_catalog, "gd_catalog_full.tsv",
           col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    EOF
  >>>

  output {
    File gd_catalog_bed  = "gd_catalog.bed"
    File gd_catalog_full = "gd_catalog_full.tsv"
  }

  runtime {
    docker: docker
    memory: "4 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
  }
}

# ------------------------------------------------------------------------
# STEP 2 (original): filter manifest to release/HC/PASS samples
# ------------------------------------------------------------------------
task PrepareManifest {
  input {
    File   manifest
    String docker
  }

  command <<<
    set -euo pipefail
    Rscript --vanilla - <<'EOF'
    library(data.table)

    ped <- as.data.frame(fread("~{manifest}"))
    ped <- ped[which(ped$release == TRUE &
                      ped$PASS_SAMPLE == TRUE &
                      ped$PASS_SAMPLE_MULTICHR == TRUE &
                      ped$PASS_SAMPLE_TERM == TRUE), ]

    rHC_samples <- ped[, c("subject_id", "CNV_ID", "snv_id", "gcnv_batch")]

    fwrite(rHC_samples, "hc_samples.tsv",
           col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    writeLines(as.character(unique(rHC_samples$gcnv_batch)), "batches.txt")
    EOF
  >>>

  output {
    File hc_samples_tsv = "hc_samples.tsv"
    File batches_list   = "batches.txt"   # consumed by the workflow via read_lines() to drive the scatter
  }

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
  }
}

# ------------------------------------------------------------------------
# STEP 3 (original): per-batch tabix extraction + dCR summary (scattered)
# ------------------------------------------------------------------------
task ProcessBatch {
  input {
    String batch_id
    File   dcr_file
    File   dcr_index
    File   gd_catalog_bed
    File   gd_catalog_full
    File   hc_samples_tsv
    Int    disk_gb
    String docker
  }

  command <<<
    set -euo pipefail

    # dcr_index must sit next to dcr_file for tabix to find it
    ln -s "~{dcr_file}" dcr.bed.gz
    ln -s "~{dcr_index}" dcr.bed.gz.tbi

    tabix -h dcr.bed.gz -R "~{gd_catalog_bed}" | gzip > GD_dCR.txt.gz

    Rscript --vanilla - <<'EOF'
    library(data.table)
    library(GenomicRanges)

    batch_id <- "~{batch_id}"

    dCR_matrix <- as.data.frame(fread("GD_dCR.txt.gz"))
    gr_dCR <- GRanges(dCR_matrix[, 1], IRanges(dCR_matrix[, 2], dCR_matrix[, 3]))

    gd_catalog  <- fread("~{gd_catalog_full}")
    gr_gd_catalog <- GRanges(gd_catalog$chr, IRanges(gd_catalog$start, gd_catalog$end))

    hc <- fread("~{hc_samples_tsv}")
    batch_samples <- hc[hc$gcnv_batch == batch_id, ]$CNV_ID

    dCR_matrix <- dCR_matrix[, names(dCR_matrix) %in% batch_samples, drop = FALSE]

    dCR_matrix[dCR_matrix > 5] <- 5
    dCR_matrix[dCR_matrix < 0] <- 0

    overlaps <- findOverlaps(gr_gd_catalog, gr_dCR)

    df_mean  <- by(dCR_matrix[subjectHits(overlaps), , drop = FALSE], queryHits(overlaps), function(x) apply(x, 2, mean))
    df_MAD   <- by(dCR_matrix[subjectHits(overlaps), , drop = FALSE], queryHits(overlaps), function(x) apply(x, 2, mad))
    df_SD    <- by(dCR_matrix[subjectHits(overlaps), , drop = FALSE], queryHits(overlaps), function(x) apply(x, 2, sd))
    df_N_int <- by(dCR_matrix[subjectHits(overlaps), , drop = FALSE], queryHits(overlaps), function(x) apply(x, 2, length))

    gd_names  <- gd_catalog$gd_id[as.numeric(names(df_mean))]
    n_samples <- ncol(dCR_matrix)
    n_gds     <- length(df_mean)

    df_dCR_GD_by_batch <- data.frame(
      sample    = rep(colnames(dCR_matrix), times = n_gds),
      batch     = batch_id,
      GD        = rep(gd_names, each = n_samples),
      dCR_mean  = unlist(df_mean),
      dCR_MAD   = unlist(df_MAD),
      dCR_SD    = unlist(df_SD),
      dCR_N_int = unlist(df_N_int),
      row.names = NULL
    )

    fwrite(df_dCR_GD_by_batch, paste0("GD_dCR_batch_", batch_id, "_summary.txt.gz"),
           col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    EOF
  >>>

  output {
    File summary_file = "GD_dCR_batch_~{batch_id}_summary.txt.gz"
  }

  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
  }
}

# ------------------------------------------------------------------------
# STEP 4 (original): rbind() all per-batch summaries into one table
# ------------------------------------------------------------------------
task MergeBatches {
  input {
    Array[File] summary_files
    String docker
  }

  command <<<
    set -euo pipefail

    # Disable pipefail briefly so zcat receiving SIGPIPE does not fail the script
    set +o pipefail
    zcat ~{summary_files[0]} | head -n 1 > header.txt
    set -o pipefail

    # Stream header + concatenated data directly to compressed output
    {
      cat header.txt
      for f in ~{sep=" " summary_files}; do
        zcat "$f" | tail -n +2
      done
    } | gzip > GD_dCR_summary_gnomAD_CNV_v4.txt.gz
  >>>

  output {
    File merged_summary = "GD_dCR_summary_gnomAD_CNV_v4.txt.gz"
  }

  runtime {
    docker: docker
    memory: "8 GB"
    cpu: 1
    disks: "local-disk 50 HDD"
  }
}

# ------------------------------------------------------------------------
# Copy the merged output to a fixed GCS destination
# ------------------------------------------------------------------------
task CopyToFinalDestination {
  input {
    File   local_file
    String destination_path
    String docker
  }

  command <<<
    set -euo pipefail
    gsutil cp "~{local_file}" "~{destination_path}"
  >>>

  output {
    File final_file = destination_path
  }

  runtime {
    docker: docker
    memory: "16 GB"
    cpu: 2
    disks: "local-disk 100 HDD"
    preemptible: 0
    maxRetries: 2
  }
}
