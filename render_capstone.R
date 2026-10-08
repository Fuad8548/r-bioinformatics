# ==============================================================================
# Module 10: Quarto Report Rendering Driver
# File: scripts/01_render_capstone.R
# Description: Programmatically renders parameterized Quarto reports into HTML
# ==============================================================================

suppressPackageStartupMessages({
    library(quarto)
})

cat("\n--- Rendering Capstone Quarto Report ---\n")

# Define target paths
qmd_file <- "reports/capstone_report.qmd"
output_file <- "capstone_final_report.html"

# Render parameterized report
quarto::quarto_render(
    input = qmd_file,
    output_file = output_file,
    execute_params = list(
        project_title = "Pan-Cancer Transcriptomic Biomarker Discovery",
        padj_cutoff   = 0.01,
        lfc_cutoff    = 1.5
    )
)

cat("\nReport successfully generated and saved to 'reports/", output_file, "'!\n", sep = "")
