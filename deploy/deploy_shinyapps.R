#!/usr/bin/env Rscript
# Deploy the app to shinyapps.io.
#
#   Rscript deploy/deploy_shinyapps.R --check     checks everything, uploads nothing
#   Rscript deploy/deploy_shinyapps.R             checks, then uploads
#
# One-time setup (R console): create a free account at https://www.shinyapps.io, open
# Account > Tokens > Show > Copy to clipboard, paste the command it gives you, it looks like
#   rsconnect::setAccountInfo(name = "<account>", token = "<token>", secret = "<secret>")
#
# Options: --allow-demo   deploy the simulated demo database (for a test deployment only)
# Environment: PHM_APP_NAME (default "phallusia-atlas"), PHM_APP_DIR (default: ../app next to this script)

args       <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
allow_demo <- "--allow-demo" %in% args
app_name   <- Sys.getenv("PHM_APP_NAME", "phallusia-atlas")

this_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])
app_dir   <- Sys.getenv("PHM_APP_DIR", file.path(dirname(normalizePath(this_file)), "..", "app"))
app_dir   <- normalizePath(app_dir, mustWork = TRUE)

fail <- function(...) { cat("\nSTOP:", ..., "\n", file = stderr()); quit(status = 1) }
ok   <- function(...) cat("ok  ", ..., "\n")

# 1. packages ---------------------------------------------------------------------------------------
need <- c("rsconnect", "shiny", "bslib", "DBI", "RSQLite", "DT", "plotly")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) fail("install these R packages first: install.packages(c(",
                       paste0('"', miss, '"', collapse = ", "), "))")
ok("packages installed (rsconnect", as.character(packageVersion("rsconnect")), ")")

# 2. database ---------------------------------------------------------------------------------------
db_real <- file.path(app_dir, "data", "phallusia.sqlite")
db_demo <- file.path(app_dir, "data", "demo.sqlite")
use_demo <- !file.exists(db_real)
if (use_demo && !allow_demo)
  fail("no real database at ", db_real, ".\n  Build it first (database/build_database.py, see README) ",
       "or add --allow-demo for a test deployment with the simulated demo data.")
db <- if (use_demo) db_demo else db_real
con <- DBI::dbConnect(RSQLite::SQLite(), db, flags = RSQLite::SQLITE_RO)
meta <- DBI::dbGetQuery(con, "SELECT key, value FROM meta")
n_genes <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM genes")$n
n_expr  <- DBI::dbGetQuery(con, "SELECT count(*) AS n FROM expression")$n
DBI::dbDisconnect(con)
is_demo <- identical(meta$value[meta$key == "is_demo"], "1")
if (is_demo && !allow_demo) fail(db, " is marked as demo data (meta is_demo = 1). Build the real database or use --allow-demo.")
ok(sprintf("database %s: %d genes, %s expression values, %.0f MB%s", basename(db), n_genes,
           format(n_expr, big.mark = ","), file.size(db) / 1e6, if (is_demo) " (DEMO DATA)" else ""))

# 3. files in the bundle ----------------------------------------------------------------------------
all_files <- list.files(app_dir, recursive = TRUE, all.files = FALSE)
all_files <- all_files[!grepl("^rsconnect/", all_files)]
drop <- if (use_demo) character() else "data/demo.sqlite"        # the real database replaces the demo one
app_files <- setdiff(all_files, drop)
if (!use_demo) app_files <- unique(c(app_files, "data/phallusia.sqlite"))
size_mb <- sum(file.size(file.path(app_dir, app_files))) / 1e6
cat(sprintf("     bundle: %d files, %.0f MB\n", length(app_files), size_mb))
for (f in app_files) cat(sprintf("       %-30s %8.1f MB\n", f, file.size(file.path(app_dir, f)) / 1e6))
if (size_mb > 900) fail("bundle is ", round(size_mb), " MB; the shinyapps.io Free and Starter plans allow 1 GB (1000 MB) per bundle.")
ok("bundle size is within the 1 GB limit of the shinyapps.io Free plan")

# 4. R packages shinyapps.io will install -----------------------------------------------------------
deps <- rsconnect::appDependencies(appDir = app_dir, appFiles = app_files)
want <- c("shiny", "bslib", "DBI", "RSQLite", "DT", "plotly")
pk <- deps[[intersect(c("Package", "package"), names(deps))[1]]]
if (!all(want %in% pk)) fail("rsconnect did not detect: ", paste(setdiff(want, pk), collapse = ", "),
                                       " (app.R must call library() for them)")
ok(length(unique(pk)), "R packages detected, including all of:", paste(want, collapse = ", "))

# 5. account -----------------------------------------------------------------------------------------
acc <- rsconnect::accounts()
acc <- acc[acc$server == "shinyapps.io", , drop = FALSE]
if (nrow(acc) == 0) {
  msg <- "no shinyapps.io account configured: run rsconnect::setAccountInfo(...) once (see the top of this file)"
  if (check_only) cat("WARN", msg, "\n") else fail(msg)
} else {
  ok("shinyapps.io account:", acc$name[1])
}

if (check_only) { cat("\nCheck finished, nothing was uploaded.\n"); quit(status = 0) }

# 6. deploy ------------------------------------------------------------------------------------------
cat("\nUploading to shinyapps.io as '", app_name, "' (this can take several minutes for the first upload) ...\n", sep = "")
rsconnect::deployApp(
  appDir      = app_dir,
  appFiles    = app_files,
  appName     = app_name,
  appTitle    = "Phallusia mammillata expression atlas",
  account     = acc$name[1],
  server      = "shinyapps.io",
  forceUpdate = TRUE,
  launch.browser = FALSE,
  lint        = FALSE
)
cat("\nDone. If the page shows an error, open the app's Logs in the shinyapps.io dashboard.\n")
