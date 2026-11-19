# Phallusia mammillata expression atlas -- Shiny app
# Reads one SQLite file made by database/build_database.py (see DB_PATH below).
# Packages: shiny, bslib, DBI, RSQLite, DT, plotly (all on CRAN). Tested versions: renv.lock in this folder.

library(shiny)
library(bslib)
library(DBI)
library(RSQLite)
library(DT)
library(plotly)

APP_VERSION <- "1.0.0"   # also in CHANGELOG.md; bump for each release

# Database: env PHALLUSIA_DB, else data/phallusia.sqlite, else the simulated demo database data/demo.sqlite
DB_PATH <- Sys.getenv("PHALLUSIA_DB", "")
if (!nzchar(DB_PATH)) DB_PATH <- if (file.exists("data/phallusia.sqlite")) "data/phallusia.sqlite" else "data/demo.sqlite"
if (!file.exists(DB_PATH)) {
  stop("Database not found: ", DB_PATH, ". Build it with database/build_database.py or set PHALLUSIA_DB.")
}
con <- dbConnect(SQLite(), DB_PATH, flags = SQLITE_RO)
onStop(function() dbDisconnect(con))

q1 <- function(sql, ...) dbGetQuery(con, sql, params = list(...))
meta <- local({ m <- dbGetQuery(con, "SELECT key, value FROM meta"); setNames(m$value, m$key) })
samples <- dbGetQuery(con, "SELECT * FROM samples ORDER BY dev_order, replicate")
stage_levels <- unique(samples$stage_label[order(samples$dev_order)])
stage_orders <- sort(unique(samples$dev_order))
IS_DEMO <- identical(unname(meta["is_demo"]), "1")

# ---------- loaded once at start-up: stage means, labels, maternal / zygotic class -----------------------------
# (about 1 s and 10 MB for the real database)
gene_tab <- dbGetQuery(con, "SELECT g.gene_id, g.description, g.biotype, e.preferred_name AS egg_name, e.description AS egg_desc
                             FROM genes g LEFT JOIN eggnog e USING(gene_id)")
gene_tab$egg_name[is.na(gene_tab$egg_name)] <- ""
gene_tab$egg_desc[is.na(gene_tab$egg_desc)] <- ""
rownames(gene_tab) <- gene_tab$gene_id

stage_mean <- local({
  d <- dbGetQuery(con, "SELECT e.gene_id, sm.dev_order, AVG(e.tpm) AS tpm FROM expression e JOIN samples sm USING(sample_idx)
                        GROUP BY e.gene_id, sm.dev_order")
  ids <- unique(d$gene_id)
  m <- matrix(0, length(ids), length(stage_orders), dimnames = list(ids, stage_levels))
  m[cbind(match(d$gene_id, ids), match(d$dev_order, stage_orders))] <- d$tpm
  m
})                                                                  # genes x stages, mean TPM of the replicates
log_mean <- log2(stage_mean + 1)
is_mt <- function(gid) startsWith(gid, "MT-")              # mitochondrial genome (MITOS2 genes)
MT_OWN_SCALE <- grepl("separate TPM scales", if (is.na(meta["tpm_note"])) "" else meta[["tpm_note"]])
MT_CLASS <- "mitochondrial (not classified)"

# Pearson correlation of stage-mean log2(TPM + 1) profiles: rows centred and scaled to length 1, so r = Z %*% z
prof_z <- local({
  z <- log_mean - rowMeans(log_mean)
  n <- sqrt(rowSums(z^2))
  z / ifelse(n > 0, n, 1)
})
prof_sd <- apply(log_mean, 1, sd)

# Maternal / zygotic class from the stage means: present in the unfertilized egg or not (explained on the Stages tab),
# plus the stage of peak mean TPM
ORIGIN_RULES <- list(min_tpm = 1)
origin_tab <- local({
  r <- ORIGIN_RULES; m <- stage_mean
  egg <- m[, 1]; peak <- apply(m, 1, max)
  cls <- ifelse(peak < r$min_tpm, "not expressed", ifelse(egg < r$min_tpm, "zygotic", "maternal"))
  peak_stage <- ifelse(peak < r$min_tpm, NA_character_, stage_levels[apply(m, 1, which.max)])
  mt <- is_mt(rownames(m)); cls[mt] <- MT_CLASS   # not transcribed from the nuclear (zygotic) genome
  data.frame(gene_id = rownames(m), origin = cls, peak_stage = peak_stage, row.names = rownames(m), stringsAsFactors = FALSE)
})
ORIGIN_LEVELS <- c("maternal", "zygotic", "not expressed", MT_CLASS)
origin_of <- function(gid) { o <- origin_tab[gid, "origin"]; ifelse(is.na(o), "no expression data", o) }

# ---------- helpers ------------------------------------------------------------------------------------------
unz <- function(b) if (is.null(b) || !length(b)) "" else rawToChar(memDecompress(b, type = "gzip"))
wrap_seq <- function(s, w = 60) {
  n <- nchar(s)
  if (is.na(n) || n == 0) return("")
  starts <- seq(1, n, by = w)
  paste(substring(s, starts, pmin(starts + w - 1, n)), collapse = "\n")
}
gc_percent <- function(s) {
  if (!nzchar(s)) return(NA_real_)
  100 * sum(strsplit(toupper(s), "")[[1]] %in% c("G", "C")) / nchar(s)
}
esc_like <- function(x) gsub("([%_\\\\])", "\\\\\\1", x)
loc_string <- function(seqid, start, end, strand) {
  paste0(seqid, ":", format(start, big.mark = ",", trim = TRUE), "-", format(end, big.mark = ",", trim = TRUE), " (", strand, ")")
}
# eggNOG text for a gene: the preferred name, else the eggNOG description ("" if none)
egg_text <- function(gid) {
  g <- gene_tab[gid, , drop = FALSE]
  ifelse(!is.na(g$egg_name) & nzchar(g$egg_name), g$egg_name, ifelse(is.na(g$egg_desc), "", g$egg_desc))
}
# Description with the eggNOG annotation added when RefSeq only says "uncharacterized"
nice_desc <- function(gid) {
  d <- gene_tab[gid, "description"]; d[is.na(d)] <- ""
  e <- egg_text(gid)
  add <- grepl("^uncharacterized", d) & nzchar(e)
  d[add] <- paste0(d[add], " [eggNOG: ", e[add], "]")
  d
}
# Short label for heatmaps / legends: gene id + (eggNOG text for uncharacterized genes, else the description)
short_label <- function(gid, width = 38) {
  d <- gene_tab[gid, "description"]; d[is.na(d)] <- ""
  e <- egg_text(gid)
  d <- ifelse(grepl("^uncharacterized", d) & nzchar(e), paste0("eggNOG: ", e), d)
  paste0(gid, "  ", ifelse(nchar(d) > width, paste0(substr(d, 1, width - 1), "…"), d))
}
# Gene id as a link that opens the gene in the app (middle-click / copy link also works: ?gene=<id>)
gene_link <- function(gid) {
  sprintf('<a href="?gene=%s" class="gene-link" data-gene="%s">%s</a>',
          utils::URLencode(gid, reserved = TRUE), htmltools::htmlEscape(gid, attribute = TRUE), htmltools::htmlEscape(gid))
}

# search: exact identifiers first, then every word must occur in the symbol / description / eggNOG / Pfam / GO text
search_genes <- function(query, max_rows = 5000) {
  query <- trimws(tolower(query))
  if (!nzchar(query)) return(NULL)
  ids <- q1("SELECT DISTINCT gene_id FROM identifiers WHERE identifier_l = ?", query)$gene_id
  terms <- unique(strsplit(query, "[[:space:],;]+")[[1]])
  terms <- terms[nzchar(terms)]
  pats <- paste0("%", esc_like(terms), "%")
  hay <- "(s.symbol_l || ' ' || s.desc_l || ' ' || s.other_l)"
  where <- paste(rep(paste0(hay, " LIKE ? ESCAPE '\\'"), length(pats)), collapse = " AND ")
  sql <- paste0(
    "SELECT g.gene_id, g.symbol, g.description, g.biotype, g.seqid, g.start, g.end, g.strand, g.peak_stage, g.peak_tpm, ",
    "e.preferred_name AS egg_name, e.pfams, s.symbol_l, s.desc_l ",
    "FROM search s JOIN genes g USING(gene_id) LEFT JOIN eggnog e USING(gene_id) WHERE ", where, " LIMIT 50000")
  d <- dbGetQuery(con, sql, params = as.list(pats))
  if (!nrow(d)) return(d)
  score <- Reduce(`+`, lapply(terms, function(t) 4 * grepl(t, d$symbol_l, fixed = TRUE) + 3 * grepl(t, d$desc_l, fixed = TRUE) + 1))
  score <- score + ifelse(d$gene_id %in% ids, 100, 0)
  d <- d[order(-score, -ifelse(is.na(d$peak_tpm), 0, d$peak_tpm)), ]
  d$origin <- origin_of(d$gene_id)
  total <- nrow(d)
  d <- head(d, max_rows)
  attr(d, "total") <- total
  d
}
get_gene <- function(gid) q1("SELECT * FROM genes WHERE gene_id = ?", gid)
get_tx <- function(gid) q1("SELECT * FROM transcripts WHERE gene_id = ? ORDER BY is_primary DESC, transcript_id", gid)
get_egg <- function(gid) q1("SELECT * FROM eggnog WHERE gene_id = ?", gid)
get_expr <- function(gid) {
  q1("SELECT sm.sample, sm.stage_label, sm.dev_order, sm.replicate, e.tpm, e.counts
      FROM expression e JOIN samples sm USING(sample_idx) WHERE e.gene_id = ? ORDER BY sm.dev_order, sm.replicate", gid)
}
resolve_ids <- function(txt) {
  tokens <- unique(trimws(unlist(strsplit(txt, "[,;[:space:]]+"))))
  tokens <- tokens[nzchar(tokens)]
  if (!length(tokens)) return(list(found = character(), missing = character()))
  r <- dbGetQuery(con, paste0("SELECT identifier_l, gene_id FROM identifiers WHERE identifier_l IN (", paste(rep("?", length(tokens)), collapse = ","), ")"),
                  params = as.list(tolower(tokens)))
  r <- r[order(match(r$identifier_l, tolower(tokens))), , drop = FALSE]   # keep the order the user typed
  list(found = unique(r$gene_id), missing = tokens[!(tolower(tokens) %in% r$identifier_l)])
}
wide_table <- function(value = c("tpm", "counts")) {
  value <- match.arg(value)
  cols <- paste(sprintf("MAX(CASE WHEN sample_idx = %d THEN %s END) AS \"%s\"", samples$sample_idx, value, samples$sample), collapse = ", ")
  dbGetQuery(con, paste0("SELECT gene_id, ", cols, " FROM expression GROUP BY gene_id"))
}
# genes with the most similar stage profile (Pearson r of log2(TPM + 1) stage means)
similar_genes <- function(gid, n = 25, min_tpm = 1) {
  if (!gid %in% rownames(prof_z)) return(NULL)
  ok <- apply(stage_mean, 1, max) >= min_tpm & prof_sd > 0
  ok <- ok & is_mt(rownames(prof_z)) == is_mt(gid)       # nuclear vs mitochondrial: separate TPM scales, separate lists
  ok[gid] <- FALSE
  r <- drop(prof_z[ok, , drop = FALSE] %*% prof_z[gid, ])
  top <- head(order(-r), n)
  data.frame(gene_id = names(r)[top], r = r[top], stringsAsFactors = FALSE)
}

# Developmental stages: short descriptions (Hotta et al. 2007 staging of Ciona, used for ascidians in general)
STAGE_INFO <- c(
  "unfertilized egg" = "Mature egg before fertilization. Everything here is maternal RNA.",
  "St3"  = "Cleavage: 4 blastomeres.",
  "St4"  = "Cleavage: 8 blastomeres; the animal-vegetal and anterior-posterior axes are set.",
  "St5b" = "Cleavage: 16 cells. Zygotic transcription of the first genes has started.",
  "St6b" = "Cleavage: 32 cells. Many cell fates (endoderm, notochord, muscle precursors) are being specified.",
  "St8"  = "Cleavage: 64 cells; the last stage before gastrulation starts.",
  "St11" = "Early gastrula: endoderm precursors invaginate.",
  "St12" = "Mid gastrula: the blastopore narrows; mesoderm moves inside.",
  "St14" = "Early neurula: the neural plate is visible.",
  "St15" = "Mid neurula: neural folds rise and the neural tube starts to close.",
  "St17" = "Initial tailbud I: the embryo starts to elongate; the tail bud appears.",
  "St19" = "Early tailbud I: the tail grows; notochord cells line up in a single row.",
  "St22" = "Mid tailbud II: the tail bends around the trunk; tissues differentiate.",
  "St24" = "Late tailbud II: sensory organs (otolith, ocellus) become pigmented.",
  "St26" = "Hatching larva: the swimming tadpole larva leaves the chorion."
)

# ---------- UI -----------------------------------------------------------------------------------------------
copy_js <- "
function copyEl(id){var t=document.getElementById(id).innerText;
 if(navigator.clipboard&&window.isSecureContext){navigator.clipboard.writeText(t);}
 else{var a=document.createElement('textarea');a.value=t;document.body.appendChild(a);a.select();document.execCommand('copy');a.remove();}}
$(document).on('keyup','#q',function(e){if(e.key==='Enter'){$('#do_search').click();}});
$(document).on('click','a.gene-link',function(e){
  if(e.ctrlKey||e.metaKey||e.shiftKey||e.button===1){return;}
  e.preventDefault(); e.stopPropagation();
  Shiny.setInputValue('open_gene', $(this).data('gene'), {priority:'event'});});
"
css <- "
.seqbox{font-family:ui-monospace,Menlo,Consolas,monospace;font-size:12.5px;white-space:pre-wrap;word-break:break-all;
        max-height:420px;overflow:auto;background:var(--bs-tertiary-bg,#f8f9fa);padding:.75rem;border:1px solid var(--bs-border-color,#dee2e6);border-radius:.375rem;}
.kv td:first-child{color:var(--bs-secondary-color,#6c757d);white-space:nowrap;padding-right:1.25rem;vertical-align:top}
.chip{display:inline-block;margin:0 .25rem .25rem 0;padding:.1rem .55rem;border-radius:1rem;font-size:.8rem;background:var(--bs-secondary-bg,#e9ecef)}
.chip-origin{background:#e3f1f5;color:#155a70}
"

ui <- page_navbar(
  id = "nav",
  fillable = FALSE,      # bslib >= 0.5 squeezes every page into the window height otherwise (overlapping controls, squashed plots)
  title = "Phallusia mammillata expression atlas",
  theme = bs_theme(version = 5, primary = "#1b6f8a", base_font = c("system-ui", "Segoe UI", "Helvetica", "Arial", "sans-serif")),
  header = tagList(
    tags$head(tags$title("Phallusia mammillata expression atlas"), tags$style(HTML(css)), tags$script(HTML(copy_js))),
    if (IS_DEMO) div(class = "alert alert-danger text-center mb-0 rounded-0 py-2",
                     tags$strong("DEMO DATA"), " - simulated sequences and expression values. This is not real biology.")
  ),
  nav_panel("Search", value = "search",
    card(
      div(class = "d-flex gap-2 align-items-start",
          div(style = "flex:1", textInput("q", NULL, width = "100%",
              placeholder = "Gene name, description, protein domain or ID: e.g. cadherin, T-box, Notch, LOC..., XP_...")),
          actionButton("do_search", "Search", class = "btn-primary")),
      div(class = "small text-muted mb-2", "Try: ", actionLink("ex1", "cadherin"), " | ", actionLink("ex2", "homeobox"), " | ",
          actionLink("ex3", "Cadherin_C_term"), " | ", actionLink("ex4", "GO:0005509"), " | ", actionLink("ex5", "mitochondrial")),
      div(class = "d-flex flex-wrap gap-3 align-items-end",
          selectInput("origin_filter", "Maternal / zygotic", c("All genes" = "", ORIGIN_LEVELS), width = "230px"),
          div(class = "small text-muted pb-3", "Click a gene ID to open it. ", actionLink("to_stages", "What do the classes mean?"))),
      uiOutput("search_msg"),
      DTOutput("results", fill = FALSE))
  ),
  nav_panel("Gene", value = "gene",
    conditionalPanel("!output.has_gene", card(card_body("No gene selected. Use the Search tab and click a gene."))),
    conditionalPanel("output.has_gene",
      div(class = "mb-2", actionLink("back", "< back to search results")),
      uiOutput("gene_header"),
      navset_card_tab(id = "gene_tabs",
        nav_panel("Expression",
          uiOutput("expr_note"),
          conditionalPanel("output.has_expr",
          layout_columns(col_widths = c(3, 4, 2, 3),
            radioButtons("expr_val", "Value", c("TPM" = "tpm", "Read counts" = "counts"), inline = TRUE),
            radioButtons("expr_type", "Plot", c("Mean + replicates" = "line", "Replicate series" = "series", "Bars (mean)" = "bar"), inline = TRUE),
            checkboxInput("expr_log", "log2(value + 1)", FALSE),
            div(class = "pt-4", downloadButton("dl_expr", "Download values", class = "btn-sm"),
                actionButton("add_cmp", "Add to Compare", class = "btn-sm btn-outline-primary"))),
          uiOutput("expr_mt_note"),
          plotlyOutput("expr_plot", height = "440px"),
          DTOutput("expr_table", fill = FALSE))),
        nav_panel("Similar profiles",
          layout_columns(col_widths = c(3, 3, 6),
            numericInput("sim_n", "Number of genes", 25, min = 5, max = 200, step = 5),
            selectInput("sim_min", "Only genes with peak TPM at least", c("1" = 1, "5" = 5, "10" = 10, "50" = 50), selected = 5),
            div(class = "pt-4", actionButton("sim_to_cmp", "Show this gene + these genes in Compare", class = "btn-sm btn-outline-primary"),
                downloadButton("dl_sim", "Download table", class = "btn-sm"))),
          uiOutput("sim_note"),
          plotlyOutput("sim_plot", height = "360px"),
          DTOutput("sim_table", fill = FALSE)),
        nav_panel("Sequence",
          layout_columns(col_widths = c(6, 6),
            selectInput("tx_sel", "Transcript / isoform", choices = NULL, width = "100%"),
            radioButtons("seq_kind", "Sequence", inline = TRUE,
                         c("Genomic (gene)" = "gene", "mRNA" = "mrna", "CDS" = "cds", "Protein" = "protein"))),
          uiOutput("seq_info"),
          div(class = "d-flex gap-2 mb-2",
              tags$button("Copy", class = "btn btn-outline-secondary btn-sm", type = "button", onclick = "copyEl('seq_text')"),
              downloadButton("dl_seq", "Download FASTA", class = "btn-sm")),
          verbatimTextOutput("seq_text", placeholder = TRUE) |> tagAppendAttributes(class = "seqbox")),
        nav_panel("Annotation", uiOutput("annotation"))
      )
    )
  ),
  nav_panel("Compare", value = "compare",
    layout_sidebar(
      sidebar = sidebar(width = 340,
        textAreaInput("cmp_genes", "Genes (IDs or symbols, one per line or comma-separated)", rows = 5, placeholder = "e.g.\nCdh2\nLOC100003"),
        textInput("cmp_term", "...and/or all genes matching a search term", placeholder = "e.g. cadherin"),
        numericInput("cmp_max", "Max genes from the search term", 40, min = 1, max = 200),
        selectInput("cmp_scale", "Colour scale", c("log2(TPM + 1)" = "log2", "Row z-score of log2(TPM + 1)" = "z", "TPM" = "tpm")),
        radioButtons("cmp_order", "Order genes",
                     c("Cluster by pattern (shape over stages)" = "pattern",
                       "Cluster by expression level" = "level",
                       "As entered" = "none"), selected = "pattern"),
        actionButton("cmp_go", "Plot heatmap", class = "btn-primary"),
        downloadButton("dl_cmp", "Download table", class = "btn-sm")),
      card(uiOutput("cmp_msg"), uiOutput("cmp_plot_ui"))
    )
  ),
  nav_panel("Stages", value = "stages",
    card(card_header("Developmental stages in this dataset"),
         p(class = "small text-muted", "Stage names and numbers follow the ascidian staging of Hotta et al. (2007, Developmental Dynamics 236:1790),",
           " defined for Ciona and used for Phallusia as well. Absolute timing depends on species and temperature; see ",
           tags$a(href = "https://www.aniseed.fr", target = "_blank", rel = "noopener", "ANISEED"), " for staged images."),
         DTOutput("stage_table", fill = FALSE)),
    card(card_header("Maternal / zygotic classes"), uiOutput("origin_text"), DTOutput("origin_counts", fill = FALSE))
  ),
  nav_panel("Data & about", value = "about",
    layout_columns(col_widths = c(6, 6),
      card(card_header("About this resource"), uiOutput("about_text")),
      card(card_header("Downloads"),
           p("Full tables for all genes (TSV)."),
           div(class = "d-flex flex-wrap gap-2",
               downloadButton("dl_tpm", "TPM matrix"), downloadButton("dl_counts", "Counts matrix"), downloadButton("dl_annot", "Gene annotation"))
      )),
    card(card_header("Samples"), DTOutput("samples_table", fill = FALSE))
  )
)

# ---------- server -------------------------------------------------------------------------------------------
server <- function(input, output, session) {
  selected <- reactiveVal(NULL)
  results <- reactiveVal(NULL)
  proxy <- dataTableProxy("results")

  run_search <- function(txt) {
    updateTextInput(session, "q", value = txt)
    results(search_genes(txt))
  }
  observeEvent(input$do_search, results(search_genes(input$q)))
  observeEvent(input$ex1, run_search("cadherin"))
  observeEvent(input$ex2, run_search("homeobox"))
  observeEvent(input$ex3, run_search("Cadherin_C_term"))
  observeEvent(input$ex4, run_search("GO:0005509"))
  observeEvent(input$ex5, run_search("mitochondrial"))
  observeEvent(input$to_stages, nav_select("nav", "stages"))

  shown_results <- reactive({
    r <- results()
    if (is.null(r) || !nrow(r) || !nzchar(input$origin_filter)) return(r)
    r2 <- r[r$origin == input$origin_filter, , drop = FALSE]
    attr(r2, "total") <- nrow(r2); attr(r2, "filtered_from") <- nrow(r)
    r2
  })
  output$search_msg <- renderUI({
    r <- shown_results()
    if (is.null(r)) return(div(class = "text-muted", sprintf("%s genes in the database. Search by name, description, eggNOG name, Pfam domain, GO term or identifier.", meta[["n_genes"]])))
    if (!nrow(r)) {
      if (!is.null(attr(r, "filtered_from"))) return(div(class = "text-warning", sprintf("None of the %d matching genes is '%s'. Choose 'All genes' to see them.", attr(r, "filtered_from"), input$origin_filter)))
      return(div(class = "text-warning", "No gene matches. Try a shorter word (e.g. 'cadher') or a different spelling."))
    }
    total <- attr(r, "total")
    filt <- if (!is.null(attr(r, "filtered_from"))) sprintf(" ('%s' only, of %d matches)", input$origin_filter, attr(r, "filtered_from")) else ""
    if (!is.null(total) && total > nrow(r))
      return(div(class = "text-muted", sprintf("%s genes match; showing the best %s. Add more words to narrow the search.",
                                               format(total, big.mark = ","), format(nrow(r), big.mark = ","))))
    div(class = "text-muted", sprintf("%d genes found%s.", nrow(r), filt))
  })
  output$results <- renderDT({
    r <- shown_results(); req(r, nrow(r) > 0)
    show <- data.frame(Gene = gene_link(r$gene_id), Description = nice_desc(r$gene_id), `eggNOG name` = r$egg_name, Type = r$biotype,
                       Location = loc_string(r$seqid, r$start, r$end, r$strand), `Peak stage` = r$peak_stage,
                       `Peak TPM` = round(r$peak_tpm, 1), `Maternal / zygotic` = r$origin, check.names = FALSE, stringsAsFactors = FALSE)
    datatable(fillContainer = FALSE, data = show, selection = "single", rownames = FALSE, class = "compact stripe hover", escape = -1,
              options = list(pageLength = 15, scrollX = TRUE, dom = "frtip"))
  })

  open_gene <- function(gid) {
    if (!nrow(get_gene(gid))) { showNotification(paste("Unknown gene:", gid), type = "error"); return() }
    selected(gid)
    updateQueryString(paste0("?gene=", utils::URLencode(gid, reserved = TRUE)), mode = "replace")
    nav_select("nav", "gene")
  }
  observeEvent(input$open_gene, open_gene(input$open_gene))
  observeEvent(input$results_rows_selected, {
    r <- shown_results(); i <- input$results_rows_selected
    if (length(i) == 1 && !is.null(r)) { open_gene(r$gene_id[i]); selectRows(proxy, NULL) }
  })
  observeEvent(input$back, { nav_select("nav", "search"); updateQueryString("?", mode = "replace") })
  observeEvent(session$clientData$url_search, {
    g <- parseQueryString(session$clientData$url_search)$gene
    if (!is.null(g) && nzchar(g)) open_gene(g)
  }, once = TRUE)

  output$has_gene <- reactive(!is.null(selected()))
  outputOptions(output, "has_gene", suspendWhenHidden = FALSE)
  output$has_expr <- reactive(!is.null(selected()) && selected() %in% rownames(stage_mean))
  outputOptions(output, "has_expr", suspendWhenHidden = FALSE)

  gene <- reactive({ req(selected()); get_gene(selected()) })
  txs <- reactive({ req(selected()); get_tx(selected()) })
  expr <- reactive({ req(selected()); get_expr(selected()) })

  output$gene_header <- renderUI({
    g <- gene(); req(nrow(g) == 1)
    links <- if (nzchar(g$ncbi_gene_id)) tags$a(href = paste0("https://www.ncbi.nlm.nih.gov/gene/", g$ncbi_gene_id), target = "_blank", rel = "noopener", "NCBI Gene") else NULL
    eg <- gene_tab[g$gene_id, ]
    egg_line <- if (nzchar(eg$egg_name) || nzchar(eg$egg_desc))
      div(class = "text-muted", "eggNOG: ", tags$strong(eg$egg_name), if (nzchar(eg$egg_desc)) paste0(if (nzchar(eg$egg_name)) " – " else "", eg$egg_desc))
    o <- origin_tab[g$gene_id, ]
    origin <- if (!is.na(o$origin) && !is_mt(g$gene_id)) tags$span(class = "chip chip-origin", title = "Maternal / zygotic class (see the Stages tab)",
                                              paste0(o$origin, if (!is.na(o$peak_stage)) paste0(", peak at ", o$peak_stage) else ""))
    div(class = "mb-3",
        h3(class = "mb-0", g$gene_id), div(class = "lead", g$description), egg_line,
        div(class = "text-muted mt-1", tags$span(class = "chip", g$biotype), tags$span(class = "chip", loc_string(g$seqid, g$start, g$end, g$strand)),
            tags$span(class = "chip", paste(g$n_transcripts, "transcript(s)")), origin, links))
  })

  # ---- expression
  output$expr_note <- renderUI({
    if (nrow(expr())) return(NULL)
    g <- gene()
    why <- if (identical(g$biotype, "pseudogene"))
      "RefSeq pseudogenes have no transcript records in the annotation, so the pipeline had nothing to quantify reads against."
    else "The gene has no transcript in the annotation that the pipeline could quantify."
    div(class = "alert alert-secondary my-2", tags$strong("No expression values for this gene. "), why)
  })
  output$expr_mt_note <- renderUI({
    gid <- selected(); req(gid)
    if (!is_mt(gid) || !MT_OWN_SCALE) return(NULL)
    div(class = "alert alert-info small my-2", tags$strong("Mitochondrial gene. "),
        "Its TPM is on the mitochondrial scale: per sample, the mitochondrial genes except the two rRNAs sum to 10^6, separately from the nuclear genes. ",
        if (startsWith(gid, "MT-RNR")) "This rRNA uses the same factor, so its values are far above 10^6 (it is much more abundant than all MT mRNAs together). ",
        "Compare it with other MT genes, not with nuclear genes. It has no maternal/zygotic class.")
  })
  expr_values <- reactive({
    d <- expr(); req(nrow(d) > 0)
    d$y <- if (input$expr_val == "tpm") d$tpm else d$counts
    if (isTRUE(input$expr_log)) d$y <- log2(d$y + 1)
    d$stage_label <- factor(d$stage_label, levels = stage_levels)
    d
  })
  output$expr_plot <- renderPlotly({
    d <- expr_values()
    tpm_lab <- if (MT_OWN_SCALE && is_mt(selected())) "TPM (mitochondrial scale)" else "TPM"
    ylab <- paste0(if (isTRUE(input$expr_log)) "log2(" else "", if (input$expr_val == "tpm") tpm_lab else "counts", if (isTRUE(input$expr_log)) " + 1)" else "")
    m <- aggregate(y ~ stage_label, d, mean)
    m <- m[order(m$stage_label), ]
    p <- plot_ly()
    if (input$expr_type == "bar") {
      p <- p |> add_bars(data = m, x = ~stage_label, y = ~y, name = "mean", marker = list(color = "#1b6f8a"),
                         hovertemplate = "%{x}<br>mean %{y:.2f}<extra></extra>")
      p <- p |> add_markers(data = d, x = ~stage_label, y = ~y, name = "replicates", text = ~sample, marker = list(color = "#222", size = 7),
                            hovertemplate = "%{text}<br>%{y:.2f}<extra></extra>")
    } else if (input$expr_type == "series") {
      cols <- c("#1b6f8a", "#e07b00", "#6a3d9a", "#33a02c")
      for (i in seq_along(sort(unique(d$replicate)))) {
        r <- sort(unique(d$replicate))[i]
        di <- d[d$replicate == r, ]; di <- di[order(di$stage_label), ]
        p <- p |> add_trace(data = di, x = ~stage_label, y = ~y, type = "scatter", mode = "lines+markers", name = paste("series", r),
                            text = ~sample, line = list(color = cols[(i - 1) %% 4 + 1], width = 2), marker = list(color = cols[(i - 1) %% 4 + 1], size = 7),
                            hovertemplate = "%{text}<br>%{y:.2f}<extra></extra>")
      }
    } else {
      p <- p |> add_trace(data = m, x = ~stage_label, y = ~y, type = "scatter", mode = "lines", name = "mean", line = list(color = "#1b6f8a", width = 3),
                          hovertemplate = "%{x}<br>mean %{y:.2f}<extra></extra>")
      p <- p |> add_markers(data = d, x = ~stage_label, y = ~y, name = "replicates", text = ~sample, marker = list(color = "#e07b00", size = 8),
                            hovertemplate = "%{text}<br>%{y:.2f}<extra></extra>")
    }
    p |> layout(xaxis = list(title = "", categoryorder = "array", categoryarray = stage_levels, tickangle = -45),
                yaxis = list(title = ylab, rangemode = "tozero"), margin = list(b = 125), legend = list(orientation = "h", y = 1.08)) |>
      config(displaylogo = FALSE)
  })
  output$expr_table <- renderDT({
    d <- expr_values()
    datatable(fillContainer = FALSE, data = data.frame(Sample = d$sample, Stage = as.character(d$stage_label), Replicate = d$replicate, TPM = round(d$tpm, 3), Counts = round(d$counts, 1)),
              rownames = FALSE, class = "compact stripe", options = list(pageLength = 10, dom = "tip"))
  })
  output$dl_expr <- downloadHandler(
    filename = function() paste0(selected(), "_expression.csv"),
    content = function(file) { d <- expr(); write.csv(d[, c("sample", "stage_label", "replicate", "tpm", "counts")], file, row.names = FALSE) })

  # ---- similar profiles (co-expression across development)
  sim <- reactive({
    req(selected())
    n <- if (isTruthy(input$sim_n) && !is.na(input$sim_n)) max(5, min(200, as.integer(input$sim_n))) else 25
    similar_genes(selected(), n = n, min_tpm = as.numeric(input$sim_min))
  })
  output$sim_note <- renderUI({
    gid <- selected(); req(gid)
    if (!gid %in% rownames(stage_mean)) return(div(class = "alert alert-secondary my-2", "No expression values for this gene, so it has no profile to compare."))
    pk <- max(stage_mean[gid, ])
    warn <- if (pk < 1) "This gene is barely expressed (peak below 1 TPM); its profile is mostly noise, so the list below is not meaningful."
            else if (prof_sd[gid] < 0.25) "This gene hardly changes over development, so many genes have a similar (flat) profile; r values are not very informative."
    tagList(
      p(class = "small text-muted mt-2", "Genes whose stage profile is most similar to this gene: Pearson correlation (r) of the mean log2(TPM + 1) per stage ",
        "(15 values per gene). Similar timing suggests, but does not prove, shared regulation. Plot: profiles as z-scores (this gene in black, top 10 in colour).",
        if (is_mt(gid)) " This is a mitochondrial gene, so only the other mitochondrial genes are compared (nuclear genes are listed for nuclear genes only)."),
      if (!is.null(warn)) div(class = "text-warning small mb-2", warn))
  })
  output$sim_plot <- renderPlotly({
    s <- sim(); gid <- selected(); req(s, nrow(s) > 0)
    zs <- function(x) { v <- log_mean[x, ]; (v - mean(v)) / ifelse(sd(v) > 0, sd(v), 1) }
    p <- plot_ly()
    for (x in head(s$gene_id, 10))
      p <- p |> add_trace(x = stage_levels, y = zs(x), type = "scatter", mode = "lines", name = short_label(x, 22), opacity = 0.75,
                          hovertemplate = paste0(x, "<br>%{x}<br>z %{y:.2f}<extra></extra>"))
    p <- p |> add_trace(x = stage_levels, y = zs(gid), type = "scatter", mode = "lines+markers", name = paste(gid, "(this gene)"),
                        line = list(color = "black", width = 4), marker = list(color = "black"), hovertemplate = "%{x}<br>z %{y:.2f}<extra></extra>")
    p |> layout(xaxis = list(title = "", categoryorder = "array", categoryarray = stage_levels, tickangle = -45),
                yaxis = list(title = "z-score of log2(TPM + 1)"), margin = list(b = 125), legend = list(font = list(size = 10))) |>
      config(displaylogo = FALSE)
  })
  sim_table_data <- reactive({
    s <- sim(); req(s, nrow(s) > 0)
    data.frame(gene_id = s$gene_id, r = round(s$r, 3), description = nice_desc(s$gene_id),
               peak_stage = stage_levels[apply(stage_mean[s$gene_id, , drop = FALSE], 1, which.max)],
               peak_tpm = round(apply(stage_mean[s$gene_id, , drop = FALSE], 1, max), 1),
               origin = origin_of(s$gene_id), stringsAsFactors = FALSE)
  })
  output$sim_table <- renderDT({
    d <- sim_table_data()
    show <- data.frame(Gene = gene_link(d$gene_id), r = d$r, Description = d$description, `Peak stage` = d$peak_stage,
                       `Peak TPM` = d$peak_tpm, `Maternal / zygotic` = d$origin, check.names = FALSE, stringsAsFactors = FALSE)
    datatable(fillContainer = FALSE, data = show, rownames = FALSE, selection = "none", class = "compact stripe hover", escape = -1,
              options = list(pageLength = 10, scrollX = TRUE, dom = "tip"))
  })
  output$dl_sim <- downloadHandler(
    filename = function() paste0(selected(), "_similar_profiles.csv"),
    content = function(file) write.csv(sim_table_data(), file, row.names = FALSE))

  # ---- sequences
  observeEvent(selected(), {
    t <- txs()
    if (!nrow(t)) { updateSelectInput(session, "tx_sel", choices = character()); return() }
    lab <- paste0(t$transcript_id, ifelse(t$is_primary == 1, " (primary)", ""), " - ",
                  ifelse(nzchar(t$protein_id), paste0(t$protein_id, ", ", t$aa_len, " aa"), ifelse(nzchar(t$biotype), t$biotype, "no protein")))
    updateSelectInput(session, "tx_sel", choices = setNames(t$transcript_id, lab), selected = t$transcript_id[which(t$is_primary == 1)[1]])
  })
  current_seq <- reactive({
    req(selected()); g <- gene(); kind <- input$seq_kind
    if (kind == "gene") {
      r <- q1("SELECT seq_zlib FROM seq_gene WHERE gene_id = ?", g$gene_id)
      req(nrow(r) == 1)
      s <- unz(r$seq_zlib[[1]])
      hdr <- sprintf(">%s %s genomic gene sequence (introns included), 5'->3' on the gene strand", g$gene_id, loc_string(g$seqid, g$start, g$end, g$strand))
      return(list(seq = s, header = hdr, kind = "genomic", note = NULL))
    }
    if (!nrow(txs())) return(list(seq = "", header = "", kind = kind, note = "This gene has no transcript record (e.g. a RefSeq pseudogene): only the genomic sequence is available."))
    req(input$tx_sel)
    r <- q1("SELECT mrna_z, cds_z, protein FROM seq_tx WHERE transcript_id = ?", input$tx_sel)
    req(nrow(r) == 1)
    s <- switch(kind, mrna = unz(r$mrna_z[[1]]), cds = unz(r$cds_z[[1]]), protein = r$protein)
    if (is.na(s)) s <- ""
    label <- c(mrna = "mRNA (spliced transcript, T instead of U)", cds = "CDS", protein = "protein")[[kind]]
    pid <- txs()$protein_id[txs()$transcript_id == input$tx_sel]
    acc <- if (kind == "protein" && length(pid) && nzchar(pid)) pid else input$tx_sel
    if (!nzchar(s)) {
      return(list(seq = "", header = "", kind = label, note = paste0("No ", if (kind == "protein") "protein" else "CDS",
        " for this transcript (a non-coding RNA such as an lncRNA, rRNA or tRNA).")))
    }
    mt <- if (startsWith(g$gene_id, "MT-") && kind %in% c("cds", "protein")) " | mitochondrial genetic code (NCBI table 13), MITOS2" else ""
    list(seq = s, header = sprintf(">%s %s | %s | %s%s%s", acc, g$gene_id, label, g$description,
                                   if (kind == "cds" && length(pid) && nzchar(pid)) paste0(" | encodes ", pid) else "", mt), kind = label, note = NULL)
  })
  output$seq_info <- renderUI({
    s <- current_seq()
    if (!is.null(s$note)) return(div(class = "text-warning mb-2", s$note))
    unit <- if (identical(s$kind, "protein")) "aa" else "bp"
    gc <- if (identical(s$kind, "protein")) "" else sprintf(", GC %.1f%%", gc_percent(s$seq))
    div(class = "text-muted small mb-2", sprintf("%s: %s %s%s", s$kind, format(nchar(s$seq), big.mark = ","), unit, gc))
  })
  fasta_text <- reactive({ s <- current_seq(); if (is.null(s$note)) paste0(s$header, "\n", wrap_seq(s$seq)) else "" })
  output$seq_text <- renderText(fasta_text())
  output$dl_seq <- downloadHandler(
    filename = function() paste0(selected(), "_", input$seq_kind, ".fasta"),
    content = function(file) writeLines(fasta_text(), file))

  # ---- annotation
  output$annotation <- renderUI({
    g <- gene(); t <- txs(); e <- get_egg(g$gene_id)
    kv <- function(...) tags$table(class = "kv table-sm mb-3", tags$tbody(...))
    row <- function(k, v) tags$tr(tags$td(k), tags$td(v))
    o <- origin_tab[g$gene_id, ]
    refseq <- kv(
      row("Gene", g$gene_id), row("Description", g$description), row("Type", g$biotype),
      row("Location", loc_string(g$seqid, g$start, g$end, g$strand)),
      row("Maternal / zygotic", if (is_mt(g$gene_id)) "not applicable (mitochondrial gene)" else if (is.na(o$origin)) "no expression data" else o$origin),
      row("Peak stage", if (is.na(o$origin)) "no expression data" else if (is.na(o$peak_stage)) "not expressed" else o$peak_stage),
      row("Transcripts", if (!nrow(t)) "none in the annotation" else tags$div(lapply(seq_len(nrow(t)), function(i)
        tags$div(t$transcript_id[i], if (nzchar(t$protein_id[i])) paste0("  ->  ", t$protein_id[i], " (", t$aa_len[i], " aa)") else "")))))
    if (startsWith(g$gene_id, "MT-")) {
      egg <- div(class = "text-muted", "Mitochondrial gene annotated with MITOS2 (ascidian mitochondrial genetic code, NCBI table 13); not run through eggNOG-mapper.")
    } else if (!nrow(e)) {
      egg <- div(class = "text-muted", "No eggNOG annotation for this gene (no orthology hit, or non-coding).")
    } else {
      chips <- function(x, link = NULL) if (!nzchar(x)) "-" else tagList(lapply(strsplit(x, ",")[[1]], function(i)
        tags$span(class = "chip", if (is.null(link)) i else tags$a(href = paste0(link, i), target = "_blank", rel = "noopener", i))))
      egg <- kv(
        row("Preferred name", if (nzchar(e$preferred_name)) e$preferred_name else "-"), row("Description", if (nzchar(e$description)) e$description else "-"),
        row("Protein domains (Pfam)", chips(e$pfams)),
        row("GO terms", chips(e$gos, "https://amigo.geneontology.org/amigo/term/")),
        row("COG category", if (nzchar(e$cog_category)) e$cog_category else "-"), row("KEGG KO", if (nzchar(e$kegg_ko)) e$kegg_ko else "-"),
        row("Orthologous groups", e$og), row("Best eggNOG hit", paste0(e$seed_ortholog, "  (e-value ", e$evalue, ")")),
        row("Annotated protein", e$protein_id))
    }
    tagList(h5("RefSeq annotation"), refseq, h5("Functional annotation (eggNOG-mapper)"), egg)
  })

  # ---- compare
  # cmp_request holds what to plot; set by the button or by "Similar profiles" (input values update one cycle late)
  cmp_request <- reactiveVal(NULL)
  observeEvent(input$cmp_go, cmp_request(list(genes = input$cmp_genes, term = input$cmp_term, max = input$cmp_max)))
  observeEvent(input$add_cmp, {
    cur <- trimws(input$cmp_genes)
    updateTextAreaInput(session, "cmp_genes", value = if (nzchar(cur)) paste0(cur, "\n", selected()) else selected())
    showNotification(paste(selected(), "added to the Compare list"), type = "message", duration = 3)
  })
  observeEvent(input$sim_to_cmp, {
    s <- sim(); req(s, nrow(s) > 0)
    ids <- c(selected(), s$gene_id)
    updateTextAreaInput(session, "cmp_genes", value = paste(ids, collapse = "\n"))
    updateTextInput(session, "cmp_term", value = "")
    cmp_request(list(genes = paste(ids, collapse = "\n"), term = "", max = input$cmp_max))
    nav_select("nav", "compare")
  })
  cmp_data <- reactive({
    rq <- cmp_request(); req(rq)
    ids <- resolve_ids(rq$genes)
    genes <- ids$found
    if (nzchar(trimws(rq$term))) {
      s <- search_genes(rq$term, max_rows = rq$max)
      if (!is.null(s) && nrow(s)) genes <- unique(c(genes, s$gene_id))
    }
    if (!length(genes)) return(list(error = "No genes to plot. Enter gene IDs or a search term.", missing = ids$missing))
    no_expr <- setdiff(genes, rownames(stage_mean))
    genes <- intersect(genes, rownames(stage_mean))
    if (!length(genes)) return(list(error = "None of these genes has expression values.", missing = ids$missing))
    list(mat = stage_mean[genes, , drop = FALSE], missing = ids$missing, no_expr = no_expr)
  })
  cmp_matrix <- reactive({
    cd <- cmp_data(); req(is.null(cd$error))
    m <- cd$mat
    v <- switch(input$cmp_scale, tpm = m, log2 = log2(m + 1), z = { l <- log2(m + 1); z <- t(scale(t(l))); z[is.nan(z)] <- 0; z })
    if (input$cmp_order != "none" && nrow(v) >= 3) {
      base <- log2(m + 1)
      if (input$cmp_order == "pattern") {
        # shape only: z-score each gene, so high and low genes with the same profile cluster together
        base <- t(scale(t(base))); base[is.nan(base)] <- 0
        ord <- hclust(dist(base))$order
      } else {
        # level and shape: distances on unscaled log2(TPM + 1), so highly expressed genes group together;
        # branches are flipped so the highest-expressed clusters come first (top of the heatmap)
        dend <- reorder(as.dendrogram(hclust(dist(base))), -rowMeans(base), agglo.FUN = mean)
        ord <- order.dendrogram(dend)
      }
      v <- v[ord, , drop = FALSE]
    }
    v
  })
  output$cmp_msg <- renderUI({
    cd <- cmp_data()
    if (!is.null(cd$error)) return(div(class = "text-warning", cd$error))
    tagList(if (length(cd$missing)) div(class = "text-warning small", "Not found: ", paste(cd$missing, collapse = ", ")),
            if (length(cd$no_expr)) div(class = "text-warning small", "No expression values (e.g. pseudogenes): ", paste(cd$no_expr, collapse = ", ")),
            if (MT_OWN_SCALE && any(is_mt(rownames(cd$mat))) && !all(is_mt(rownames(cd$mat))))
              div(class = "text-warning small", "Mitochondrial (MT-) and nuclear genes are on separate TPM scales; compare their shapes (z-score), not their levels."),
            div(class = "text-muted small", sprintf("%d genes; mean of the two replicates per stage. Uncharacterized genes are labelled with their eggNOG name.", nrow(cd$mat))))
  })
  output$cmp_plot_ui <- renderUI({
    cd <- cmp_data(); req(is.null(cd$error))
    plotlyOutput("cmp_plot", height = paste0(max(380, 20 * nrow(cd$mat) + 200), "px"))
  })
  output$cmp_plot <- renderPlotly({
    v <- cmp_matrix()
    lab <- short_label(rownames(v))
    zl <- switch(input$cmp_scale, tpm = "TPM", log2 = "log2(TPM+1)", z = "z-score")
    p <- plot_ly(x = colnames(v), y = lab, z = unname(v), type = "heatmap",
                 colorscale = if (input$cmp_scale == "z") "RdBu" else "Viridis", reversescale = input$cmp_scale == "z",
                 hovertemplate = "%{y}<br>%{x}<br>%{z:.2f}<extra></extra>", colorbar = list(title = zl))
    p |> layout(xaxis = list(title = "", tickangle = -45, categoryorder = "array", categoryarray = colnames(v)),
                yaxis = list(title = "", autorange = "reversed", type = "category"), margin = list(l = 260, b = 150)) |> config(displaylogo = FALSE)
  })
  output$dl_cmp <- downloadHandler(
    filename = function() "compare_mean_tpm.csv",
    content = function(file) { cd <- cmp_data(); req(is.null(cd$error)); write.csv(round(cd$mat, 4), file) })

  # ---- stages
  output$stage_table <- renderDT({
    st <- unique(samples[, c("dev_order", "stage", "hotta_stage", "stage_label")])
    st <- st[order(st$dev_order), ]
    key <- ifelse(is.na(st$hotta_stage) | !nzchar(st$hotta_stage), st$stage, st$hotta_stage)
    desc <- unname(STAGE_INFO[key]); desc[is.na(desc)] <- ""
    runs <- vapply(st$dev_order, function(k) paste(samples$run_accession[samples$dev_order == k], collapse = ", "), "")
    datatable(fillContainer = FALSE, data = data.frame(`#` = st$dev_order, Stage = st$stage, `Hotta stage` = st$hotta_stage, Description = desc, `ENA runs` = runs,
                         check.names = FALSE), rownames = FALSE, class = "compact stripe", options = list(dom = "t", pageLength = 50, ordering = FALSE))
  })
  output$origin_text <- renderUI({
    r <- ORIGIN_RULES
    tagList(
      p("Each gene is classed by whether its RNA is already in the unfertilized egg (mean TPM of the two replicates). ",
        "A maternal gene can also be transcribed by the embryo later: polyA RNA-seq cannot tell new zygotic transcripts from maternal ones ",
        "of the same gene (that needs e.g. intron reads or metabolic labelling). The peak stage is the stage with the highest mean TPM."),
      tags$ul(
        tags$li(tags$strong("maternal"), sprintf(": at least %g TPM in the unfertilized egg.", r$min_tpm)),
        tags$li(tags$strong("zygotic"), sprintf(": below %g TPM in the egg, and at least %g TPM at some later stage.", r$min_tpm, r$min_tpm)),
        tags$li(tags$strong("not expressed"), sprintf(": below %g TPM at every stage.", r$min_tpm)),
        tags$li(tags$strong(MT_CLASS), ": genes of the mitochondrial genome (MT-). They are transcribed in the mitochondria, not from the zygotic nuclear genome, so the classes do not apply.")),
      p(class = "small text-muted", "Use the filter on the Search tab to list one class. Counts in this database:"))
  })
  output$origin_counts <- renderDT({
    tab <- table(factor(origin_tab$origin, levels = ORIGIN_LEVELS))
    bt <- gene_tab[origin_tab$gene_id, "biotype"]
    pc <- table(factor(origin_tab$origin[bt == "protein_coding"], levels = ORIGIN_LEVELS))
    datatable(fillContainer = FALSE, data = data.frame(Class = names(tab), `All genes` = as.integer(tab), `Protein-coding` = as.integer(pc), check.names = FALSE),
              rownames = FALSE, class = "compact", options = list(dom = "t", ordering = FALSE))
  })

  # ---- about / downloads
  output$about_text <- renderUI({
    tagList(
      tags$table(class = "kv table-sm mb-3", tags$tbody(
        tags$tr(tags$td("App version"), tags$td(APP_VERSION, " (", tags$a(href = "https://github.com/AlenaKizenko/phmamm-expression-atlas",
                                                                          target = "_blank", rel = "noopener", "GitHub"), ")")),
        lapply(c("assembly", "annotation", "reads", "pipeline", "eggnog", "built", "n_genes"), function(k)
        tags$tr(tags$td(c(assembly = "Genome assembly", annotation = "Gene annotation", reads = "RNA-seq data", pipeline = "Quantification",
                          eggnog = "Functional annotation", built = "Database built", n_genes = "Genes")[[k]]), tags$td(meta[[k]]))))),
      p(class = "small text-muted", "Expression is given as TPM (transcripts per million) and read counts per gene, from 15 developmental stages (two replicates each). ",
        "Gene names are RefSeq names; many genes are 'uncharacterized' - use the eggNOG annotation and protein domains to find candidates. ",
        "Mitochondrial genes (MT-...) were annotated with MITOS2.",
        if (!is.na(meta["tpm_note"])) paste0(" TPM: ", meta[["tpm_note"]], ". The mitochondrial transcripts take 15-30% of Salmon's TPM, depending on stage, so on one scale they would shift every nuclear gene between stages.")))
  })
  output$samples_table <- renderDT(datatable(fillContainer = FALSE, data = samples[, c("sample", "run_accession", "ena_sample", "stage_label", "replicate")], rownames = FALSE,
                                             class = "compact stripe", options = list(pageLength = 10, dom = "tip")))
  tsv <- function(f) downloadHandler(filename = function() f$name, content = function(file) write.table(f$data(), file, sep = "\t", quote = FALSE, row.names = FALSE))
  output$dl_tpm <- tsv(list(name = "phallusia_gene_tpm.tsv", data = function() wide_table("tpm")))
  output$dl_counts <- tsv(list(name = "phallusia_gene_counts.tsv", data = function() wide_table("counts")))
  output$dl_annot <- tsv(list(name = "phallusia_gene_annotation.tsv", data = function() {
    d <- dbGetQuery(con, "SELECT g.gene_id, g.description, g.biotype, g.seqid, g.start, g.end, g.strand, g.ncbi_gene_id, e.preferred_name, e.description AS eggnog_description, e.pfams, e.gos, e.kegg_ko, e.og FROM genes g LEFT JOIN eggnog e USING(gene_id)")
    d$maternal_zygotic <- origin_of(d$gene_id)
    d$peak_stage <- origin_tab[d$gene_id, "peak_stage"]
    d
  }))
}

shinyApp(ui, server)
