# ── PhenoSuite Provenance Tracker ────────────────────────────────────────────
# Lightweight R5 reference class that captures session provenance and generates
# a JSON sidecar + headless replay script alongside analysis outputs.
#
# Dependencies (already in Docker image): jsonlite, digest, tools
# ─────────────────────────────────────────────────────────────────────────────

ProvenanceTracker <- setRefClass("ProvenanceTracker",
  fields = list(
    app_name        = "character",
    session_id      = "character",
    session_start   = "POSIXct",
    analysis_start  = "ANY",       # POSIXct or NULL
    analysis_end    = "ANY",       # POSIXct or NULL
    input_files     = "list",
    output_dir      = "character",
    parameters      = "list",
    seeds           = "list",
    session_info    = "ANY",       # sessionInfo() output
    docker_digest   = "character",
    git_sha         = "character",
    custom_metadata = "list",
    replay_template = "ANY",       # function or NULL
    .input_ids      = "character"  # track registered fileInput IDs for filtering
  ),

  methods = list(

    # ── Constructor ──────────────────────────────────────────────────────────
    initialize = function(app_name, session = NULL, output_dir = tempdir(),
                          replay_template = NULL) {
      app_name        <<- app_name
      session_id      <<- if (!is.null(session)) session$token else "standalone"
      session_start   <<- Sys.time()
      analysis_start  <<- NULL
      analysis_end    <<- NULL
      input_files     <<- list()
      output_dir      <<- output_dir
      parameters      <<- list()
      seeds           <<- list()
      custom_metadata <<- list()
      .input_ids      <<- character(0)

      session_info    <<- utils::sessionInfo()
      docker_digest   <<- Sys.getenv("PHENOSUITE_IMAGE_DIGEST", "unknown")
      git_sha         <<- Sys.getenv("PHENOSUITE_GIT_SHA", "unknown")
      replay_template <<- replay_template
    },

    # ── Register uploaded input files ────────────────────────────────────────
    # file_info: the data.frame from input$fileN (columns: name, size, type, datapath)
    # input_id: the Shiny inputId string (e.g. "file1") used to filter from parameters
    register_input = function(file_info, input_id = NULL) {
      if (!is.null(input_id)) {
        .input_ids <<- unique(c(.input_ids, input_id))
      }
      if (is.data.frame(file_info)) {
        for (i in seq_len(nrow(file_info))) {
          hash <- digest::digest(file = file_info$datapath[i], algo = "sha256")
          input_files[[length(input_files) + 1L]] <<- list(
            original_name = file_info$name[i],
            size_bytes    = file_info$size[i],
            content_type  = if (!is.null(file_info$type)) file_info$type[i] else NA_character_,
            sha256        = hash
          )
        }
      }
    },

    # ── Capture all UI parameters ────────────────────────────────────────────
    # input: the Shiny input object (or a plain list for testing)
    capture_parameters = function(input) {
      if (inherits(input, "reactivevalues")) {
        params <- shiny::reactiveValuesToList(input)
      } else {
        params <- as.list(input)
      }
      # Remove Shiny internals (dot-prefixed keys)
      params <- params[!grepl("^\\.", names(params))]
      # Remove fileInput entries (they contain temp paths, not meaningful params)
      if (length(.input_ids) > 0) {
        params <- params[!names(params) %in% .input_ids]
      }
      # Convert non-JSON-safe values
      params <- lapply(params, function(x) {
        if (is.null(x)) return("__NULL__")
        if (is.function(x)) return("__FUNCTION__")
        if (inherits(x, "data.frame")) return(as.list(x))
        x
      })
      parameters <<- params
    },

    # ── Seed tracking ────────────────────────────────────────────────────────
    set_seed = function(value, label = "default") {
      base::set.seed(value)
      seeds[[label]] <<- value
    },

    # Record known seeds used by PhenoSuite shared utilities.
    # Call this at tracker init to document the hardcoded seeds without
    # modifying the utility functions themselves.
    record_known_seeds = function() {
      seeds[["clustering_leiden"]]  <<- 99L
      seeds[["downsampling"]]      <<- 99L
      seeds[["scatter_shuffle"]]   <<- 99L
    },

    # ── Timestamps ───────────────────────────────────────────────────────────
    analysis_started = function() {
      analysis_start <<- Sys.time()
    },

    analysis_completed = function() {
      analysis_end <<- Sys.time()
      write_provenance()
    },

    # ── Register a single output file ────────────────────────────────────────
    register_output = function(file_path) {
      # Outputs are hashed at write time; this is for explicit tracking
      custom_metadata[["registered_outputs"]] <<- c(
        custom_metadata[["registered_outputs"]], file_path
      )
    },

    # ── Write provenance.json + replay.R ─────────────────────────────────────
    write_provenance = function() {
      dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

      # ── Hash output files ──────────────────────────────────────────────────
      output_entries <- list()
      all_files <- list.files(output_dir, recursive = TRUE, full.names = TRUE)
      # Exclude provenance files themselves
      all_files <- all_files[!basename(all_files) %in% c("provenance.json", "replay.R")]
      for (f in all_files) {
        finfo <- file.info(f)
        hash <- tryCatch(
          digest::digest(file = f, algo = "sha256"),
          error = function(e) paste0("ERROR: ", e$message)
        )
        output_entries[[length(output_entries) + 1L]] <- list(
          relative_path = sub(paste0("^", gsub("([.+*?^${}()|\\[\\]])", "\\\\\\1", output_dir), "/?"), "", f),
          sha256        = hash,
          size_bytes    = as.numeric(finfo$size)
        )
      }

      # ── Parse sessionInfo ──────────────────────────────────────────────────
      si <- session_info
      pkg_list <- list()
      if (!is.null(si$otherPkgs)) {
        for (p in si$otherPkgs) {
          pkg_list[[p$Package]] <- p$Version
        }
      }
      if (!is.null(si$loadedOnly)) {
        for (p in si$loadedOnly) {
          pkg_list[[p$Package]] <- p$Version
        }
      }

      # ── Build sidecar ─────────────────────────────────────────────────────
      fmt_time <- function(t) {
        if (is.null(t)) return(jsonlite::unbox(NA))
        jsonlite::unbox(format(t, "%Y-%m-%dT%H:%M:%S%z"))
      }

      sidecar <- list(
        phenosuite_version = jsonlite::unbox("1.0.0"),
        app_name           = jsonlite::unbox(app_name),
        session_id         = jsonlite::unbox(session_id),
        timestamps = list(
          session_start  = fmt_time(session_start),
          analysis_start = fmt_time(analysis_start),
          analysis_end   = fmt_time(analysis_end)
        ),
        environment = list(
          r_version            = jsonlite::unbox(si$R.version$version.string),
          platform             = jsonlite::unbox(si$R.version$platform),
          docker_image_digest  = jsonlite::unbox(docker_digest),
          git_sha              = jsonlite::unbox(git_sha),
          packages             = pkg_list
        ),
        inputs     = input_files,
        parameters = parameters,
        seeds      = seeds,
        outputs    = output_entries
      )
      if (length(custom_metadata) > 0) {
        sidecar$custom_metadata <- custom_metadata
      }

      json_path <- file.path(output_dir, "provenance.json")
      writeLines(
        jsonlite::toJSON(sidecar, pretty = TRUE, auto_unbox = FALSE, null = "null"),
        json_path
      )

      # ── Generate replay script ─────────────────────────────────────────────
      replay_path <- file.path(output_dir, "replay.R")
      writeLines(.build_replay_script(), replay_path)

      message(sprintf("[provenance] Written %s and %s", json_path, replay_path))
    },

    # ── Internal: build replay.R content ─────────────────────────────────────
    .build_replay_script = function() {
      buf <- new.env(parent = emptyenv())
      buf$lines <- character(0)
      add <- function(...) buf$lines <- c(buf$lines, ...)

      add("#!/usr/bin/env Rscript")
      add(sprintf("# Auto-generated replay script by PhenoSuite Provenance System"))
      add(sprintf("# Generated: %s", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
      add(sprintf("# App: %s", app_name))
      add(sprintf("# Original session: %s", session_id))
      add("")

      # ── Environment info ───────────────────────────────────────────────────
      add("# ── Environment ─────────────────────────────────────────────────")
      add(sprintf("# R version: %s", session_info$R.version$version.string))
      add(sprintf("# Docker image digest: %s", docker_digest))
      add(sprintf("# Git SHA: %s", git_sha))
      add("")

      # ── Library loads ──────────────────────────────────────────────────────
      add("# ── Libraries ───────────────────────────────────────────────────")
      if (!is.null(session_info$otherPkgs)) {
        for (p in session_info$otherPkgs) {
          add(sprintf("library(%s)", p$Package))
        }
      }
      add("")

      # ── Seeds ──────────────────────────────────────────────────────────────
      if (length(seeds) > 0) {
        add("# ── Random seeds ────────────────────────────────────────────────")
        add("seeds <- list(")
        seed_lines <- sapply(names(seeds), function(nm) {
          sprintf("  %s = %s", nm, paste(deparse(seeds[[nm]]), collapse = " "))
        }, USE.NAMES = FALSE)
        add(paste(seed_lines, collapse = ",\n"))
        add(")")
        add("")
      }

      # ── Parameters ─────────────────────────────────────────────────────────
      add("# ── Parameters ──────────────────────────────────────────────────")
      add("params <- list(")
      if (length(parameters) > 0) {
        param_lines <- sapply(names(parameters), function(nm) {
          val <- parameters[[nm]]
          # deparse() can return multiple lines for complex objects; collapse them
          deparsed <- paste(deparse(val, width.cutoff = 500L), collapse = " ")
          sprintf("  `%s` = %s", nm, deparsed)
        }, USE.NAMES = FALSE)
        add(paste(param_lines, collapse = ",\n"))
      }
      add(")")
      add("")

      # ── Input files ────────────────────────────────────────────────────────
      add("# ── Input files ─────────────────────────────────────────────────")
      add("# Place your input files in the same directory as this script.")
      for (entry in input_files) {
        add(sprintf("# Expected: %s (SHA-256: %s)", entry$original_name, entry$sha256))
      }
      if (length(input_files) > 0) {
        fnames <- vapply(input_files, function(e) e$original_name, character(1))
        add(sprintf("input_files <- c(%s)",
                     paste(sprintf('"%s"', fnames), collapse = ", ")))
      } else {
        add("input_files <- character(0)")
      }
      add("")

      # ── Verify input hashes ────────────────────────────────────────────────
      add("# ── Verify input integrity ──────────────────────────────────────")
      add("if (requireNamespace('digest', quietly = TRUE)) {")
      add("  expected_hashes <- c(")
      if (length(input_files) > 0) {
        hash_lines <- vapply(input_files, function(e) {
          sprintf('    "%s"', e$sha256)
        }, character(1))
        add(paste(hash_lines, collapse = ",\n"))
      }
      add("  )")
      add("  for (i in seq_along(input_files)) {")
      add("    if (file.exists(input_files[i])) {")
      add("      actual <- digest::digest(file = input_files[i], algo = 'sha256')")
      add("      if (actual != expected_hashes[i]) {")
      add("        warning(sprintf('Hash mismatch for %s: expected %s, got %s',")
      add("                        input_files[i], expected_hashes[i], actual))")
      add("      } else {")
      add("        message(sprintf('Verified: %s', input_files[i]))")
      add("      }")
      add("    } else {")
      add("      warning(sprintf('Input file not found: %s', input_files[i]))")
      add("    }")
      add("  }")
      add("}")
      add("")

      # ── Output directory ───────────────────────────────────────────────────
      add("# ── Output ──────────────────────────────────────────────────────")
      add('out_dir <- file.path(getwd(), "replay_output")')
      add("dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)")
      add("")

      # ── App-specific body or fallback ──────────────────────────────────────
      if (is.function(replay_template)) {
        add("# ── Analysis (app-specific) ──────────────────────────────────")
        body_lines <- tryCatch(
          replay_template(parameters, input_files, "out_dir"),
          error = function(e) {
            sprintf("# ERROR generating replay body: %s", e$message)
          }
        )
        add(body_lines)
      } else {
        add("# ── Analysis ────────────────────────────────────────────────────")
        add(sprintf("# This replay script was generated for the '%s' app.", app_name))
        add("# No app-specific replay template was provided.")
        add("# To replay, use the parameters and input files above with the")
        add("# appropriate PhenoSuite utility functions.")
        add("#")
        add("# Typical usage pattern:")
        add("#   source('/srv/shiny-server/phenomenalist/utils/RunPhenomenalist-shiny/RunPhenomenalist-shiny.R')")
        add("#   # Then call the relevant analysis function with params and input_files")
      }
      add("")
      add("message('Replay complete.')")

      buf$lines
    }
  )
)
