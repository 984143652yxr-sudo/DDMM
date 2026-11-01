# Small helpers shared by run.R, report.R and the tools/ scripts.
# Nothing in this file affects simulated data or fitted results.

`%||%` <- function(a, b) if (is.null(a) || !length(a)) b else a

# ---------------------------------------------------------------------------
# Package location and paths
# ---------------------------------------------------------------------------

# User-supplied paths are relative to the shell's working directory.
absolute_path <- function(path, base) {
  path <- path.expand(path)
  if (!grepl('^/', path)) path <- file.path(base, path)
  normalizePath(path, mustWork = FALSE)
}

# Show paths inside the package relative to it (shorter log lines).
pretty_path <- function(path, root = getwd()) {
  path <- normalizePath(path, mustWork = FALSE)
  root <- normalizePath(root)
  if (startsWith(path, paste0(root, '/'))) substring(path, nchar(root) + 2) else path
}

# ---------------------------------------------------------------------------
# Command-line parsing:  positional arguments plus --key value / --flag
# ---------------------------------------------------------------------------

parse_command_line <- function(args, flags = character(), options = character()) {
  positional <- character()
  values <- list()
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (startsWith(a, '--')) {
      key <- sub('^--', '', a)
      if (grepl('=', key)) {
        values[[sub('=.*$', '', key)]] <- sub('^[^=]*=', '', key)
      } else if (key %in% flags) {
        values[[key]] <- TRUE
      } else if (key %in% options) {
        if (i == length(args)) stop('Option --', key, ' needs a value.')
        values[[key]] <- args[i + 1]
        i <- i + 1
      } else {
        stop('Unknown option --', key, '. Known: ',
             paste0('--', c(options, flags), collapse = ' '))
      }
    } else {
      positional <- c(positional, a)
    }
    i <- i + 1
  }
  unknown <- setdiff(names(values), c(flags, options))
  if (length(unknown)) stop('Unknown option --', unknown[1])
  list(positional = positional, options = values)
}

# "1:50" -> 1..50 ;  "51:100" ;  "7" ;  "1,3,5" ;  "1:3,10"
parse_replicates <- function(spec) {
  parts <- strsplit(gsub('\\s', '', spec), ',', fixed = TRUE)[[1]]
  ids <- unlist(lapply(parts, function(part) {
    if (grepl('^[0-9]+:[0-9]+$', part)) {
      ends <- as.integer(strsplit(part, ':', fixed = TRUE)[[1]])
      seq.int(ends[1], ends[2])
    } else if (grepl('^[0-9]+$', part)) {
      as.integer(part)
    } else stop('Cannot read replicate specification "', spec,
                '". Use e.g. 1:50, 51:100, 7 or 1,3,5.')
  }))
  ids <- sort(unique(ids))
  if (!length(ids) || any(ids < 1) || any(ids > 999))
    stop('Replicate IDs must lie in 1..999 (seeds reserve three digits).')
  ids
}

parse_number_list <- function(spec) as.numeric(strsplit(spec, ',', fixed = TRUE)[[1]])
parse_name_list <- function(spec) strsplit(gsub('\\s', '', spec), ',', fixed = TRUE)[[1]]

# ---------------------------------------------------------------------------
# Safe file writing
# ---------------------------------------------------------------------------

# Write to a temporary name and rename: an interrupted job never leaves a
# half-written checkpoint behind.
save_rds_atomic <- function(object, path) {
  temporary <- paste0(path, '.tmp', Sys.getpid())
  saveRDS(object, temporary)
  if (!file.rename(temporary, path)) stop('Could not write ', path)
  invisible(path)
}

write_csv_atomic <- function(x, path) {
  temporary <- paste0(path, '.tmp', Sys.getpid())
  utils::write.csv(x, temporary, row.names = FALSE)
  if (!file.rename(temporary, path)) stop('Could not write ', path)
  invisible(path)
}

write_lines_atomic <- function(text, path) {
  temporary <- paste0(path, '.tmp', Sys.getpid())
  writeLines(text, temporary)
  if (!file.rename(temporary, path)) stop('Could not write ', path)
  invisible(path)
}

# ---------------------------------------------------------------------------
# Locks
# A lock is a directory (creation is atomic) holding owner.txt with host/pid.
# A lock whose owner process no longer exists on this host is stale and is
# removed automatically, so a killed job never blocks the next one.
# ---------------------------------------------------------------------------

lock_owner <- function(path) {
  file <- file.path(path, 'owner.txt')
  if (!file.exists(file)) return(NULL)
  fields <- readLines(file, warn = FALSE)
  values <- sub('^[^=]*=', '', fields)
  names(values) <- sub('=.*$', '', fields)
  as.list(values)
}

lock_is_stale <- function(path) {
  owner <- lock_owner(path)
  if (is.null(owner)) {
    # Owner file not written yet (a lock being created right now) or a lock
    # left by a crash in that instant; treat as stale only after a minute.
    age <- as.numeric(difftime(Sys.time(), file.info(path)$mtime, units = 'secs'))
    return(is.finite(age) && age > 60)
  }
  same_host <- identical(owner$host, Sys.info()[['nodename']])
  same_host && !isTRUE(tools::pskill(as.integer(owner$pid), 0L))
}

acquire_lock <- function(path) {
  for (attempt in 1:2) {
    if (dir.create(path, showWarnings = FALSE)) {
      writeLines(c(paste0('host=', Sys.info()[['nodename']]),
                   paste0('pid=', Sys.getpid()),
                   paste0('started=', format(Sys.time(), '%Y-%m-%d %H:%M:%S'))),
                 file.path(path, 'owner.txt'))
      return(TRUE)
    }
    if (lock_is_stale(path)) unlink(path, recursive = TRUE) else return(FALSE)
  }
  FALSE
}

release_lock <- function(path) unlink(path, recursive = TRUE)

describe_lock <- function(path) {
  owner <- lock_owner(path)
  if (is.null(owner)) return('unknown owner')
  sprintf('pid %s on %s since %s', owner$pid, owner$host, owner$started)
}

# ---------------------------------------------------------------------------
# Logging: every message goes to the console and to the run's run.log.
# ---------------------------------------------------------------------------

.log <- new.env()
.log$file <- NULL

log_to <- function(path) .log$file <- path

say <- function(...) {
  text <- paste0(format(Sys.time(), '%H:%M:%S'), '  ', paste0(...))
  cat(text, '\n', sep = '')
  if (!is.null(.log$file)) cat(text, '\n', sep = '', file = .log$file, append = TRUE)
  invisible(text)
}

# ---------------------------------------------------------------------------
# External programs (snpknock2 for HMM knockoffs, BOLT-LMM)
# Order: environment variable, then PATH, then no fallback.
# ---------------------------------------------------------------------------

find_program <- function(variable, command) {
  candidates <- c(Sys.getenv(variable), Sys.which(command))
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    if (file.exists(candidate)) return(normalizePath(candidate))
  }
  ''
}

snpknock2_path <- function() find_program(
  'SNPKNOCK_BIN', 'snpknock2')

bolt_path <- function() find_program(
  'BOLT_BIN', 'bolt')

require_programs <- function(methods) {
  missing <- character()
  if ('Knockoff' %in% methods && !nzchar(snpknock2_path()))
    missing <- c(missing, 'snpknock2 (set SNPKNOCK_BIN=/path/to/snpknock2)')
  if ('BOLT' %in% methods && !nzchar(bolt_path()))
    missing <- c(missing, 'BOLT-LMM (set BOLT_BIN=/path/to/bolt)')
  if (length(missing)) stop('Cannot find ', paste(missing, collapse = ' and '), '.')
  invisible(TRUE)
}

require_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Install the R package(s): ', paste(missing, collapse = ', '))
}
