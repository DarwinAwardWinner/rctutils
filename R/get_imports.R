## Inspired by https://stackoverflow.com/a/79835313/125921.


#' Find the parent environment containing a name.
#'
#' @param what A character string, the name to search for
#' @param envir The environment at which to start the search. If the name is not
#'   found in this environment, then its parent (enclosing) environments are
#'   searched. This can also be a function, in which case the function's
#'   environment is used as the starting point.
#' @param stop_after An environment after which the search should be terminated.
#'   The parent and further ancestors of this environment will not be searched.
#' @param namespace_only If TRUE, then only return the name of the found
#'   environment if it is a namespace. Otherwise return NA.
#' @inheritParams utils::find
#'
#' @returns The name of the environment in which the specified name is found or
#'   NA if the environment is unnamed, or the name is not found in any
#'   environment.
#' @export
#'
#' @examples
#' find_in_environment("median", envir = stats::median)
#' find_in_environment("stop", envir = stats::median)
find_in_environment <- function(what, envir, mode = "any", stop_after = .BaseNamespaceEnv, namespace_only = TRUE) {
    stopifnot(is.character(what))
    if (length(what) > 1L) {
        warning("elements of 'what' after the first will be ignored")
        what <- what[1L]
    }
    if (is.function(envir)) {
        envir <- environment(envir)
    }
    envir <- as.environment(envir)
    stop_after <- as.environment(stop_after)
    current_env <- envir
    found <- FALSE
    while (!found) {
        if (exists(what, current_env, mode = mode, inherits = FALSE)) {
            found <- TRUE
        } else {
            if (identical(current_env, stop_after)) {
                break
            }
            current_env <- parent.env(current_env)
            if (is.null(current_env) || identical(current_env, emptyenv())) {
                break
            }
        }
    }
    if (found) {
        if (namespace_only && !isNamespace(current_env)) {
            return(NA_character_)
        }
        env_name <- environmentName(current_env)
        if (is.null(env_name) || is.na(env_name) || env_name == "") {
            return(NA_character_)
        } else {
            return(env_name[1])
        }
    } else {
        return(NA_character_)
    }
}


#' Get a table of all imported functions in a function body
#'
#' @param fun The function to search for imported functions
#' @param ignore_env_names List of package names to ignore. By default, "base"
#'   is ignored since base functions do not need to be declared in packages. The
#'   global and empty environments are also ignored, since functions in those
#'   environments are not importable.
#'
#' @returns A data frame listing all functions imported by `fun`, along with
#'   various ways to import those functions into your package.
#' @export
#'
#' @examples
#'
#' get_imports(stats::xtabs)
#'
get_imports <- function(fun, envir = environment(fun), ignore_env_names = c("base", "R_GlobalEnv", "R_EmptyEnv", NA)) {
    req_ns("codetools")
    stopifnot(is.function(fun))
    globals <- codetools::findGlobals(fun)
    pks <- vapply(globals, find_in_environment, character(1), namespace_only = FALSE, envir = envir, mode = "function")
    pks <- vapply(pks, function(x) c(x, NA_character_)[1], character(1))
    data.frame(Name = globals, Env_Name = pks) |>
        mutate(
            Package = sub("^package:", "", .data$Env_Name),
            NAMESPACE_Import = ifelse(
                is.na(.data$Package),
                NA,
                paste0("importFrom(", .data$Package, ",", .data$Name, ")")
            ),
            Roxygen2_Import = ifelse(
                is.na(.data$Package),
                NA,
                paste0("#' @importFrom ", .data$Package, " ", .data$Name)
            ),
            usethis_Import = ifelse(
                is.na(.data$Package),
                NA,
                paste0("usethis::use_import_from(\"", .data$Package, "\", \"", .data$Name, "\")")
            )
        ) |>
        filter(!.data$Package %in% ignore_env_names)
}
