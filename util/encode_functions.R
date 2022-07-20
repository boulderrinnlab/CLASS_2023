library(tidyverse)
library(httr)
library(janitor)

contstruct_query <- function(experiment_accession,
                             base_url = "https://www.encodeproject.org/report.tsv?",
                             file_format = "fastq",
                             type = "File",
                             status = "released",
                             fields = c("accession", "read_count", "md5sum",
                                        "controlled_by", "paired_end",
                                        "paired_with", "replicate", "target")) {
  query <- paste(list(paste0("type=", type),
                      paste0("status=", status),
                      paste0("file_format=", file_format),
                      paste0("dataset=%2Fexperiments%2F", experiment_accession, "%2F"),
                      map_chr(fields, ~paste0("field=", .))) %>%
                   flatten(),
                 collapse = "&")
  url <- paste0(base_url, query)
  return(url)
}
encode_file_info <- function(experiment_accession,
                             base_url = "https://www.encodeproject.org/report.tsv?",
                             file_format = "fastq",
                             type = "File",
                             status = "released",
                             fields = c("accession", "read_count", "md5sum",
                                        "controlled_by", "paired_end",
                                        "paired_with", "replicate", "target")) {
  path <- "report.tsv?"
  base_url <- modify_url("https://www.encodeproject.org/", path = path)
  url <- contstruct_query(experiment_accession,
                          base_url = base_url,
                          file_format,
                          type,
                          status,
                          fields)
  resp <- GET(url)
  if (http_error(resp)) {
    error_message <- content(resp, type = "text/html", encoding = "UTF-8") %>%
      xml_find_all("//p") %>%
      xml_text() %>%
      first()
    stop(
      sprintf(
        "ENCODE API request failed [%s]\n%s",
        status_code(resp),
        error_message
      ),
      call. = FALSE
    )
  }
  if (http_type(resp) != "text/tsv") {
    stop("API did not return text/tsv", call. = FALSE)
  }
  body <- read_tsv(content(resp, "text"), skip = 1) %>%
    clean_names()
  return(body)
}