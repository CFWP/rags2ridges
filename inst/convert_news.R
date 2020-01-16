
library("stringr")
news <- readLines("inst/NEWS.Rd")

cat(head(news, n = 100), sep="\n")


news_md <-
  news %>%
  tail(n = -2) %>%
  str_replace("\\\\code\\{(.+?)\\}", "`\\1`") %>%
  str_replace("\\\\code\\{(.+?)\\}", "`\\1`") %>%
  str_replace("\\\\cpkg\\{(.+?)\\}", "[\\1](https://cran.r-project.org/package=\\1)") %>%
  str_replace("\\\\cpkg\\{(.+?)\\}", "[\\1](https://cran.r-project.org/package=\\1)") %>%
  str_replace("\\\\cpkg\\{(.+?)\\}", "[\\1](https://cran.r-project.org/package=\\1)") %>%
  str_replace("^[[:space:]]*\\\\title\\{", "# ") %>%
  str_replace("^[[:space:]]*\\\\section\\{", "## ") %>%
  str_replace("^[[:space:]]*\\\\subsection\\{", "### ") %>%
  str_replace("^[[:space:]]*\\\\itemize\\{", "") %>%
  str_replace("^      \\\\item", "  *") %>%  # Level 1 items
  str_replace("^          \\\\item", "    -") %>%  # Level 2 items
  str_replace("[[:space:]]*\\}\\{*$", "") %>%
  str_replace("^        ([^[:space:]].*)", "    \\1") %>%  # Fix indentation
  str_replace("^            ([^[:space:]].*)", "      \\1")  # Fix indentation level 2


# Remove all blank lines
news_md <- news_md[news_md != ""]

# Add blank lines before each header
i <- 1
while (i <= length(news_md)) {
  is_header <- grepl("^#{2,3} ", news_md[i])
  if (is_header) {
    n_blanks <- 4 - str_count(news_md[i], "#")
    news_md <- c(news_md[seq_len(i - 1)], rep("", n_blanks), news_md[i:length(news_md)])
    i <- i + n_blanks + 1  # more  as we have now interjected a blank line
  } else {
    i <- i + 1
  }
}

news_md %>%
  # head(n = 100) %>%
  cat(sep = "\n")

# news %>%
#   head(n = 18) %>%
#   str_view_all("^\\\\title\\{(.*)\\}$")

writeLines(news_md, con = "NEWS.md")
