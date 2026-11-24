# lintr over the R scripts run by the pipeline. Style rules are relaxed to
# what the existing scripts follow; the checks that remain catch real errors.
linters <- lintr::linters_with_defaults(
  line_length_linter = NULL,
  object_name_linter = NULL,
  indentation_linter = NULL,
  commented_code_linter = NULL,
  object_usage_linter = NULL,
  brace_linter = NULL,
  quotes_linter = NULL,
  infix_spaces_linter = NULL,
  commas_linter = NULL,
  spaces_left_parentheses_linter = NULL,
  object_length_linter = NULL,
  trailing_whitespace_linter = NULL,
  return_linter = NULL,
  whitespace_linter = NULL,
  assignment_linter = NULL,
  paren_body_linter = NULL,
  seq_linter = NULL,
  vector_logic_linter = NULL
)
files <- list.files("workflow/scripts", pattern = "\\.R$", full.names = TRUE)
results <- unlist(lapply(files, lintr::lint, linters = linters), recursive = FALSE)
for (lint in results) print(lint)
if (length(results) > 0) quit(status = 1)
