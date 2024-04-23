# Contribution Guidelines

## Issues

If you encounter a problem with ostRich, please try the following before creating a new issue:

1. Read the documentation thoroughly, and check that the input matches what is expected (in most cases, this is a named list)
2. Restart R and rerun the offending code
3. Update ostRich to the latest version

If none of these solutions work, please create a new issue that includes:

1. A clear statement of the problem in the title
2. A small reproducible example
3. Additional detailed explanation, as needed


## Pull Requests

All contributed code must adhere to the tidyverse style guide: https://style.tidyverse.org/index.html. When in doubt, refer to the existing codebase.

If adding new functionality, please include proper unit tests.

Verify that `devtools::check(document = TRUE)` runs without errors, warnings, or notes before submitting a pull request.
