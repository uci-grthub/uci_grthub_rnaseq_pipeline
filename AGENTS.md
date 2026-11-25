# Notes for coding agents

- Software comes from pixi (`pixi.toml`): run commands with `pixi run ...`.
  Pick the environment: `default` for snakemake and Python, `-e tools` for
  aligners/QC, `-e r` for R. Add a tool with
  `pixi add --feature tools -c bioconda <pkg>==<version>`.
- Rules get their tools through `params.env=pixi_env(...)` plus `envmodules:`
  (workflow/rules/common.smk). A new rule needs both.
- Scripts run by rules live in `workflow/scripts/` and are called as
  `{SCRIPTS}/name`. Do not use paths relative to the working directory: runs
  happen from a project directory, not the repo.
- Don't put project data, configs or outputs in the repo.
- Any change to results must update `tests/expected/` in the same PR, with the
  reason given in the PR description; it is a major version bump.
- Before pushing: `pixi run lint`, then `pixi run tests/run_test.sh`.
- Formatting (`snakefmt`, pinned) reflows rule code, and Snakemake reruns
  rules whose code changed. Never reformat in the same PR as a logic change.
- PR titles use Conventional Commits (`feat:`, `fix:`, `feat!:`).
