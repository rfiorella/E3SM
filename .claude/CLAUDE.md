<!-- markdownlint-disable -->
# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.
@../AGENTS.md

## claude-workflows integration

Project initialized via `/init-project`. Configuration at
`.claude/project-config.yml`. Specs live under `specs/`; campaign
manifests under `campaigns/`. Active campaign:
`campaigns/wiso.campaign.md` (EAMxx water tracers + isotope ratios).

Workflow commands:
- `/new-spec model-e3sm` — start a new spec (interview).
- `/resolve-spec specs/<path>.md` — validate a spec against the
  rubric in the `spec-schema` skill.
- `/run-spec specs/<path>.md` — execute a validated spec.
- `/new-campaign` / `/run-campaign campaigns/<path>.md` — stacked-PR
  chain mode.
- `/restack-campaign campaigns/<path>.md` — rebase the chain when
  `master` advances or a PR merges.

Skills loaded for this project (per `.claude/project-config.yml`):
`spec-schema`, `e3sm-build-and-test`, `e3sm-platforms`,
`eamxx-cpp-conventions`, `tracer-conservation`, `regression-baseline`,
`scientific-modeling-conventions`.
