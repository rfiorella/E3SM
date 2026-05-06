# EAMXX-wiso

Water-isotope-enabled fork of E3SM, focused on the EAMxx (SCREAM) atmosphere
component. Active development happens against EAMxx as a **standalone CMake
build**, not full coupled compsets.

## claude-workflows

This project uses the claude-workflows framework. Configuration lives at
`.claude/project-config.yml`. Specs are stored under `specs/` in **immutable**
mode — they are not git-tracked, and revisions are written as new dated files
rather than edits.

Common entry points:
- `/new-spec model-e3sm` — start a new spec for an EAMxx code change
- `/run-spec specs/<file>.md` — execute an approved spec
- `/resolve-spec specs/<file>.md` — address validator gaps

## Source layout

EAMxx source: `components/eamxx/`
- `src/` — C++ source (physics, dynamics, control, share)
- `tests/` — unit and integration tests
- `scripts/test-all-eamxx` — primary test driver
- `cmake/` — build configuration

## Platform

Default platform: `docker-local`.

Path mapping (host → container):
- `/Users/rfiorella/code/E3SM/EAMXX-wiso` → `/code/E3SM/EAMXX-wiso`
- `/data/...` exists only inside the container (baselines, run dirs).

Container image: `rfiorella/model-containers:e3sm-openmpi-dev-latest`.

## Build & test

Standalone EAMxx, run from inside the container:

```
cd /code/E3SM/EAMXX-wiso/components/eamxx
./scripts/test-all-eamxx ...
```

Target resolution for isotope-development iteration is `ne4` (lowest available)
to keep the inner loop fast.
