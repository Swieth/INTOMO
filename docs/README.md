# INTOMO Documentation

This folder contains the technical documentation for **INTOMO v1.0**, a MATLAB-based GNSS troposphere tomography library.

## Documents

| File | Contents |
|------|----------|
| [01-overview.md](01-overview.md) | What INTOMO is, the two-component architecture, and a high-level pipeline diagram |
| [02-setup.md](02-setup.md) | Prerequisites, installation paths, and required folder structure |
| [03-configuration.md](03-configuration.md) | Every setting in `conf.m` and the `switches` structure |
| [04-inputs.md](04-inputs.md) | Required input files with full format specifications |
| [05-workflow.md](05-workflow.md) | The nine processing stages mapped to MATLAB functions |
| [06-outputs.md](06-outputs.md) | Every produced `.mat` artefact with its on-disk path and contents |
| [07-data-structures.md](07-data-structures.md) | Reference for key MATLAB structures (`model`, `station`, `values`, etc.) |
| [08-folder-map.md](08-folder-map.md) | Annotated directory tree for the whole repository |
| [09-issues.md](09-issues.md) | Identified logical bugs, portability problems, robustness gaps, and code-style issues |

## Quick start

1. Read [02-setup.md](02-setup.md) to understand the prerequisites and folder layout.
2. Read [03-configuration.md](03-configuration.md) and edit `Tomography/conf.m` for your project.
3. Place input data as described in [04-inputs.md](04-inputs.md).
4. Run `RUNINTOMO.m` from the MATLAB command window.
5. Retrieve results as described in [06-outputs.md](06-outputs.md).

For a deeper understanding of what happens internally, see [05-workflow.md](05-workflow.md).

## Source references

The documentation was derived from:

- [`Tomography/RUNINTOMO.m`](../Tomography/RUNINTOMO.m) — main driver script
- [`Tomography/conf.m`](../Tomography/conf.m) — configuration file
- [`Tomography/ENGINE/`](../Tomography/ENGINE) — core processing functions
- [`Tomography/READ/`](../Tomography/READ) — data readers
- [`Tomography/PPROCESS/`](../Tomography/PPROCESS) — pre-processing utilities
- [`RayTracer/`](../RayTracer) — 3-D ray tracing routines
- [`INTOMO_variables.md`](../INTOMO_variables.md) — original variable reference
- [`Instruction.md`](../Instruction.md) — original instruction document
