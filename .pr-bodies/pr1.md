**Summary**

Ports pyMDMix to Python 3 and Biskit 3 so that it installs from `environment.yml` with current numpy, scipy and AmberTools.

**Changes**

- Normalization of leading tabs, automatic 2to3 conversion and manual syntax and semantic fixes.
- Port to Biskit 3 and to the Python 3 standard library (`subprocess`, `multiprocessing`, `pickle`, text and binary modes).
- `GridData`: avoid elementwise array-versus-string comparison (numpy >= 1.25).
- Projects, replicas, systems and the solvent database pickled by Python 2 can be loaded: latin1 decoding, and `Biskit 2` class references (module `Biskit.PDBModel`, ...) are translated to their Biskit 3 modules and subpackages (`biskit.pdbModel`, `biskit.core.localpath`, ...).
- Private temporary folder in `tools.tempDir()` (the tests no longer remove the shared system temporary folder), endian-aware DCD reader, ASCII and NetCDF trajectories.
- Python 3.12+ and numpy 2: `distutils.dir_util.create_tree` (removed in 3.12) is replaced by `tools.create_tree`; `if m.chainBreaks():` no longer relies on the truth value of an array.
- Biskit 3.0.1 writes integer PDB fields with `repr(value)`; with numpy >= 2, `repr(numpy.int32(1))` is `np.int32(1)` and every PDB written by `PDBModel.writePdb` is corrupted. `tools.py` makes the `repr` used by Biskit's `FortranFormat` convert numpy scalars to Python numbers (limited to that module, no global numpy option is changed).
- Packaging (`setup.py`, `requirements.txt`, `environment.yml`, Dockerfile) and documentation for Python 3.

**Compatibility**

Python 2 is no longer supported.

**Validation**

- `python pyMDMix/test.py all` on Python 3.11.17 (numpy < 2) and Python 3.13.16 (numpy 2.5.4): 36 tests, 0 failures, 1 skipped (see PR 4).
- Projects created with Python 2.7.15, pyMDMix 0.2.8 and Biskit 2.4.3 (an empty project and a project with two replicas, `pep` system, TIP3P) open with this branch: system, force fields and replicas are listed by `mdmix info project`.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_016M7mpPRFVJcSvhJWhCpc6p
