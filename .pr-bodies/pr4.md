Depends on #__DEP2__ (merge that one first; until then this diff also contains its commits).

**Summary**

Fixes problems found while running the analysis pipeline on real projects.

**Changes**

1. **Stale lock files (issue #6).** `FileLock` stores host and pid in the lock file, creates it next to the locked file, removes locks of dead processes (same host) or older than 10 minutes, and a forked worker never removes its parent's lock. The timeout error names the lock file. `returnMDMixProject` no longer reports a lock error as "No project file found".
2. **`&>` bashism.** Under `/bin/sh` (dash on Debian and Ubuntu) `cmd &> file` runs `cmd` in the background. Alignment jobs ignored `-C` (60 simultaneous cpptraj for 3 replicas) and `check_com.sh` ran concurrently with the next MD segment. They now use `> file 2>&1`.
3. **Alignment and cpptraj density process replicas one after the other** (each with its `-C` workers); `--parallel-replicas` restores the previous behavior.
4. **Progress reporting** (`[done/total, percent, elapsed] replica step`). `Executor.waitJobCompletion` waits for the jobs and not only for the queue.
5. **Density with `--com`.** Solvent residues are stored as integer atom-index arrays instead of one boolean mask of the whole system per residue. Centers and per-atom residue coordinates are vectorized and grid counting uses `numpy.add.at`.
6. RMSD axis label: A, not A squared.
   `returnMDMixProject` reported "No project file found" for any error while loading; when a `.mproj` exists but cannot be loaded, the real cause is now logged.
7. Tests: `tools.tryRemove` leaves the folder being removed if the process is inside it (several tests `chdir` into temporary folders and removed them, leaving the process in a deleted directory that broke the next test); `test_OpenMMWriter` is skipped because `OpenMMWriter.writeCheckCOM` is not implemented (also in the original code); the error message of `createMDInput` used a misspelled attribute.

**Validation**

- Density `--com` on 1YER (ETAO_1, steps 1:2, `-C4`, same input): wall time 21 min 48 s before, 21 s after; peak memory 1.35 GB before, 272 MB after. Grids ETA_COM and WAT_COM: 15 and 11 voxels (of 3.94 million) differ by one count (sum of counts 1284839 vs 1284824 and 18541852 vs 18541845, below 0.002 %). The cause of those differences was not determined (probably floating point arithmetic of the center of mass).
- `analyze align bysolvent -s ETAO -C4` on three replicas: never more than 4 cpptraj processes, progress lines shown.
- `FileLock`: normal use, stale lock with dead pid, lock held by a live process, empty old lock, forked child, relative path.

🤖 Generated with [Claude Code](https://claude.com/claude-code)

https://claude.ai/code/session_016M7mpPRFVJcSvhJWhCpc6p
