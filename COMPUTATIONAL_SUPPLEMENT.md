# Restored computational supplement

The complete working supplement is in [ancillary/](ancillary/README.md).
It includes the recovered exact rational closure routines, new standalone
certificate verifiers, fresh GAP 4.14.0 / CTblLib 1.3.9 exports, tests, and
regenerated JSON outputs.

Start with:

```sh
.venv/bin/python ancillary/reproduce.py
```

A fresh-export reproduction also passed using `--gap .tools/gap-4.14.0/gap`.
See the ancillary README for setup on another machine.

**The current manuscript has two known discrepancies:** the Fi23 decomposition
must be M_91(Q) + M_3(Q) + Q^4, and PSL(2,29) is missing from the rational-class
example list. Read [the detailed findings](ancillary/PAPER_DISCREPANCIES.md).
The Monster certificate reproduces unchanged. The PDF is not certified release
ready by this restoration.
