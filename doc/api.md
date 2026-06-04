# API reference

```{note}
This page is a scaffold placeholder. The narrative API guide (the lifecycle
walk-through, the mu / normalisation convention, error handling, and worked
examples) is written under bead **lfunctions-rrt**, which lands on this same
branch. TODO(rrt): replace this note and the prose below; keep or extend the
Breathe directive at the bottom so the reference stays generated from the
header.
```

The public interface is declared in
[`include/glfunc.h`](https://github.com/edgarcosta/lfunctions/blob/main/include/glfunc.h).
All work flows through one lifecycle: initialise an `Lfunc_t`, supply the Euler
factors, call `Lfunc_compute`, query the results, then clear.

## Full header reference

The declarations below are generated from `include/glfunc.h` by Doxygen and
rendered through Breathe, so the header is the single source of truth. As the
header's comments are migrated to Doxygen style, their descriptions appear here
automatically.

```{doxygenfile} glfunc.h
:project: liblfun
```
