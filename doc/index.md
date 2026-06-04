# liblfun

`liblfun` is a C library that numerically computes data for motivic
L-functions of degree 2 through 9: the root number, the zeros on the critical
line, the analytic rank, the leading Taylor coefficient at the central point,
special values, and plot samples. The numerical methods are described in
*Computing motivic L-functions* by Bober, Booker, Costa, Lee, Platt, and
Sutherland (in preparation).

The library exposes an opaque-pointer C interface in
[`include/glfunc.h`](https://github.com/edgarcosta/lfunctions/blob/main/include/glfunc.h):
create an `Lfunc_t`, supply the Euler factors, call `Lfunc_compute`, then query
the results. See the [API reference](api.md) for the full lifecycle and every
entry point.

## Building

The build is driven by [Autoconf](https://www.gnu.org/software/autoconf/);
`./configure` is committed, so a checkout needs no autotools. See the
[README](https://github.com/edgarcosta/lfunctions/blob/main/README.md) for
dependencies, `./configure` flags, and the make targets.

To build these docs locally, from the `doc/` directory:

```sh
make html
```

`make html` runs Doxygen over `include/glfunc.h` and then Sphinx, writing HTML
to `_build/html`. It needs `sphinx-build` (with `myst-parser` and `breathe`)
and the `doxygen` binary installed; see the repository README for a one-time
virtualenv setup.

<!--
Glob toctree: every top-level page (api.md, plus topic pages a sibling PR adds,
e.g. rational.md / sympow.md) is auto-included in alphabetical order, so adding
doc/<name>.md needs no edit here and never triggers a "document isn't included
in any toctree" warning under -W. `api` sorts first among the planned pages and
is the entry point; if a page must precede it, list that page explicitly above
the glob. The glob already matches api.md, so it is never empty.
-->
```{toctree}
:maxdepth: 2
:caption: Contents
:glob:

*
```

## Indices

- {ref}`genindex`
- {ref}`search`
