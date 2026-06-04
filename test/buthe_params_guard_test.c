// init_buthe must RETURN a fatal Lerror_t (ERR_BUTHE_PARAMS) when (b,h) violate
// h < 2*pi*b/5; never call exit(), which would kill the host process. The
// guard is unreachable via the public API for degree <= 9 (deg-9 h=8<2*pi*b/5),
// so we drive init_buthe directly with a tiny L->B in a forked child: broken
// tree exit(0)s inside init_buthe (parent sees status 0); fixed tree returns the
// flag and the child _exit(42)s. A positive loop then confirms the shipped
// params (B=512/r) never trip the guard. Asserts on the child exit code /
// returned flag. Exit 0 = pass.
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <flint/arb.h>
#include "glfunc_internals.h"

#define CHILD_OK  42   // fixed tree: init_buthe returned a fatal ERR_BUTHE_PARAMS
#define CHILD_BAD 43   // fixed tree: returned, but wrong or non-fatal flag

static Lfunc *mk(double B, uint64_t degree) {
  Lfunc *L = (Lfunc *)calloc(1, sizeof(Lfunc));
  assert(L);
  L->degree = degree;
  arb_init(L->pi); arb_const_pi(L->pi, 200);
  arb_init(L->B);  arb_set_d(L->B, B);
  return L;
}
static void release(Lfunc *L) {     // clear mk's fields + the 6 init_buthe inits
  arb_clear(L->buthe_Wf); arb_clear(L->buthe_Winf); arb_clear(L->buthe_Ws);
  arb_clear(L->buthe_b);  arb_clear(L->buthe_C);    arb_clear(L->buthe_h);
  arb_clear(L->pi); arb_clear(L->B); free(L);
}

int main(void) {
  // ---- negative: tiny B trips the guard; init_buthe must return, not exit ----
  pid_t pid = fork();
  assert(pid >= 0);
  if (pid == 0) {                                  // child
    Lfunc *L = mk(1.0, 2);                          // b=B/8=0.125 -> 2*pi*b/5 << 8
    Lerror_t e = init_buthe(L, 200);                // broken tree: exit(0) here
    int code = ((e & ERR_BUTHE_PARAMS) && fatal_error(e)) ? CHILD_OK : CHILD_BAD;
    release(L);
    _exit(code);
  }
  int st = 0; assert(waitpid(pid, &st, 0) == pid);
  assert(WIFEXITED(st));                            // not killed by a signal
  if (WEXITSTATUS(st) == 0) {
    fprintf(stderr, "FAIL: init_buthe terminated the process (exit(0) not removed)\n");
    return 1;                                       // <-- fails on the broken tree
  }
  assert(WEXITSTATUS(st) == CHILD_OK);              // fixed tree, correct fatal flag

  // ---- positive: shipped params (B = 512/r) never trip, all degrees 2..9 ----
  for (uint64_t r = 2; r <= MAX_DEGREE; r++) {
    Lfunc *L = mk(512.0 / (double)r, r);
    assert(init_buthe(L, 200) == ERR_SUCCESS);
    release(L);
  }
  printf("PASS: init_buthe returns ERR_BUTHE_PARAMS (no exit); shipped params OK.\n");
  return 0;
}
