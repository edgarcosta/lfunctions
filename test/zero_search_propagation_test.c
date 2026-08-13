#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>

#define sign zero_search_test_sign
#define direction zero_search_test_direction
#define isolate_zero zero_search_test_isolate_zero
#define stat_point zero_search_test_stat_point
#define find_zeros zero_search_test_find_zeros
#define upsample_stride zero_search_test_upsample_stride
#include "../src/zeros.c"
#undef upsample_stride
#undef find_zeros
#undef stat_point
#undef isolate_zero
#undef direction
#undef sign

typedef struct {
  bool success;
  long value;
} script_step_t;

static const script_step_t *active_script;
static size_t active_script_length;
static size_t active_script_position;

static void set_script(const script_step_t *script, size_t length)
{
  active_script = script;
  active_script_length = length;
  active_script_position = 0;
}

bool zero_search_test_upsample_stride(arb_ptr result, arb_ptr t0,
                                      Lfunc *L, uint64_t side,
                                      uint64_t prec)
{
  (void)t0;
  (void)L;
  (void)side;
  (void)prec;
  assert(active_script_position < active_script_length);
  script_step_t step = active_script[active_script_position++];
  if (!step.success)
    return false;
  arb_set_si(result, step.value);
  return true;
}

typedef struct {
  Lfunc L;
  arb_t samples[4];
  arb_t zeros[MAX_ZEROS];
} zero_fixture_t;

static void fixture_init(zero_fixture_t *fixture)
{
  fixture->L.degree = 1;
  fixture->L.fft_NN = 16;
  fixture->L.wprec = DEFAULT_TARGET_PREC;
  arb_init(fixture->L.one_over_A);
  arb_one(fixture->L.one_over_A);
  arb_init(fixture->L.pi);
  arb_const_pi(fixture->L.pi, DEFAULT_TARGET_PREC);
  for (size_t i = 0; i < 4; i++)
    arb_init(fixture->samples[i]);
  for (size_t i = 0; i < MAX_ZEROS; i++)
    arb_init(fixture->zeros[i]);
  fixture->L.u_values_off[0] = fixture->samples;
  fixture->L.zeros[0] = fixture->zeros;
}

static void fixture_clear(zero_fixture_t *fixture)
{
  for (size_t i = 0; i < MAX_ZEROS; i++)
    arb_clear(fixture->zeros[i]);
  for (size_t i = 0; i < 4; i++)
    arb_clear(fixture->samples[i]);
  arb_clear(fixture->L.pi);
  arb_clear(fixture->L.one_over_A);
}

static void set_samples(zero_fixture_t *fixture, long a, long b,
                        long c, long d)
{
  arb_set_si(fixture->samples[0], a);
  arb_set_si(fixture->samples[1], b);
  arb_set_si(fixture->samples[2], c);
  arb_set_si(fixture->samples[3], d);
}

static void test_second_midpoint_failure(void)
{
  static const script_step_t script[] = {{true, 2}, {false, 0}};
  zero_fixture_t fixture = {0};
  arb_t z1, z2;
  fixture_init(&fixture);
  set_samples(&fixture, 1, 2, 3, 4);
  arb_init(z1);
  arb_init(z2);
  set_script(script, sizeof(script) / sizeof(script[0]));

  Lerror_t status = zero_search_test_stat_point(
      z1, z2, 1, &fixture.L, 0, DEFAULT_TARGET_PREC, true);
  assert(status == ERR_STAT_POINT);
  assert(active_script_position == active_script_length);

  arb_clear(z2);
  arb_clear(z1);
  fixture_clear(&fixture);
}

static void test_second_isolation_failure(void)
{
  static const script_step_t script[] = {
      {true, 2}, {true, -1}, {true, 0}, {false, 0}};
  zero_fixture_t fixture = {0};
  arb_t z1, z2;
  fixture_init(&fixture);
  set_samples(&fixture, 1, 2, 3, 4);
  arb_init(z1);
  arb_init(z2);
  set_script(script, sizeof(script) / sizeof(script[0]));

  Lerror_t status = zero_search_test_stat_point(
      z1, z2, 1, &fixture.L, 0, DEFAULT_TARGET_PREC, true);
  assert(status == ERR_UPSAMPLE);
  assert(active_script_position == active_script_length);

  arb_clear(z2);
  arb_clear(z1);
  fixture_clear(&fixture);
}

static void test_warning_does_not_commit_stale_outputs(void)
{
  static const script_step_t seed_script[] = {
      {true, -1}, {true, 0}, {true, 0}};
  static const script_step_t warning_script[] = {{true, 0}};
  zero_fixture_t fixture = {0};
  fixture_init(&fixture);
  set_samples(&fixture, 3, 1, 3, 4);

  set_script(seed_script,
             sizeof(seed_script) / sizeof(seed_script[0]));
  assert(zero_search_test_find_zeros(&fixture.L, 0) == ERR_SUCCESS);
  assert(active_script_position == active_script_length);
  assert(!arb_is_zero(fixture.zeros[0]));
  assert(!arb_is_zero(fixture.zeros[1]));

  set_script(warning_script,
             sizeof(warning_script) / sizeof(warning_script[0]));
  assert(zero_search_test_find_zeros(&fixture.L, 0) == ERR_DBL_ZERO);
  assert(active_script_position == active_script_length);
  assert(arb_is_zero(fixture.zeros[0]));
  assert(arb_is_zero(fixture.zeros[1]));

  fixture_clear(&fixture);
}

int main(int argc, char **argv)
{
  if (argc == 1 || strcmp(argv[1], "midpoint") == 0)
    test_second_midpoint_failure();
  if (argc == 1 || strcmp(argv[1], "isolation") == 0)
    test_second_isolation_failure();
  if (argc == 1 || strcmp(argv[1], "stale") == 0)
    test_warning_does_not_commit_stale_outputs();
  return 0;
}
