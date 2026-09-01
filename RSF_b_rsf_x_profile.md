# Piecewise-linear `b_rsf(x)` profile (one-phase VS–VW–VS)

Re-apply this on a clean branch. It is independent of the y–z multigrid work.

## Why

Herrendörfer-style VS–VW–VS used to be built from several RSF phases (constant `b` plus two-point linear ramps). Phase ratios at the contacts interrupt RSF. This change puts the whole along-strike `b(x)` profile on **one phase**: piecewise-linear through a knot list, clamped outside.

A two-point input is unchanged (old linear ramp is the `n = 2` case).

## Input (no length parameter)

Do **not** add `n_b_rsf`. The parser counts entries (max 16). `b_rsf_val` and `b_rsf_x` must have the same length, at least 2, and `b_rsf_x` must be strictly increasing.

```
b_rsf     = 0.001
b_rsf_val = 0.001 0.017 0.017 0.001
b_rsf_x   = -47000.0 -43000.0 33000.0 37000.0
```

`b_rsf` remains the fallback if the vectors are omitted. `b_rsf_x` is a length (scaled by `unit_length`); `b_rsf_val` is dimensionless.

Evaluation at `x`:

| interval | `b` |
|---|---|
| `x ≤ x[0]` | `val[0]` (clamp) |
| `x[i-1] < x ≤ x[i]` | linear from `val[i-1]` to `val[i]` |
| `x ≥ x[n-1]` | `val[n-1]` (clamp) |

With the four knots above: VS (`0.001`) outside, linear ramps in the 4 km bands, VW (`0.017`) in the middle.

The phase **box must cover the whole profile**. Knots outside the box never apply. For the example above use something like:

```
bounds = -75600.0 75600.0 ...
```

not a left-half box that ends at `x = 0`.

## Files to change

| file | what |
|---|---|
| `src/LaMEM.h` | `#define _max_b_rsf_knots_ 16` |
| `src/phase.h` | `Material_t` arrays of size `_max_b_rsf_knots_`, plus `nb_rsf` |
| `src/parsing.h` / `src/parsing.cpp` | `getScalarParamUpTo` (read 1..max entries, return count) |
| `src/phase.cpp` | parse vectors, check lengths / monotonic `x`, set `b_rsf_trans` |
| `src/constEq.cpp` | `getBRsf` piecewise linear + clamp |

`getBRsf` call sites stay the same (`getBRsf(mat, ctx->x_coor)`).

## 1. `src/LaMEM.h`

Next to the other `_max_*` defines:

```c
#define _max_b_rsf_knots_ 16
```

## 2. `src/phase.h` — `Material_t`

Replace the old 2-point fields:

```c
PetscScalar  b_rsf;             // used if no x-profile
PetscScalar  b_rsf_val[_max_b_rsf_knots_];
PetscScalar  b_rsf_x  [_max_b_rsf_knots_];
PetscInt     nb_rsf;            // knot count (0 if unused)
PetscInt     b_rsf_trans;       // 1 if profile is active
```

Old code had `b_rsf_val[2]`, `b_rsf_x[2]` only.

## 3. Parser: `getScalarParamUpTo`

`getScalarParam` requires **exactly** `num` values, so a 4-knot line fails. Add a sibling that accepts `1..maxnum` and returns the actual count in `nread`.

Declaration in `src/parsing.h`:

```c
PetscErrorCode getScalarParamUpTo(
		FB          *fb,
		ParamType    ptype,
		const char  *key,
		PetscScalar *val,
		PetscInt     maxnum,
		PetscInt    *nread,
		PetscScalar  scal);
```

Implementation: copy `getScalarParam`, but

- start with `*nread = 0`
- pass `maxnum` into `PetscOptionsGetScalarArray` / `FBGetScalarArray`
- if optional and not found, return 0 with `*nread = 0` (do **not** require `nval == maxnum`)
- scale only the `nval` entries that were read
- set `*nread = nval`

See current `src/parsing.cpp` function `getScalarParamUpTo`.

## 4. `src/phase.cpp` — `DBMatReadPhase`

After reading scalar `b_rsf`, replace the old 2-point block (sentinel `PETSC_MAX_REAL` on `b_rsf_x[0]`) with:

1. `b_rsf_trans = 0`, `nb_rsf = 0`
2. `getScalarParamUpTo` for `b_rsf_val` (scale `1.0`) and `b_rsf_x` (scale `scal->length`)
3. If either count is nonzero:
   - error if counts differ
   - error if count `< 2`
   - error if `b_rsf_x` is not strictly increasing
   - `nb_rsf = n_x`, `b_rsf_trans = 1`

## 5. `src/constEq.cpp` — `getBRsf`

```c
static inline PetscScalar getBRsf(Material_t *mat, PetscScalar x)
{
	PetscInt    i, n;
	PetscScalar x0, x1, t;

	if(!mat->b_rsf_trans) return mat->b_rsf;

	n = mat->nb_rsf;

	if(x <= mat->b_rsf_x[0])   return mat->b_rsf_val[0];
	if(x >= mat->b_rsf_x[n-1]) return mat->b_rsf_val[n-1];

	for(i = 1; i < n; i++)
	{
		if(x <= mat->b_rsf_x[i])
		{
			x0 = mat->b_rsf_x[i-1];
			x1 = mat->b_rsf_x[i];
			t  = (x - x0) / (x1 - x0);
			return mat->b_rsf_val[i-1] + t*(mat->b_rsf_val[i] - mat->b_rsf_val[i-1]);
		}
	}

	return mat->b_rsf_val[n-1];
}
```

For `n = 2` this is the old ramp (linear between the two knots, clamp outside).

## Checks when porting

- Two-point files (`b_rsf_val = a b`, `b_rsf_x = x0 x1`) still parse and give the same `b(x)`.
- Four-knot VS–VW–VS on **one** phase; box spans all knots.
- Mismatched vector lengths abort.
- Unsorted / repeated `b_rsf_x` abort.
- Omit both vectors → constant `b_rsf` as before.
- `unit_length` scaling: knots in metres in the `.dat` file, stored nondimensional.

## Not part of this change

- y–z geometric multigrid
- `munmap_chunk` / `asprintf` / incomplete `FDSTAGCreate` probe
- z-dependent `a_rsf` (still two-point linear)
