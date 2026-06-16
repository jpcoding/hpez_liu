# Isoline (zero-crossing) error control — SZ3 vs QoZ/HPEZ vs SPERR how-to

How to build and run an **isoline-preserving** compression case on three
compressors, preserving the single isovalue **0** on the hurricane field
`Pf48.bin.f32` (dims `100×500×500`, float32, 25,000,000 values).

The "isoline" QoI (`qoi=4`) shrinks the local error bound near an isovalue so
that the decompressed data never crosses it:
`interpret_eb(x) = min(global_eb, |x − isovalue|)`, enforced by
`check_compliance: (x − iso)·(x̂ − iso) ≥ 0`. The isoline source is the same
algorithm in all three (`Isoline.hpp`); only the surrounding framework differs:

- **SZ3** (`szcompressor/SZ3@qoi_error_control`) and **QoZ** (`HPEZ-QoZ2.0`)
  enforce it through their QoI quantizer.


- **SPERR** (`jpcoding/SPERR@qoi_isoline`) — isoline ported in for this work —
  enforces it through SPERR's existing per-point **outlier-correction** path
  (the QoI is given a `<T>`-only interface, no spatial `N`).

---

## 1. Setting isovalues

### QoZ / HPEZ (this repo) — explicit `isovalues=` config key

QoZ now supports an explicit, comma-separated `isovalues` key under
`[QoISettings]`. This is the recommended way:

```ini
[QoISettings]
qoi=4
isovalues=0              ; single zero-crossing isoline
; isovalues=-100,0,250   ; multiple isolines, all preserved
```

If `isovalues` is omitted it defaults to `{0}`. Implementation: the key is
parsed in `Config::loadcfg` (`include/QoZ/utils/Config.hpp`, fills
`conf.isovalues` and sets `qoiIsoNum`), and consumed in `GetQOI` `case 4:`
(`include/QoZ/qoi/QoIInfo.hpp`). The stock SZ3 auto-generation path (evenly
spaced isovalues from `qoiIsoNum`) is commented out in QoZ, so the explicit
list is the mechanism.

### SPERR (`jpcoding/SPERR@qoi_isoline`) — explicit `--qoi_isovalues` CLI flag

The isoline port adds a CLI flag (no config file): pass `--qoi_id 4` and a
comma-separated `--qoi_isovalues`. `--qoi_tol` only needs to be > 0 to engage
the QoI path (the isoline ignores its value).

```bash
--qoi_id 4 --qoi_tol 1e-2 --qoi_isovalues 0            # single zero-crossing isoline
--qoi_id 4 --qoi_tol 1e-2 --qoi_isovalues -100,0,250   # multiple isolines
```

If `--qoi_isovalues` is omitted it defaults to `{0}`. Implementation: parsed in
`utilities/sperr3d.cpp` into `QoIMeta.isovalues`, consumed in `GetQOI` `case 4:`
(`include/qoi/QoIInfo.hpp`) which builds `QoI_Isoline<T>` from
`include/qoi/Isoline.hpp`.

### SZ3 (`qoi_error_control`) — still needs a code edit

Stock SZ3 only exposes `qoi=4` + `qoiIsoNum=N` and **auto-generates** `N`
evenly-spaced isovalues across `[min, max]` — there's no explicit isovalue
config key. To preserve a specific value (e.g. 0), edit `GetQOI` `case 4:`
(`include/SZ3/qoi/QoIInfo.hpp`):

```cpp
case 4:{
    // FORCED: preserve the single isovalue 0 (zero-crossing)
    std::vector<T> values;
    values.push_back((T)0);
    return std::make_shared<SZ::QoI_Isoline<T, N>>(conf.dims, values, conf.absErrorBound);
}
```

---

## 2. Build

### SZ3 (`qoi_error_control` branch)

```bash
git clone --depth 1 -b qoi_error_control https://github.com/szcompressor/SZ3.git SZ3_qoi
cd SZ3_qoi
# (apply the case-4 edit above to include/SZ3/qoi/QoIInfo.hpp)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j4
# -> build/test/sz
```

**macOS Clang fix (if it fails to compile):** `std::max(double, float)` type
mismatch in `include/SZ3/api/impl/SZInterp.hpp:72` and
`include/SZ3/api/impl/SZLorenzoReg.hpp:153`. Cast the second arg:
```cpp
max_abs_eb = std::max(max_abs_eb, (double)fabs(sampling_data[i] - samples[i]));
```

### QoZ / HPEZ (this repo)

QoZ requires **SymEngine** (+ FLINT/GMP/MPFR) for its symbolic QoIs — not used
by the isoline path, but the build links them.

```bash
brew install symengine            # pulls in flint, gmp, mpfr
# (apply the case-4 edit above to include/QoZ/qoi/QoIInfo.hpp)
GMP=$(brew --prefix gmp); MPFR=$(brew --prefix mpfr)
SE=$(brew --prefix symengine);  FLINT=$(brew --prefix flint)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 \
  -DCMAKE_PREFIX_PATH="$SE" \
  -DCMAKE_EXE_LINKER_FLAGS="-L$GMP/lib -L$MPFR/lib -L$SE/lib -L$FLINT/lib -lflint -lmpfr -lgmp"
cmake --build build -j4
# -> build/test/hpez
```

Notes: `-DCMAKE_POLICY_VERSION_MINIMUM=3.5` is needed for SymEngine's old CMake
config; the `-lflint` link flag resolves `__fmpz_clear_mpz` (inlined from
SymEngine headers).

### SPERR (`jpcoding/SPERR@qoi_isoline`)

SPERR also pulls in SymEngine for its QoI framework, so it needs the same deps
plus `libmpc` and `zstd`. On macOS the include/lib dirs must be passed
explicitly (and `zstd` must be linked, since SPERR links it `INTERFACE` and
macOS shared libs require all symbols resolved):

```bash
git clone -b qoi_isoline git@github.com:jpcoding/SPERR.git SPERR_qoi && cd SPERR_qoi
brew install symengine libmpc zstd     # + flint, gmp, mpfr (symengine deps)
GMP=$(brew --prefix gmp); MPFR=$(brew --prefix mpfr); MPC=$(brew --prefix libmpc)
SE=$(brew --prefix symengine); FLINT=$(brew --prefix flint); ZSTD=$(brew --prefix zstd)
INC="-I$GMP/include -I$MPFR/include -I$MPC/include -I$FLINT/include -I$SE/include -I$ZSTD/include"
LNK="-L$GMP/lib -L$MPFR/lib -L$MPC/lib -L$SE/lib -L$FLINT/lib -L$ZSTD/lib -lflint -lmpc -lmpfr -lgmp -lzstd"
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_UNIT_TESTS=OFF \
  -DCMAKE_POLICY_VERSION_MINIMUM=3.5 -DCMAKE_PREFIX_PATH="$SE" \
  -DCMAKE_CXX_FLAGS="$INC" \
  -DCMAKE_EXE_LINKER_FLAGS="$LNK" -DCMAKE_SHARED_LINKER_FLAGS="$LNK"
cmake --build build -j4
# -> build/bin/sperr3d
```

---

## 3. Config files

Both use an SZ2/3-style `.config`. Isoline = `qoi=4`.

`iso_sz3.config` (SZ3 — value comes from the code edit above; `qoiIsoNum` only
controls how many auto-isovalues, ignored once you hard-code one):
```ini
[QoISettings]
qoi=4
qoiEB=1e-2
qoiIsoNum=1
```

`iso_qoz.config` (QoZ — set the isovalue(s) explicitly):
```ini
[QoISettings]
qoi=4
qoiEB=1e-2
isovalues=
```

SPERR uses no config file — everything is on the command line (see §4).

---

## 4. Run (compress + decompress + stats in one shot)

`-i` input, `-o` decompressed output, `-3 nx ny nz`, `-c` config,
`-M REL 0.001` value-range relative bound (cap; the isoline tightens it near 0),
`-a` print stats. For this field the value range is 6636.14, so
**REL 0.001 ≡ abs bound ≈ 6.6361**.

```bash
DATA=/Users/pjiao/data/hurricane_100x500x500/Pf48.bin.f32

# SZ3
SZ3_qoi/build/test/sz   -f -i $DATA -o Pf48.sz3.out -3 500 500 100 -c iso_sz3.config -M REL 0.001 -a

# QoZ/HPEZ  (-q 0 = "SZ3.1, no QoZ features", for a fair baseline comparison)
build/test/hpez         -f -i $DATA -o Pf48.qoz.out -3 500 500 100 -c iso_qoz.config -M REL 0.001 -q 0 -a
```

SPERR has no relative mode — only `--pwe` (absolute point-wise error). To match
`REL 0.001` you compute the bound externally as `0.001 × value_range`
(`= 6.6361` here) and pass it to `--pwe`. Isovalues go through `--qoi_id 4
--qoi_isovalues 0`; `--qoi_tol` just needs to be > 0 to engage the QoI path
(the isoline ignores its value):

```bash
# SPERR (--pwe = 0.001 * 6636.14 = 6.6361, computed externally)
SPERR_qoi/build/bin/sperr3d -c --ftype 32 --dims 500 500 100 \
    --pwe 6.6361 --qoi_id 4 --qoi_tol 1e-2 --qoi_isovalues -100,0,250 \
    --bitstream Pf48.sperr.stream --decomp_f Pf48.sperr.out --print_stats $DATA
```

Look for `isovalues: 0.0000…` in the output to confirm the forced isovalue.

`-q` levels for `hpez`: 0=SZ3.1, 1=QoZ1, 2=HPEZ-L2, 3=HPEZ-L3, 4=HPEZ-L4.

---

## 5. Verify zero-isoline preservation (repo-independent)

The tools' own `-a`/`--print_stats` checks can be misleading: **SZ3's "different
cells for isovalue …" reports against an auto-computed comparison isovalue
(≈ −93.67), NOT 0** — ignore it. QoZ and SPERR both print a correct isoline
check (`Max qoi error = 0`). The authoritative, tool-independent check counts
points whose sign relative to 0 flips (exactly the `check_compliance` violation
`(x)(x̂) < 0`):

```python
import numpy as np
o = np.fromfile('/Users/pjiao/data/hurricane_100x500x500/Pf48.bin.f32', dtype=np.float32)
for name, f in [('SZ3','Pf48.sz3.out'), ('QoZ','Pf48.qoz.out'), ('SPERR','Pf48.sperr.out')]:
    d = np.fromfile(f, dtype=np.float32)
    viol = int(np.sum((o>0)&(d<0)) + np.sum((o<0)&(d>0)))   # opposite sides of 0
    rng = float(o.max()-o.min()); mse = float(np.mean((o-d)**2))
    psnr = 20*np.log10(rng) - 10*np.log10(mse)
    print(f'{name}: maxAbsErr={np.max(np.abs(o-d)):.3f} PSNR={psnr:.2f} zero-crossing violations={viol}')
```

---

## 6. Results (REL 0.001 ≡ abs bound 6.6361)

| Compressor | Compression ratio | PSNR | Max abs err | Zero-crossing violations |
|---|---|---|---|---|
| SZ3 `qoi_error_control` | 65.5 | 72.6 | 4.0  | **0 / 25,000,000** |
| QoZ `HPEZ-QoZ2.0` (`-q 0`) | 95.5 | 69.2 | 6.6  | **0 / 25,000,000** |
| SPERR `qoi_isoline` (`--pwe 6.6361`) | 87.2 | 76.7 | 6.6  | **0 / 25,000,000** |

### QoZ level sweep (`-q 0…4`, REL 0.001)

| `-q` | Level | CR | PSNR | Max abs err | Zero-crossing violations |
|---|---|---|---|---|---|
| 0 | SZ3.1   | 95.5 | 69.2 | 6.64 | **0 / 25,000,000** |
| 1 | QoZ1    | 71.2 | 70.8 | 6.64 | **0 / 25,000,000** |
| 2 | HPEZ-L2 | 69.1 | 71.1 | 6.64 | **0 / 25,000,000** |
| 3 | HPEZ-L3 | 69.1 | 71.1 | 6.64 | **0 / 25,000,000** |
| 4 | HPEZ-L4 | 90.5 | 71.2 | 6.64 | **0 / 25,000,000** |

---

## 7. Verify QoI with the built-in `qoi_val` tool (QoZ)

`test/qoi_val.cpp` builds to `build/test/qoi_val`. It re-evaluates the QoI on
both the original and the decompressed file and reports the QoI error — for the
isoline (`qoi=4`) it confirms the isovalue is preserved (`Max qoi error = 0`).
Unlike the `-a` stat printed during compression, this is a standalone check you
can run on any pair of files.

Give it the **same config** used for compression (so it picks up `qoi=4` and
`isovalues=0`), the data type (`-f`), and the dims (`-3 nx ny nz`). `-i` is the
original input, `-o` is the decompressed output:

```bash
DATA=/Users/pjiao/data/hurricane_100x500x500/Pf48.bin.f32

build/test/qoi_val -f -3 500 500 100 -c iso_qoz.config -i $DATA -o Pf48.qoz.out
```

Look for `Max qoi error = 0` in the `QoI error info:` block to confirm the
zero-isoline was preserved (this matches the zero-crossing-violations count from
§5). The tool also prints QoI PSNR/NRMSE and, for 3D data, QoI Laplacian and
gradient-length errors.
