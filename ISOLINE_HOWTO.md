# Isoline (zero-crossing) error control — SZ3 vs QoZ/HPEZ how-to

How to build and run an **isoline-preserving** compression case on both
compressors, preserving the single isovalue **0** on the hurricane field
`Pf48.bin.f32` (dims `100×500×500`, float32, 25,000,000 values).

The "isoline" QoI (`qoi=4`) shrinks the local error bound near an isovalue so
that the decompressed data never crosses it:
`interpret_eb(x) = min(global_eb, |x − isovalue|)`, enforced by
`check_compliance: (x − iso)·(x̂ − iso) ≥ 0`. The isoline source is
functionally identical in both repos (`Isoline.hpp`); only the surrounding
framework differs.

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
isovalues=0
```

---

## 4. Run (compress + decompress + stats in one shot)

`-i` input, `-o` decompressed output, `-3 nx ny nz`, `-c` config,
`-M ABS 30` global absolute bound (cap; the isoline tightens it near 0),
`-a` print stats.

```bash
DATA=/Users/pjiao/data/hurricane_100x500x500/Pf48.bin.f32

# SZ3
SZ3_qoi/build/test/sz   -f -i $DATA -o Pf48.sz3.out -3 500 500 100 -c iso_sz3.config -M ABS 30 -a

# QoZ/HPEZ  (-q 0 = "SZ3.1, no QoZ features", for a fair baseline comparison)
build/test/hpez         -f -i $DATA -o Pf48.qoz.out -3 500 500 100 -c iso_qoz.config -M ABS 30 -q 0 -a
```

Look for `isovalues: 0.0000…` in the output to confirm the forced isovalue.

`-q` levels for `hpez`: 0=SZ3.1, 1=QoZ1, 2=HPEZ-L2, 3=HPEZ-L3, 4=HPEZ-L4.

---

## 5. Verify zero-isoline preservation (repo-independent)

The tools' own `-a` checks can be misleading: **SZ3's "different cells for
isovalue …" reports against an auto-computed comparison isovalue (≈ −93.67),
NOT 0** — ignore it. QoZ does print a correct isoline check (`Max qoi error =
0`). The authoritative check counts points whose sign relative to 0 flips
(exactly the `check_compliance` violation `(x)(x̂) < 0`):

```python
import numpy as np
o = np.fromfile('/Users/pjiao/data/hurricane_100x500x500/Pf48.bin.f32', dtype=np.float32)
for name, f in [('SZ3','Pf48.sz3.out'), ('QoZ','Pf48.qoz.out')]:
    d = np.fromfile(f, dtype=np.float32)
    viol = int(np.sum((o>0)&(d<0)) + np.sum((o<0)&(d>0)))   # opposite sides of 0
    rng = float(o.max()-o.min()); mse = float(np.mean((o-d)**2))
    psnr = 20*np.log10(rng) - 10*np.log10(mse)
    print(f'{name}: maxAbsErr={np.max(np.abs(o-d)):.3f} PSNR={psnr:.2f} zero-crossing violations={viol}')
```

---

## 6. Results (global ABS = 30)

| Compressor | Compression ratio | PSNR (dB) | Max abs err | Zero-crossing violations |
|---|---|---|---|---|
| SZ3 `qoi_error_control` | 80.5  | 67.8 | 8.0  | **0 / 25,000,000** |
| QoZ `HPEZ-QoZ2.0` (`-q 0`) | 137.4 | 60.0 | 28.2 | **0 / 25,000,000** |

Both preserve the zero isoline perfectly. The CR/PSNR difference comes from the
surrounding framework (QoZ's eb auto-tuning settled on a looser effective bound
→ higher CR, lower PSNR), **not** from the isoline logic, which is identical.

Data composition: 92.6% > 0, 7.0% < 0, 0.4% exactly 0 — so the zero isoline is
physically meaningful for this field.
