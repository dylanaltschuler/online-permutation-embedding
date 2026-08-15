# Numerical certificates for the scaling constants

The repository concerns the finite interval-arithmetic certificates in the manuscript. The manuscript uses these certified computations together with separate analytic arguments to obtain rigorous estimates on some scaling limits. 

Contained files:

```text
intervals.json   candidate interval enclosures for gamma, beta_plus, beta_minus
certify.py       rigorous FLINT/Arb checker
slacks.json      recorded certification slacks and finite bounds
generate_candidates.py   optional non-rigorous interval generator
```

## Requirement

Use Python 3 with `python-flint` installed:

```bash
python -m pip install python-flint
```

On Windows/PowerShell, replace `$py` by the Python executable for your environment:

```powershell
$py = "C:\path\to\envs\cert\python.exe"
& $py -m pip install python-flint
```

## Recurrences checked

All three dynamic programs (defining beta_k^+, beta_k^-, and gamma_k) use the same scalar kernel

```text
G(p,q) = min_{0 <= a < b <= 1}
         (1 + p log(b/a) + q log((1-a)/(1-b))) / (b-a).
```

The file `intervals.json` contains interval enclosures for the three sequences

```text
gamma[0] = 0,
gamma[k] = (1/k) sum_{i=1}^k G(gamma[i-1], gamma[k-i]),

beta_plus[0] = 0,
beta_plus[k] = max_{1 <= i <= k} G(beta_plus[i-1], beta_plus[k-i]),

beta_minus[0] = 0,
beta_minus[k] = min_{1 <= i <= k} G(beta_minus[i-1], beta_minus[k-i]).
```

## Quick check

The following command checks only the first 20 values from `intervals.json` and writes a small local slack file. It is a smoke test for the installation, not the full proof.

### Bash/macOS/Linux

```bash
python certify.py intervals.json --max-k 20 --slacks slacks_K20.json
```

### PowerShell/Windows

```powershell
& $py certify.py intervals.json --max-k 20 --slacks slacks_K20.json
```

The command should end with `wrote slacks_K20.json` and each checked quantity should print `PASS`.

## Full checks

To regenerate one combined `slacks.json`, run:

```bash
python certify.py intervals.json \
  --prec 256 \
  --jobs 8 \
  --manuscript-goals \
  --slacks slacks.json
```

The option `--jobs 8` uses eight worker processes. It changes wall-clock time only; it does not change the inequalities checked. Adjust `8` to the number of cores you want to use, or omit the option for a serial run.

For targeted reruns, check one quantity at a time:

```bash
python certify.py intervals.json --only gamma \
  --prec 256 \
  --jobs 8 \
  --manuscript-goals \
  --slacks slacks_gamma.json

python certify.py intervals.json --only beta_plus \
  --prec 256 \
  --jobs 8 \
  --manuscript-goals \
  --slacks slacks_beta_plus.json

python certify.py intervals.json --only beta_minus \
  --prec 256 \
  --jobs 8 \
  --manuscript-goals \
  --slacks slacks_beta_minus.json
```

PowerShell equivalents replace `python` by `& $py` and use backticks for line continuation if desired.

## What `certify.py` proves

Suppose `intervals.json` gives intervals `[L_k,U_k]` for a sequence.

For `gamma`, the checker verifies with outward-rounded Arb arithmetic that

```text
L_k <= (1/k) sum_i G(L_{i-1}, L_{k-i}),
U_k >= (1/k) sum_i G(U_{i-1}, U_{k-i}).
```

For `beta_plus`, it verifies

```text
L_k <= max_i G(L_{i-1}, L_{k-i}),
U_k >= max_i G(U_{i-1}, U_{k-i}).
```

For `beta_minus`, it verifies

```text
L_k <= min_i G(L_{i-1}, L_{k-i}),
U_k >= min_i G(U_{i-1}, U_{k-i}).
```

The lower `beta_minus` inclusion requires every split to be lower-bounded. To avoid unnecessary expensive scalar root certificates, `certify.py` first tries the analytic lower bound

```text
G(p,q) >= (sqrt(p) + sqrt(q))^2,
```

and, on boundary inputs, the stronger bound

```text
G(0,q) >= (sqrt(q) + 1/2)^2.
```

Splits not certified by these lower bounds are checked by the same Arb scalar certificate used for the other quantities. The upper `beta_minus` inclusion only needs one split whose certified upper endpoint is at most `U_k`; the script uses the split `i=2` first and falls back to a full scan if needed.

Since `G` is monotone in both variables, these interval inclusions imply by induction that the true dynamic-program values lie in the proposed intervals.

After certifying the `beta_minus` intervals, the checker also evaluates the
finite lower bound from the manuscript

```text
c_typ >= (2 / ((K+1)(K+2)) * sum_{j=1}^K sqrt(beta_minus[j]))^2.
```

For this estimate the checker uses `K=312` and the certified lower endpoint of
each `beta_minus[j]` interval. This cutoff is independent of the `K=5000` used
for the other finite estimates.

## Finite estimates recorded in `slacks.json`

The included `slacks.json` consolidates successful checks at 256-bit Arb
precision. It records:

```text
0.4893404633455635... <= c_typ <= 0.4996662905439521... < 0.49967,
0.5054731461429961... <= c_+     <= 0.5056753576203792...,
1/4 <= c_- <= 0.4886678999987247... < 0.48867.
```

The lower bound on `c_typ` uses the certified `beta_minus` intervals through
`K=312`; the remaining numerical bounds use `K=5000`. The lower bound
`c_- >= 1/4` is analytic, as in the manuscript; the numerical certificate
supplies the upper bound on `c_-`.

## Generating candidate intervals

The generator uses ordinary floating point arithmetic and is not part of the proof. The proof is the successful run of `certify.py` on the candidate intervals.

To regenerate all candidate intervals:

```bash
python generate_candidates.py --K 5000 --out intervals.json
```

By default, the generator uses the radius

```text
abs_pad + rel_pad*k^2
```

with `abs_pad=1e-12` and `rel_pad=1e-9`. The included `intervals.json` uses a wider `beta_minus` padding (`1e-8*k^2`) so the upper recurrence inclusion has comfortable slack while still proving `c_- < 0.48867`.

The `beta_minus` generator uses a fast extrapolation of the split `i=2` after checking a prefix by floating point full scan. This is only a way to propose candidates; `certify.py` verifies the full recurrence rigorously.

## Running time

The direct certificate is quadratic in `K` because the recurrences contain all splits `i=1,...,k`. The `beta_minus` check is faster in practice because of the analytic lower-bound screen described above. The `gamma` and `beta_plus` checks are the expensive parts of the full certificate.

## File guide

- `README.md`: this file.
- `intervals.json`: the single manuscript-aligned candidate interval file for `gamma`, `beta_plus`, and `beta_minus`.
- `certify.py`: the single rigorous checker. It reads `intervals.json` and writes `slacks.json` by default.
- `slacks.json`: recorded recurrence slacks, finite estimates, and manuscript-goal slacks.
- `generate_candidates.py`: optional non-rigorous generator for new candidate intervals.
