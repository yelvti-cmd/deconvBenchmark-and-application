# Third-party method dependencies

The repository contains study-specific scripts and parameters, not copies of
third-party method source code. Install or obtain each method from its upstream
repository and cite the corresponding publication.

| Method | Upstream source | Recorded version | Repository treatment |
|---|---|---:|---|
| BayesPrism | <https://github.com/Danko-Lab/BayesPrism> | 2.2.2 | Called by the included representative wrapper |
| bMIND/MIND | <https://github.com/randel/MIND> | MIND 0.3.3 | Install the upstream R package; do not treat `prior_function_new.R` as study-authored code |
| swCAM | <https://github.com/Lululuella/swCAM> | Record the exact upstream commit used | Obtain `sCAMfastNonNeg.R` from upstream; do not vendor it as study-authored code |
| Scissor | <https://github.com/sunduanchen/Scissor> | 2.0.0 | Install the upstream R package; do not treat `scissor_function.R` as study-authored code |

For a frozen release, record the exact release tag or commit hash used for any
GitHub dependency that does not provide a stable package version. General
software behavior and installation requirements are governed by the upstream
projects; this repository documents the study-specific inputs and parameters.
