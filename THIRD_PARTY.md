# Third-party code

## `third_party/smooth2a.m`

* **What:** `smooth2a` - 2-D matrix smoothing with a (2*Nr+1) x (2*Nc+1) moving-average kernel.
* **Author:** Greg Reeves, Division of Biology, Caltech, 2009 (the author's name and the original
  copyright header are kept inside the file).
* **Licence:** none stated in the file. It is not a licence that allows redistribution without
  conditions, so **check before publishing**: either keep it with this attribution (the usual
  courtesy, and the file is only a plotting helper), or replace it with an own few lines /
  `conv2` / `imgaussfilt`-style smoothing, or drop it.
* **Used by:** `src/meshare.m` only, i.e. the optional field plots (strain and martensite
  fractions over the beam). The simulations themselves do not need it.
