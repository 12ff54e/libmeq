# **M**agnetic **EQ**uilibrium library

Given poloidal magnetic flux $\psi_{p}$ (analytic or numeric), calculate its Boozer coordinate representation.

2D field can be represented as BSpline or Zernike series, controlled by macro `MEQ_ZERNIKE_SERIES_`

# GEQDSK coordinate conventions

`GFileRawData` keeps flux, profiles, field and current in the file's units and
coordinate convention. Its `cocos` member records that convention. Reading a
file selects it in this order:

1. An API override: `data.read(stream, Cocos(11))`.
2. A `COCOS=11` or `COCOS_OUT=11` tag in the first line, including extra metadata.
3. A CHEASE producer name, matched without regard to case: COCOS 2.
4. EFIT, also used when the producer is unknown: COCOS 1, 3, 5 or 7, inferred
   from the signs of current, outward flux change, toroidal field and q.
   Missing signs retain the COCOS 1 defaults.

Producer detection is a guess. CHEASE can export other conventions, and EFIT
variants use different flux and q signs. Use an explicit index when known.
The reader rejects malformed COCOS tags and indices outside 1–8 and 11–18.
For data filled by hand, set `data.cocos = Cocos(2)` if it uses CHEASE's
convention; the default is COCOS 1.

`Cocos::factors_to` gives the flux, flux-derivative, toroidal and safety-factor
conversion factors. It preserves the physical field and current directions.
The flux-derivative factor applies to both `p_prime` and `f_f_prime`.
Indices 11–18 use full poloidal flux; indices 1–8 use flux divided by 2π.

`MagneticEquilibrium` converts its callback values to COCOS 2: phi runs
clockwise from above, theta runs counterclockwise in the (R,Z) plane, and
poloidal flux uses Wb/rad. It subtracts the axis flux. Contour searches use an
outward-increasing working flux, but returned flux, its derivatives, and the
Jacobian keep their signs. Profiles remain ordered from axis to edge.
With `use_si=false`, lengths use R_axis and fields use the signed F_axis/R_axis
in COCOS 2; pressure uses B_axis²/μ₀. Thus normalized flux and `psi_delta()`
can be negative. The source data are unchanged.
With `use_si=true`, pressure is in pascals, including its axis value.

The sign rules follow [Sauter and Medvedev, COCOS, Table I and Appendix C](https://crppwww.epfl.ch/~sauter/COCOS/Sauter_COORD_CONVENTIONS_COCOS_2012_updated_after_reprint_for_Appendices_and_refs.pdf).
See also the [EFIT g-file conventions](https://efit-ai.gitlab.io/efit/files.html)
and [CHEASE documentation](https://crppwww.epfl.ch/~sauter/chease/).

# Note

This repository is separated from [mct](https://github.com/12ff54e/mct), which depends this library to construct magneitc equlibrium.
