# M2 fidelity review — 14-item traps checklist vs fesom_ice_maevp.c

Walked 2026-06-11 against the final M2 diff (plan
`docs/plans/20260611-mevp-sea-ice-rheology.md`, Fidelity-traps checklist).
Fortran reference: `port2/fesom2/src/ice_maEVP.F90`, `EVPdynamics_m` :429-882.
C: `src/fesom_ice_maevp.c` (line numbers at review time).

| # | Trap | C evidence | Verdict |
|---|------|-----------|---------|
| 1 | `rdt = FULL ice_dt`; drag carries rdt; `drag·u_w` outside the `rdt·(…)` group | `rdt = ice->ice_dt` (:112, no `/evp_rheol_steps`); `drag = rdt*cd_oce_ice*umod*DENSITY_0*inv_thickness` (:290); `rhsu = u_ice + drag*u_w + rdt*(inv_thickness*stress_x + u_rhs) + beta_evp*u_aux` (:294-296) | ✅ |
| 2 | 0.5 in sigma11/22 update only — NOT in pressure_fac, NOT in sigma12 | `pressure_fac = det2*pstar*msum*exp(…)` no 0.5 (:203-204); `sigma12 = det1*sigma12 + pressure*eps12*vale` no 0.5 (:258); `sigma11/22 = det1*… + 0.5*pressure*(…)` (:259-260) | ✅ |
| 3 | NO theta_io rotation | no `theta_io`/`ax`/`ay` in the file; drag term uses `u_w` directly (:294, :297) | ✅ |
| 4 | Element mask `msum>0.01` (3-vertex MEAN of m_ice); node mask `a_ice>=0.01` | `msum = (m0+m1+m2)*val3; if (msum > 0.01)` (:199-201); `if (a_ice[i] >= 0.01)` (:174). No per-vertex m/a>0 test (that is std-EVP's mask) | ✅ |
| 5 | Mass regularization verbatim `M/((1+M²)·area)` | `mass[i] = ms / ((1.0 + ms*ms) * area[i*nl+0])` (:179) | ✅ |
| 6 | ssh2rhs assembly UNGUARDED (all 3 nodes); stress-divergence assembly GUARDED `< myDim_nod2D` | unguarded writes rhs_a/rhs_m at n0/n1/n2 (:158-160); guarded `if (eln[k] < mesh->myDim_nod2D)` (:266) | ✅ |
| 7 | rhs_a/rhs_m zeroed over myDim ONLY; `/area` scaling inside the ice branch only | zero loop `row < myDim_nod2D` (:142-145); `rhs_a[i] /= area` inside `if (a_ice>=0.01)` (:182-183) | ✅ |
| 8 | aux init AND final copy cover myDim+eDim; no exchange after final copy | both loops `row < N` with `N = myDim+eDim` (:117-120, :353-356); nothing after the final copy but the UF dump | ✅ |
| 9 | `bc_index_nod2D` multiplies det; build = owned edges, `myList_edge2D > edge2D_in` convention | `det = bc_index_nod2D[i] / (obd²+rf²)` (:304); M1 rebuild in fesom_ice.c: loop `ed < myDim_edge2D`, `myList_edge2D[ed] <= edge2D_in → continue` (the std-EVP BC convention; replaces the multi-rank-unsafe `edge_tri[2nd]<0` test) | ✅ |
| 10 | `delta_min` ADDITIVE | `pressure = pressure_fac/(delta + delta_min)` (:254). (std-EVP's `max(delta, delta_min)` correctly NOT copied) | ✅ |
| 11 | sigma NOT zeroed on entry — persists across calls | no sigma zeroing anywhere in fesom_ice_maevp.c; ice->work.sigma* are calloc-zeroed once at init (= Fortran alloc-time zeros); in an mEVP run no other code writes them | ✅ |
| 12 | -r8 literals double; index translations | all literals double (`1.0/3.0`, `0.5`, `0.01`, `9.0`, `4.0` — no float suffixes); `area(1,i)` → `area[i*nl+0]` (:179); `gradient_sca(1:3/4:6)` → `gs[0..2]/gs[3..5]` (:153, :222); `meancos = val3*metric_factor[el]` (:224) | ✅ |
| 13 | Non-ice nodes SKIPPED, not zeroed (NO else) | node solve `if (ice_nod[i]) { … }` with no else (:283-310); a non-ice node's u_ice_aux keeps its entry value through the final copy | ✅ |
| 14 | `uice_old`/`vice_old` never touched | `u_old`/`uice_old` do not appear in fesom_ice_maevp.c | ✅ |

Extra checks (not in the 14 but verified):
- eps11/eps22/eps12 formulas with metric terms match :702-705 term-by-term
  (eps11 `− Σv·meancos`; eps12 `0.5·(Σ(dy·u+dx·v) + Σu·meancos)`).
- delta = `sqrt(eps1² + vale·(eps2² + 4·eps12²))` matches :712.
- Velocity solve det/numerators match :814-817 sign-for-sign
  (`+rf·rhsv` / `−rf·rhsu`); umod uses AUX velocity (:806).
- rhsu reads `u_ice[i]` (entry value, NOT aux) per :810.
- Edge BC zeroes the AUX arrays at both endpoints, owned edges only
  (:824-838); cavity block (use_cavity) not ported, same exclusion as std-EVP.
- Exchange order: two blocking `fesom_exchange_nod2D(aux)` then the rhs zero
  loop — result-identical collapse of Fortran's begin/zero/end overlap
  (:861-872); the zero loop touches only u_rhs/v_rhs.
- Main-loop element body entirely inside `if (ice_el)`; non-ice elements keep
  stale eps/sigma (never read — same as Fortran).
- Element/node precomputes zero-then-cavity-cycle order preserved (:626-631,
  :653-658).
- Dump points mirror the Fortran M0 patch exactly: Q(a,m,msnow,elev),
  U0(u,v), F(stress_x,stress_y,u_w,v_w), P node(inv_thickness,mass,rhs_a,
  rhs_m), P elem(pressure_fac), it1/it2/it60/it{steps}(u_aux,v_aux +
  sigma11/12/22, after edge-BC before exchange), UF(u,v) — same filenames,
  call-counter gate ≤4 on both sides.
- Init-order safety: the bc_index_nod2D rebuild reads partit->myList_edge2D;
  fesom_partit population (fesom_main.c:~274-283, serial synth identity)
  precedes fesom_ice_init (fesom_main.c:382) on both paths.
