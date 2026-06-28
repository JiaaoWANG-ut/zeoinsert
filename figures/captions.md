# Figure Captions (Nature-style, English)

Full captions for Figures 1–4. Source-level summaries also appear in each
`fig*.py` module docstring.

---

## Figure 1 | Accessibility-aware packing workflow and pore map in FAU

**Accessibility-aware Monte Carlo packing workflow and steric pore map in FAU zeolite.**

(a) Four-step workflow: (1) host framework and guest library; (2) periodic
flood-fill on a voxel grid to classify accessible vs sterically blocked voids
(red sodalite cages); (3) simulated-annealing Monte Carlo insertion with
blocked-pore penalty; (4) validated host–guest model with guests confined to
accessible channels. (b) Steric accessibility map of FAU showing eight
inaccessible sodalite β-cages (red) that are topologically connected in a
perfect crystal but closed to guest-sized probes (probe radius 2.5 Å). These
regions are invisible to geometry-only packers. (c) Two-dimensional voxel
cross-section (z-slice) of the unit cell: grey, solid framework; blue,
accessible void; red, blocked void. (d) Volume partition of the FAU unit cell
(grid 64³): 67% solid, 31% accessible, 2% sterically blocked.

---

## Figure 2 | Probe physics, ablation, and convergence

**Guest-size-dependent pore closure and ablation of accessibility-aware packing.**

(a) Accessible cell volume as a function of probe radius for five frameworks
(FAU, LTL, ERI, OFF, MAZ). Vertical dotted lines mark effective radii of Li⁺,
H₂O, CO₂, N₂, and EC. (b) Fraction of void volume that becomes inaccessible
(in closed cages) as probe radius increases; ERI and OFF show sharp closure
above ~3 Å. (c) Fraction of CO₂ guests placed in closed cages vs loading in FAU;
accessibility-aware packing (green) remains at 0% across loadings 8–48, whereas
accessibility-blind packing (red) yields 4–12.5% misplacement (mean ± s.d., five
seeds). (d) OVITO render of accessibility-aware packing: 32 CO₂ molecules in
accessible supercages/channels only. (e) Accessibility-blind packing: misplaced
guests highlighted in red trapped in sodalite cages. (f) Simulated-annealing
convergence for 70 CO₂ molecules in FAU (five seeds): overlap penalty energy
(blue, symlog scale) and temperature schedule (grey dashed, right axis).

---

## Figure 3 | Application cases vs Packmol

**Electrolyte, CO₂ adsorption, and CO₂/H₂O separation compared with Packmol.**

(a) Accessibility-aware packing in FAU for three cases (top row): Li-salt
electrolyte (EC/DMC/LiPF₆, 20 molecules), CO₂ adsorption (32 CO₂), and
CO₂/H₂O competitive loading (16+16). (b) Packmol baseline (bottom row); guests
violating steric accessibility or minimum-distance cutoffs highlighted in red.
(c) Number of guests whose centres fall in closed cages: ours = 0 for all cases;
Packmol = 3 (electrolyte), 4 (CO₂), 5 (separation). (d) Stacked physical
violations for Packmol: closed-cage placement (red), framework clash (orange),
guest–guest overlap (gold); green diamonds mark zero total violations for our
method. (e) Distribution of minimum guest–framework distances across all cases;
vertical dotted line at 2.2 Å cutoff. Packmol shows a tail below the cutoff;
ours does not.

---

## Figure 4 | Cross-framework generalization

**Generalization across zeolites, MOFs, and COFs.**

(a–c) OVITO renders of accessibility-aware CO₂ packing in representative hosts:
MOF-5 (MOF), UiO-66 (MOF), and COF-5 (COF). (d) Heat map of achieved loading
(molecules per unit cell) for six frameworks (FAU, LTL, ERI, MOF-5, UiO-66,
COF-5) and four guests (CO₂, H₂O, EC, LiPF₆); annotation `0v` denotes zero
physical violations at guest-size-scaled loading; non-zero values report
violation counts. (e) Accessible cell volume fraction (probe 1.8 Å) for each
framework, coloured by material class (blue, zeolite; purple, MOF; orange,
COF). (f) Closed-cage misplacement by Packmol on cubic orthogonal hosts (FAU,
MOF-5) only; our method achieves zero misplacement (green diamonds). Packmol
places 1–3 guests in closed cages for EC/LiPF₆ in FAU. Packmol supports
orthogonal cells only; non-orthogonal MOF/COF hosts are shown for our method in
(a–c) and (d).

---

## Notes for manuscript

- All distance cutoffs: 2.2 Å (guest–framework and guest–guest).
- Steric blocked regions: periodic flood-fill with probe radius matched to guest
  size (default 1.8–2.5 Å depending on analysis).
- Packmol version 21.2.1; tolerance 2.2 Å; fixed framework; `inside box` constraint.
- MC engine: simulated annealing, Metropolis acceptance, conflict-driven moves.
