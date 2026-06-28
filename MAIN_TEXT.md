# Accessibility-aware Monte Carlo packing of guest molecules in nanoporous frameworks

## Introduction

Constructing physically valid host–guest configurations is a prerequisite for atomistic studies of adsorption, ion transport and separation in nanoporous materials such as zeolites, metal–organic frameworks (MOFs) and covalent organic frameworks (COFs). Conventional molecular packers place guests by enforcing only minimum interatomic distances, treating the host as a static obstacle field. This geometry-only view ignores a key physical constraint: a cage may be topologically void in the perfect crystal yet sterically closed to a molecule of finite size, because its connecting windows are narrower than the guest. Guests artificially seeded inside such closed cages cannot enter or leave during a real experiment, biasing loading, dynamics and free-energy estimates.

Here we introduce an accessibility-aware Monte Carlo (MC) packing workflow that couples a guest-size-dependent steric pore map directly into the insertion energy. We show that the method eliminates closed-cage placements and framework clashes across zeolites, MOFs and COFs, whereas a state-of-the-art geometry-only packer (Packmol) systematically traps guests in inaccessible pores.

## Results

### An accessibility-aware packing workflow

Our workflow proceeds in four stages (**Fig. 1a**): (1) the host framework and a guest molecule library are loaded; (2) a periodic flood-fill on a voxel grid partitions the unit cell into solid framework, accessible void and sterically blocked void; (3) simulated-annealing MC inserts guests under a combined overlap and blocked-pore penalty; and (4) the validated host–guest model confines all guests to accessible channels. The central physical object is the steric accessibility map. For FAU zeolite probed at a guest-sized radius of 2.5 Å, eight sodalite β-cages are classified as inaccessible (**Fig. 1b**, red), even though they remain topologically connected in the ideal crystal. These regions are precisely the volumes that geometry-only packers cannot detect.

A two-dimensional voxel cross-section makes the partition explicit (**Fig. 1c**): grey voxels are solid framework, blue voxels are accessible void, and red voxels are blocked void. Integrating over the full grid (64³), the FAU unit cell decomposes into 67% solid, 31% accessible and 2% sterically blocked volume (**Fig. 1d**). The blocked fraction is small but, as shown below, decisive for the physical validity of the packed model.

### Probe physics and the cost of accessibility blindness

Pore accessibility is intrinsically guest-size dependent. Scanning the probe radius reveals that the accessible cell volume decreases in discrete steps as the probe grows (**Fig. 2a**), each step corresponding to a window closing to molecules of a given size; vertical markers indicate the effective radii of Li⁺, H₂O, CO₂, N₂ and ethylene carbonate (EC). The complementary quantity—the fraction of void volume locked inside closed cages—rises sharply for cage-type frameworks such as ERI and OFF above ~3 Å (**Fig. 2b**), confirming that accessibility cannot be inferred from a single static pore network.

To quantify the consequence for packing, we measured the fraction of CO₂ guests placed in closed cages as a function of loading in FAU (**Fig. 2c**). Accessibility-aware packing remains at 0% misplacement across loadings of 8–48 molecules, whereas accessibility-blind packing traps 4–12.5% of guests (mean ± s.d., five seeds). The structural origin is visualized in **Fig. 2d**, where 32 CO₂ molecules occupy only accessible supercages and channels, versus **Fig. 2e**, where accessibility-blind insertion strands several guests (red) inside sodalite cages. The MC engine converges reliably: for 70 CO₂ molecules in FAU, the overlap penalty energy decays to zero under the linear annealing schedule across all five seeds (**Fig. 2f**).

### Application cases outperform a geometry-only baseline

We next benchmarked three representative loading scenarios in FAU against Packmol (**Fig. 3a,b**): a Li-salt electrolyte (EC/DMC/LiPF₆, 20 molecules), CO₂ adsorption (32 CO₂), and competitive CO₂/H₂O loading (16+16). In every case our method places all guests in accessible space (**Fig. 3a**), whereas Packmol seeds guests that violate steric accessibility or minimum-distance cutoffs (**Fig. 3b**, red).

Counting guests whose centres fall in closed cages, our method yields zero across all three cases, while Packmol traps 3 (electrolyte), 4 (CO₂) and 5 (separation) guests (**Fig. 3c**). Decomposing the Packmol failures by type—closed-cage placement, framework clash and guest–guest overlap—shows that all three violation modes occur, whereas our method registers zero total violations in every case (**Fig. 3d**, green diamonds). The distribution of minimum guest–framework distances summarizes the difference (**Fig. 3e**): Packmol exhibits a tail below the 2.2 Å cutoff, indicating physical clashes, while our distribution lies entirely above it.

### Generalization across zeolites, MOFs and COFs

Finally, we tested transferability beyond zeolites. Accessibility-aware packing produces clash-free CO₂ configurations in MOF-5, UiO-66 and COF-5 (**Fig. 4a–c**). Across a matrix of six frameworks (FAU, LTL, ERI, MOF-5, UiO-66, COF-5) and four guests (CO₂, H₂O, EC, LiPF₆), the method achieves guest-size-scaled loadings with zero physical violations in nearly all cells (**Fig. 4d**, annotated `0v`). The underlying accessible volume fraction varies systematically with material class (**Fig. 4e**): zeolites, MOFs and COFs span a wide range of porosity, yet the packing remains valid throughout.

On the orthogonal hosts where Packmol is applicable (FAU, MOF-5), the geometry-only baseline again places 1–3 guests in closed cages for the larger EC/LiPF₆ species, whereas our method achieves zero misplacement (**Fig. 4f**, green diamonds). Because Packmol supports orthogonal cells only, the non-orthogonal MOF and COF hosts are reported for our method alone (**Fig. 4a–c,d**).

## Discussion

Across all tested systems, coupling a guest-size-dependent steric accessibility map into the MC insertion energy is sufficient to eliminate the closed-cage placements and framework clashes that geometry-only packers produce. The accessible-volume analysis (**Fig. 2a,b**; **Fig. 4e**) shows why a single distance cutoff is inadequate: accessibility is a function of guest size, not of the host alone, and it changes discontinuously as windows close. By making this constraint explicit, the workflow yields host–guest models that are physically realizable starting points for molecular dynamics and free-energy calculations, and it transfers without modification from zeolites to MOFs and COFs (**Fig. 4**). We anticipate that accessibility-aware packing will be particularly valuable for high-throughput screening, where manual inspection of every packed structure is infeasible and silently trapped guests would otherwise corrupt downstream property predictions.
