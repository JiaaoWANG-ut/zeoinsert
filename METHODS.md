# Methods

## Overview

We developed an **accessibility-aware Monte Carlo (MC) packing** workflow to construct physically valid host–guest structures for periodic porous frameworks (zeolites, MOFs, COFs). The method combines (i) a guest-size-dependent steric pore map obtained by periodic flood-fill on a voxel grid, and (ii) simulated-annealing MC insertion with overlap penalties and a blocked-pore penalty. Unlike geometry-only packers (e.g. Packmol), our approach prevents guest molecules from being placed in sterically closed cages that are topologically void but inaccessible to molecules of a given size.

The workflow comprises four steps: (1) load host framework and guest molecular geometries; (2) build a steric accessibility grid; (3) perform MC packing with accessibility constraints; (4) validate the final structure against unified physical-violation metrics.

---

## Periodic geometry and minimum-image convention

The host unit cell is defined by a $3 \times 3$ matrix $\mathbf{H}$ whose rows are the lattice vectors $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$. Cartesian coordinates $\mathbf{r}$ and fractional coordinates $\mathbf{f}$ are related by

$$
\mathbf{r} = \mathbf{f}\,\mathbf{H}, \qquad \mathbf{f} = \mathbf{r}\,\mathbf{H}^{-1}.
$$

All interatomic distances are evaluated under **minimum-image convention** (MIC). For displacement vectors $\{\Delta\mathbf{r}_k\}$, the MIC-adjusted displacement is

$$
\Delta\mathbf{r}_k^{\mathrm{MIC}} = \Delta\mathbf{r}_k - \mathrm{round}(\Delta\mathbf{f}_k)\,\mathbf{H},
\qquad \Delta\mathbf{f}_k = \Delta\mathbf{r}_k\,\mathbf{H}^{-1},
$$

and the minimum distance between point sets $A$ and $B$ is

$$
d_{\min}(A,B) = \min_{i,j} \left\|\Delta\mathbf{r}_{ij}^{\mathrm{MIC}}\right\|.
$$

Guest centres are wrapped into the unit cell via $\mathbf{f} \leftarrow \mathbf{f} - \lfloor \mathbf{f} \rfloor$.

---

## Steric pore accessibility analysis

### Voxelization

The unit cell is discretized into an $N_g \times N_g \times N_g$ grid ($N_g = 48$–$64$ in this work). Voxel centres in fractional coordinates are

$$
\mathbf{f}_{ijk} = \frac{1}{N_g}\left(i+\tfrac{1}{2},\, j+\tfrac{1}{2},\, k+\tfrac{1}{2}\right),
\qquad i,j,k \in \{0,\ldots,N_g-1\}.
$$

For each voxel centre $\mathbf{r}_{ijk} = \mathbf{f}_{ijk}\mathbf{H}$, the distance to the nearest framework atom is computed (via a $3\times3\times3$ replicated KD-tree for efficiency). A voxel is classified as **solid** if

$$
d_{\min}(\mathbf{r}_{ijk},\,\text{framework}) < r_{\mathrm{probe}},
$$

where $r_{\mathrm{probe}}$ is the probe radius, set to match the effective radius of the guest species (typically 1.2–2.5 Å). Remaining voxels constitute the **void** set $\mathcal{V}$.

### Topological vs steric accessibility

Two probe modes are distinguished:

| Mode | $r_{\mathrm{probe}}$ | Physical meaning |
|------|---------------------|----------------|
| Topological | 1.2 Å | Small probe; detects only truly enclosed voids in a perfect crystal |
| Steric | $\approx r_{\mathrm{guest}}$ | Guest-sized probe; marks cages/windows too narrow for the target molecule |

For zeolite FAU with $r_{\mathrm{probe}} = 2.5$ Å, eight sodalite $\beta$-cages (644 voxels each at $N_g=64$) are sterically blocked despite being topologically connected in the perfect structure.

### Periodic flood-fill

Void voxels are partitioned into connected components using 26-neighbour connectivity on a non-periodic label field, followed by **union–find merging** of labels that touch across opposite unit-cell faces (to recover periodicity).

A void component $C \subseteq \mathcal{V}$ is **accessible** if it touches at least one face of the unit cell:

$$
C \cap \partial\Omega \neq \varnothing,
$$

where $\partial\Omega$ denotes the union of the six cell-boundary voxel layers. The accessible void set is

$$
\mathcal{A} = \bigcup \left\{ C \subseteq \mathcal{V} \;\middle|\; C \cap \partial\Omega \neq \varnothing \right\}.
$$

**Blocked** void voxels are initially $\mathcal{B}_{\mathrm{raw}} = \mathcal{V} \setminus \mathcal{A}$. Isolated speckle clusters smaller than $n_{\min}$ voxels (default $n_{\min}=4$) are reassigned to $\mathcal{A}$ to suppress grid discretization artefacts. Manual blocked regions can be appended for user-defined confinement.

The accessible volume fraction is

$$
\phi_{\mathrm{acc}} = \frac{|\mathcal{A}|}{N_g^3},
\qquad
\phi_{\mathrm{blocked}} = \frac{|\mathcal{B}|}{|\mathcal{V}|}.
$$

### Probe-radius sweep

To quantify guest-size-dependent pore closure, $r_{\mathrm{probe}}$ is scanned over $[0.6,\,3.4]$ Å in steps of 0.2 Å. For each radius, $\phi_{\mathrm{acc}}(r_{\mathrm{probe}})$ and $\phi_{\mathrm{blocked}}(r_{\mathrm{probe}})$ are recorded, revealing step-like reductions in accessible volume when cage windows close to probes of increasing size.

---

## Guest molecular representation

Each guest species $s$ is provided as an atomic coordinate file (XYZ). Atoms are centred at the origin:

$$
\mathbf{p}_a^{(s)} \leftarrow \mathbf{p}_a^{(s)} - \frac{1}{N_s}\sum_{a=1}^{N_s}\mathbf{p}_a^{(s)}.
$$

Guest molecule $i$ of species $s_i$ is placed by a **centre** $\mathbf{f}_i$ (fractional) and a **rotation** $\mathbf{R}_i \in \mathrm{SO}(3)$:

$$
\mathbf{x}_{i,a} = \mathbf{R}_i\,\mathbf{p}_a^{(s_i)} + \mathbf{f}_i\,\mathbf{H},
\qquad a = 1,\ldots,N_{s_i}.
$$

Rotations are drawn uniformly on $\mathrm{SO}(3)$ via the quaternion method; small rotations for local MC moves are generated by the axis–angle (Rodrigues) formula.

For multi-species loading, species counts are assigned either explicitly or by stoichiometric ratio with largest-remainder apportionment.

---

## Overlap energy function

The total overlap energy is a sum of three pairwise penalty terms:

$$
E = E_{\mathrm{fw}} + E_{\mathrm{gg}} + E_{\mathrm{block}}.
$$

### Guest–framework overlap

For guest $i$ with framework cutoff $d_{\mathrm{fw}}^{(i)}$ (default 2.2 Å):

$$
E_{\mathrm{fw}} = \sum_{i=1}^{N_{\mathrm{mol}}} \max\!\left(0,\; d_{\mathrm{fw}}^{(i)} - d_{\min}^{(i,\mathrm{fw})} \right)^2,
$$

where $d_{\min}^{(i,\mathrm{fw})} = d_{\min}(\text{atoms of guest } i,\,\text{framework})$ under MIC.

### Guest–guest overlap

For each unordered pair $(i,j)$, $i<j$, with cutoff $d_{\mathrm{mm}}^{(ij)} = \max(d_{\mathrm{mm}}^{(i)}, d_{\mathrm{mm}}^{(j)})$:

$$
E_{\mathrm{gg}} = \sum_{i<j} \max\!\left(0,\; d_{\mathrm{mm}}^{(ij)} - d_{\min}^{(ij)} \right)^2.
$$

### Blocked-pore penalty

When accessibility-aware packing is enabled, a flat penalty is applied if the centre of guest $i$ falls in a blocked voxel:

$$
E_{\mathrm{block}} = \sum_{i=1}^{N_{\mathrm{mol}}} \lambda_{\mathrm{block}}\,\mathbb{1}\!\left[\mathbf{f}_i \in \mathcal{B}\right],
\qquad \lambda_{\mathrm{block}} = 100\;\text{a.u.}
$$

Voxel membership is determined by $\lfloor \mathbf{f}_i \cdot N_g \rfloor \mod N_g$. When accessibility blocking is disabled (ablation baseline), $E_{\mathrm{block}} \equiv 0$ but misplacement is still measured post hoc against the steric grid.

---

## Simulated-annealing Monte Carlo packing

### Initialization

Each guest centre $\mathbf{f}_i$ is drawn uniformly from accessible voxel centres $\{\mathbf{f} : \mathbf{f} \in \mathcal{A}\}$ (or uniformly in the unit cell if no grid is used). Rotations $\mathbf{R}_i$ are random. The initial energy $E$ is computed from the overlap function above.

### Move types

At each MC step, the guest with the largest per-molecule conflict score

$$
c_i = e_{\mathrm{fw}}^{(i)} + e_{\mathrm{block}}^{(i)} + \sum_{j \neq i} e_{\mathrm{gg}}^{(ij)}
$$

is selected for trial displacement. Three move types are applied with probabilities $(p_{\mathrm{small}}, p_{\mathrm{jump}}, p_{\mathrm{rot}}) = (0.70,\,0.10,\,0.20)$:

1. **Small move** ($p_{\mathrm{small}}$): translate centre by $\Delta\mathbf{r} = \delta_t \,\boldsymbol{\xi}$ ($\boldsymbol{\xi}$ unit random vector, $\delta_t = 0.8$ Å) and apply a small random rotation (maximum angle $\pm 20°$).
2. **Big jump** ($p_{\mathrm{jump}}$): re-sample $\mathbf{f}_i$ from $\mathcal{A}$ and draw a new uniform $\mathbf{R}_i$.
3. **Large rotation** ($p_{\mathrm{rot}}$): keep centre fixed, apply random rotation up to $\pm 60°$.

### Metropolis acceptance with linear cooling

The temperature schedule is linear:

$$
T(n) = \max\!\left(T_{\min},\; T_0\left(1 - \frac{n}{N_{\mathrm{MC}}}\right)\right),
\qquad T_0 = 1.0,\; T_{\min} = 0.02,
$$

where $n$ is the MC step and $N_{\mathrm{MC}}$ is the maximum iteration count (typically $8\times10^3$–$3\times10^4$).

A trial move changing the energy by $\Delta E$ is accepted with Metropolis probability

$$
P_{\mathrm{acc}} = \min\!\left(1,\; \exp\!\left(-\frac{\Delta E}{T(n)}\right)\right).
$$

Moves with $\Delta E \le 0$ are always accepted. Energy differences are evaluated incrementally (only terms involving the moved molecule are recomputed).

### Convergence criterion

The simulation terminates when $E < \varepsilon$ with $\varepsilon = 10^{-6}$ a.u., or when $N_{\mathrm{MC}}$ steps are reached.

---

## Physical-violation metrics (validation)

All packing methods (ours, Packmol, accessibility-blind ablation) are scored by the same post hoc judge:

| Metric | Criterion | Definition |
|--------|-----------|------------|
| Inaccessible placement | Centre in $\mathcal{B}$ | $\mathbb{1}[\mathbf{f}_i \in \mathcal{B}]$ |
| Framework clash | $d_{\min}^{(i,\mathrm{fw})} < d_{\mathrm{fw}}$ | Per-molecule minimum guest–framework distance |
| Guest overlap | $d_{\min}^{(ij)} < d_{\mathrm{mm}}$ | Per-molecule pair count |

The total violation count for $N_{\mathrm{mol}}$ guests is

$$
N_{\mathrm{viol}} = N_{\mathrm{inacc}} + N_{\mathrm{gf}} + N_{\mathrm{gg}},
$$

where $N_{\mathrm{inacc}} = \sum_i \mathbb{1}[\mathbf{f}_i \in \mathcal{B}]$, $N_{\mathrm{gf}}$ counts molecules with framework clash, and $N_{\mathrm{gg}}$ counts molecules involved in at least one guest–guest overlap.

---

## Baseline: Packmol

As an external geometry-only baseline, we used Packmol v21.2.1. The framework atoms are supplied as a fixed `structure` with zero degrees of freedom; guest molecules are inserted with `inside box` constraints and a global `tolerance` of 2.2 Å (matching our overlap cutoffs). Packmol enforces minimum interatomic distances but has **no knowledge of pore accessibility** and therefore can place guests in sterically closed cages. Packmol supports orthogonal unit cells only; non-orthogonal MOF/COF hosts are evaluated with our method alone.

---

## Target loading for cross-framework generalization

For each framework–guest pair, the target number of molecules per unit cell is scaled by accessible volume and inversely by guest size (relative to CO₂, $r_{\mathrm{CO}_{2}} = 1.65$ Å):

$$
N_{\mathrm{target}} = \min\!\left(N_{\mathrm{load,max}},\; \max\!\left(4,\; \left\lfloor \frac{|\mathcal{A}|}{700}\left(\frac{r_{\mathrm{CO}_{2}}}{r_{\mathrm{guest}}}\right)^2 \right\rfloor \right)\right),
\qquad N_{\mathrm{load,max}} = 24.
$$

The steric grid for each pair uses $r_{\mathrm{probe}} = \max(r_{\mathrm{guest}},\,1.8\;\text{Å})$.

---

## Computational details

| Parameter | Value |
|-----------|-------|
| Grid resolution $N_g$ | 48 (generalization) / 64 (FAU analysis) |
| Connectivity | 26-neighbour |
| Min. blocked cluster size | 4 voxels |
| Framework cutoff $d_{\mathrm{fw}}$ | 2.2 Å |
| Guest–guest cutoff $d_{\mathrm{mm}}$ | 2.2 Å |
| Blocked-pore penalty $\lambda_{\mathrm{block}}$ | 100 a.u. |
| Translation step $\delta_t$ | 0.8 Å |
| Small rotation limit | $\pm 20°$ |
| Large rotation limit | $\pm 60°$ |
| $T_0 / T_{\min}$ | 1.0 / 0.02 |
| MC steps $N_{\mathrm{MC}}$ | $8\times10^3$–$3\times10^4$ |

Structures are read/written in CIF/XYZ format via OVITO 3.15. Pore grids are cached as compressed NumPy archives (`.npz`). High-quality structure figures are rendered with OVITO TachyonRenderer (ambient occlusion, CPK colouring).

All calculations were performed on a single CPU (Python 3.12, NumPy 2.4, SciPy 1.18). Typical wall times: steric grid build $\sim$1.5 s per framework ($N_g=64$); MC packing of 32 CO₂ in FAU $\sim$1–3 s; Packmol baseline $\sim$1 s.

---

## Software availability

The implementation is organized as follows:

| Module | Function |
|--------|----------|
| `pore_accessibility.py` | Voxel grid, periodic flood-fill, probe sweep |
| `mc_engine.py` | Simulated-annealing MC `pack()` |
| `error_metrics.py` | Unified violation scoring |
| `baselines/run_packmol.py` | Packmol wrapper |
| `render_ovito.py` | Tachyon structure rendering |

Source code is available at https://github.com/JiaaoWANG-ut/zeoinsert.

---

## Statistical analysis

Ablation data (accessibility-aware vs blind packing) report mean $\pm$ s.d. over five random seeds ($s = 0,1,2,3,4$) at each loading level. Convergence traces show individual seed trajectories and the seed-averaged overlap energy. Probe-radius sweeps are deterministic (no stochastic component).
