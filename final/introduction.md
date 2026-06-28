# Introduction

Crystalline nanoporous materials—zeolites, metal–organic frameworks (MOFs) and covalent organic frameworks (COFs)—provide high-surface-area, chemically tunable scaffolds that underpin gas storage, molecular separation, ion transport and heterogeneous catalysis<sup>1–4</sup>. The principle of reticular chemistry, in which molecular building blocks are assembled into predictable extended networks, has expanded this design space to a practically unlimited number of candidate frameworks<sup>1</sup>. This explosion of structures has, in turn, motivated large curated repositories such as the IZA Database of Zeolite Structures<sup>5</sup>, the Computation-Ready, Experimental MOF (CoRE MOF) database<sup>6</sup>, hypothetical-MOF libraries<sup>7</sup>, and the CURATED COFs collection<sup>8</sup>, which together place hundreds of thousands of frameworks within reach of computational study.

Molecular simulation has become the workhorse for navigating this space. High-throughput computational screening, typically built on grand canonical Monte Carlo (GCMC) sampling of adsorption, now routinely ranks thousands to millions of frameworks for a target application and guides experimental synthesis<sup>3,4,9,10</sup>. Every such calculation begins from an atomistic host–guest configuration: molecules must be placed inside the pore network in a way that respects both thermodynamic and kinetic accessibility before any dynamics, free-energy or uptake estimate can be trusted<sup>11,12</sup>. The fidelity of this initial placement therefore propagates directly into the reliability of downstream predictions.

Two broad strategies are used to generate these starting structures. General-purpose packing tools such as Packmol build initial configurations by solving a constrained optimization that enforces a minimum interatomic separation between all atoms, treating the framework as a static obstacle field<sup>13,14</sup>. Packmol further restricts guest placement to orthogonal bounding boxes (`inside box` constraints) and therefore cannot be applied directly to unit cells with non-orthogonal lattice vectors—a limitation that excludes many MOFs and COFs from automated packing workflows. Random-insertion moves in GCMC similarly propose trial positions that are accepted or rejected on the basis of interaction energy and Boltzmann statistics<sup>12</sup>. Both views are adequate for open, continuously connected pore spaces, but they share a blind spot: a cavity may be topologically void in the perfect crystal yet sterically closed to a molecule of finite size, because the windows connecting it to the percolating network are narrower than the guest.

The geometric characterization of pore space is itself a mature field. Tools such as Zeo++, which performs a Voronoi decomposition of the void space, and the grid-based PoreBlazer compute pore-size distributions, accessible surface areas and accessible volumes for a spherical probe of a chosen radius, and can identify guest-inaccessible regions<sup>15,16</sup>. In the specific context of adsorption simulations, the standard remedy for inaccessible pockets is "pore blocking": a flood-fill or energy-grid segmentation isolates disconnected cavities, which are then excluded from sampling either by inserting repulsive blocking spheres or by setting their grid energies to prohibitively large values<sup>17–19</sup>. Neglecting this step is known to cause systematic overestimation of uptake, and its importance has been documented repeatedly across zeolites and MOFs<sup>17,18</sup>. Crucially, however, these accessibility analyses have been developed to *count* adsorption correctly within GCMC, and they remain decoupled from the general-purpose packers used to *build* host–guest structures for molecular dynamics: the packer enforces distances, while accessibility is treated, if at all, as a separate post hoc correction.

This gap has tangible consequences. The molecular-sieving effect—long established in zeolite chemistry—dictates precisely which species can enter the internal pore space of a framework<sup>20</sup>. When a packer or insertion algorithm seeds a guest inside a cage whose windows are too small for it to traverse, the resulting configuration is kinetically trapped: the molecule can neither enter nor leave during a physical experiment, and any trajectory launched from such a state samples a non-equilibrium ensemble, biasing predicted loadings, distorting diffusion pathways and skewing free-energy estimates. A geometry-only packer such as Packmol will therefore generate configurations that appear clash-free yet are physically unrealizable<sup>13,14</sup>.

Here we introduce an accessibility-aware Monte Carlo packing workflow that couples a guest-size-dependent steric pore map directly into the insertion energy. By analysing the three-dimensional connectivity of the pore network with a probe radius matched to the guest, the algorithm distinguishes open, accessible cages from sterically closed voids *before* any molecule is placed, and penalizes blocked-cage occupancy during simulated annealing rather than correcting for it afterwards. The workflow operates on general periodic unit cells and does not require orthogonal box constraints. We show that the method eliminates closed-cage placements and framework clashes across a diverse set of nanoporous materials—zeolites, MOFs and COFs—whereas a state-of-the-art geometry-only packer (Packmol) systematically traps guests in inaccessible pores on the subset of hosts where it is applicable.

---

## References

1. Yaghi, O. M.; O’Keeffe, M.; Ockwig, N. W.; Chae, H. K.; Eddaoudi, M.; Kim, J. Reticular synthesis and the design of new materials. *Nature* **2003**, *423*, 705–714. DOI: 10.1038/nature01650.

2. Furukawa, H.; Cordova, K. E.; O’Keeffe, M.; Yaghi, O. M. The chemistry and applications of metal–organic frameworks. *Science* **2013**, *341*, 1230444. DOI: 10.1126/science.1230444.

3. Smit, B.; Maesen, T. L. M. Molecular simulations of zeolites: Adsorption, diffusion, and shape selectivity. *Chem. Rev.* **2008**, *108*, 4125–4184. DOI: 10.1021/cr8002642.

4. Düren, T.; Bae, Y.-S.; Snurr, R. Q. Using molecular simulation to characterise metal–organic frameworks for adsorption applications. *Chem. Soc. Rev.* **2009**, *38*, 1237–1247. DOI: 10.1039/b803498m.

5. Baerlocher, C.; McCusker, L. B. Database of Zeolite Structures. International Zeolite Association, http://www.iza-structure.org/databases/.

6. Chung, Y. G.; et al. Advances, updates, and analytics for the computation-ready, experimental metal–organic framework database: CoRE MOF 2019. *J. Chem. Eng. Data* **2019**, *64*, 5985–5998. DOI: 10.1021/acs.jced.9b00835.

7. Wilmer, C. E.; Leaf, M.; Lee, C. Y.; Farha, O. K.; Hauser, B. G.; Hupp, J. T.; Snurr, R. Q. Large-scale screening of hypothetical metal–organic frameworks. *Nat. Chem.* **2012**, *4*, 83–89. DOI: 10.1038/nchem.1192.

8. Ongari, D.; Yakutovich, A. V.; Talirz, L.; Smit, B. Building a consistent and reproducible database for adsorption evaluation in covalent–organic frameworks. *ACS Cent. Sci.* **2019**, *5*, 1663–1675. DOI: 10.1021/acscentsci.9b00619.

9. Colón, Y. J.; Snurr, R. Q. High-throughput computational screening of metal–organic frameworks. *Chem. Soc. Rev.* **2014**, *43*, 5735–5749. DOI: 10.1039/C4CS00070F.

10. Boyd, P. G.; Lee, Y.; Smit, B. Computational development of the nanoporous materials genome. *Nat. Rev. Mater.* **2017**, *2*, 17037. DOI: 10.1038/natrevmats.2017.37.

11. Ernst, M.; Evans, J. D.; Gryn’ova, G. Host–guest interactions in framework materials: Insight from modeling. *Chem. Phys. Rev.* **2023**, *4*, 041303. DOI: 10.1063/5.0144827.

12. Frenkel, D.; Smit, B. *Understanding Molecular Simulation: From Algorithms to Applications*, 2nd ed.; Academic Press: San Diego, 2002.

13. Martínez, J. M.; Martínez, L. Packing optimization for automated generation of complex system’s initial configurations for molecular dynamics and docking. *J. Comput. Chem.* **2003**, *24*, 819–825. DOI: 10.1002/jcc.10216.

14. Martínez, L.; Andrade, R.; Birgin, E. G.; Martínez, J. M. Packmol: A package for building initial configurations for molecular dynamics simulations. *J. Comput. Chem.* **2009**, *30*, 2157–2164. DOI: 10.1002/jcc.21224.

15. Willems, T. F.; Rycroft, C. H.; Kazi, M.; Meza, J. C.; Haranczyk, M. Algorithms and tools for high-throughput geometry-based analysis of crystalline porous materials. *Microporous Mesoporous Mater.* **2012**, *149*, 134–141. DOI: 10.1016/j.micromeso.2011.08.020.

16. Sarkisov, L.; Harrison, A. Computational structure characterisation tools in application to ordered and disordered porous materials. *Mol. Simul.* **2011**, *37*, 1248–1257. DOI: 10.1080/08927022.2011.592832.

17. Gómez-Álvarez, P.; Ruiz-Salvador, A. R.; Hamad, S.; Calero, S. Importance of blocking inaccessible voids on modeling zeolite adsorption: Revisited. *J. Phys. Chem. C* **2017**, *121*, 4462–4470. DOI: 10.1021/acs.jpcc.7b00031.

18. Zhang, K.; Nalaparaju, A.; Chen, Y.; Jiang, J. Crucial role of blocking inaccessible cages in the simulation of gas adsorption in a paddle-wheel metal–organic framework. *RSC Adv.* **2013**, *3*, 16152–16158. DOI: 10.1039/c3ra42213e.

19. Kim, J.; Martin, R. L.; Rübel, O.; Haranczyk, M.; Smit, B. High-throughput characterization of porous materials using graphics processing units. *J. Chem. Theory Comput.* **2012**, *8*, 1684–1693. DOI: 10.1021/ct200787v.

20. Breck, D. W. *Zeolite Molecular Sieves: Structure, Chemistry and Use*; Wiley: New York, 1974.
