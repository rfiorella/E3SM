# EAMXX water isotope ratio implementation order

This provides a list of steps to implement water isotope ratios and generic tracers 
in EAMXX - each of these steps should become a spec, run in a campaign, with each spec 
having its own PR as described in the skill.

### Design principles:
- Water tracer/isotope ratio code should be kept modular and in its own namespace so that it can be applied cleanly
if parameterizations change.
- Water tracers should be implemented as a way of tracking a feature or a subset of the hydrological cycle. Isotope ratios are a special case of water tracer where proportionality across the phase changed is not exactly maintained (e.g., alpha != 1)
- As much of the following implementation as is practical should occur in a new scream::water_tracers namespace.
- Existing tests must never be changed. Tests building on unit tests already available for specific parameterizations should test I/O with the expanded water variable arrays.


### References 
I have created a references folder in this repo. It contains:
- an implementation of water tracers and isotope ratios in CAM that can be used for inspiration
- some implementation specs that were from prior attempts to implement water tracers.
- The publications and proposals listed below in the list of expected PRs.

### List of expected PRs:
1. Implement the concept of water species, including bulk water (H2O),
H216O, HDO, H218O, H217O, and HTO. Allow additional tracers to be specified in CMakeLists for water_tracers namespace using add_tracer.
2. Extend the water vapor arrays (qv) to have an extra dimension that
allows setting for multiple water tracers. The first index of this dimension is defined
as the currently implemented qv, and the values in this slice of the array must remain
unchanged throughout development. Add tests examining ability of model to have this dimension for qv > 1 (e.g., H2O and HDO, and H2O, H216O, HDO, and H218O).
3. Extend the cloud liquid and ice arrays to have an extra tracer dimension, 
similar to what was just done for qv. Add tests examining ability of model to have this dimension for these arrays > 1 (e.g., H2O and HDO, and H2O, H216O, HDO, and H218O).
4. Extend the liquid and solid precipitation arrays to have an extra tracer
dimension, following what has been previously done for qv, cloud liquid, and cloud ice. Add tests examining ability of model to have this dimension for these arrays > 1 (e.g., H2O and HDO, and H2O, H216O, HDO, and H218O).
5. Generate a list of phase changes throughout all of the parameterizations in the model.
Group them first by type (evaporation, condensation, mixed-phase, those occurring under
supersaturation conditions, deposition, sublimation) and then by parameterization (SHOC, P3, etc.)
6. Implement equilibrium fractionation functions for water phase changes. Use Horita and Wesolowski 1994
for liquid->vapor equilibrium (with Barkan and Luz for H217O relationship to H218O) and provide a few options for ice-vapor equilibrium: a) Merlivat and Nief 1967,  Lamb et al. 2017, and Ellehoj et al 2013 for HDO and 
b) Majoube 1971 and Ellehoj et al. 2013 for H218O. These functions should return a fractionation factor. These options should be cmake selectable. 
Test each at a range of Earth-relevent temperatures and ensure they match values that would be 
obtained from equations in these publications.
7. Implement functions to return alpha_diff (the fractionation factor due to molecular diffusivity) for HDO, H218O and H217O relative to H216O. Use Merlivat 1978 for HDO and Hellman and Harvey 2020 for H218O and H217O. Tests should resolve as follows: a) if water type is H2O or H216O = 1; b) if water type is H218O = 1.0283, c) if water type is H217O = 1.0146; d) if water type is hdo = 0.9755.
8. Implement functions to determine the net fractionation factor for non-equilibrium processes, blending the equilibrium and non-equilibrium processes as required. For ocean evaporation, use the Brutsaert evaporation model coupled to the Craig-Gordon model. For hydrometeor partial evaporation, use the Stewart 1975 model. For mixed-phase preciptation, use Ciais and Jouzel 1994. Tests should be designed based on original publications. In all cases, tests should require that fractionation factors for H2O and H216O water types equal 1.
9. Develop a utility function to take the ratio of two tracers to maintain proportionality or implement phase changes - same as wtrc_ratio in CAM.
10. Implement hooks to calculate tracer/isotopologue partitioning during ocean evaporation. Write tests occurring for an ocean surface with 
water isotope ratio equivalent to SMOW, a surface temperature of 28°C, and RH of 80%.
11. Implement hooks to calculate tracer/isotopologue partitioning during any phase changes in SHOC. 
12. Implement hooks to calculate tracer/isotopologue partitioning during any phase changes in P3.
13. Implement hooks to calculate tracer/isotopologue partitioning during any phase changes in Zhang-McFarlane deep convection scheme.
14. Implement hooks to calculate tracer/isotopologue partitioning during any phase changes in other parameterizations not captured here.
15. Implement methane oxidation by OH radical and potential impacts on HDO concentrations. Allow to be turned on or off by cmake.
16. For the next few PRs we will be considering sub-types of tracers within a water type. Use Fiorella et al. 2021 and the proposal as a guide for these specs 16-20. To start, develop a method for allowing 
tracking of a particular predefined region of the evaporation field with the region of interest being specified in CMake if possible. The region of interest should be able to be specified by a lat/lon box or a shapefile.
17. Implement a second method of tagging evaporation fields by using spherical harmonic basis functions to decompose the evaporation field or converting to cartesian coordinates and back to lat/lon after the model runs. Additionally, implement a method using needlets for regional refinement in particular receptor regions. Again these should be addable tracers in CMake.
18. Implement a method of specifying a function or scalar multiplicative factor for these evaporation field tracers - for example, allowing the tracer to be multiplied by the current time or surface temperature to track the mean time or temperature of evaporation.
19. Implement a method of tracking condensation properties - this should allow preserving a quantity associated with vapor to liquid or ice fluxes, such as temperature or pressure of the layer. The tracer should be addable with CMake, with the quantity to use defined in an add_tracer CMake Call
20. Implement a method of tracking parcel-integrated quantities, e.g., Q*dt, so that the residence time of water vapor in the atmosphere can be estimated. As before this should be definable in CMake.
21. Develop a comprehensive guide and test cases using the doubly periodic EAMxx configuration that can be used for testing.




