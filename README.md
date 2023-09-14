

# Description
This code is used to address symmetry breaking bifurcations using modified squared-operator iteration method.

# Abstract:
We address symmetry breaking bifurcations (SBBs) in the ground-state (GS) and dipole-mode (DM)
solitons of the system of one-dimensional linearly coupled nonlinear Schrodinger (NLS) equations, modeling the
propagation of light in a dual-core planar waveguide with the Kerr nonlinearity and two different types of PT -
symmetric potential (localized and delocalized). The PT -symmetric potential is employed to obtained different
types of solutions. A supercritical pitchfork bifurcation occurs in families of symmetric solutions of both the GS
and DM types, i.e., spatially even and odd ones. A novel feature of the system is interplay between breakings of
the PT and inter-core symmetries. Stability of symmetric GS and DM modes and their asymmetric counterparts,
produced by SBBs of both types, is explored by means of the linear-stability analysis and direct simulations.
It is found that the instability of PT -symmetric solutions takes place prior to the inter-core symmetry breaking.
Surprisingly, stable inter-core-symmetric GS solutions may remain stable while the PT symmetry is broken. Fully
asymmetric GS and DM solitons are only partially stable. In addition to the systems with the PT -symmetric
potentials, we construct symmetric and asymmetric GS solitons under the action of a pure imaginary localized
potential, for which the SBB is subcritical. These results exhibit that stable solitons can still be found in dissipative
systems. Finally, excitations of symmetric and asymmetric GS solitons are investigated by making the potential’s
parameters or the system’s coupling constant functions of the propagation distance, showing that GS solitons can
be converted from an asymmetric shape onto a symmetric one under certain conditions. These results may pave
the way for the study of linear and nonlinear phenomena in a dual-core planar waveguide with PT -symmetric
potential and the related physical experimental designs.


# Dual-core_code

Program pt_pu.m is used to calculate PT symmetry breaking. And Program eigg.m is to calculate eigenvalues.

Program MSOM_asym.m is to compute the numerical solutions including both symmetric and asymmetric ones.

Programs mu_Psym.m and mu_p_asym.m are to obtain the relationship between power P and propagation constant.

Program ssb_W0_j.m is to obtain the boundary of the breaking of the inter-core symmetry for symmetric and asymmetric modes.

Program ssb_theta.m is to obtain the asymmetry characteristic θ.

Programs evolution_j1.m and evolution_v0.m are used to obtain the numerically simulated evolution.
