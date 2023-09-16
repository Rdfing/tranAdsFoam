# tranAdsFoam [![Github All Releases](https://img.shields.io/github/downloads/Rdfing/interAdsFoam/total.svg)]()

Macro-scale transport adsorption solver. The solver is developed to simulate chemical adsorption or particle filtration through a fixed bed reactor.

## Description
This solver is developed to simulate turbulence mixing and chemical adsorption in urban water treatment systems. The solver integrates VOF with chemical adsorption.
- Non-equilibrium (Langmuir-like) adsorption model 
- First-order time splitting for transport and reaction (for now)
- Stiff ODE solver for reaction
- Examples of inverse modeling of adsorption and headloss are provided in the reference

Here is an example 
![Alt Text](https://github.com/Rdfing/tranAdsFoam/blob/main/Case1_column_test_SLR_40/solution_sample.gif)

## Prerequisites
OpenFOAM-v2212

## Authors
Haochen Li, PhD <br />
Department of Civil and Environmental Engineering <br />
University of Tennessee, Knoxville

## Reference
Li, H., Sansalone, J., 2022. InterAdsFoam: An Open-Source CFD Model for Granular Media–Adsorption Systems with Dynamic Reaction Zones Subject to Uncontrolled Urban Water Fluxes. J. Environ. Eng. 148, 04022049. https://doi.org/10.1061/(ASCE)EE.1943-7870.0002027

## Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and distributor of the OpenFOAM software and owner of the OPENFOAM®  and OpenCFD®  trademarks.
