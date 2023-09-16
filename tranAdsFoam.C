/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    tranAdsLumpFoam

Group
    grpBasicSolvers

Description

    Chemical transport and adsorption solver based on LFD and 
    langmuir isotherm.

    Frist-order time splitting with stiff ODE solver, 

    Personally used for fixed bed simulation with fast reaction.

    Author: Haochen Li

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fvOptions.H"
#include "simpleControl.H"

#include "adsorptionODE.H" // define the adsorption ODE

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Passive scalar transport equation solver."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating scalar transport\n" << endl;

    #include "CourantNo.H"

    // constant parameter to simplify the expression
    const dimensionedScalar alpha = ((1-porosity)/porosity)*rhop;

    // Create the selected ODE system solver
    // ODE system parameter
    const scalar alpha_1 = alpha.value();

    const scalar kf_1 = kf.value(); 
    const scalar kr_1 = kr.value();

    // langmuir isotherm parameters, pay attention to the unit!
    const scalar qMax_1 = Langmuir_qMax.value();

    // Initial condition
    //const scalar y0 = 2e-20;

    // Create the ODE system
    myODE ode(alpha_1, kf_1, qMax_1, kr_1);

    // construct the solver
    word ODESolverName ("rodas23");
    dictionary dict;
    dict.add("solver", ODESolverName);
    autoPtr<ODESolver> odeSolver = ODESolver::New(ode, dict);

    // Required to store initial condition 
    scalarField yStart(ode.nEqns());
    // Required to store dydx
    scalarField dyStart(ode.nEqns());

    // some solving controls
    
    //odeSolver->maxSteps() = 20000; //10000 is the default
    //odeSolver->absTol() = 1e-12; //
    odeSolver->relTol() = 1e-4; //1e-4; is the default



    // First-order splitting scheme (see the JEE paper SI)
    // transport equation time loop
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        // extract the time
        const dimensionedScalar deltaT(mesh.time().deltaT()); 

        // initialize the varibles that pass the data to ODE
        volScalarField q_ODE = q; //current q
        volScalarField Y_ODE = Y; //current Y

        volScalarField dYdt = Y/deltaT; //current Y
        
        scalar ODEdeltaT = 0; // time step size for the ODE solver
        //volScalarField Y0 = Y; //current Y


        // Step 1 solve the ODE over deltaT
        Info << "Step 1 ODE: Solving for reaction" << endl;
        ODEdeltaT = deltaT.value();
        #include "solveAdsorption.H"

        forAll(dYdt,celli)
        {
            dYdt[celli] = (max(Y_ODE[celli],0.0)-max(Y[celli],0.0))/ODEdeltaT;
        }


        // Step 2 solve the transport equation over deltaT
        Info << "Step 2 PDE: Solving for transport" << endl;
        while (simple.correctNonOrthogonal())
        {
            fvScalarMatrix YEqn
            (
                fvm::ddt(Y)
              + fvm::div(phi/porosity, Y)
              - fvm::laplacian(Dh, Y)
              ==
                dYdt
            );

            YEqn.relax();
            fvOptions.constrain(YEqn);
            YEqn.solve();
            fvOptions.correct(Y);

            Y.max(0.0); //clip off the negtive values
        }  

        runTime.write();
        
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
