Info<< "Reading field Y\n" << endl;


volScalarField Y
(
    IOobject
    (
        "Y",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

// q is the adsorbed phases concentration [kg adsorbed/kg adsorbent]
// so is dimensionless
volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


Info<< "Reading field U\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);


// volScalarField qStar // equilibirum q
// (
//     IOobject
//     (
//         "qStar",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );

// volScalarField Omega // reaction rate
// (
//     IOobject
//     (
//         "Omega",
//         runTime.timeName(),
//         mesh,
//         IOobject::NO_READ,
//         IOobject::NO_WRITE
//     ),
//     mesh
// );



Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
(
    IOobject
    (
        "transportProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);


Info<< "Reading transportProperties\n" << endl;

// Switch for Strang split or simple split
//const Switch Strang	(transportProperties.lookup("Strang"));
// const word timeSplitMethod (transportProperties.lookup("timeSplitMethod"));
// Info << "timeSplitMethod = " << timeSplitMethod << " \n" << endl;

// number of iteration in the ODE solving
//dimensionedScalar nODEstep_("nODEstep", dimless, transportProperties);

// molecular diffusivity
const dimensionedScalar D("D", dimViscosity, transportProperties);

// adsorbent density
const dimensionedScalar rhop("rhop", dimDensity, transportProperties);
// porosity
const dimensionedScalar porosity("porosity", dimless, transportProperties);


// Langmuir isotherm parameters
// the isotherm plot must be in the SI unit, C (kg/m3) and Qe (g/g)
const dimensionedScalar Langmuir_qMax("Langmuir_qMax", dimless, transportProperties);
//const dimensionedScalar Langmuir_b("Langmuir_b", dimless, transportProperties);

// effective diffusivity based on Millington-Quirk correlation
// const dimensionedScalar Dm = D*porosity/(Foam::pow(porosity,-1.0/3.0));

// hydrodynamics dispersion coefficient 
const dimensionedScalar Dh("Dh", dimViscosity, transportProperties);
// const dimensionedScalar Dh = Dm + Dm; // low reynold number asumption and ignore the transverse direction


// Forward and backward reaction rate, both has unit of 1/s
const dimensionedScalar kf("kf", dimensionSet(0,0,-1,0,0,0,0), transportProperties);
const dimensionedScalar kr("kr", dimensionSet(0,0,-1,0,0,0,0), transportProperties);


#include "createPhi.H"

#include "createFvOptions.H"
