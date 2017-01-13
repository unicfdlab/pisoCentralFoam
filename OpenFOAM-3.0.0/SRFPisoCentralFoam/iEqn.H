fvScalarMatrix iEqn
(
    fvm::ddt(rho,i)
    +
    fvm::div(phiPos,i)
    +
    fvm::div(phiNeg,i)
    -
    fvm::laplacian(turbulence->alphaEff(), i)
    ==
    dpdt
    +
    dissip
    +
    fvOptions(rho,i)
);

fvOptions.constrain(iEqn);

iEqn.solve();

// Calculate enthalpy out of rothalpy

h = i + 0.5*magSqr(Urot) - 0.5*magSqr(Urel);
h.correctBoundaryConditions();

// Bound the enthalpy using TMin and TMax
//volScalarField Cp = thermo.Cp();
//h = Foam::min(h, TMax*Cp);
//h = Foam::max(h, TMin*Cp);
//h.correctBoundaryConditions();

// Re-initialise rothalpy based on limited enthalpy
i = h - 0.5*magSqr(Urot) + 0.5*magSqr(Urel);

fvOptions.correct(h);

thermo.correct();

Info << "max/min T:" << max(T).value() << "/" << min(T).value() << endl;
