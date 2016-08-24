
volVectorField HbyA ("HbyA", Urel);
volScalarField rAUrel ("rAUrel", 1.0 / UrelEqn.A());
HbyA = UrelEqn.H() * rAUrel;

psi_pos = fvc::interpolate(psi, pos, "reconstruct(psi)");
psi_neg = fvc::interpolate(psi, neg, "reconstruct(psi)");

psiUrel_pos= fvc::interpolate(psi*HbyA, pos, "reconstruct(Urel)");
psiUrel_neg= fvc::interpolate(psi*HbyA, neg, "reconstruct(Urel)");

phiv_pos= (psiUrel_pos / psi_pos) & mesh.Sf();
phiv_neg= (psiUrel_neg / psi_neg) & mesh.Sf();

MRF.makeRelative(phiv_pos);
MRF.makeRelative(phiv_neg);

c = sqrt(thermo.Cp()/thermo.Cv() / psi);

cSf_pos = fvc::interpolate(c, pos, "reconstruct(psi)")*mesh.magSf();
cSf_neg = fvc::interpolate(c, neg, "reconstruct(psi)")*mesh.magSf();

ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);

a_pos = ap/(ap - am);
aSf = am*a_pos;
a_neg = 1.0 - a_pos;

phiv_pos *= a_pos;
phiv_neg *= a_neg;
aphiv_pos = phiv_pos - aSf;
aphiv_neg = phiv_neg + aSf;
phid_pos = aphiv_pos * psi_pos;
phid_neg = aphiv_neg * psi_neg;

surfaceScalarField Dp_pos
(
    "Dp_pos",
    fvc::interpolate(rho*rAUrel, pos, "reconstruct(Dp)")
);

surfaceScalarField Dp_neg
(
    "Dp_neg",
    fvc::interpolate(rho*rAUrel, neg, "reconstruct(Dp)")
);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn_pos
    (
	fvm::div(phid_pos,p) - fvm::laplacian(Dp_pos*a_pos,p)
    );

    fvScalarMatrix pEqn_neg
    (
	fvm::div(phid_neg,p) - fvm::laplacian(Dp_neg*a_neg,p)
    );

    fvScalarMatrix pEqn
    (
        fvm::ddt(psi,p)
        +
        pEqn_pos
        +
        pEqn_neg
        ==
        fvOptions(psi,p,rho.name())
    );

    fvOptions.constrain(pEqn);

    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        phiPos = pEqn_pos.flux();
        phiNeg = pEqn_neg.flux();
        phi = phiPos + phiNeg;
    }
}

p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
p_neg = fvc::interpolate(p, neg, "reconstruct(p)");

gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());
#include "limitPressureGradient.H"
Urel = HbyA - rAUrel * gradp;
Urel.correctBoundaryConditions();
fvOptions.correct(Urel);

Info << "max(Urel): " << max(Urel).value() << endl;

Info << "max/min of p: " << max(p).value() << "/" << min(p).value() << endl;

solve
(
    fvm::ddt(rho) + fvc::div(phi)
);

#include "compressibleContinuityErrs.H"

rho = thermo.rho();
