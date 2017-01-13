
volVectorField HbyA ("HbyA", U);
volScalarField rAU ("rAU", 1.0 / UEqn.A());
HbyA = UEqn.H() * rAU;

psi_own = fvc::interpolate(psi, own, "reconstruct(psi)");
psi_nei = fvc::interpolate(psi, nei, "reconstruct(psi)");

psiU_own= fvc::interpolate(psi*HbyA, own, "reconstruct(U)");
psiU_nei= fvc::interpolate(psi*HbyA, nei, "reconstruct(U)");

phiv_own= (psiU_own / psi_own) & mesh.Sf();
phiv_nei= (psiU_nei / psi_nei) & mesh.Sf();

MRF.makeRelative(phiv_own);
MRF.makeRelative(phiv_nei);

c = sqrt(thermo.Cp()/thermo.Cv() / psi);

cSf_own = fvc::interpolate(c, own, "reconstruct(psi)")*mesh.magSf();
cSf_nei = fvc::interpolate(c, nei, "reconstruct(psi)")*mesh.magSf();

ap = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
am = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);

alpha_own = ap/(ap - am);
aSf = am*alpha_own;
alpha_nei = 1.0 - alpha_own;

phiv_own *= alpha_own;
phiv_nei *= alpha_nei;

aphiv_own = phiv_own - aSf;
aphiv_nei = phiv_nei + aSf;

phid_own  = aphiv_own * psi_own;
phid_nei  = aphiv_nei * psi_nei;

Dp_own    = fvc::interpolate(rho*rAU, own, "reconstruct(Dp)");
Dp_nei    = fvc::interpolate(rho*rAU, nei, "reconstruct(Dp)");

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pEqn_own
    (
	fvm::div(phid_own,p) - fvm::laplacian(Dp_own*alpha_own,p)
    );
    
    fvScalarMatrix pEqn_nei
    (
	fvm::div(phid_nei,p) - fvm::laplacian(Dp_nei*alpha_nei,p)
    );
    
    fvScalarMatrix pEqn
    (
        fvm::ddt(psi,p)
        +
        pEqn_own
        +
        pEqn_nei
        ==
        fvOptions(psi,p,rho.name())
    );
    
    fvOptions.constrain(pEqn);
    
    pEqn.solve(mesh.solver(p.select(pimple.finalInnerIter())));
    
    if (pimple.finalNonOrthogonalIter())
    {
        phi_own = pEqn_own.flux();
        phi_nei = pEqn_nei.flux();
        phi = phi_own + phi_nei;
    }
}

p_own = fvc::interpolate(p, own, "reconstruct(p)");
p_nei = fvc::interpolate(p, nei, "reconstruct(p)");

rho_own = p_own*psi_own;
rho_nei = p_nei*psi_nei;

gradp = fvc::div((alpha_own*p_own + alpha_nei*p_nei)*mesh.Sf());
#include "limitPressureGradient.H"
U = HbyA - rAU * gradp;
U.correctBoundaryConditions();
fvOptions.correct(U);

Info << "max(U): " << max(U).value() << endl;

Info << "max/min p: " << max(p).value() << "/" << min(p).value() << endl;

#include "massEqn.H"

#include "compressibleContinuityErrs.H"

rho = thermo.rho();
