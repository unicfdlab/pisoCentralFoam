volVectorField HbyA ("HbyA", Urel);
volScalarField rAUrel ("rAUrel", 1.0 / UrelEqn.A());
HbyA = UrelEqn.H() * rAUrel;

psi_pos = fvc::interpolate(psi, pos, "reconstruct(psi)");
psi_neg = fvc::interpolate(psi, neg, "reconstruct(psi)");

psiUrel_pos= fvc::interpolate(psi*HbyA, pos, "reconstruct(Urel)");
psiUrel_neg= fvc::interpolate(psi*HbyA, neg, "reconstruct(Urel)");

phiv_pos= (psiUrel_pos / psi_pos) & mesh.Sf();
phiv_neg= (psiUrel_neg / psi_neg) & mesh.Sf();
        
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

for (label nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
{
    fvScalarMatrix pEqn_pos
    (
	fvm::div(phid_pos,p) - fvm::laplacian(Dp_pos*a_pos,p)
    );
    
    fvScalarMatrix pEqn_neg
    (
	fvm::div(phid_neg,p) - fvm::laplacian(Dp_neg*a_neg,p)
    );
    
    solve
    (
	fvm::ddt(psi,p)
	+
	pEqn_pos
	+
	pEqn_neg
    );
    
    if (nonOrth == nNonOrthCorr)
    {
	phiPos = pEqn_pos.flux();
	phiNeg = pEqn_neg.flux();
	phi    = phiPos + phiNeg;
    }
}

p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
p_neg = fvc::interpolate(p, neg, "reconstruct(p)");

gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());
#include "limitPressureGradient.H"
Urel = HbyA - rAUrel * gradp;
Urel.correctBoundaryConditions();

Info << "max(Urel): " << max(Urel).value() << endl;

solve
(
    fvm::ddt(rho) + fvc::div(phi)
);
#include "compressibleContinuityErrs.H"

// Bound the pressure
if (min(p) < pMin || max(p) > pMax)
{
	p.max(pMin);
	p.min(pMax);
	p.correctBoundaryConditions();
}

/*// Bound the velocity
volScalarField magUrel = mag(Urel);

if (max(magUrel) > UrelMax)
{
	volScalarField Urellimiter =
	    pos(magUrel - UrelMax)*UrelMax/(magUrel + smallUrel)
	    + neg(magUrel - UrelMax);
	Urellimiter.max(scalar(0));
	Urellimiter.min(scalar(1));

	Urel *= Urellimiter;
	Urel.correctBoundaryConditions();
}*/

rho = thermo.rho();

// Bound rho
volScalarField R = thermo.Cp() - thermo.Cv();

volScalarField rhoMin = pMin/(R*TMax);
volScalarField rhoMax = pMax/(R*TMin);

rho = Foam::min(rho, rhoMax);
rho = Foam::max(rho, rhoMin);
rho.relax();
