psi_own   = fvc::interpolate(psi, own, "reconstruct(psi)");
psi_nei   = fvc::interpolate(psi, nei, "reconstruct(psi)");

p_own     = fvc::interpolate(p, own, "reconstruct(p)");
p_nei     = fvc::interpolate(p, nei, "reconstruct(p)");

psiU_own  = fvc::interpolate(psi*U, own, "reconstruct(U)");
psiU_nei  = fvc::interpolate(psi*U, nei, "reconstruct(U)");

phiv_own  = (psiU_own / psi_own) & mesh.Sf();
phiv_nei  = (psiU_nei / psi_nei) & mesh.Sf();

MRF.makeRelative(phiv_own);
MRF.makeRelative(phiv_nei);

c         = sqrt(thermo.Cp()/thermo.Cv() / psi);

cSf_own   = fvc::interpolate(c, own, "reconstruct(psi)")*mesh.magSf();
cSf_nei   = fvc::interpolate(c, nei, "reconstruct(psi)")*mesh.magSf();

ap        = max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero);
am        = min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero);

alpha_own = ap/(ap - am);
aSf       = am*alpha_own;
alpha_nei = 1.0 - alpha_own;

phiv_own *= alpha_own;
phiv_nei *= alpha_nei;
aphiv_own = phiv_own + aSf;
aphiv_nei = phiv_nei - aSf;

phi_own   = aphiv_own * p_own * psi_own;
phi_nei   = aphiv_nei * p_nei * psi_nei;
phi = phi_own + phi_nei;

#include "updateKappa.H"

