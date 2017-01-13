aphiv_own = phi_own / (p_own * psi_own);
aphiv_nei = phi_nei / (p_nei * psi_nei);

dimensionedScalar amaxSmall ("amaxSmall", amaxSf.dimensions(), SMALL * min(mesh.magSf()).value());
amaxSf = max(mag(aphiv_own), mag(aphiv_nei)) + amaxSmall;

surfaceScalarField amaxSfbyDelta
(
    mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
);

surfaceScalarField Maf
(
    mag(phi) / (psi_own*p_own*alpha_own + psi_nei*p_nei*alpha_nei)
    / (cSf_own*alpha_own + cSf_nei*alpha_nei)
);

Info << "max/min Maf: " << max(Maf).value() << "/" << min(Maf).value() << endl;

dimensionedScalar cDeltaT = runTime.deltaT();

#ifdef PISOCENTRALFOAM_LTS
if (LTS)
{
    cDeltaT.value() = 1.0 / gMax(trDeltaT().internalField());
}
#endif

kappa = 
    max
    (
        min
        (
            Maf / (amaxSfbyDelta/mesh.magSf() * cDeltaT),
            scalar(1.0)
        ),
        0.0
    );

forAll(kappa.boundaryField(), iPatch)
{
    fvsPatchField<scalar>& kappapf = kappa.boundaryField()[iPatch];
    if (isA<coupledFvsPatchField<scalar> > (kappapf))
    {
        forAll (kappapf, iFace)
        {
            kappapf[iFace] = 0.0;
        }
    }
}

Info << "max / min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
phi_nei += (1.0 - kappa) * phi_own;
phi_own *= kappa;



