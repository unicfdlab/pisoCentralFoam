fvVectorMatrix UrelEqn
(
    fvm::ddt(rho,Urel)
    +
    fvm::div(phiPos,Urel)
    +
    fvm::div(phiNeg,Urel)
    +
    rho*SRF->Su()
    +
    turbulence->divDevRhoReff(Urel)
    ==
    fvOptions(rho, Urel)
);

fvOptions.constrain(UrelEqn);

if (pimple.momentumPredictor())
{
    solve
    (
	UrelEqn == -gradp
    );

    fvOptions.correct(Urel);
}
