    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    surfaceScalarField own
    (
        IOobject
        (
            "own",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("own", dimless, 1.0)
    );
    
    surfaceScalarField nei
    (
        IOobject
        (
            "nei",
            runTime.timeName(),
            mesh
        ),
        mesh,
        dimensionedScalar("neg", dimless, -1.0)
    );
    
    surfaceScalarField phi
    (
        IOobject
        (
            "phi",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        linearInterpolate(rho*U) & mesh.Sf()
    );
    
    surfaceScalarField psi_own
    (
        "psi_own",
        fvc::interpolate(psi, own, "reconstruct(psi)")
    );
    surfaceScalarField psi_nei
    (
        "psi_nei",
        fvc::interpolate(psi, nei, "reconstruct(psi)")
    );
    
    surfaceVectorField psiU_own
    (
        "psiU_own",
        fvc::interpolate(psi*U, nei, "reconstruct(U)")
    );
    surfaceVectorField psiU_nei
    (
        "psiU_nei",
        fvc::interpolate(psi*U, nei, "reconstruct(U)")
    );
    
    surfaceScalarField p_own
    (
        "p_own",
        fvc::interpolate(p, own, "reconstruct(p)")
    );
    surfaceScalarField p_nei
    (
        "p_nei",
        fvc::interpolate(p, nei, "reconstruct(p)")
    );
    
    surfaceScalarField rho_own
    (
        "rho_own",
        fvc::interpolate(rho, own, "reconstruct(rho)")
    );
    surfaceScalarField rho_nei
    (
        "rho_nei",
        fvc::interpolate(rho, nei, "reconstruct(rho)")
    );
    
    surfaceScalarField phiv_own
    (
        "phiv_own",
        (psiU_own / psi_own)  & mesh.Sf()
    );
    surfaceScalarField phiv_nei
    (
        "phiv_nei",
        (psiU_nei / psi_nei) & mesh.Sf()
    );

    surfaceScalarField cSf_own
    (
        "cSf_own",
        fvc::interpolate(c, own, "reconstruct(psi)")*mesh.magSf()
    );
    surfaceScalarField cSf_nei
    (
        "cSf_nei",
        fvc::interpolate(c, nei, "reconstruct(psi)")*mesh.magSf()
    );

    surfaceScalarField ap
    (
        "ap",
        max(max(phiv_own + cSf_own, phiv_nei + cSf_nei), v_zero)
    );
    surfaceScalarField am
    (
        "am",
        min(min(phiv_own - cSf_own, phiv_nei - cSf_nei), v_zero)
    );

    surfaceScalarField alpha_own
    (
        "alpha_own",
        ap/(ap - am)
    );
    surfaceScalarField amaxSf
    (
        "amaxSf",
        max(mag(am), mag(ap))
    );
    surfaceScalarField aSf
    (
        "aSf",
        am*alpha_own
    );
    surfaceScalarField alpha_nei
    (
        "a_nei",
        1.0 - alpha_own
    );

    phiv_own *= alpha_own;
    phiv_nei *= alpha_nei;

    surfaceScalarField aphiv_own
    (
        "aphiv_own",
        phiv_own + aSf
    );
    surfaceScalarField aphiv_nei
    (
        "aphiv_nei",
        phiv_nei - aSf
    );

    surfaceScalarField phid_own
    (
        "phid_own",
        aphiv_own * psi_own
    );
    surfaceScalarField phid_nei
    (
        "phid_nei",
        aphiv_nei * psi_nei
    );

    surfaceScalarField phi_own
    (
        "phi_own",
        phi * 0.0
    );

    surfaceScalarField phi_nei
    (
        "phi_nei",
        phi
    );

    surfaceScalarField kappa
    (
        IOobject
        (
            "kappa",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimless
    );
    forAll(kappa, iFace)
    {
        kappa[iFace] = 0.0;
    }
    forAll(kappa.boundaryField(), iPatch)
    {
        kappa.boundaryField()[iPatch] = 0.0;
    }

    surfaceScalarField Dp_own
    (
        "Dp_own",
        fvc::interpolate(rho*dimensionedScalar("zero", dimTime*dimVolume/dimMass, 0.0), own, "reconstruct(Dp)")
    );
    
    surfaceScalarField Dp_nei
    (
        "Dp_nei",
        fvc::interpolate(rho*dimensionedScalar("zero", dimTime*dimVolume/dimMass, 0.0), nei, "reconstruct(Dp)")
    );
    
