Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();
thermo.correct();

volScalarField& p = thermo.p();
volScalarField& h = thermo.he();
volScalarField& psi = const_cast<volScalarField&>(thermo.psi());
const volScalarField& T = thermo.T();

Info<< "Reading field Urel\n" << endl;

volVectorField Urel
(
    IOobject
    (
        "Urel",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating SRF model\n" << endl;
autoPtr<SRF::SRFModel> SRF
(
    SRF::SRFModel::New(Urel)
);

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    Urel + SRF->U()
);

volVectorField Urot("Urot", SRF->U());

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField dpdt
(
    "dpdt",
    p / runTime.deltaT() * 0.0
);

volScalarField dissip
(
    "dissip",
    dpdt
);

volVectorField gradp
(
    "gradp",
    fvc::grad(p)
);
volVectorField limGradP
(
    "limGradP",
    fvc::grad(p)
);

DynamicList<label> badQualityCells(0);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
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
    linearInterpolate(rho*Urel) & mesh.Sf()
);


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        Urel,
        phi,
        thermo
    )
);

Info<< "Creating rotalpy energy i\n" << endl;

    volScalarField i
    (
        IOobject
        (
            "i",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        h
    );

i += 0.5*magSqr(Urel);
i -= 0.5*magSqr(Urot);
