    surfaceScalarField meshPhi
    (
        "volMeshPhi",
        phiv_own * 0.0
    );
    
    surfaceScalarField rel_phiv_own
    (
       "rel_phiv_own",
       phiv_own * 0.0
    );
    
    surfaceScalarField rel_phiv_nei
    (
       "rel_phiv_nei",
       phiv_own * 0.0
    );
    
//    surfaceScalarField meshPhid_own
//    (
//        "meshPhid_own",
//        alpha_own*psi_own*meshPhi
//    );
    
//    surfaceScalarField meshPhid_nei
//    (
//        "meshPhid_nei",
//        alpha_nei*psi_nei*meshPhi
//    );
    
    surfaceScalarField meshPhi_own
    (
        "meshPhi_own",
        phi_own * 0.0
    );
    
    surfaceScalarField meshPhi_nei
    (
        "meshPhi_nei",
        phi_nei * 0.0
    );
    