    K = 0.5*magSqr(U);
    KChange = fvc::ddt(rho,K) + fvc::div(phi_own,K) + fvc::div(phi_nei,K);
    KSource = fvc::ddt(p) + fvc::div
            (
                (
                    linearInterpolate
                    (
                        (-turbulence->devRhoReff()) & U
                    ) & mesh.Sf()
                )()
            );

