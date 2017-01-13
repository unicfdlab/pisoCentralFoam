/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    pisoCentralDyMFoam

Description
    Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor with support of moving meshes

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "cellQuality.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "createFields.H"
    #include "createMRF.H"
    #include "createTimeControls.H"
    #include "createFvOptions.H"
    bool checkMeshCourantNo =
            readBool(pimple.dict().lookup("checkMeshCourantNo"));
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    #include "initContinuityErrs.H"
    #include "readCourantType.H"
    Info<< "\nStarting time loop\n" << endl;
    
    #include "createSurfaceFields.H"
    Info<< "Creating turbulence model\n" << endl;
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New
        (
            rho,
            U,
            phi,
            thermo
        )
    );
    #include "markBadQualityCells.H"
    #include "createSurfaceFieldsDyM.H"
    #include "initKappaField.H"

    while (runTime.run())
    {
        #include "acousticCourantNo.H"
        #include "compressibleCourantNo.H"
        #include "readTimeControls.H"
        #include "setDeltaT.H"
        
        runTime++;
        
        psi.oldTime();
        rho.oldTime();
        p.oldTime();
        U.oldTime();
        h.oldTime();
        Ek.oldTime();
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // --- Move mesh and update fluxes
        {
            // Do any mesh changes
            mesh.update();
           
           //make fluxes relative
            if (mesh.changing())
            {
                meshPhi     = fvc::meshPhi(rho,U)();
                
                meshPhi_own = kappa*alpha_own*p_own*psi_own*meshPhi;
                meshPhi_nei = ((1.0 - kappa)*alpha_own*p_own*psi_own + alpha_nei*p_nei*psi_nei)*meshPhi;
                
                phi_own -= meshPhi_own;
                phi_nei -= meshPhi_nei;
                phi = phi_own + phi_nei;
                
                if (checkMeshCourantNo)
                {
                    #include "customMeshCourantNo.H"
                }
                
                #include "markBadQualityCells.H"
            }
        }
        
        // --- Solve density
        #include "massEqn.H"
        
        // --- Solve momentum
        #include "UEqn.H"
        
        // --- Solve energy
        #include "hEqn.H"
        
        // --- Solve pressure (PISO)
        {
            while (pimple.correct())
            {
                #include "pEqnDyM.H"
            }
        }
        
        // --- Solve turbulence
        turbulence->correct();

        dpdt = fvc::ddt(p) - fvc::div(meshPhi,p);
        EkChange = fvc::ddt(rho,Ek) + fvc::div(phi_own,Ek) + fvc::div(phi_nei,Ek);

        //make fluxes absolute
        phi_own += meshPhi_own;
        phi_nei += meshPhi_nei;
        phi = phi_own + phi_nei;

        //update kappa
        #include "updateKappa.H"
    
        //dpdt = fvc::ddt(p) - fvc::div(meshPhi*alpha_own*p_own) - fvc::div(meshPhi*alpha_nei*p_nei);
//        EkChange = fvc::ddt(rho,Ek) + fvc::div(phi_own,Ek) + fvc::div(phi_nei,Ek)
//                -fvc::div((meshPhi_nei + (1.0 - kappa)*meshPhi_own),Ek) - fvc::div(kappa*meshPhi_own,Ek);


        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
