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
    pisoCentralFoam

Description
    Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "psiThermo.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "cellQuality.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    
    #include "createTime.H"
    #include "createMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    scalar initialDeltaT = -VGREAT;
    #include "createFields.H"
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
    
    #include "createFvOptions.H"
    #include "createMRF.H"
    #include "initContinuityErrs.H"
    #include "readCourantType.H"
    
    #include "markBadQualityCells.H"
    
    #include "initKappaField.H"
    
    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
        #include "readAdditionalPimpleControl.H"
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "acousticCourantNo.H"
            #include "compressibleCourantNo.H"
            #include "readTimeControls.H"
            #include "setDeltaT.H"
        }
        
        runTime++;
        
        psi.oldTime();
        rho.oldTime();
        p.oldTime();
        U.oldTime();
        h.oldTime();
        K.oldTime();
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // --- Solve density
        #include "massEqn.H"
        
        // --- Solve turbulence
        turbulence->correct();
        
        // --- Solve momentum
        #include "UEqn.H"
        
        // --- Solve energy
        if (!updateEnergyInPISO)
        {
            #include "hEqn.H"
        }
        
        // --- Solve pressure (PISO)
        {
            while (pimple.correct())
            {
                if (updateEnergyInPISO) //update each iteration before pressure
                {
                    #include "hEqn.H"
                }
                #include "pEqn.H"
                if (updateEnergyInPISO)
                {
                    #define PISOCENTRALFOAM_LTS
                    #include "updateKappa.H"
                    #include "updateMechanicalFields.H"
                }
            }
            
            if (!updateEnergyInPISO)
            {
                #include "updateKappa.H"
                #include "updateMechanicalFields.H"
            }
        }
        

        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
