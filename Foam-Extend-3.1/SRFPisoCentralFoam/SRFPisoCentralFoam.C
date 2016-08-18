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
#include "basicPsiThermo.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "cellQuality.H"

#include "RASModel.H"
#include "SRFModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readCourantType.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;
    
    #include "createSurfaceFields.H"
    #include "markBadQualityCells.H"
    
    while (runTime.run())
    {
	#include "readPISOControls.H"
	#include "acousticCourantNo.H"
	#include "compressibleCourantNo.H"
	#include "readTimeControls.H"
	#include "setDeltaT.H"
	#include "readFieldBounds.H"
	
	runTime++;
        /*
	psi.oldTime();
	rho.oldTime();
	p.oldTime();
	Urel.oldTime();
	h.oldTime();
	i.oldTime();*/
	psi.storePrevIter();
	rho.storePrevIter();
	p.storePrevIter();
	Urel.storePrevIter();
	h.storePrevIter();
	i.storePrevIter();

        Info<< "Time = " << runTime.timeName() << nl << endl;

	// --- Solve density
	solve
	(
	    fvm::ddt(rho) + fvc::div(phi)
	);
	
	// --- Solve momentum
	#include "UrelEqn.H"
	
	// --- Solve energy
	#include "iEqn.H"
	
	// --- Solve pressure (PISO)
	{
	    for (label corr = 1; corr <= nCorr; corr++)
	    {
		#include "pEqn.H"
	    }
	    #include "updateKappa.H"
	}
	
	// --- Solve turbulence
	turbulence->correct();

	U = Urel + SRF->U();
	
	dpdt = fvc::ddt(p);
	
	runTime.write();

	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
