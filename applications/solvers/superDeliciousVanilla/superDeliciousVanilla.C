/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    superDeliciousVanilla (the code name for the new unified solver)

Description
    Transient solver for wind-energy and atmospheric boundary layer flows.
    The incompressible equations are solved along with a Boussinesq buoyancy
    term. There is flexibility in turbulence modeling. Although the solver
    is primarily meant for and tested for large-eddy simulations, it could
    also be used with a RANS turbulence model.

    The code includes planetary Coriolis forces, pressure gradient forces
    driving the wind, ability to use mesoscale forcings, coupling with
    atmospheric-style lower rough-wall boundary conditions, and coupling
    with actuator turbine models.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "turbulentTransportModel.H"
#include "wallDist.H"
#include "fixedFluxPressureFvPatchScalarField.H"
#include "processorFvPatch.H"
#include "timeVaryingMappedInletOutletFvPatchField.H"
#include "fvOptions.H"
#include "pimpleControl.H"
#include "ABL.H"
#include "adjustPhiExtended.H"
#include "horizontalAxisWindTurbinesALMOpenFAST.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    #include "computeDivergence.H"



    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    // This does an initial pressure solve assuming hydrostatic conditions.  It solves
    // dp/dx = buoyancy.	
    #include "initializePressureField.H"


    // Time stepping loop.
    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "setDeltaT.H"

        word timeNameOld = runTime.timeName();

        runTime++;

        Info << "Time Step = " << runTime.timeIndex() << ", ";
        Info << "Time = " << timeNameOld << " to " << runTime.timeName() << " s";
        Info << nl << endl;
      
        // Test to see if simulation begins without t^(n-1) (*_0) fields.  
        // Without those fields, the backward d/dt scheme reverts to Euler 
        // implicit in outer iteration 0, but uses full backward in 
        // iterations >0, and the U,T^(n-1) data is incorrect.  Outer
        // iterations >0 are only required to achieve 2nd-order accuracy
        // in time, but since that is not possible in the first time step 
        // without t^(n-1) data, we really only need 1 outer iteration, 
        // anyway.
        bool limitOuterLoop = ((U.nOldTimes() == 0) || (T.nOldTimes() == 0)) ? true : false;

        // Update the once-per-time-step source terms.
        // - perturbation zone forcing.
        momentumPerturbationZones.update();
        temperaturePerturbationZones.update();
        actuatorTurbineArray.update();
        hurricaneMesoTerm.update();


        // Outer-iteration loop.
        while (pimple.loop())
        {
            // Update the once-per-outer-iteration source terms.
            // - Geostrophic/mesoscale forcing.
            momentumGeoMesoTerm.update(pimple.finalPimpleIter());
            temperatureGeoMesoTerm.update(pimple.finalPimpleIter());

            // - Rayleigh/viscious damping forcing.
            momentumSpongeLayers.update();

            // - Coriolis forcing.
            Coriolis.update();

            // - Buoyancy forcing.
            Boussinesq.update();


            // Predictor step.
            Info << "   Predictor" << endl;

            #include "UEqn.H"
            #include "turbulenceCorrect.H"
            #include "TEqn.H"

            // Corrector steps.
            int corrIter = 1;
            while (pimple.correct())
            {
                Info << "   Corrector Step " << corrIter << endl;

                // Update the once-per-inner-iteration source terms.
                // - Buoyancy forcing.
              //Boussinesq.update();

                #include "pEqn.H"
              //#include "turbulenceCorrect.H"
              //#include "TEqn.H"
                corrIter++;
            }

            // If starting without the t^(n-1) data, advance through the remaining
            // outer iterations without doing anything by calling pimple.loop().
            if (limitOuterLoop)
            {
                while (pimple.loop());
                Info << endl;
                break;
            }

            Info << endl;
        }

        // Compute the continuity errors.
        #include "computeDivergence.H"

        // Update timeVaryingMappedInletOutlet parameters
        #include "updateFixesValue.H"

        // Write the solution if at write time.
        runTime.write();

        // Report timing.
        Info << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
             << "  ClockTime = " << runTime.elapsedClockTime() << " s"
             << nl << nl << nl << nl << endl;
    }

    actuatorTurbineArray.end();

    Info << "Ending the simulation" << endl;

    return 0;
}


// ************************************************************************* //
