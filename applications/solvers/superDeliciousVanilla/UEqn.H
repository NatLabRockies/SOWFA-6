    // Solve the momentum equation
    
    // Set up the solution matrix
    fvVectorMatrix UEqn
    (
        fvm::ddt(U) 
      + fvm::div(phi, U)
      + turbulence->divDevReff(U)
    ==
      - fvc::div(Rwall)
      + Coriolis.force()
      + momentumGeoMesoTerm.force()
      + hurricaneMesoTerm.source()
      + momentumDampingLayers.source()
      + momentumPerturbationZones.source()
      + (1.0 - scalar(bodyForceStabilization)) * actuatorTurbineArray.source()
      + fvOptions(U)
    );

    UEqn.relax();

    fvOptions.constrain(UEqn);

    if (pimple.momentumPredictor())
    {
        solve
        (
            UEqn
         ==
            fvc::reconstruct
            (
                (
                  - fvc::snGrad(p_rgh)
                  + Boussinesq.source()
                  + scalar(bodyForceStabilization) * fvc::flux(actuatorTurbineArray.source())/mesh.magSf()
                )*mesh.magSf()
            )
        );

        fvOptions.correct(U);
    }

  //if (pimple.finalPimpleIter() && pimple.finalPISOIter())
    if (pimple.finalPimpleIter())
    {
        momentumPerturbationZones.perturb(U);
    }
