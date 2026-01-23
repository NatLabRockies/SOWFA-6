{
    kappat = turbulence->nut()/Prt;
    kappat.correctBoundaryConditions();

    volScalarField kappaEff("kappaEff", turbulence->nu()/Pr + kappat);

    fvScalarMatrix TEqn
    (
        fvm::ddt(T)
      + fvm::div(phi, T)
      - fvm::laplacian(kappaEff, T)
     ==
        fvc::div(qwall)
      + temperatureGeoMesoTerm.force()
      + temperatureDampingLayers.source()
      + temperaturePerturbationZones.source()
      + fvOptions(T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

  //if (pimple.finalPimpleIter() && pimple.finalPISOIter())
    if (pimple.finalPimpleIter())
    {
        Info << "perturbing temperature field..." << endl;
        temperaturePerturbationZones.perturb(T);
    }
}
