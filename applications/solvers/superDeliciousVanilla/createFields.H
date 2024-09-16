// TAKE CARE OF MISCELLANEOUS INITIAL WORK
    // - Get access to the ABLProperties dictionary
    IOdictionary ABLProperties
    (
        IOobject
        (
            "ABLProperties",
            runTime.time().constant(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    
    
    //- Automatically update fixesValue for timeVaryingMappedInletOutlet BC?
    bool updateFixesValue(ABLProperties.lookupOrDefault<bool>("updateFixesValue",false));

   
    //- When using outflow pressure boundary conditions that are Neumann, the
    //  adjustPhi function is called to match outflow mass flow rate to that of
    //  the inflow boundaries.  In computing inflow and outflow mass flow rate, 
    //  adjustPhi does not include coupled boundaries, including those of type
    //  cyclic, cyclicAMI, cyclicACMI.  If you want to run a setup with a concurrent
    //  precursor, you will need to count one of the precursor boundary faces as
    //  an inflow face.  Therefore, SOWFA includes a modified adjustPhi function called
    //  adjustPhiExtended that will account for selected coupled boundaries as
    //  inflow, so read the list of those here.
    List<word> coupledInflowPatches;
    if (ABLProperties.found("coupledInflowPatches"))
    {
        List<word> coupledInflowPatches_(ABLProperties.lookup("coupledInflowPatches"));
        coupledInflowPatches = coupledInflowPatches_;
        Info << "Specifying the following coupled inflow patches: " << coupledInflowPatches << endl;
    }


    //- Indicate if this is a restart, or try to guess.
    bool isRestart = (runTime.timeIndex() == 0) ? false : true;
    isRestart = runTime.controlDict().lookupOrDefault<bool>("isRestart",isRestart);
 
    if (isRestart)
    {
        Info << "Treating this case as a restart of a previous case" << endl;
    }
    else
    {
        Info << "Treating this case as a cold start with no restart data" << endl;
    }

    if (isRestart && (runTime.timeIndex() == 0))
    {
        Warning << "User specified to start the solver in 'restart' mode with no " << endl
                << "no indication of previous solution history.  Beware that 'restart' "
                << "mode assumes all boundary conditions are updated already and will "
                << "not be updated before entering the time loop." << endl;
    }


    //- The initial pressure field is often set to zero, which is not a great initial
    //  condition.  A better option is to solve the pressure field for hydrostatic
    //  conditions.  If the user wants this, set this flag to true.
    bool initPressureField = isRestart ? false : true;
    initPressureField = runTime.controlDict().lookupOrDefault<bool>("initHydroStatic",initPressureField);


    //- Should the velocity flux non-coupled BCs be updated before entering the time loop.
    bool updatePhiUncoupledBoundariesOnRestart(ABLProperties.lookupOrDefault<bool>("updatePhiUncoupledBoundariesOnRestart",false));

    if (isRestart && updatePhiUncoupledBoundariesOnRestart)
    {
        Warning << "User specified to start solver in 'restart' mode and to also " << endl
                << "update the phi non-coupled boundary values, which will make the " << endl
                << "restart an 'imperfect' restart.  The values of the solution after " << endl
                << "start up will not be identical to the non-restarted case.  This " << endl
                << "configuration is meant for restarts where non-coupled boundary " << endl
                << "conditions have changed (i.e., periodic to inflow/outflow)." << endl;
    }






// CREATE AND READ THE RESOLVED-SCALE OR MEAN FIELDS
    // - Create and read the gravitational acceleration uniform field
     Info<< "Creating and reading the gravitational acceleration field, g" << endl;
     uniformDimensionedVectorField g
     (
         IOobject
         (
             "g",
             runTime.constant(),
             mesh,
             IOobject::MUST_READ,
             IOobject::NO_WRITE
         )
     );


    // - Create and read the modified pressure field [m^2/s^2].
    Info << "Creating and reading the density-normalized modified pressure field, p_rgh" << endl;
    volScalarField p_rgh
    (
        IOobject
        (
            "p_rgh",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    
    // - Create and read the density-normalized pressure field [m^2/s^2].
    Info << "Creating and reading the density-normalized pressure field, p" << endl;
    volScalarField p
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        p_rgh
    );



    // The quantity p_rgh is required on cell faces (flux) to project the velocity flux given the
    // updated pressure field, so add it to the list kept internally in OpenFOAM.
    mesh.setFluxRequired(p_rgh.name());
    
    
    // - Create and read the potential temperature field [K].
    Info << "Creating and reading potential temperature field, T" << endl;
    volScalarField T
    (
        IOobject
        (
            "T",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    
    // - Create and read the velocity field [m/s].
    Info << "Creating and reading velocity field, U" << endl;
    volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    
    // - Create and calculate the velocity flux [m^3/s].
    Info << "Creating and calculating velocity flux field, phi" << endl;
    IOobject phiFieldHeader
    (
        "phi",
         runTime.timeName(),
         mesh,
         IOobject::MUST_READ
     );

     bool phiPresent = (phiFieldHeader.typeHeaderOk<surfaceScalarField>(true)) ? true : false;
     if (phiPresent)
     {
         Info << "  -phi field file present, so using existing phi, but will update uncoupled boundaries if desired" << endl;
     }
     else
     {
         Info << "  -phi field file not present, so using simple interpolation of velocity to faces (i.e., fvc::flux(U), no Rhie-Chow interpolation)" << endl;
     }
    
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
         fvc::flux(U)
     );

  
     
     // If this is not a restart, but it is a cold start, then update the
     // U and T boundary conditions.  This is called twice in case mapped
     // patches are used to make sure that all patches are in sync.  For
     // example a mapped patch may be updated first and map from an un-
     // updated BC, so it will not effectively update in the first pass,
     // but it will in the second pass.
     if (!isRestart)
     {
         for (int i = 0; i < 2; i++)
         {
           //Info << "U.correctBoundaryConditions()" << endl;
             U.correctBoundaryConditions();
             phi = fvc::flux(U);

           //Info << "T.correctBoundaryConditions()" << endl;
             T.correctBoundaryConditions();
         }
     }


     // if the phi file was present, and the user wants the uncoupled boundaries
     // updated, do that here; otherwise recompute phi through interpolation of
     // velocity field.
     if (phiPresent)
     {
         if (updatePhiUncoupledBoundariesOnRestart)
         {
             Info << "Updating phi.boundaryField()" << endl;
             forAll(phi.boundaryField(), patchi)
             {
                 const fvPatchVectorField& UPatch = U.boundaryField()[patchi];
                 if (!phi.boundaryField()[patchi].coupled())
                 {
                     surfaceScalarField::Boundary& bphi = phi.boundaryFieldRef();
                     bphi[patchi] = UPatch & mesh.boundary()[patchi].Sf();
                 }
             }
         }
       //else
       //{
       //    Info << "Computing full phi field as phi = fvc::flux(U)" << endl;
       //    phi = fvc::flux(U);
       //}
     }




// SET UP PRESSURE REFERENCE LOCATION/HANDLING
    // This switch dictates whether or not the pressure reference cell is set explicitly
    // in the p_rgh solve.  If true, pressure is constrained at the pressure reference
    // cell by manipulating the matrix row for that cell to remove the null space.  This
    // assures that the pressure level is constrained, but it also means the continuity
    // equation is no longer solved at that cell, so the divergence error can be significant
    // enough there to cause localized oscillations.  Iterative solvers can still converge
    // even with a null space, but more iterations may be needed, so this switch can be
    // set to false.
    bool activatePressureRefCell(pimple.dict().lookupOrDefault<bool>("activatePressureRefCell", true));
    label pRefCell = 0;
    scalar pRefValue = 0.0;
    vector hRef_ = vector::zero;
    dimensionedVector hRef("hRef",dimLength,hRef_);

    if (p_rgh.needReference())
    {
        Info << "Pressure reference needed in this case" << endl;


        if (activatePressureRefCell)
        {
             Info << "Pressure reference cell matrix row modification enabled" << endl;
        }
        else
        {
             Info << "Pressure reference cell matrix row modification NOT enabled" << endl;
        }
        
        
        // Deal with reference pressure level. (When calling setRefCell, do NOT set the 
        // last argument to true!!!)
        Info << "Setting reference pressure level" << endl;
        setRefCell
        (
            p,
            p_rgh,
            pimple.dict(),
            pRefCell,
            pRefValue
      //    true
        );
        
        
        // Find the height of the pressure reference cell
        label pRefCellGlobal = 0;
        label pRefCellProc = 0;
        if (pRefCell != -1)
        {
            hRef_ = mesh.C()[pRefCell];
            pRefCellGlobal = pRefCell;
            pRefCellProc = Pstream::myProcNo();
        }
        reduce(hRef_,sumOp<vector>());
        reduce(pRefCellGlobal,sumOp<label>());
        reduce(pRefCellProc,sumOp<label>());
        hRef.value() = hRef_;
        Info << "Pressure reference location is " << hRef.value() << " m at cell " << pRefCellGlobal << " on processor " << pRefCellProc << endl;
    }

    else 
    {
        Info << "Pressure reference NOT needed in this case" << endl;
    }





 

// SET THE VISCOSITY MODEL, REFERENCE TEMPERATURE, AND PRANDTL NUMBER
    // - Set up the laminar (molecular) transport model.
    Info << "Creating and reading laminar transport properties" << endl;
    singlePhaseTransportModel laminarTransport(U, phi);
    
    
    // - Read in reference temperature [K] and laminar and turblent Prandtl numbers.
    //   Note that turbulent Prandtl number may be overriden by a stability-
    //   dependent value based on the choice of SGS turbulence model.
    dimensionedScalar TRef("TRef", dimTemperature, laminarTransport);
    dimensionedScalar Pr("Pr", dimless, laminarTransport);
    dimensionedScalar Prt("Prt", dimless, laminarTransport);




// INSTANTIATE OBJECTS OF THE SOURCE TERM CLASSES
    // - Coriolis force term [m/s^2].
    Info << "Creating Coriolis force" << endl;
    CoriolisForce Coriolis(U);
    
    
    // - Geostrophic/mesoscale source terms [<field units>/s].
    Info << "Creating geostrophic/mesoscale source term(s)" << endl;
    drivingForce<vector> momentumGeoMesoTerm(ABLProperties, "momentum", Coriolis, U);
    drivingForce<scalar> temperatureGeoMesoTerm(ABLProperties, "temperature", Coriolis, T);


    // - Hurricane source term [m/s^2].
    Info << "Creating G. Bryan's 'simple' hurricane source term" << endl;
    hurricaneSource hurricaneMesoTerm(ABLProperties, U);
    
    
    // - Buoyancy term [m/s^2].
    Info << "Creating buoyancy term" << endl;
    buoyancyModel Boussinesq(T, TRef, hRef, g, ABLProperties);
    #include "adjustPressureLevel.H"
    
    
    // - Rayleigh/viscous damping layer source terms [<field units>/s].
    Info << "Creating Rayleigh damping layer(s)" << endl;
    spongeLayer momentumSpongeLayers("spongeList",U);
    
    
    // - Flow perturbation zone [<field units>/s].
    Info << "Creating perturbation zone(s)" << endl;
    perturbationZone<vector> momentumPerturbationZones(ABLProperties,U);
    perturbationZone<scalar> temperaturePerturbationZones(ABLProperties,T);
    momentumPerturbationZones.perturb(U);
    temperaturePerturbationZones.perturb(T);


    // - Wind turbine actuator line/disk array [m/s^2].
    Info << "Creating actuator disk/line wind turbine array" << endl;
    turbineModels::horizontalAxisWindTurbinesALMOpenFAST actuatorTurbineArray(U);
    
    // Allows the user to apply body force with Rhie-Chow type interpolation for stabilization. See:
    //
    //    Rethore, P-E. M., & Sørensen, N. N., "Actuator Disc Model Using a Modified Rhie-Chow/SIMPLE Pressure
    //    Correction Algorithm. Comparison with Analytical Solutions." In Conference proceedings (online) European
    //    Wind Energy Association (EWEA), 2008.
    //
    bool bodyForceStabilization = mesh.schemesDict().lookupOrDefault<bool>("bodyForceStabilization",true);
    
    
    // - User set fvOptions, such as additional source terms.
    #include "createFvOptions.H"






// CREATE THE TURBULENCE QUANTITIES
    // - Set up the turbulence model.
    Info << "Creating and reading turbulent transport properties" << endl;
    autoPtr<incompressible::turbulenceModel> turbulence
    (
        incompressible::turbulenceModel::New(U, phi, laminarTransport)
    );
    
    
    // - Create and read the kinematic turbulent thermal conductivity [m2/s].
    Info << "Creating and reading turbulent thermal conductivity field, kappat" << endl;
    volScalarField kappat
    (
        IOobject
        (
            "kappat",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    
    // - Create and read the wall shear stress field [m^2/s^2].
    Info << "Creating and reading the wall shear stress field, Rwall" << endl;
    volSymmTensorField Rwall
    (
        IOobject
        (
            "Rwall",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
    
    
    // - Create and read the wall temperature flux field [K-m/s].
    Info << "Creating and reading the wall temperature flux field, qwall" << endl;
    volVectorField qwall
    (
        IOobject
        (
            "qwall",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
   
    Info << "Creating SGS temperature flux field, qSgs..." << endl;
    volVectorField qSgs
    (
        IOobject
        (
            "qSgs",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("qSgs",dimensionSet(0, 1, -1, 1, 0, 0, 0),Zero)
    );

   
    // Create the deviatoric SGS stress field
    Info << "Creating SGS stress field, RdevSgs..." << endl;
    volSymmTensorField RdevSgs
    (
        IOobject
        (
            "RdevSgs",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedSymmTensor("RdevSgs",dimensionSet(0, 2, -2, 0, 0, 0, 0),Zero)
    ); 


    // - if this is not a restart, update turbulence BCs and do an initial field solve.
    if (!isRestart)
    {
        kappat = turbulence->nut()/Prt;
        kappat.correctBoundaryConditions();    
        qwall.correctBoundaryConditions();    
        Rwall.correctBoundaryConditions();    
        turbulence->validate();
    }
