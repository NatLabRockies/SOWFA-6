// This does an initial pressure solve assuming hydrostatic conditions.  It solves
// dp/dx = buoyancy.
if (initPressureField)
{
    Info << "Performing initialization of pressure field to hydrostatic conditions..." << endl;

    // Get the buoyancy term.
    Boussinesq.update();
    surfaceScalarField buoyancyTerm = Boussinesq.source();

    // If using fixedFluxPressure boundary conditions, set the gradient
    // at the boundary to the bouyancy term.
    setSnGrad<fixedFluxPressureFvPatchScalarField>
    (
        p_rgh.boundaryFieldRef(),
        buoyancyTerm.boundaryField()
    );

    // Start with dp_rgh/dx = bouyancy.  Take the gradient of both sides to get
    // laplacian(p_rgh) = gradient(buoyancy).  Solve that for p_rgh.
    fvScalarMatrix p_rghEqn
    (
        fvm::laplacian(p_rgh) == fvc::div(buoyancyTerm*mesh.magSf())
    );
    p_rghEqn.solve(mesh.solver(p_rgh.select(pimple.finalInnerIter())));

    // Adjust the overall pressure level.
    #include "adjustPressureLevel.H"
}
