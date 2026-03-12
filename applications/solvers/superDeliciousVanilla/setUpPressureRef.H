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
    point pRefPoint = Zero;
    label pRefProc = 0;
    scalar pRefValue = 0.0;
    vector hRef_ = Zero;
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
        
        // Lookup the pressure reference point and value from the PIMPLE dictionary.
        // Check to see if the user specified the pressure reference cell on the master.
        if (pimple.dict().found("pRefCell"))
        {
	    pRefProc = 0;
            if (Pstream::master())
            {
                pRefCell = readLabel(pimple.dict().lookup("pRefCell"));
                if (pRefCell < 0 || pRefCell >= mesh.nCells())
                {
                    FatalIOErrorInFunction
                    (
                        pimple.dict()
                    )   << "Illegal master cellID " << pRefCell
                        << ". Should be 0.." << mesh.nCells()
                        << exit(FatalIOError);
                }
            }
            else
            {
                pRefCell = -1;
            }
        }
        // Check to see if the user specified the pressure reference point.
        else if (pimple.dict().found("pRefPoint"))
        {
            point pRefPoint(pimple.dict().lookup("pRefPoint"));
            simpleMeshSearch::findCellSimple(mesh,pRefCell,pRefPoint,pRefProc);
        }
        // If nothing was specified, exit with an error message.
        else
        {
            FatalIOErrorInFunction
            (
                pimple.dict()
            )   << "Unable to set reference cell for pressure field "
                << nl
                << "    Please supply either pRefCell" 
                << " or pRefPoint" << nl << exit(FatalIOError);
        }

        pRefValue = readScalar(pimple.dict().lookup("pRefValue"));

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
