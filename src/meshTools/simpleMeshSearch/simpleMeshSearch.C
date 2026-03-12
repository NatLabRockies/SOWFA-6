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

\*---------------------------------------------------------------------------*/

#include "simpleMeshSearch.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

bool Foam::simpleMeshSearch::findCellSimple
(
    fvMesh& mesh,
    label& refCell,
    point& refPoint,
    label& refProc
)
{
    // Loop through all grid cells and compute distance to reference cell.
    // The cell with minimum distance is the one that might contain the
    // reference point.  Add a small random perturbation to each grid cell
    // center location to break ties between two cells.
    Random rndGen(123456);
    scalar minDis = VGREAT;
    label minDisCell = -1;

    point pert = 1.0E-4*(2.0*rndGen.sample01<vector>()-vector::one);
    forAll(mesh.C(),i)
    {
        point pert = 1.0E-4*(2.0*rndGen.sample01<vector>()-vector::one);
        point p = mesh.C()[i] + pert;
        scalar dis = Foam::mag(p - refPoint);
	if (dis < minDis)
        {
            minDis = dis;
            minDisCell = i;
	}
    }

    // Have all processors compare the minimum distance found, and choose
    // the smallest of those as the possible reference cell of all found
    // by the processors.
    scalar minDisGlobal = returnReduce<scalar>(minDis, minOp<scalar>());
    refCell = (minDis == minDisGlobal ? minDisCell : -1);
    refProc = (refCell == -1 ? 0 : Pstream::myProcNo());
    reduce(refProc, sumOp<label>());

    // Just because a minimum distance was found does not mean that the 
    // reference point actually lies within the found cell.  Check that
    // and return true only if the reference point actually lies within
    // the found cell. Here I use the pointInCell function to see if the
    // reference point is really inside the cell, and I use the FACE_PLANES
    // option, which seems most robust.  That code takes the dot product
    // of the face normal and the vector from face center to reference 
    // point for all cell faces.  If the dot product is positive (the 
    // reference point is pointed into the cell from the face) for all
    // faces, one can deduce that the reference point must be inside
    // the cell. 
    bool isInside = false;
    if (Pstream::myProcNo() == refProc)
    {
        isInside = mesh.pointInCell(refPoint,refCell,polyMesh::FACE_PLANES);
    }
    reduce(isInside, sumOp<bool>());

    return isInside;
}

// ************************************************************************* //
