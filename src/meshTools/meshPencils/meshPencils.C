/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "meshPencils.H"
#include "ListListOps.H"
#include "PstreamCombineReduceOps.H"
#include "mergePoints.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(meshPencils, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::meshPencils::findPencils_()
{
    // On each processor, do a search through all cells and make a list of 
    // all distinct pencil locations. For illustration, take a pencil that is
    // along z.  The pencil "location" is it's (x,y) values.  It's "position"
    // is it's z value.

    // This is the processor's own list of distinct pencil locations. The list
    // is initialized as empty.
    pointField pencilLocations(0);

    // Loop over all cells.
    forAll(mesh_.cells(),cellI)
    {

        // Get the cell's "location" and "position".
        const point& cellCenter = mesh_.C()[cellI];
        point cellLoc(Zero);
        cellLoc[0] = cellCenter[orthogonalDirectionIndices_[0]];
        cellLoc[1] = cellCenter[orthogonalDirectionIndices_[1]];
        scalar cellPos = cellCenter[directionIndex_];

        // If the cell is within the position extents along the pencil,
        // then proceed
        if ( (cellPos >= extentLowerLim_) && 
             (cellPos <= extentUpperLim_) )
        {
            // If the pencil location list is empty, add this cell's location
            // to the pencil location list; otherwise...
            bool addToList = true;
            if (!pencilLocations.empty())
            {
                // ...compare this cell's location against all pencil locations
                // in the list, and if it is in common, don't add to the list.
                // In fact, upon the first find that it is in common, break out
                // of the for loop and move on.
                forAll(pencilLocations,pencilI)
                {
                    const point pencilLoc = pencilLocations[pencilI];
                    if( (mag(cellLoc[0] - pencilLoc[0]) < tol_) &&
                        (mag(cellLoc[1] - pencilLoc[1]) < tol_) )
                    {
                        addToList = false;
                        break;
                    }
                }
            }

            // If the cell location passes the test, add it to the pencil location list.
            if (addToList)
            {
                pencilLocations.append(cellLoc);
            }
        }
    }


    // Combine the plane location lists from each processor onto the master
    combineReduce(pencilLocations,ListPlusEqOp<pointField>());


    // Remove duplicate locations using the built-in mergePoints functionality.
    labelList pMap;
    Foam::mergePoints
    (
        pencilLocations,
        tol_,
        false,
        pMap,
        pencilLocations_
    );
    numberOfPencils_ = pencilLocations_.size();



    // Make a list of lists of cell ID labels.  Each list within the list contains
    // the cell ID labels corresponding to a specific pencil on a processor.  The overall list
    // should contain as many cell ID labels lists as there are distinct pencils.
    // Also sum up the volume of cells in each pencil on each processor.
    forAll(pencilLocations_,pencilI)
    {
        const point pencilLoc = pencilLocations_[pencilI];

        numCellPerPencil_.append(0);
        totVolPerPencil_.append(0.0);

        DynamicList<label> pencilCells;

        forAll(mesh_.cells(),cellI)
        {
            const point& cellCenter = mesh_.C()[cellI];
            point cellLoc;
            cellLoc[0] = cellCenter[orthogonalDirectionIndices_[0]];
            cellLoc[1] = cellCenter[orthogonalDirectionIndices_[1]];
            scalar cellPos = cellCenter[directionIndex_];
            if ( (cellPos >= extentLowerLim_) && 
                 (cellPos <= extentUpperLim_) )
            {
                if( (mag(cellLoc[0] - pencilLoc[0]) < tol_) &&
                    (mag(cellLoc[1] - pencilLoc[1]) < tol_) )
                {
                    pencilCells.append(cellI);
                    numCellPerPencil_[pencilI]++;
                    totVolPerPencil_[pencilI] += mesh_.V()[cellI];
                }
            }
        }
        pencilCellList_.append(pencilCells);
    }


    // Gather the volumes per plane to get a global value.
    reduce(numCellPerPencil_,sumOp<List<label> >());
    reduce(totVolPerPencil_,sumOp<List<scalar> >());
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::meshPencils::meshPencils
(
    const fvMesh& mesh,
    const scalar extentLowerLim,
    const scalar extentUpperLim,
    const word direction
)
:
    mesh_(mesh),
    extentLowerLim_(extentLowerLim),
    extentUpperLim_(extentUpperLim)
{

    if (direction == "x")
    {
        directionIndex_ = 0;
        orthogonalDirectionIndices_[0] = 1;
        orthogonalDirectionIndices_[1] = 2;
    }
    else if (direction == "y")
    {
        directionIndex_ = 1;
        orthogonalDirectionIndices_[0] = 0;
        orthogonalDirectionIndices_[1] = 2;
    }
    else if (direction == "z")
    {
        directionIndex_ = 2;
        orthogonalDirectionIndices_[0] = 0;
        orthogonalDirectionIndices_[1] = 1;
    }

    // Set initialization flag
    initialized_ = false;

    // Initialize
    initialize();
}


Foam::meshPencils::meshPencils
(
    const fvMesh& mesh
)
:
    mesh_(mesh),
    extentLowerLim_(-1.0E9),
    extentUpperLim_( 1.0E9)
{
    // Default to z-pencils
    directionIndex_ = 2;
    orthogonalDirectionIndices_[0] = 0;
    orthogonalDirectionIndices_[1] = 1;

    // Set initialization flag
    initialized_ = false;

    // Initialize
    initialize();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::meshPencils::~meshPencils()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::meshPencils::initialize()
{
    if (not initialized_)
    {
        tol_ = 1.0E-8 * mesh_.bounds().mag();
        findPencils_();
        initialized_ = true;
    }
}



void Foam::meshPencils::findMax
(
    const volScalarField& vField,
    scalarField vMax,
    scalarField maxLoc
)
{
    
    Info << " " << endl;    
}


/*
template<class Type>
List<Type> Foam::meshPencils::average
(
    const GeometricField<Type, fvPatchField, volMesh>& vField
)
{
    List<Type> vFieldMean(numberOfPlanes_,Type::zero);

    forAll(planeLocationValues_,planeI)
    {
        for (label i = 0; i < numCellPerPlane_[planeI]; i++)
        {
            label cellI = planesCellList_[planeI][i];
            vFieldMean[planeI] += vField[cellI] * mesh_.V()[cellI];
        }
    }

    reduce(vFieldMean,sumOp<List<Type> >());

    forAll(planeLocationValues_,planeI)
    {
        vFieldMean[planeI] /= totVolPerPlane_[planeI];
    }

    return vFieldMean;
}


template<class Type>
Type Foam::meshPencils::average
(
    const GeometricField<Type, fvPatchField, volMesh>& vField,
    label planeI
)
{
    Type vFieldMean(Type::zero);

    for (label i = 0; i < numCellPerPlane_[planeI]; i++)
    {
        label cellI = planesCellList_[planeI][i];
        vFieldMean += vField[cellI] * mesh_.V()[cellI];
    }

    reduce(vFieldMean,sumOp<Type>());

    vFieldMean /= totVolPerPlane_[planeI];

    return vFieldMean;
}


// Explicit instantiation of template function average<Type>
template vector Foam::meshPencils::average<vector>(const volVectorField&, label);
template List<vector> Foam::meshPencils::average<vector>(const volVectorField&);
template symmTensor Foam::meshPencils::average<symmTensor>(const volSymmTensorField&, label);
template List<symmTensor> Foam::meshPencils::average<symmTensor>(const volSymmTensorField&);

// Specialization for average<scalar>
// Reason: scalar::zero does not exist, so initialization must use 0.0 instead
// Note that template specialization must be within Foam namespace,
// and providing the namespace as part of the type does not work in this case.
namespace Foam
{
    template<>
    List<scalar> meshPencils::average<scalar>
    (
        const GeometricField<scalar, fvPatchField, volMesh>& vField
    )
    {
        List<scalar> vFieldMean(numberOfPlanes_,0.0);
    
        forAll(planeLocationValues_,planeI)
        {
            for (label i = 0; i < numCellPerPlane_[planeI]; i++)
            {
                label cellI = planesCellList_[planeI][i];
                vFieldMean[planeI] += vField[cellI] * mesh_.V()[cellI];
            }
        }
    
        reduce(vFieldMean,sumOp<List<scalar> >());
    
        forAll(planeLocationValues_,planeI)
        {
            vFieldMean[planeI] /= totVolPerPlane_[planeI];
        }
    
        return vFieldMean;
    }


    template<>
    scalar meshPencils::average<scalar>
    (
        const GeometricField<scalar, fvPatchField, volMesh>& vField,
        label planeI
    )
    {
        scalar vFieldMean(0.0);
    
        for (label i = 0; i < numCellPerPlane_[planeI]; i++)
        {
            label cellI = planesCellList_[planeI][i];
            vFieldMean += vField[cellI] * mesh_.V()[cellI];
        }
    
        reduce(vFieldMean,sumOp<scalar>());
    
        vFieldMean /= totVolPerPlane_[planeI];
    
        return vFieldMean;
    }
}
*/


// ************************************************************************* //
