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

#include "dampingLayer.H"
#include "geometricTransformations.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


//- Initialize this class
template<class Type>
void Foam::dampingLayer<Type>::initialize()
{
    t_ = runTime_.value();

    subDict_ = dict_.subOrEmptyDict("dampingLayers." & field_.name());
    subDictList_ = subDict_.toc();

    // Initial messages.
    nLayers_ = subDictList_.size();
    if (nLayers_ > 0)
    {
        Info << "  -using " << nLayers_ << " damping layers for field "
             << field_.name() << "..." << endl;
    }
    else
    {
        Info << "  -no damping layers specified for field "
             << field_.name() << ". Skipping..." << endl;
    }

    // Size the lists.
    adjacentBoundary_.setSize(nLayers_);
    dampingFunctionType_.setSize(nLayers_);
    useWallDist_.setSize(nLayers_);
    layerThickness_.setSize(nLayers_);
    dampingTimeScale_.setSize(nLayers_);
    boundaryNormal_.setSize(nLayers_);
    boundaryPoint_.setSize(nLayers_);
    gridCellList_.setSize(nLayers_);
    distanceFromBoundary_.setSize(nLayers_);

    
    // Call functions required to initialize this damping layers.
    readSubDict();
    inputChecks();
    defineWallDist();
    findBoundaryNormal();
    findDistanceFromBoundary();
}



//- Read the sub-dictionaries
template<class Type>
void Foam::dampingLayer<Type>::readSubDict()
{
    for (int m = 0; m < nLayers_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));
        
        adjacentBoundary_[m] = subSubDict.lookupOrDefault<word>("associatedBoundary","none");
    
        useWallDist_[m] = subSubDict.lookupOrDefault<bool>("useWallDistance",0);
    
        layerThickness_[m] = subSubDict.lookupOrDefault<scalar>("thickness", Zero);

        dampingTimeScale_[m] = subSubDict.lookupOrDefault<scalar>("timeScale", Foam::VGREAT);

        dampingFunctionType_[m] = subSubDict.lookupOrDefault<word>("function","none");
    }
}



//- Check for necessary inputs.  If they are not there, throw an error.
template<class Type>
void Foam::dampingLayer<Type>::inputChecks()
{
    for (int m = 0; m < nLayers_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));

        // Check to see that associated boundary name is provided.  
        if (!subSubDict.found("associatedBoundary"))
        {
            FatalErrorInFunction << "Must specify  'associatedBoundary'."
                                 << abort(FatalError);
        }
    }
}



//- Check to see if wall distance will be used.  If so, compute it.
template<class Type>
void Foam::dampingLayer<Type>::defineWallDist()
{
    // The wall distance function can be expensive, so only make an object of the wall
    // distance function if necessary.

    useWallDistAny_ = false;
    forAll(useWallDist_,m)
    {
        if (useWallDist_[m])
        {
            useWallDistAny_ = 1;
        }
    }

    if (useWallDistAny_)
    {
        wallDist wallDistance(mesh_);
        zAgl_ = wallDistance.y();
    }
}



//- Find boundary normal.
template<class Type>
void Foam::dampingLayer<Type>::findBoundaryNormal()
{
    for (int m = 0; m < nLayers_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));
    
        // Find the patch number for the associated patch.  All processors
        // seem to know about non-processor patches even if they do not contain
        // any patch faces.  If, for some reason, the processor has no knowledge
        // of the patch, the patch number will remain -1.
        label patchNum = -1;    
        forAll(mesh_.boundary(),i)
        {
            const word patchName = mesh_.boundary()[i].name();
            patchNum = (adjacentBoundary_[m] == patchName) ? i : patchNum;
        }
       
        // Check to see how many patch points of the associated lateral boundary
        // that this processor contains.  The localPoints() function gets just the
        // points that this processor uses, as opposed to the points() function that
        // does some parallel reduce.
        const vectorField& boundaryPointsLocal = mesh_.boundaryMesh()[patchNum].localPoints();
        label patchSize = boundaryPointsLocal.size();
    
        // If this processor contains at least some of the patch faces, set
        // the hasPatch variable to 1.
        label hasPatch = (patchSize > 0) ? 1 : 0;
    
        // Get the normal unit vector of the first patch face.  We assume that the patch is
        // planar so all faces should have the normal of the patch.  Only do this if
        // this processor actually has faces on the patch, but then parallel communicate
        // the result. In OpenFOAM, the boundary normal points out of the domain, so we want
        // the negative of it that points into the domain.
        vector boundaryNormal = vector::zero;
        if (hasPatch == 1)
        {
            boundaryNormal = -mesh_.Sf().boundaryField()[patchNum][0];
            boundaryNormal /= mag(boundaryNormal);
        }
        reduce(boundaryNormal, sumOp<vector>());
        reduce(hasPatch, sumOp<label>());
        boundaryNormal /= scalar(hasPatch);

        boundaryNormal_[m] = boundaryNormal;

        // Get a point on the patch. Do this by averaging all points on the patch.
        const vectorField& boundaryPoints = mesh_.boundaryMesh()[patchNum].points();
        scalar nBoundaryPoints = scalar(boundaryPoints.size());
        point boundaryPoint = vector::zero;
        forAll (boundaryPoints, i)
        {
            boundaryPoint += boundaryPoints[i];
        }
        boundaryPoint /= nBoundaryPoints;
        boundaryPoint_[m] = boundaryPoint;
    }
}



//- Find normal distance from damping start plane for each cell.
template<class Type>
void Foam::dampingLayer<Type>::findDistanceFromBoundary()
{
    for (int m = 0; m < nLayers_; m++)
    {
        DynamicList<label> gridCellList;
        DynamicList<scalar>distanceFromBoundary;
        forAll(mesh_.C(),j)
        {
            point meshPoint = mesh_.C()[j];
            if (useWallDist_[m])
            {
                meshPoint.z() = zAgl_[j];
            }
            scalar distance = (meshPoint - boundaryPoint_[m]) & boundaryNormal_[m];

            if (distance <= layerThickness_[m])
            {
                gridCellList.append(j);
                distanceFromBoundary.append(distance);
            }
        }
        gridCellList_[m] = gridCellList;
        distanceFromBoundary_[m] = distanceFromBoundary;
    }
}



//- Update the damping source term.
template<class Type>
void Foam::dampingLayer<Type>::update()
{
    for (int m = 0; m < nLayers_; m++)
    {
        forAll(gridCellList_[m], j)
        {
            Type source = Zero;

            if ((dampingFunctionType_[m] == "sineSquared") ||
                (dampingFunctionType_[m] == "sinSquared"))
            {
                source = (1.0/dampingTimeScale_[m]) * 
                         Foam::sqr(Foam::sin(distanceFromBoundary_[m][j])) * 
                         field_[gridCellList_[m][j]];
            }
            else if (dampingFunctionType_[m] == "linear")
            {
            }
            else if ((dampingFunctionType_[m] == "cosine") ||
                     (dampingFunctionType_[m] == "cos"))
            {
            }

            source_[j] = source;
        }
        
    }
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::dampingLayer<Type>::dampingLayer 
(
    const IOdictionary& dict,
    volFieldType& field
)
:
    // Set the pointer to the input dictionary
    dict_(dict),

    // Set the pointer to runTime
    runTime_(field.time()),

    // Set the pointer to the mesh
    mesh_(field.mesh()),

    // Set the pointer to the field being perturbed
    field_(field),

    // Initially, set the height above ground as absolute height.
    zAgl_(mesh_.C() & vector(0,0,1)),

    // Initialize the perturbation source field
    source_
    (
        IOobject
        (
            "dampingSource." & field_.name(),
            runTime_.timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<Type>
        (
            "zero",
            dimensionSet(field_.dimensions()/dimTime),
            Zero
        )
    )


{
    // Initialize the object of the class.
    initialize();
    update();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::dampingLayer<Type>::~dampingLayer()
{}




// ************************************************************************* //
