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
    locationType_.setSize(nLayers_);
    adjacentBoundary_.setSize(nLayers_);
    layerOrigin_.setSize(nLayers_);
    layerVec_i_.setSize(nLayers_);
    layerVec_j_.setSize(nLayers_);
    layerVec_k_.setSize(nLayers_);

    
    // Call functions required to initialize this damping layers.
    readSubDict();
    inputChecks();
    defineWallDist();
}



template<class Type>
void Foam::dampingLayer<Type>::readSubDict()
{
    for (int m = 0; m < nLayers_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));
        
        locationType_[m] = subSubDict.lookupOrDefault<word>("locationType","none");
    
        adjacentBoundary_[m] = subSubDict.lookupOrDefault<word>("associatedBoundary","none");
    
        useWallDist_[m] = subSubDict.lookupOrDefault<bool>("useWallDistance",0);
    
        layerOrigin_[m] = subSubDict.lookupOrDefault<vector>("origin", vector::zero);
        layerVec_i_[m] = subSubDict.lookupOrDefault<vector>("i", vector::zero);
        layerVec_j_[m] = subSubDict.lookupOrDefault<vector>("j", vector::zero);
        layerVec_k_[m] = subSubDict.lookupOrDefault<vector>("k", vector::zero);
    }
}



template<class Type>
void Foam::dampingLayer<Type>::inputChecks()
{
    for (int m = 0; m < nLayers_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));

        // Check to see that boundary is provided for "lateralBoundary" option.  
        if (locationType_[m] == "lateralBoundary")
        {
            if (!subSubDict.found("associatedBoundary"))
            {
                FatalErrorInFunction << "Must specify  'associatedBoundary' when "
                                     << "'locationType' is 'adjacentBoundary'."
                                     << abort(FatalError);
            }
        }
       
        // Check to see that origin, i-, j-, and k-vectors are provided for
        // "arbitraryBox" option.
        else if (locationType_[m] == "arbitraryBox")
        {
            if (!subSubDict.found("origin") ||
                !subSubDict.found("i") || 
                !subSubDict.found("j") || 
                !subSubDict.found("k")) 
            {
                FatalErrorInFunction << "Must specify  'origin', 'i', 'j', "
                                     << "'k' when 'locationType' is 'arbitraryBox'."
                                     << abort(FatalError);
            }
        }
      
        // If neither option is selected, throw an error.
        else
        {
            FatalErrorInFunction << "Must specify that 'locationType' is either "
                                 << "'lateralBoundary' or 'arbitraryBox'."
                                 << abort(FatalError);
        }
    }
}



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



template<class Type>
void Foam::dampingLayer<Type>::setDampingStrength()
{
    for (int m = 0; m < nLayers_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));
    
        // First, we need to get the origin, orientation, and size of the perturbation
        // zone.  These parameters are given for the "arbitraryBox" specification, but
        // if the zone is adjacent to a specified boundary, we need to calculate these
        // parameters.
        scalar xLength = 0.0;
        scalar yLength = 0.0;
        scalar zLength = 0.0;
    
        // if the perturbation zone is associated with a lateral boundary, it 
        // means that it is a zone adjacent to and touching the given boundary.  
        // This assumes that the boundary is lateral (not an upper or lower) and
        // a planar surface, but that plane need not be Cartesian oriented.  The
        // zone will also conform to complex terrain.
        if (locationType_[m] == "lateralBoundary")
        {
    
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
    
            // This is the meat of this part of the code.  Here, we set the cell perturbation zone's
            // defining i-, j-, k-vectors and origin in the coordinate system of the domain.  We assume
            // that the k-vector is up, and the i-vector is normal to the boundary and pointing inward
            // to the domain.  The j-vector completes this with the right-hand rule.  The vector
            // magnitudes are the x-, y-, and z-lengths of the cell perturbation zone.  Because not
            // all processors can see the patch, nor can each processor see the full patch, some 
            // parallel reduces happen in finding the bounding box and boundary normal.
    
            // Get the bounding box of the patch.  With the second argument of the boundaryBox
            // constructor set to "true", a parallel reduce happens to get the bouding box
            // of the entire patch, not just this processor's portion of the patch.
            boundBox boundaryBounds(boundaryPointsLocal,true);
    
            // The dimensions of the perturbation zone are given for depth normal to the 
            // lateral boundary (xLength), and the height (zLength), but the width of
            // the zone is the width of the boundary patch which is computed from the
            // bounding box width.
            xLength = layerThickness_[m];
            yLength = sqrt(sqr(boundaryBounds.max().x() - boundaryBounds.min().x()) + 
                           sqr(boundaryBounds.max().y() - boundaryBounds.min().y()));
            zLength = layerHeight_[m];
    
            // Get the normal unit vector of the first patch face.  We assume that the patch is
            // planar so all faces should have the normal of the patch.  Only do this if
            // this processor actually has faces on the patch, but then parallel 
            // communicate the result.
            vector boundaryNormal = vector::zero;
            if (hasPatch == 1)
            {
                boundaryNormal = mesh_.Sf().boundaryField()[patchNum][0];
                boundaryNormal /= mag(boundaryNormal);
            }
            reduce(boundaryNormal, sumOp<vector>());
            reduce(hasPatch, sumOp<label>());
            boundaryNormal /= scalar(hasPatch);
    
            // Define the i-, j-, and k-vectors of the cell perturbation zone.  The i-vector
            // is the opposite of the patch normal (the patch normal points out; the i-
            // vector points in.  The k-vector is up.  The j-vector is orthogonal to the
            // others following the right-hand rule.
            boxVec_i_[m] = -boundaryNormal;
            boxVec_k_[m] =  vector(0.0,0.0,1.0);
            boxVec_j_[m] =  boxVec_k_[m] ^ boxVec_i_[m];
    
            boxVec_i_[m] *= xLength;
            boxVec_j_[m] *= yLength;
            boxVec_k_[m] *= zLength;
    
            // To find the origin in the local coordinate system.
            vectorField boundaryPointsLocalP = transformGlobalCartToLocalCart(boundaryPointsLocal,boxVec_i_[m],boxVec_j_[m],boxVec_k_[m]);
            boundBox boundaryBoundsP(boundaryPointsLocalP,true);
            boxOrigin_[m] = boundaryBoundsP.min();
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
