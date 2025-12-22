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
    blendingFunctionType_.setSize(nLayers_);
    blendingFraction_.setSize(nLayers_);
    useWallDist_.setSize(nLayers_);
    layerThickness_.setSize(nLayers_);
    dampingTimeScale_.setSize(nLayers_);
    boundaryNormal_.setSize(nLayers_);
    boundaryPoint_.setSize(nLayers_);
    gridCellList_.setSize(nLayers_);
    distanceFromBoundary_.setSize(nLayers_);
    dampingStrength_.setSize(nLayers_);
    referenceValue_.setSize(nLayers_);
    dampedComponents_.setSize(nLayers_);

    
    // Call functions required to initialize this damping layers.
    readSubDict();
    inputChecks();
    defineWallDist();
    findBoundaryNormal();
    findDistanceFromBoundary();
    setDampingStrength();
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

        blendingFunctionType_[m] = subSubDict.lookupOrDefault<word>("blendingFunction","none");

        blendingFraction_[m] = subSubDict.lookupOrDefault<scalar>("blendingFraction",1.0);

        referenceValue_[m] = subSubDict.lookupOrDefault<Type>("referenceValue", Zero);

        dampedComponents_[m] = subSubDict.lookupOrDefault<Type>("dampedComponents", Zero);
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
          //Pout << "adjacentBoundary_[" << m << "] = " << adjacentBoundary_[m] << tab << "patchName = " << patchName << endl;
        }
      //Pout << "Patch Number: " << patchNum << endl;
       
        // Check to see how many patch points of the associated lateral boundary
        // that this processor contains.  The localPoints() function gets just the
        // points that this processor uses, as opposed to the points() function that
        // does some parallel reduce.
      //const vectorField& boundaryPointsLocal = mesh_.boundaryMesh()[patchNum].localPoints();
        const vectorField& boundaryPointsLocal = mesh_.Cf().boundaryField()[patchNum];
        label patchSize = boundaryPointsLocal.size();
    
        // If this processor contains at least some of the patch faces, set
        // the hasPatch variable to 1.
        bool hasPatch = (patchSize > 0) ? true : false;
        label procsHavingPatch = (patchSize > 0) ? 1 : 0;
      //Pout << "PatchSize (local) = " << patchSize << tab << "hasPatch = " << hasPatch << endl;

        // Get the normal unit vector of the first patch face.  We assume that the patch is
        // planar so all faces should have the normal of the patch.  Only do this if
        // this processor actually has faces on the patch, but then parallel communicate
        // the result. In OpenFOAM, the boundary normal points out of the domain, so we want
        // the negative of it that points into the domain.
        vector boundaryNormal = vector::zero;
        scalar boundaryArea = 0.0;
        if (hasPatch)
        {
            forAll(mesh_.Sf().boundaryField()[patchNum], i)
            {
                vector boundaryNormal_ = -mesh_.Sf().boundaryField()[patchNum][i];
                scalar boundaryArea_ = mag(boundaryNormal_);
                if (i == 0)
                {
                    boundaryNormal = boundaryNormal_ / mag(boundaryNormal_);
                }
                boundaryArea += boundaryArea_;
            }
          //Pout << "boundaryNormal (before reduce) = " << boundaryNormal << endl;
          //Pout << "boundaryArea (before reduce) = " << boundaryArea << endl;
        }
        reduce(boundaryNormal, sumOp<vector>());
        reduce(boundaryArea, sumOp<scalar>());
        reduce(procsHavingPatch, sumOp<label>());
        boundaryNormal /= scalar(procsHavingPatch);
        boundaryNormal_[m] = boundaryNormal;
      //Pout << "boundaryNormal (after reduce) = " << boundaryNormal << endl;
      //Pout << "boundaryArea (after reduce) = " << boundaryArea << endl;

        // Get a point on the patch. Do this by averaging all points on the patch. 
        // Each processor that has part of the patch will find it's patch area-weighted
        // average location.  Those averages then will be parallel reduced and
        // averaged together to find the final average patch location.
        point boundaryPoint = vector::zero;
        if (hasPatch)
        {
            forAll (boundaryPointsLocal, i)
            {
                scalar faceArea = mag(mesh_.Sf().boundaryField()[patchNum][i]);
                point a = mesh_.Cf().boundaryField()[patchNum][i];
                boundaryPoint += boundaryPointsLocal[i] * faceArea;
              //Pout << "boundaryPoint = " << boundaryPointsLocal[i] << tab << "a = " << a << endl;
            }
        }
        reduce(boundaryPoint, sumOp<vector>());
        boundaryPoint /= boundaryArea;
        boundaryPoint_[m] = boundaryPoint;
      //Pout << "boundaryPoint = " << boundaryPoint << endl;
      //Info << "processors having patch = " << procsHavingPatch << endl;
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
            // Distance to wall is difference between field point and a point on the 
            // boundary dotted with the boundary normal direction
            scalar distance = (meshPoint - boundaryPoint_[m]) & boundaryNormal_[m];

            // If a grid cell is within the damping layer thickness, store its
            // index.
            if ((distance >= 0.0) && 
                (distance <= layerThickness_[m]))
            {
                gridCellList.append(j);
                distanceFromBoundary.append(distance);
            }
        }

        gridCellList_[m] = gridCellList;
        distanceFromBoundary_[m] = distanceFromBoundary;
    }
}





//- Compute the damping strength.
template<class Type>
void Foam::dampingLayer<Type>::setDampingStrength()
{
    // Loop over damping layers
    for (int m = 0; m < nLayers_; m++)
    {
        DynamicList<scalar> strength_;

        // Loop over any cell within the damping zone
        forAll(gridCellList_[m], j)
        {
            scalar strength = 0.0;
            scalar depthFull = (1.0 - blendingFraction_[m]) * layerThickness_[m];
            scalar depthBlended = blendingFraction_[m] * layerThickness_[m];
            scalar distanceFraction = (distanceFromBoundary_[m][j] - depthFull) / max(1.0E-6,depthBlended);

            if (distanceFromBoundary_[m][j] < depthFull)
            {
                strength = 1.0;
            }
            else
            {
                // Sine squared strength function
                if ((blendingFunctionType_[m] == "sineSquared") ||
                    (blendingFunctionType_[m] == "sinSquared"))
                {
                    strength = Foam::sqr(
                                           Foam::sin((Foam::constant::mathematical::pi/2.0) * 
                                           (1.0 - distanceFraction))
                                        );
                }

                // Linear strength function
                else if (blendingFunctionType_[m] == "linear")
                {
                    strength = 1.0 - distanceFraction;
                }

                // Quadratic strength function
                else if ((blendingFunctionType_[m] == "parabolic") ||
                         (blendingFunctionType_[m] == "quadratic"))
                {
                    strength = Foam::sqr(1.0 - distanceFraction);
                }

                // Cosine strength function
                else if ((blendingFunctionType_[m] == "cosine") ||
                         (blendingFunctionType_[m] == "cos"))
                {
                    strength = 0.5 * (1.0 + Foam::cos(Foam::constant::mathematical::pi * distanceFraction));
                }
            }

            strength_.append(strength);

            label cellID = gridCellList_[m][j];
            if (strength > sourceStrength_[cellID])
            {
                sourceStrength_[cellID] = strength;
            }
        }
        dampingStrength_[m] = strength_;
    }
}





//- Update the damping source term.
template<class Type>
void Foam::dampingLayer<Type>::update()
{
    // Zero the damping source term.
    source_ *= 0.0;

    // Loop over damping layers
    for (int m = 0; m < nLayers_; m++)
    {
        // Loop over any cell within the damping zone
        forAll(gridCellList_[m], j)
        {
            label cellID = gridCellList_[m][j];
            Type diff = referenceValue_[m] - field_[cellID];

            // The final source is s = (1/tau) * strength * componentMask * (u_ref - u)
            Type source = (1.0 / dampingTimeScale_[m]) * dampingStrength_[m][j] * cmptMultiply(diff,dampedComponents_[m]);
            
            // Deal with overlapping damping layer sources.  Choose the largest source, component-wise.
            setSource(source,source_[cellID]);
        }      
    }

    // Linearly interpolate the source to processor boundaries.
    source_.correctBoundaryConditions();
}





// Function to compare sources from multiple damping layers where they overlap and choose
// the maximum.  This is done component-wise.
void Foam::dampingLayer<scalar>::setSource(scalar& sourceIn, scalar& sourceOut)
{
    sourceOut = Foam::mag(sourceIn) > Foam::mag(sourceOut) ? sourceIn : sourceOut;
}

void Foam::dampingLayer<vector>::setSource(vector& sourceIn, vector& sourceOut)
{
    sourceOut.x() = Foam::mag(sourceIn.x()) > Foam::mag(sourceOut.x()) ? sourceIn.x() : sourceOut.x();
    sourceOut.y() = Foam::mag(sourceIn.y()) > Foam::mag(sourceOut.y()) ? sourceIn.y() : sourceOut.y();
    sourceOut.z() = Foam::mag(sourceIn.z()) > Foam::mag(sourceOut.z()) ? sourceIn.z() : sourceOut.z();
}

void Foam::dampingLayer<symmTensor>::setSource(symmTensor& sourceIn, symmTensor& sourceOut)
{
    sourceOut.xx() = Foam::mag(sourceIn.xx()) > Foam::mag(sourceOut.xx()) ? sourceIn.xx() : sourceOut.xx();
    sourceOut.xy() = Foam::mag(sourceIn.xy()) > Foam::mag(sourceOut.xy()) ? sourceIn.xy() : sourceOut.xy();
    sourceOut.xz() = Foam::mag(sourceIn.xz()) > Foam::mag(sourceOut.xz()) ? sourceIn.xz() : sourceOut.xz();
    sourceOut.yy() = Foam::mag(sourceIn.yy()) > Foam::mag(sourceOut.yy()) ? sourceIn.yy() : sourceOut.yy();
    sourceOut.yz() = Foam::mag(sourceIn.yz()) > Foam::mag(sourceOut.yz()) ? sourceIn.yz() : sourceOut.yz();
    sourceOut.zz() = Foam::mag(sourceIn.zz()) > Foam::mag(sourceOut.zz()) ? sourceIn.zz() : sourceOut.zz();
}

void Foam::dampingLayer<tensor>::setSource(tensor& sourceIn, tensor& sourceOut)
{
    sourceOut.xx() = Foam::mag(sourceIn.xx()) > Foam::mag(sourceOut.xx()) ? sourceIn.xx() : sourceOut.xx();
    sourceOut.xy() = Foam::mag(sourceIn.xy()) > Foam::mag(sourceOut.xy()) ? sourceIn.xy() : sourceOut.xy();
    sourceOut.xz() = Foam::mag(sourceIn.xz()) > Foam::mag(sourceOut.xz()) ? sourceIn.xz() : sourceOut.xz();
    sourceOut.yx() = Foam::mag(sourceIn.yx()) > Foam::mag(sourceOut.yx()) ? sourceIn.yx() : sourceOut.yx();
    sourceOut.yy() = Foam::mag(sourceIn.yy()) > Foam::mag(sourceOut.yy()) ? sourceIn.yy() : sourceOut.yy();
    sourceOut.yz() = Foam::mag(sourceIn.yz()) > Foam::mag(sourceOut.yz()) ? sourceIn.yz() : sourceOut.yz();
    sourceOut.zx() = Foam::mag(sourceIn.zx()) > Foam::mag(sourceOut.zx()) ? sourceIn.zx() : sourceOut.zx();
    sourceOut.zy() = Foam::mag(sourceIn.zy()) > Foam::mag(sourceOut.zy()) ? sourceIn.zy() : sourceOut.zy();
    sourceOut.zz() = Foam::mag(sourceIn.zz()) > Foam::mag(sourceOut.zz()) ? sourceIn.zz() : sourceOut.zz();
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

    // Initially, set the height above ground as absolute height
    zAgl_(mesh_.C() & vector(0,0,1)),

    // Number of component in the field to be damped, which is
    // also the number of components in the damping force field
    nComponents_(pTraits<Type>::nComponents),

    // Initialize the perturbation source field
    source_
    (
        IOobject
        (
            "dampingSource." & field_.name(),
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensioned<Type>
        (
            "zero",
            dimensionSet(field_.dimensions()/dimTime),
            Zero
        )
    ),

    sourceStrength_
    (
        IOobject
        (
            "DDD." & field_.name(),
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar
        (
            "zero",
            dimensionSet(dimless),
            0.0
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
