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

#include "perturbationZone.H"
#include "geometricTransformations.H"
#include "wallDist.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::perturbationZone<Type>::initialize()
{
    t_ = runTime_.value();

    subDict_ = dict_.subOrEmptyDict("flowPerturbationZones." & field_.name());
    subDictList_ = subDict_.toc();

    // Initial messages.
    nZones_ = subDictList_.size();
    if (nZones_ > 0)
    {
        Info << "  -using " << nZones_ << " perturbation zones for field " << field_.name() << "..." << endl;
    }
    else
    {
        Info << "  -no perturbation zones specified for field " << field_.name() << ". Skipping..." << endl;
    }

    // Size the lists.
    locationType_.setSize(nZones_);
    adjacentBoundary_.setSize(nZones_);
    boundaryOffset_.setSize(nZones_);
    layerThickness_.setSize(nZones_);
    layerHeight_.setSize(nZones_);
    useWallDist_.setSize(nZones_);
    boxOrigin_.setSize(nZones_);
    boxVec_i_.setSize(nZones_);
    boxVec_j_.setSize(nZones_);
    boxVec_k_.setSize(nZones_);
    res_.setSize(nZones_);
    dims_.setSize(nZones_);
    points_.setSize(nZones_);
    zoneBox_.setSize(nZones_);
    cellBox_.setSize(nZones_);
    gridCellsInCellBox_.setSize(nZones_);
    fluctuationMagMode_.setSize(nZones_);
    fluctuations_.setSize(nZones_);
    fluctuationMagnitude_.setSize(nZones_);
    EckertNumber_.setSize(nZones_);
    updateMode_.setSize(nZones_);
    updatePeriod_.setSize(nZones_);
    updatePeriodSlab_.setSize(nZones_);
    updatePeriodScalar_.setSize(nZones_);
    PBLHeight_.setSize(nZones_);
    applicationMode_.setSize(nZones_);
    clipAtTwoThirdsPBLHeight_.setSize(nZones_);
    activationTime_.setSize(nZones_);
    deactivationTime_.setSize(nZones_);
    lastUpdateTime_.setSize(nZones_,runTime_.value());
    lastUpdateTimeSlab_.setSize(nZones_);

    // Set the last update time far in the past so that an update will happen 
    // on the first time step.
    for (int m = 0; m < nZones_; m++)
    {
        lastUpdateTime_[m] -= VGREAT;
    }

    readSubDict();
    inputChecks();
    defineWallDist();
    createPerturbationCells();
    identifyGridCellsInCellBox();
}


template<class Type>
void Foam::perturbationZone<Type>::readSubDict()
{
    for (int m = 0; m < nZones_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));
        
        locationType_[m] = subSubDict.lookupOrDefault<word>("locationType","none");
    
        adjacentBoundary_[m] = subSubDict.lookupOrDefault<word>("associatedBoundary","none");

	boundaryOffset_[m] = subSubDict.lookupOrDefault<scalar>("boundaryOffset",0.0);
    
        layerThickness_[m] = subSubDict.lookupOrDefault<scalar>("thickness",0.0);
        layerHeight_[m] = subSubDict.lookupOrDefault<scalar>("height",0.0);
    
        useWallDist_[m] = subSubDict.lookupOrDefault<bool>("useWallDistance",0);
    
        boxOrigin_[m] = subSubDict.lookupOrDefault<vector>("origin", vector::zero);
        boxVec_i_[m] = subSubDict.lookupOrDefault<vector>("i", vector::zero);
        boxVec_j_[m] = subSubDict.lookupOrDefault<vector>("j", vector::zero);
        boxVec_k_[m] = subSubDict.lookupOrDefault<vector>("k", vector::zero);
    
        res_[m] = subSubDict.lookupOrDefault<vector>("resolution", vector::zero);
    
        vector dimsVec = subSubDict.lookupOrDefault<vector>("dimensions", vector::zero);
        List<label> dims;
        for (int i = 0; i < 3; i++)
        {
            dims.append(label(dimsVec[i]));
        }
        dims_[m] = dims;

        fluctuationMagMode_[m] = subSubDict.lookupOrDefault<word>("fluctuationMagnitudeMode","");
    
        fluctuationMagnitude_[m] = subSubDict.lookupOrDefault<Type>("fluctuationScale",Zero);

        EckertNumber_[m] = subSubDict.lookupOrDefault<scalar>("EckertNumber",0.2);

        updateMode_[m] = subSubDict.lookupOrDefault<word>("updateMode","fixedFrequency");

        updatePeriodScalar_[m] = subSubDict.lookupOrDefault<scalar>("updatePeriodScalar",1.0);

        if (subSubDict.found("PBLHeight"))
        {
            PBLHeight_[m] = Function1<scalar>::New("PBLHeight",subSubDict);
        }

        updatePeriod_[m] = subSubDict.lookupOrDefault<scalar>("updatePeriod",10.0);

        applicationMode_[m] = subSubDict.lookupOrDefault<word>("applicationMode","sourceTerm");

        clipAtTwoThirdsPBLHeight_[m] = subSubDict.lookupOrDefault<bool>("clipAtTwoThirdsPBLHeight",false);

        activationTime_[m] = subSubDict.lookupOrDefault<scalar>("activationTime",-1.0E6);

        deactivationTime_[m] = subSubDict.lookupOrDefault<scalar>("deactivationTIme",1.0E6);
    }
}



template<class Type>
void Foam::perturbationZone<Type>::inputChecks()
{
    for (int m = 0; m < nZones_; m++)
    {
        dictionary subSubDict(subDict_.subDict(subDictList_[m]));

        // Check to see that PBLHeight is provided for the "clipAtTwoThirdsPBLHeight"
        // option.
        if (clipAtTwoThirdsPBLHeight_[m])
        {
            if (!subSubDict.found("PBLHeight"))
            {
                FatalErrorInFunction << "Must specify 'PBLHeight' when "
                                     << "'clipAtTwoThirdsPBLHeight' is true."
                                     << abort(FatalError);
            }
        }

        // Check to see that updatePeriod is provided for the "fixedFrequency"
        // updateMode option.
        if (updateMode_[m] == "fixedFrequency")
        {
            if (!subSubDict.found("updatePeriod"))
            {
                FatalErrorInFunction << "Must specify 'updatePeriod' when "
                                     << "'updateMode' is 'fixedFrequency'."
                                     << abort(FatalError);
            }
        }

        // Check to see that PBLHeight is provided for the "windAtPBLHeight"
        // updateMode option.
        if (updateMode_[m] == "windAtPBLHeight")
        {
            if (!subSubDict.found("PBLHeight"))
            {
                FatalErrorInFunction << "Must specify 'PBLHeight' when "
                                     << "'updateMode' is 'windAtPBLHeight'."
                                     << abort(FatalError);
            }
        }
 
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
    
        // Make sure that either dimensions or resolution is given.
        if (!subSubDict.found("dimensions") && !subSubDict.found("resolution"))
        {
            FatalErrorInFunction << "Must specify either perturbation zone dimensions or resolution."
                                 << abort(FatalError);
        }

        // Do checks on fluctuation magnitude type and associated arguments.
        if ((fluctuationMagMode_[m] == "manual") && (!subSubDict.found("fluctuationScale")))
        {
           FatalErrorInFunction << "Must specify fluctuationScale in 'manual' fluctuation magnitude mode."
                                << abort(FatalError);
        }

        if ((fluctuationMagMode_[m] == "Eckert") && (!subSubDict.found("EckertNumber")))
        {
           FatalErrorInFunction << "Must specify EckertNumber in 'Eckert' fluctuation magnitude mode."
                                << abort(FatalError);
        }

        if (activationTime_[m] > deactivationTime_[m])
        {
           FatalErrorInFunction << "Activation time must be earlier than deactivation time."
                                << abort(FatalError);
        }
    }
}




template<class Type>
void Foam::perturbationZone<Type>::defineWallDist()
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
void Foam::perturbationZone<Type>::createPerturbationCells()
{
    for (int m = 0; m < nZones_; m++)
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
    
            // To find the origin in the local coordinate system, transform the boundary points
            // into the patch-local coordinate system, create a global bounding box (i.e., one
            // that is parallel aware and collects points from all processes), and find the minimum
            // of that bounding box (i.e., the "lower-left" point in this box).
            vectorField boundaryPointsLocalP = transformGlobalCartToLocalCart(boundaryPointsLocal,boxVec_i_[m],boxVec_j_[m],boxVec_k_[m]);
            boundBox boundaryBoundsP(boundaryPointsLocalP,true);
            boxOrigin_[m] = boundaryBoundsP.min();

            // If the perturbation zone is offset normal to the boundary, apply that offset to the
            // box origin in the i direction.
            boxOrigin_[m].x() += boundaryOffset_[m];

	    // If over complex terrain, the box origin z coordinate will be the minimum terrain height
	    // on this boundary.  If we want height to be height above ground (terrain conforming
	    // perturbation zone, we need to zero the box origin.
	    if (useWallDist_[m])
            {
                boxOrigin_[m].z() = 0.0;
            }
        }
    
      
    
        // Get the resolution and dimensions of the perturbation zone cells (not the grid
        // cells, but the perturbation cells).  Note, that the user should specify either
        // dimensions or resolution, not both.  If both are defined, the code will throw
        // an error.  If resolution is defined, the nearest resolution that will fit within
        // the perturbation zone will be used.
        xLength = mag(boxVec_i_[m]);
        yLength = mag(boxVec_j_[m]);
        zLength = mag(boxVec_k_[m]);
    
        if (!subSubDict.found("dimensions") && subSubDict.found("resolution"))
        {
            dims_[m][0] = label(round(xLength/res_[m].x()));
            dims_[m][1] = label(round(yLength/res_[m].y()));
            dims_[m][2] = label(round(zLength/res_[m].z()));
        }
        res_[m].x() = xLength/dims_[m][0];
        res_[m].y() = yLength/dims_[m][1];
        res_[m].z() = zLength/dims_[m][2];
    
        
    
        // Build out the list of perturbation cell centers in the perturbation zone local
        // coordinate system.  Ultimately, create the bounding box for each cell.
        label nPoints = dims_[m][0] * dims_[m][1] * dims_[m][2];
    
        List<vector> points(nPoints,vector::zero);
        List<boundBox> cellBox(nPoints);
        List<Type> fluctuations(nPoints,Zero);
    
        vector dby2(0.5*res_[m].x(), 0.5*res_[m].y(), 0.5*res_[m].z());
    
        int ii = 0;
        
        for (int k = 0; k < dims_[m][2]; k++)
        {
            for (int j = 0; j < dims_[m][1]; j++)
            {
                for (int i = 0; i < dims_[m][0]; i++)
                {
                    vector d(i*res_[m].x(), j*res_[m].y(), k*res_[m].z());
                    points[ii] = boxOrigin_[m] + dby2 + d;
    
                    List<vector> vertices(8);
                    vertices[0] = points[ii] + vector(-dby2.x(),-dby2.y(),-dby2.z());
                    vertices[1] = points[ii] + vector( dby2.x(),-dby2.y(),-dby2.z());
                    vertices[2] = points[ii] + vector(-dby2.x(), dby2.y(),-dby2.z());
                    vertices[3] = points[ii] + vector( dby2.x(), dby2.y(),-dby2.z());
                    vertices[4] = points[ii] + vector(-dby2.x(),-dby2.y(), dby2.z());
                    vertices[5] = points[ii] + vector( dby2.x(),-dby2.y(), dby2.z());
                    vertices[6] = points[ii] + vector(-dby2.x(), dby2.y(), dby2.z());
                    vertices[7] = points[ii] + vector( dby2.x(), dby2.y(), dby2.z());
                    cellBox[ii] = boundBox(vertices,false);
    
                    ii++;
                }
            }
        }
    
        points_[m] = points;
    
        cellBox_[m] = cellBox;

        fluctuations_[m] = fluctuations;


        // Build the list of bounding boxes for each zone.
        List<vector> vertices(8);
        vertices[0] =  boxOrigin_[m] + vector(0.0*xLength,0.0*yLength,0.0*zLength); 
        vertices[1] =  boxOrigin_[m] + vector(1.0*xLength,0.0*yLength,0.0*zLength); 
        vertices[2] =  boxOrigin_[m] + vector(1.0*xLength,1.0*yLength,0.0*zLength); 
        vertices[3] =  boxOrigin_[m] + vector(0.0*xLength,1.0*yLength,0.0*zLength); 
        vertices[4] =  boxOrigin_[m] + vector(0.0*xLength,0.0*yLength,1.0*zLength); 
        vertices[5] =  boxOrigin_[m] + vector(1.0*xLength,0.0*yLength,1.0*zLength); 
        vertices[6] =  boxOrigin_[m] + vector(1.0*xLength,1.0*yLength,1.0*zLength); 
        vertices[7] =  boxOrigin_[m] + vector(0.0*xLength,1.0*yLength,1.0*zLength); 
        zoneBox_[m] = boundBox(vertices,false);

        // In case perturbation horizontal slabs are updated independently in time,
        // size the updatePeriodSlab_ variable accordingly.
        updatePeriodSlab_[m].setSize(dims_[m][2],0.0);
        lastUpdateTimeSlab_[m].setSize(dims_[m][2],t_-VGREAT);
    }
}



template<class Type>
void Foam::perturbationZone<Type>::identifyGridCellsInZone(int m, List<label>& gridCellList)
{
    // The perturbations cells are defined in a local coordinate system, so rotate 
    // the mesh cell points into this local coordinate system.  Also, if desired, make 
    // the mesh cell height, height above ground.
    forAll(mesh_.C(),j)
    {
        vector meshPoint = mesh_.C()[j];
        if (useWallDist_[m])
        {
            meshPoint.z() = zAgl_[j];
        }
        meshPoint = transformGlobalCartToLocalCart(meshPoint,boxVec_i_[m],boxVec_j_[m],boxVec_k_[m]);
        if (zoneBox_[m].contains(meshPoint))
        {
            gridCellList.append(j);
        }
    }
}



template<class Type>
void Foam::perturbationZone<Type>::identifyGridCellsInCellBox()
{

    // Find the mesh cells contained within the perturbation cells using the 
    // bounding box functionality.
    for (int m = 0; m < nZones_; m++)
    {
        // To speed things up, start with a coarse search.  Get the collection of all
        // grid cells that lie within this perturbation zone.
        Info << "zoneBox_[" << m << "] = " << zoneBox_[m] << endl;
        List<label> gridCellsInZone;
        identifyGridCellsInZone(m, gridCellsInZone);
        
        // Now do a finer search to find the grid cells contained in each perturbation box.
        List<List<label> > gridCellsInCellBox_zone;
        forAll(points_[m],i)
        {
            List<label> gridCellsInCellBox_cell;
            forAll(gridCellsInZone,j)
            {
                // Get the grid cell center location and replace the z-value with that held in the zAgl
                // variable if height above ground, versus absolute height, is to be used.
                label jj = gridCellsInZone[j];
                vector meshPoint = mesh_.C()[jj];
                if (useWallDist_[m])
                {
                    meshPoint.z() = zAgl_[jj];
                }
                meshPoint = transformGlobalCartToLocalCart(meshPoint,boxVec_i_[m],boxVec_j_[m],boxVec_k_[m]);
                
                if (cellBox_[m][i].contains(meshPoint))
                {
                    gridCellsInCellBox_cell.append(jj);
                }
            }
            gridCellsInCellBox_zone.append(gridCellsInCellBox_cell);
        }
    
        gridCellsInCellBox_[m] = gridCellsInCellBox_zone;
    }
}



template<class Type>
void Foam::perturbationZone<Type>::identifyGridCellsAtHeight(int m, scalar h, List<label>& gridCellList)
{
    // Get all the grid cells in this zone.
    List<label> gridCellsInZone;
    identifyGridCellsInZone(m, gridCellsInZone);

    // Find the grid cells within a certain height band (the height can be absolute or 
    // relative to the ground).  The first guess at a band width that will capture cell
    // centers is some estimate of grid cell height, which we will approximate as the
    // cube-root of the smallest cell volume in the list of cells.  We successively
    // double the bandwidth until grid cells are captured.

    // Set the initial bandwidth.
    scalar bandWidth = great;
    forAll(gridCellsInZone, j)
    {
        label jj = gridCellsInZone[j];
        scalar l = Foam::cbrt(mesh_.V()[jj]);
        if (l < bandWidth)
        {
            bandWidth = l;
        }
    }
    Info << "gridCellsInZone.size() = " << gridCellsInZone.size() << endl;
    Info << "bandWidth = " << bandWidth << endl;
    Info << "h = " << h << endl;

    // Search for grid cells within the height band.
    int i = 0;
    if (gridCellsInZone.size() > 0)
    {
        while (gridCellList.size() == 0)
        {
            i++;
            forAll(gridCellsInZone, j)
            {
                label jj = gridCellsInZone[j];
                vector meshPoint = mesh_.C()[jj];
                if (useWallDist_[m])
                {
                    meshPoint.z() = zAgl_[jj];
                }
                if ((meshPoint.z() < h + 0.5*bandWidth) && (meshPoint.z() > h - 0.5*bandWidth))
                {
                    gridCellList.append(jj);
                }
            }
            bandWidth *= 2.0;
        }
    }
}



template<class Type>
void Foam::perturbationZone<Type>::getVelocityInZone(int m,
                                                     vector& velAvg,
                                                     vector& velMin,
                                                     vector& velMax)
{
  //Info << "In getVelocityInZone..." << endl;

    // Get the current PBL height.
    scalar PBLHeightCurrent = (clipAtTwoThirdsPBLHeight_[m]) ? PBLHeight_[m]->value(t_) : 1.0E15;


    // Do the volume-weighted average and the minimum and maximum magnitude
    // velocities over the identified cells.  Parallel reduce these because the cells
    // at height are usually spread across processors.
    scalar totalVol = 0.0;
    velAvg = vector::zero;
    velMin = great*vector::one;
    velMax = vector::zero;

    // Processor local quantities
    forAll(gridCellsInCellBox_[m],i)
    {
        forAll(gridCellsInCellBox_[m][i],j)
        {
            label jj = gridCellsInCellBox_[m][i][j];

            scalar h = mesh_.C()[jj].z();

            if (h <= (2.0/3.0)*PBLHeightCurrent)
            {
                vector Uh = U_[jj];
                Uh.z() = 0.0;

                velAvg += Uh*mesh_.V()[jj];
                totalVol += mesh_.V()[jj];

                velMin = (Foam::magSqr(Uh) < Foam::magSqr(velMin)) ? Uh : velMin;

                velMax = (Foam::magSqr(Uh) > Foam::magSqr(velMax)) ? Uh : velMax;
            }
        }
    }

    // Parallel reduce.
    reduce(velAvg, sumOp<vector>());
    reduce(totalVol, sumOp<scalar>());
    velAvg /= totalVol;

    reduce(velMin, minMagSqrOp<vector>());

    reduce(velMax, maxMagSqrOp<vector>());

  //Info << "velAvg = " << velAvg << endl;
  //Info << "velMin = " << velMin << endl;
  //Info << "velMax = " << velMax << endl;
}



template<class Type>
void Foam::perturbationZone<Type>::getVelocityAtHeight(int m, 
                                                       scalar h, 
                                                       vector& velAvg, 
                                                       vector& velMin, 
                                                       vector& velMax)
{
  //Info << "In getVelocityAtHeight..." << endl;
    // Get the list of grid cells at the desired height.
    List<label> gridCellList;
    identifyGridCellsAtHeight(m,h,gridCellList);

    // Do the volume-weighted average and the minimum and maximum magnitude
    // velocities over the identified cells.  Parallel reduce these because the cells
    // at height are usually spread across processors.
    scalar totalVol = 0.0;
    velAvg = vector::zero;
    velMin = great*vector::one;
    velMax = vector::zero;

    // Processor local quantities
    forAll(gridCellList,j)
    {
        label jj = gridCellList[j];

        vector Uh = U_[jj];
        Uh.z() = 0.0;

        velAvg += Uh*mesh_.V()[jj];
        totalVol += mesh_.V()[jj];
        
        velMin = (Foam::magSqr(Uh) < Foam::magSqr(velMin)) ? Uh : velMin;  

        velMax = (Foam::magSqr(Uh) > Foam::magSqr(velMax)) ? Uh : velMax;
    }

    // Parallel reduce.
    reduce(velAvg, sumOp<vector>());
    reduce(totalVol, sumOp<scalar>());
    velAvg /= totalVol;

    reduce(velMin, minMagSqrOp<vector>());

    reduce(velMax, maxMagSqrOp<vector>());

  //Info << "height = " << h << endl;
  //Info << "velAvg = " << velAvg << endl;
  //Info << "velMin = " << velMin << endl;
  //Info << "velMax = " << velMax << endl;
}



template<class Type>
void Foam::perturbationZone<Type>::getVelocityOverSlabs(int m, 
                                                        List<vector>& velAvg, 
                                                        List<vector>& velMin, 
                                                        List<vector>& velMax)
{
  //Info << "in getVelocityOverSlabs()..." << endl;

    List<scalar> totalVol(dims_[m][2],0.0);

    // Processor local quantities.
    int ii = 0;
    for (int k = 0; k < dims_[m][2]; k++)
    {
        velAvg.append(vector::zero);
        velMin.append(great*vector::one);
        velMax.append(vector::zero);
        for (int j = 0; j < dims_[m][1]; j++)
        {
            for (int i = 0; i < dims_[m][0]; i++)
            {
                forAll(gridCellsInCellBox_[m][ii],n)
                {
                    label jj = gridCellsInCellBox_[m][ii][n];

                    vector Uh = U_[jj];
                    Uh.z() = 0.0;

                    velAvg[k] += Uh*mesh_.V()[jj];
                    totalVol[k] += mesh_.V()[jj];

                    velMin[k] = (Foam::magSqr(Uh) < Foam::magSqr(velMin[k])) ? Uh : velMin[k];

                    velMax[k] = (Foam::magSqr(Uh) > Foam::magSqr(velMax[k])) ? Uh : velMax[k];
                }
                ii++;
            }
        }
    }

    // Parallel reduce.
    reduce(velAvg, sumOp<List<vector> >());
    reduce(totalVol, sumOp<List<scalar> >());
    for (int k = 0; k < dims_[m][2]; k++)
    {
        velAvg[k] /= totalVol[k];
        reduce(velMin[k], minMagSqrOp<vector>());
        reduce(velMax[k], maxMagSqrOp<vector>());
    }
  //Info << "velAvg = " << velAvg;
  //Info << "velMin = " << velMin;
  //Info << "velMax = " << velMax;
}




template<class Type>
void Foam::perturbationZone<Type>::updateCellFluctuations(int m, int k)
{
    // Populate all perturbation cells with random fluctuation values. The
    // fluctuations are of Type, so they will have fluctuations on all
    // n dimensions of Type.  The perturbations are created with the 
    // built-in random number generator.

    // We may need the PBL height.
    scalar PBLHeightCurrent = (fluctuationMagMode_[m] == "Eckert") ?  PBLHeight_[m]->value(t_) : 1.0E15;
  //Info << "t = " << t_ << tab << "PBLHeightCurrent = " <<  PBLHeightCurrent << endl;

    // If the perturbation magnitude is not directly specified,  update it
    // here.
    if (fluctuationMagMode_[m] == "Eckert")
    {
        // Sample the mean, min, and max velocity within the perturbation slab
        // at 1.1 times the PBL height.
        vector velAvg;
        vector velMin;
        vector velMax;
        getVelocityAtHeight(m,1.1*PBLHeightCurrent,velAvg,velMin,velMax);
        scalar rho = 1.225; 
        scalar cp = 1004.0;
        
        fluctuationMagnitude_[m] = magSqr(velAvg)/(rho*cp*EckertNumber_[m]) * pTraits<Type>::one;

        Info << "Eckert Fluctuation Magnitude = " << fluctuationMagnitude_[m] << endl;
    }

    // Initialize a list of fluctuations to zero for every perturbation 
    // cell.  First test to see if the entire perturbation field is to be
    // updated (k == -1) or if a specific slab is to be updated (k >= 0).
    label nPoints = dims_[m][0] * dims_[m][1];
    nPoints = (k == -1) ? nPoints * dims_[m][2] : nPoints;
    List<Type> fluctuations(nPoints,Zero);

    // To assure that all processors have a consistent fluctuation list,
    // only the master calls the random number generator, and then the
    // result is parallel communicated.  The fluctuation is created by
    // creating a random number from a set uniformly distributed between
    // 0 and 1, multiplying it by a fluctuation magnitude, and then shifting
    // the number down so that the set is centered upon zero.
    if (Pstream::master())
    {
        forAll(fluctuations,ii)
        {
            fluctuations[ii] = cmptMultiply(fluctuationMagnitude_[m],
                                          (randGen_.sample01<Type>() - 0.5*pTraits<Type>::one));
        }
    }

    reduce(fluctuations, sumOp<List<Type> >());  

    if (k == -1)
    {
        fluctuations_[m] = fluctuations;
    }
    else
    {
        int iii = 0;
        for (int ii = k*nPoints; ii < (k+1)*nPoints; ii++)
        {
            scalar h = points_[m][ii].z();
            scalar mask = 1.0;
            if (clipAtTwoThirdsPBLHeight_[m])
            {
                mask = (h <= (2.0*3.0)*PBLHeightCurrent) ? 1.0 : 0.0;
            }
            fluctuations_[m][ii] = fluctuations[iii] * mask;
            iii++;
        }
    }
}



template<class Type>
void Foam::perturbationZone<Type>::updateSourceTerm(int m, int k)
{
    // Figure out which perturbation cells are to be modified.
    label nPoints = dims_[m][0] * dims_[m][1];
    label kStart = (k == -1) ? 0 : k;
    label kEnd = (k == -1) ? dims_[m][2] : k+1;

    // Populate the source term field with values such that if the source term
    // is placed on the RHS of the governing equation, the field will update
    // by the perturbed value over one time step.
    scalar dt = runTime_.deltaT().value();
    for (int i = kStart*nPoints; i < kEnd*nPoints; i++)
    {
        forAll(gridCellsInCellBox_[m][i],j)
        {
            label kk = gridCellsInCellBox_[m][i][j];
            source_[kk] = fluctuations_[m][i] / dt;
        }
    }
    source_.correctBoundaryConditions();
}



template<class Type>
void Foam::perturbationZone<Type>::zeroSourceTerm(int m, int k)
{
    // Figure out which perturbation cells are to be modified.
    label nPoints = dims_[m][0] * dims_[m][1];
    label kStart = (k == -1) ? 0 : k;
    label kEnd = (k == -1) ? dims_[m][2] : k+1;

    // Zero out the source term.
    for (int i = kStart*nPoints; i < kEnd*nPoints; i++)
    {
        forAll(gridCellsInCellBox_[m][i],j)
        {
            label kk = gridCellsInCellBox_[m][i][j];
            source_[kk] = Zero;
        }
    }
    source_.correctBoundaryConditions();
}



template<class Type>
void Foam::perturbationZone<Type>::updatePerturbationField(int m, int k)
{
    // Figure out which perturbation cells are to be modified.
    label nPoints = dims_[m][0] * dims_[m][1];
    label kStart = (k == -1) ? 0 : k;
    label kEnd = (k == -1) ? dims_[m][2] : k+1;

    // Populate the perturbation field that can be directly added to the field.
    for (int i = kStart*nPoints; i < kEnd*nPoints; i++)
    {
        forAll(gridCellsInCellBox_[m][i],j)
        {
            label kk = gridCellsInCellBox_[m][i][j];
            fieldPerturbation_[kk] = fluctuations_[m][i]; 
        }
    }
}



template<class Type>
void Foam::perturbationZone<Type>::zeroPerturbationField(int m, int k)
{
    // Figure out which perturbation cells are to be modified.
    label nPoints = dims_[m][0] * dims_[m][1];
    label kStart = (k == -1) ? 0 : k;
    label kEnd = (k == -1) ? dims_[m][2] : k+1;

    // Zero out the perturbation field.
    for (int i = kStart*nPoints; i < kEnd*nPoints; i++)
    {
        forAll(gridCellsInCellBox_[m][i],j)
        {
            label kk = gridCellsInCellBox_[m][i][j];
            fieldPerturbation_[kk] = Zero;
        }
    }
}



template<class Type>
void Foam::perturbationZone<Type>::update()
{
    // Get current time and time since last perturbation update.
    t_ = runTime_.value();


    // Loop over perturbation zones.
    for (int m = 0; m < nZones_; m++)
    {
        // Only do perturbations if within the time perturbations are active.
        if ((t_ >= activationTime_[m]) && (t_ < deactivationTime_[m]))
        {

            // Get the time since the last perturbation update.
            scalar timeSinceUpdate = t_ - lastUpdateTime_[m];
    
            // For now, the only mode is fixed frequency update.  Each
            // zone can update at its own frequency, though.
            if (updateMode_[m] != "slabLocalWind")
            {
                if (timeSinceUpdate >= updatePeriod_[m])
                {
                    if (updateMode_[m] != "fixedFrequency")
                    {
                        vector velAvg;
                        vector velMin; 
                        vector velMax;
                        vector vel;
                        if (updateMode_[m] == "maxWind")
                        {
                            getVelocityInZone(m,velAvg,velMin,velMax);
                            vel = velMax;
                        }
                        else if (updateMode_[m] == "windAtPBLHeight")
                        {
                            scalar PBLHeightCurrent = PBLHeight_[m]->value(t_);
                            getVelocityAtHeight(m,PBLHeightCurrent,velAvg,velMin,velMax);
                            vel = velAvg;
                        }
    
                        // Get the flow-through time across the width of the perturbation zone. This becomes the update time.
                        vector d = boxVec_i_[m];
                        scalar mag_d = max(mag(d),1.0E-6);
                        vector n = d/mag_d;
    
                        updatePeriod_[m] = updatePeriodScalar_[m] * (mag_d/(max(mag(vel & n),1.0E-6)*sign(vel & n)));
    
                      //Info << "vel = " << vel << tab << "d = " << d << tab << "mag(d) = " << mag_d << tab << "(vel & n) = " << (vel & n) << endl;
                      //Info << "updatePeriod = " << updatePeriod_[m] << endl;
                    }
    
                    if (updatePeriod_[m] >= 0.0)
                    {
                        updateCellFluctuations(m);
                    
                        // Update the perturbation field if the field is to be directly
                        // updated; otherwise leave it zero.
                        if (applicationMode_[m] == "direct")
                        {
                            updatePerturbationField(m);
                        }
    
                        // If indirectly perturbing the field through source terms, update
                        // the source term.
                        else if (applicationMode_[m] == "sourceTerm")
                        {
                            (timeSinceUpdate > GREAT) ? updatePerturbationField(m) : updateSourceTerm(m);
                        }
                    }
    
                    // Mark this as the time of last update of perturbations.
                    updatePeriod_[m] = mag(updatePeriod_[m]);
                    lastUpdateTime_[m] = t_;
                       
                }
                else
                {
                    zeroSourceTerm(m);
                    zeroPerturbationField(m);
                }
                Info << updatePeriod_[m] << endl;
            }
            else if (updateMode_[m] == "slabLocalWind")
            {
                // Get the slab local velocities
                List<vector> velAvg;
                List<vector> velMin;
                List<vector> velMax;
                getVelocityOverSlabs(m,velAvg,velMin,velMax);
    
                for (int k = 0; k < dims_[m][2]; k++)
                {
                    // Get the time since the last perturbation update.
                    scalar timeSinceUpdate = t_ - lastUpdateTimeSlab_[m][k];
    
                    if (timeSinceUpdate >= updatePeriodSlab_[m][k])
                    {
                        // Get the flow-through time across the width of the perturbation zone. This becomes the update time.
                        vector d = boxVec_i_[m];
                        scalar mag_d = max(mag(d),1.0E-6);
                        vector n = d/mag_d;
    
                        vector vel = velAvg[k];
    
                        updatePeriodSlab_[m][k] = updatePeriodScalar_[m] * (mag_d/(max(mag(vel & n),1.0E-6)*sign(vel & n)));
    
                      //Info << "vel = " << vel << tab << "d = " << d << tab << "mag(d) = " << mag_d << tab << "(vel & n) = " << (vel & n) << endl;
                      //Info << "updatePeriod = " << updatePeriodSlab_[m][k] << endl;       
    
                        if (updatePeriodSlab_[m][k] >= 0.0)
                        {
                            updateCellFluctuations(m,k);
    
                            // Update the perturbation field if the field is to be directly
                            // updated; otherwise leave it zero.
                            if (applicationMode_[m] == "direct")
                            {
                                updatePerturbationField(m,k);
                            }
    
                            // If indirectly perturbing the field through source terms, update
                            // the source term.
                            else if (applicationMode_[m] == "sourceTerm")
                            {
                                (timeSinceUpdate > GREAT) ? updatePerturbationField(m,k) : updateSourceTerm(m,k);
                            }
                        }
    
                        // Mark this as the time of last update of perturbations.
                        updatePeriodSlab_[m][k] = mag(updatePeriodSlab_[m][k]);
                        lastUpdateTimeSlab_[m][k] = t_;
                    }
                    else
                    {
                        zeroSourceTerm(m,k);
                        zeroPerturbationField(m,k);
                    }
                }
              //Info << updatePeriodSlab_[m] << endl;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::perturbationZone<Type>::perturbationZone 
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

    // Set the pointer to the velocity field
    U_(field.db().objectRegistry::template lookupObject<volVectorField>("U")),

    // Set the heigh above ground as first absolute height.
    zAgl_(mesh_.C() & vector(0,0,1)),

    // Initialize the random number generator.
    randGen_(label(0)),

    // Initialize the perturbation source field
    source_
    (
        IOobject
        (
            "perturbationSource." & field_.name(),
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
    ),

    // Initialize the fluctuation field to zero.
    fieldPerturbation_
    (
        IOobject
        (
            "perturbationField." & field_.name(),
            runTime_.timeName(),
            mesh_
        ),
        mesh_,
        dimensioned<Type>
        (
            "zero",
            dimensionSet(field_.dimensions()),
            Zero
        )
    )

{
   // Initialize the object of this class and perform the first perturbation
   // cell update and application to the field.
   initialize();
   update();
}




// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::perturbationZone<Type>::~perturbationZone()
{}




// ************************************************************************* //
