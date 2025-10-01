/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

Class
    Foam::dampingLayer<Type>

Description
    Add general class to add a damping/relaxation source, for example a 
    Rayleigh damping, to the flow field.  This is a template class, so it can
    be applied to any field.

Usage

SourceFiles
    dampingLayer.C

\*---------------------------------------------------------------------------*/

#ifndef dampingLayer_H
#define dampingLayer_H

#include "scalar.H"
#include "fvCFD.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class dampingLayer Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class dampingLayer
{
    // private data:

        typedef GeometricField<Type, fvPatchField, volMesh> volFieldType;

        //- Constants
            const IOdictionary& dict_;

            //- Runtime pointer
            const Time& runTime_;

            //- Mesh pointer
            const fvMesh& mesh_;

            //- Field pointer (field cannot be modified)
            const volFieldType& field_;


        //- Damping source field
        volFieldType source_;

        //- Damping subdictionary within the overall dictionary.
        dictionary subDict_;

        //- List of the names of the damping layers.
        wordList subDictList_;

        //- Number of damping layers.
        label nLayers_;

        //- Current time.
        scalar t_;

        //- Field for storing height above the ground
        volScalarField zAgl_;

        //- Defines how the layer's location is defined, choices are:
        //     -boundary: a layer adjacent to and touching a given boundary.
        //     -arbitrary: a layer that is stand-alone in the domain.
        List<word> locationType_;

        //- For zones where the locationType is "boundary", then use these 
        //  location definitions:

            //-- name of boundary to which the zone will be adjacent.
            List<word> adjacentBoundary_;

            //-- use height relative to the ground instead of absolute height
            List<bool> useWallDist_;
            bool useWallDistAny_;

        //- For zones where the locationType is "arbitraryBox", then use 
        //  these location definitions:

            //-- origin of the layer
            List<vector> layerOrigin_;

            //-- vectors defining the arbitrary layer 
            List<vector> layerVec_i_;
            List<vector> layerVec_j_;
            List<vector> layerVec_k_;




    // Private Member Functions

        //- Initialize this class
        void initialize();

        //- Read the sub-dictionaries
        void readSubDict();

        //- Check for necessary inputs.  If they are not there, throw an error.
        void inputChecks();

        //- Check to see if wall distance will be used.  If so, compute it.
        void defineWallDist();




public:

    //- Runtime type information
    ClassName("dampingLayer");


    // Constructors
    dampingLayer 
    (
        const IOdictionary& dict,
        volFieldType& field
    );


    // Destructor
    virtual ~dampingLayer();

    // Overall update function that applications can call
    void update();

    // Return the source term field.
    volFieldType& source()    
    {
        return source_;
    }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "dampingLayer.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
