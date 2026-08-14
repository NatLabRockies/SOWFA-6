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
    Foam::meshPencils

Description
    Data structure containing a set of mesh "pencils" in a given direction.
    A "pencil" is a line or column of grid cells.

SourceFiles
    meshPencils.

\*---------------------------------------------------------------------------*/

#ifndef meshPencils_H
#define meshPencils_H

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

/*---------------------------------------------------------------------------*\
                           Class meshPencils Declaration
\*---------------------------------------------------------------------------*/

class meshPencils
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Lower and upper limit of along-pencil coordinate to consider
        //  within the domain.
        const scalar extentLowerLim_;
        const scalar extentUpperLim_;

        //- Index corresponding to the direction along the pencil.
        label directionIndex_;

        //- Index pair corresponding to the two directions not along the pencil.
        Pair<label> orthogonalDirectionIndices_;

        //- Total number of pencils
        label numberOfPencils_;

        //- Total volume per pencil (global value)
        List<scalar> totVolPerPencil_;

        //- List of lists of cell ID labels corresponding to a specific pencil on a processor
        List<List<label> > pencilCellList_;

        //- Number of cells per pencil on a processor
        labelList numCellPerPencil_;

        //- List of plane locations
        List<Pair<scalar>> pencilLocations_;

        //- Initialization flag
        bool initialized_;


    // Private Member Functions

        //- Find cells per pencil 
        void findPencils_();


public:

    //- Declare name of the class and its debug switch
    ClassName("meshPencils");


    // Static member data

        //- tolerance on plane locations
        static scalar tol_;


    // Constructors

        //- Construct from components
        meshPencils
        (
            const fvMesh& mesh,
            const word& direction,
            const scalar extentLowerLim,
            const scalar extentUpperLim
        );

        //- Construct default (z-pencils)
        meshPencils
        (
            const fvMesh& mesh
        );


    // Destructor
    virtual ~meshPencils();


    // Member functions

        //- Initialize
        void initialize();

        //- Find the location of the maximum of a field (on all or a
        //  particular pencil).
        template<class Type>
        List<scalar> maxLocation
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        );

        template<class Type>
        scalar maxLocation
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField,
            label pencilI
        );

        //- Find the value of the maximum of a field (on all or a
        //  particular pencil).
        template<class Type>
        List<Type> maxValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        );

        template<class Type>
        Type maxValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField,
            label pencilI
        );

        //- Find the location of the minimum of a field (on all or a
        //  particular pencil).
        template<class Type>
        List<scalar> minLocation
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        );

        template<class Type>
        scalar minLocation
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField,
            label pencilI
        );

        //- Find the value of the minimum of a field (on all or a
        //  particular pencil).
        template<class Type>
        List<Type> minValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        );

        template<class Type>
        Type minValue
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField,
            label pencilI
        );

        //- Compute pencil average (on all or a particular pencil).
        template<class Type>
        List<Type> average
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField
        );

        template<class Type>
        Type average
        (
            const GeometricField<Type, fvPatchField, volMesh>& vField,
            label pencilI
        );

        //- Access
        
            label numberOfPencils() const
            {
                return numberOfPencils_;
            }
    
            List<scalar>& totVolPerPencil()
            {
                return totVolPerPencil_;
            }
    
            List<List<label> >& pencilCellList()
            {
                return pencilCellList_;
            }
    
            labelList& numCellPerPencil()
            {
                return numCellPerPencil_;
            }
    
            List<Pair<scalar> >& pencilLocations()
            {
                return pencilLocations_;
            }


};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Loop across all pencils in \a meshPencils
// \par Usage
// \code
// forAllPencils(meshPencils, i)
// {
//      statements;
// }
// \endcode
// \sa forAllPencils
#define forAllPencils(meshPencils, i) \
    for (Foam::label i=0; i<(meshPencils).numberOfPencils(); i++)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
