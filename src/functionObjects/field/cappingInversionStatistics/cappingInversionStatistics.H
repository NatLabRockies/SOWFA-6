/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2020 OpenCFD Ltd.
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
    Foam::functionObjects::cappingInversionStatistics

Group
    grpFieldFunctionObjects

Description
    Computes the capping inversion height.



See also
    - Foam::functionObject
    - Foam::functionObjects::fvMeshFunctionObject
    - Foam::functionObjects::logFiles

SourceFiles
    cappingInversionStatistics.C

\*---------------------------------------------------------------------------*/

#ifndef cappingInversionStatistics_H
#define cappingInversionStatistics_H

#include "fvMeshFunctionObject.H"
#include "logFiles.H"
#include "volFields.H"
#include "uniformDimensionedFields.H"
#include "meshPencils.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                       Class cappingInversionStatistics Declaration
\*---------------------------------------------------------------------------*/

class cappingInversionStatistics
:
    public fvMeshFunctionObject,
    public logFiles
{
    // Private Member Data

        //- Mesh pencils
        meshPencils zPencils_;

        //- Capping inversion height
        uniformDimensionedScalarField zi_;

        



    // Private Member Functions

      cappingInversionStatistics(const cappingInversionStatistics&);

      void operator=(const cappingInversionStatistics&);


protected:

    // Protected Member Functions

        //- Write the file header
        virtual void writeFileHeader(const label i);


public:

    //- Runtime type information
    TypeName("cappingInversionStatistics");


    // Constructors

        //- Construct from Time and dictionary
        cappingInversionStatistics
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cappingInversionStatistics();


    // Member Functions

        //- Read the cappingInversionStatistics data
        virtual bool read(const dictionary&);

        //- Compute the capping inversion height.
        virtual bool execute();

        //- Write the cappingInversionStatistics
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
