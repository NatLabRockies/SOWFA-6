/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "cappingInversionStatistics.H"
#include "addToRunTimeSelectionTable.H"
#include "uniformDimensionedFields.H"
//#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cappingInversionStatistics, 0);

    addToRunTimeSelectionTable
    (
        functionObject, 
        cappingInversionStatistics, 
        dictionary
    );
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::cappingInversionStatistics::writeFileHeader(const label i)
{
    // Add headers to output data
    writeHeader(file(), "Capping inversion height");
    writeCommented(file(), "Time (s)");
    writeTabbed(file(), "zi");
    file() << endl;
}




// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::cappingInversionStatistics::cappingInversionStatistics
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    zPencils_(mesh_),
    zi_
    (
        IOobject
        (
            "zi",
            runTime.timeName(),
            obr_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        scalar(750.0)
    )


{
    // Read the dictionary input for this function object.
    read(dict);

    resetName(typeName);    

    
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::cappingInversionStatistics::~cappingInversionStatistics()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool Foam::functionObjects::cappingInversionStatistics::read(const dictionary& dict)
{
    fvMeshFunctionObject::read(dict);

    return false;
}



bool Foam::functionObjects::cappingInversionStatistics::execute()
{

    // Get the temperature field
    bool foundT = obr_.foundObject<volScalarField>("T");
    if (!foundT)
    {
        FatalErrorInFunction
            << "Unable to find temperature field (T) in the database. "
            << exit(FatalError);
    }
    const volScalarField& T = obr_.lookupObject<volScalarField>("T");
    
    // Compute the gradient of the temperature field.
    volVectorField Tgrad = fvc::grad(T);


    return true;
}



bool Foam::functionObjects::cappingInversionStatistics::write()
{
    Log << type() << " " << name() << "write: " << nl;

    logFiles:: write();

    if (Pstream::master())
    {
        file()
            << mesh_.time().value()
            << tab << zi_.value()
            << endl;
    }

    Log << "Capping inversion height: " << endl;

    return true;
}


// ************************************************************************* //
