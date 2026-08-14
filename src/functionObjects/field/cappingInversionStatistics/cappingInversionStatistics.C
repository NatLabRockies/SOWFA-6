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
#include "zeroGradientFvPatchFields.H"
#include "fvcGrad.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(cappingInversionStatistics, 0);
    addToRunTimeSelectionTable(functionObject, cappingInversionStatistics, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


void Foam::functionObjects::cappingInversionStatistics::writeFileHeader(const label i)
{
    writeHeader(file(), "Stabilization blending factor");
    writeCommented(file(), "Time");
    writeTabbed(file(), "Scheme1");
    writeTabbed(file(), "Scheme2");
    writeTabbed(file(), "Blended");
    file()  << endl;
}


bool Foam::functionObjects::cappingInversionStatistics::calc()
{
  //read(dict_);
   // init(false);
    return true;
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
        scalar(0.0)
    )
{

    // Read the dictionary input for this function object.
    read(dict);

    // Call the initialization function for the first time.
    init(true);

}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::cappingInversionStatistics::read
(
    const dictionary& dict
)
{
    return false;
}


bool Foam::functionObjects::cappingInversionStatistics::execute()
{
    return true;
}


bool Foam::functionObjects::cappingInversionStatistics::write()
{
    if (Pstream::master())
    {
        writeTime(file());
        file()
          //  << tab << nCellsScheme1
          //<< tab << nCellsScheme2
          //<< tab << nCellsBlended
            << endl;
    }

    return true;
}


// ************************************************************************* //
