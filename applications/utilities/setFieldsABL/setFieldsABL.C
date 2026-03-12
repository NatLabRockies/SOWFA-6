/*---------------------------------------------------------------------------*\
This file was modified or created at the National Renewable Energy
Laboratory (NREL) on January 6, 2012 in creating the SOWFA (Simulator for
Offshore Wind Farm Applications) package of wind plant modeling tools that
are based on the OpenFOAM software. Access to and use of SOWFA imposes
obligations on the user, as set forth in the NWTC Design Codes DATA USE
DISCLAIMER AGREEMENT that can be found at
<http://wind.nrel.gov/designcodes/disclaimer.html>.
\*---------------------------------------------------------------------------*/

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

Application
    setFieldsABL

Description
    Initializes the flow field for turbulent atmospheric boundary layer LES.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "singlePhaseTransportModel.H"
#include "interpolateXY.H"
#include "interpolateSplineXY.H"
#include "Random.H"
#include "wallDist.H"
#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addDictOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// Read in the existing solution files.   
Info << "Reading field U" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volVectorField U0
(
    IOobject
    (
        "U_0",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("U_0",dimVelocity,Zero)
);

Info << "Reading field T" << endl;
volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);

volScalarField T0
(
    IOobject
    (
        "T_0",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("T_0",dimTemperature,Zero)
);

Info << "Reading field p_rgh" << endl;
volScalarField p_rgh
(
    IOobject
    (
        "p_rgh",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


// Compute the velocity flux at the faces.  This is needed
// by the laminar transport model.
Info<< "Creating/Calculating face flux field, phi..." << endl;
surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);


// Read the gravitational acceleration.  This is needed 
// for calculating dp/dn on boundaries.
Info << "Reading gravitational acceleration..." << endl;
uniformDimensionedVectorField g
(
    IOobject
    (
        "g",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);


// Read the value of TRef in the transportProperties file.
singlePhaseTransportModel laminarTransport(U, phi);
dimensionedScalar TRef(laminarTransport.lookup("TRef"));


// Use Tref and the T field to compute rhok, which is needed
// to calculate dp/dn on boundaries.
Info<< "Creating the kinematic density field, rhok..." << endl;
volScalarField rhok
(
    IOobject
    (
        "rhok",
        runTime.timeName(),
        mesh
    ),
    1.0 - (T - TRef)/TRef
);


// Get access to the input dictionary.
const word dictName("setFieldsABLDict");
#include "setSystemMeshDictionaryIO.H"
Info << "Reading " << dictPath << endl;
IOdictionary dict(dictIO);
/*
IOdictionary dict
(
    IOobject
    (
        "dict",
        runTime.time().system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);
*/



// Read in the dict entries.
word velocityInitType(dict.lookup("velocityInitType"));
word temperatureInitType(dict.lookup("temperatureInitType"));
word tableInterpTypeU(dict.lookupOrDefault<word>("tableInterpTypeU","linear"));
word tableInterpTypeT(dict.lookupOrDefault<word>("tableInterpTypeT","linear"));
scalar deltaU(dict.lookupOrDefault<scalar>("deltaU",1.0));
scalar deltaV(dict.lookupOrDefault<scalar>("deltaV",1.0));
scalar zPeak(dict.lookupOrDefault<scalar>("zPeak",0.03));
scalar Uperiods(round(dict.lookupOrDefault<scalar>("Uperiods",4)));
scalar Vperiods(round(dict.lookupOrDefault<scalar>("Vperiods",4)));
scalar xMin(dict.lookupOrDefault<scalar>("xMin",0.0));
scalar yMin(dict.lookupOrDefault<scalar>("yMin",0.0));
scalar zMin(dict.lookupOrDefault<scalar>("zMin",0.0));
scalar xMax(dict.lookupOrDefault<scalar>("xMax",3000.0));
scalar yMax(dict.lookupOrDefault<scalar>("yMax",3000.0));
scalar zMax(dict.lookupOrDefault<scalar>("zMax",1000.0));
scalar zRef(dict.lookupOrDefault<scalar>("zRef",600.0));
bool useWallDistZ(dict.lookupOrDefault<bool>("useWallDistZ",false));
bool scaleVelocityWithHeight(dict.lookupOrDefault<bool>("scaleVelocityWithHeight",false));
scalar zInversion(dict.lookupOrDefault<scalar>("zInversion",600.0));
scalar Ug(dict.lookupOrDefault<scalar>("Ug",15.0));
scalar UgDir(dict.lookupOrDefault<scalar>("UgDir",270.0));
scalar Tbottom(dict.lookupOrDefault<scalar>("Tbottom",300.0));
scalar Ttop(dict.lookupOrDefault<scalar>("Ttop",304.0));
scalar dTdz(dict.lookupOrDefault<scalar>("dTdz",0.003));
scalar widthInversion(dict.lookupOrDefault<scalar>("widthInversion",80.0));
scalar TPrimeScale(dict.lookupOrDefault<scalar>("TPrimeScale",0.0));
scalar z0(dict.lookupOrDefault<scalar>("z0",0.016));
scalar kappa(dict.lookupOrDefault<scalar>("kappa",0.40));
List<List<scalar> > profileTable(dict.lookup("profileTable"));
bool updateInternalFields(dict.lookupOrDefault<bool>("updateInternalFields",true));
bool updateBoundaryFields(dict.lookupOrDefault<bool>("updateBoundaryFields",true));

// If velocity is only to be scaled, set wallDist and scaleVelocityWithHeight to true
if ( velocityInitType == "scaleGiven" )
{
    useWallDistZ = true;
    scaleVelocityWithHeight = true;
    temperatureInitType = "none";
}

// Change the table profiles from scalar lists to scalar fields
scalarField zProfile(profileTable.size(),0.0);
scalarField UProfile(profileTable.size(),0.0);
scalarField VProfile(profileTable.size(),0.0);
scalarField TProfile(profileTable.size(),0.0);
forAll(zProfile,i)
{
   zProfile[i] = profileTable[i][0];
   UProfile[i] = profileTable[i][1];
   VProfile[i] = profileTable[i][2];
   TProfile[i] = profileTable[i][3];
}

// Get distance from the wall only if required.
vector up = vector::zero;
up.z() = 1.0;
volScalarField d(mesh.C() & up);    
if (useWallDistZ)
{
    Info << "Calculating wall distance..." << endl;
    wallDist dWall(mesh);
    d = dWall.y();
}

// Now calculate the field quantities.
scalar uStar = (kappa*Ug)/(Foam::log(zInversion/z0));
Info << "u* = " << uStar << " m/s" << endl;

// Calculate the wind vector direction.
Info << "Calculating wind vector..." << endl;
if (UgDir > 180.0)
{
    UgDir = UgDir - 180.0;
}
else
{
    UgDir = UgDir + 180.0;
}
UgDir = 90.0 - UgDir;
if (UgDir < 0.0)
{
    UgDir = UgDir + 360.0;
}
UgDir = UgDir * ((Foam::constant::mathematical::pi)/180.0);

// Calculate the wind vector.
vector UgToVector;
UgToVector.x() = Foam::cos(UgDir);
UgToVector.y() = Foam::sin(UgDir);
UgToVector.z() = 0.0;
vector UgVec = Ug * UgToVector;

// Compute the domain extents.
scalar xExtent = xMax - xMin;
scalar yExtent = yMax - yMin;
scalar zExtent = zMax - zMin;

// Update the interior fields.
if (updateInternalFields)
{
    // Velocity.
    Info << "Updating internal U field..." << endl;
    forAll(U,cellI)
    {
        scalar x = mesh.C()[cellI].x() - xMin;
        scalar y = mesh.C()[cellI].y() - yMin;
        scalar z = 0.0;
        scalar zAbsolute = 0.0;
        scalar zSurface = 0.0;
        scalar zColumn = 0.0;
        scalar velScalar = 1.0;
        if (useWallDistZ)
        {
            z = d[cellI];
            zAbsolute = mesh.C()[cellI].z();
            zSurface = zAbsolute - z;
            zColumn = zMax - zSurface;
            if (scaleVelocityWithHeight)
            {
                velScalar = zRef/zColumn;
            }
        }
        else
        {
            z = mesh.C()[cellI].z() - zMin;
        }
        vector UPrime = vector::zero;
        UPrime.x() = deltaU * Foam::exp(0.5) * Foam::cos(Uperiods * 2.0 * Foam::constant::mathematical::pi * y/yExtent) * 
                    (z/(zPeak*zExtent)) * Foam::exp(-0.5*Foam::pow((z/(zPeak*zExtent)),2));
        UPrime.y() = deltaV * Foam::exp(0.5) * Foam::sin(Vperiods * 2.0 * Foam::constant::mathematical::pi * x/xExtent) * 
                    (z/(zPeak*zExtent)) * Foam::exp(-0.5*Foam::pow((z/(zPeak*zExtent)),2));
        UPrime.z() = 0.0;

        if ((z <= zInversion) && (velocityInitType == "log"))
        {
            U[cellI] = velScalar*(uStar/kappa)*Foam::log(z/z0)*UgToVector;

        }
        else if (((z > zInversion) && (velocityInitType == "log")) || (velocityInitType == "geostrophic"))
        {
            U[cellI] = velScalar*UgVec;
        }
        else if (velocityInitType == "table")
        {
            if (tableInterpTypeU == "cubic")
            {
                U[cellI].x() = velScalar*interpolateSplineXY(z,zProfile,UProfile);
                U[cellI].y() = velScalar*interpolateSplineXY(z,zProfile,VProfile);
            }
            else
            {
                U[cellI].x() = velScalar*interpolateXY(z,zProfile,UProfile);
                U[cellI].y() = velScalar*interpolateXY(z,zProfile,VProfile);
            }
        }
        else if (velocityInitType == "scaleGiven")
        {
            U[cellI].x() = velScalar*U[cellI].x();
            U[cellI].y() = velScalar*U[cellI].y();
        }

        if( velocityInitType != "scaleGiven")
        {
            U[cellI] += UPrime;
        }
    }

    // Potential temperature.
    Info << "Updating internal T field..." << endl;
    Random TRandom(label(0));
    if ( temperatureInitType != "none" )
    {
        forAll(T,cellI)
        {
            scalar TPrime = TPrimeScale * (TRandom.sample01<scalar>() - 0.5);
            scalar z = 0.0;
            if (useWallDistZ)
            {
                z = d[cellI];
            }
            else
            {
                z = mesh.C()[cellI].z() - zMin;
            }

            T[cellI] = Tbottom;
            if ((z >= zInversion - 0.5*widthInversion) && (z <= zInversion + 0.5*widthInversion) && (temperatureInitType == "simple"))
            {
                T[cellI] = Tbottom + ((Ttop - Tbottom)/Foam::max(widthInversion,1.0E-10)) * (z -(zInversion - 0.5*widthInversion));
            }
            else if ((z > zInversion + 0.5 * widthInversion) && (temperatureInitType == "simple"))
            {
                T[cellI] = Ttop + dTdz * (z - (zInversion + 0.5*widthInversion));
            }
            else if (temperatureInitType == "table")
            {
                if (tableInterpTypeT == "cubic")
                {
                    T[cellI] = interpolateSplineXY(z,zProfile,TProfile);
                }
                else
                {
                    T[cellI] = interpolateXY(z,zProfile,TProfile);
                }
            }

            if (z/zMax < zPeak)
            {
                T[cellI] += TPrime;
            }
        }
    }

    // Modified pressure.
    Info << "Updating internal p_rgh field..." << endl;
    if ( velocityInitType != "scaleGiven" )
    {
        forAll(p_rgh,cellI)
        {
            p_rgh[cellI] = 0.0;
        }
    }
}


// Update the boundary field.
if (updateBoundaryFields)
{
    Info << "Updating boundaries..." << endl;
    U.correctBoundaryConditions();
    T.correctBoundaryConditions();
    p_rgh.correctBoundaryConditions();
}


// Copy the intial fields for U and T to the old time.
// This allows the 2nd-order backward time differencing
// to not revert to 1st-order in the first time step.
U0 = U;
T0 = T;


// Write out the updated fields.
Info<< "Writing field U" << endl;
U.write();
U0.write();
Info<< "Writing field T" << endl;
T.write(); 
T0.write();
Info<< "Writing field p_rgh" << endl;
p_rgh.write();


return 0;
}


// ************************************************************************* //

