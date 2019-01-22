/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = ef195594f0485c87e797fc51fa01cacf5b498bee
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void ground_ef195594f0485c87e797fc51fa01cacf5b498bee(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchScalarField,
    groundFixedValueFvPatchScalarField
);


const char* const groundFixedValueFvPatchScalarField::SHA1sum =
    "ef195594f0485c87e797fc51fa01cacf5b498bee";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

groundFixedValueFvPatchScalarField::
groundFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct ground sha1: ef195594f0485c87e797fc51fa01cacf5b498bee"
            " from patch/DimensionedField\n";
    }
}


groundFixedValueFvPatchScalarField::
groundFixedValueFvPatchScalarField
(
    const groundFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct ground sha1: ef195594f0485c87e797fc51fa01cacf5b498bee"
            " from patch/DimensionedField/mapper\n";
    }
}


groundFixedValueFvPatchScalarField::
groundFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct ground sha1: ef195594f0485c87e797fc51fa01cacf5b498bee"
            " from patch/dictionary\n";
    }
}


groundFixedValueFvPatchScalarField::
groundFixedValueFvPatchScalarField
(
    const groundFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct ground sha1: ef195594f0485c87e797fc51fa01cacf5b498bee"
            " as copy\n";
    }
}


groundFixedValueFvPatchScalarField::
groundFixedValueFvPatchScalarField
(
    const groundFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct ground sha1: ef195594f0485c87e797fc51fa01cacf5b498bee "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

groundFixedValueFvPatchScalarField::
~groundFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy ground sha1: ef195594f0485c87e797fc51fa01cacf5b498bee\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void groundFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs ground sha1: ef195594f0485c87e797fc51fa01cacf5b498bee\n";
    }

//{{{ begin code
    #line 34 "/home/statisdisc/OpenFOAM/statisdisc/run/partitionedShallowWater/surfaceHeatingSingleFluid/0/theta.boundaryField.ground"
vector dir = vector(1,1,0);
            scalarField var = (patch().Cf() & dir);
            scalarField value = 300 + var/1000.;
            operator==(value);
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

