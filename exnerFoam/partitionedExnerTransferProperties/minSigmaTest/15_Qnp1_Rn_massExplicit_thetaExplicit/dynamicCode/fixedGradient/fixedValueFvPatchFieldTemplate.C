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
    // SHA1 = 813a50c402ef2c3ee93508a8ff923908489e9239
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void fixedGradient_813a50c402ef2c3ee93508a8ff923908489e9239(bool load)
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
    fixedGradientFixedValueFvPatchScalarField
);


const char* const fixedGradientFixedValueFvPatchScalarField::SHA1sum =
    "813a50c402ef2c3ee93508a8ff923908489e9239";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

fixedGradientFixedValueFvPatchScalarField::
fixedGradientFixedValueFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF)
{
    if (false)
    {
        Info<<"construct fixedGradient sha1: 813a50c402ef2c3ee93508a8ff923908489e9239"
            " from patch/DimensionedField\n";
    }
}


fixedGradientFixedValueFvPatchScalarField::
fixedGradientFixedValueFvPatchScalarField
(
    const fixedGradientFixedValueFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct fixedGradient sha1: 813a50c402ef2c3ee93508a8ff923908489e9239"
            " from patch/DimensionedField/mapper\n";
    }
}


fixedGradientFixedValueFvPatchScalarField::
fixedGradientFixedValueFvPatchScalarField
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
        Info<<"construct fixedGradient sha1: 813a50c402ef2c3ee93508a8ff923908489e9239"
            " from patch/dictionary\n";
    }
}


fixedGradientFixedValueFvPatchScalarField::
fixedGradientFixedValueFvPatchScalarField
(
    const fixedGradientFixedValueFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf)
{
    if (false)
    {
        Info<<"construct fixedGradient sha1: 813a50c402ef2c3ee93508a8ff923908489e9239"
            " as copy\n";
    }
}


fixedGradientFixedValueFvPatchScalarField::
fixedGradientFixedValueFvPatchScalarField
(
    const fixedGradientFixedValueFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF)
{
    if (false)
    {
        Info<<"construct fixedGradient sha1: 813a50c402ef2c3ee93508a8ff923908489e9239 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

fixedGradientFixedValueFvPatchScalarField::
~fixedGradientFixedValueFvPatchScalarField()
{
    if (false)
    {
        Info<<"destroy fixedGradient sha1: 813a50c402ef2c3ee93508a8ff923908489e9239\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void fixedGradientFixedValueFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs fixedGradient sha1: 813a50c402ef2c3ee93508a8ff923908489e9239\n";
    }

//{{{ begin code
    #line 32 "/home/statisdisc/OpenFOAM/statisdisc/run/partitionedShallowWater/surfaceHeatingSingleFluid/0/theta.boundaryField.ground"
const scalar t = this->db().time().value();
            const scalar radius = 3.5e3;
            vector dir = vector(1,0,0);
            scalarField var = (patch().Cf() & dir);
            scalarField value = min(scalar(radius),max(scalar(-radius),var));
            value = -(value-radius)*(value+radius)*1e5 + 1;
            Info << min(value) << " " << max(value) << endl;
            operator==(value);
//}}} end code

    this->fixedValueFvPatchField<scalar>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

