/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2021 OpenCFD Ltd.
    Copyright (C) YEAR AUTHOR, AFFILIATION
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
#include "PatchFunction1.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// dynamicCode:
// SHA1 = 7c66ca5fd1f6b26464b202f883aea45a035e62ca
//
// unique function name that can be checked if the correct library version
// has been loaded
extern "C" void vaporizationFlowRate_7c66ca5fd1f6b26464b202f883aea45a035e62ca(bool load)
{
    if (load)
    {
        // Code that can be explicitly executed after loading
    }
    else
    {
        // Code that can be explicitly executed before unloading
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    vaporizationFlowRateFixedValueFvPatchVectorField
);

} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::
vaporizationFlowRateFixedValueFvPatchVectorField::
vaporizationFlowRateFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(p, iF)
{
    if (false)
    {
        printMessage("Construct vaporizationFlowRate : patch/DimensionedField");
    }
}


Foam::
vaporizationFlowRateFixedValueFvPatchVectorField::
vaporizationFlowRateFixedValueFvPatchVectorField
(
    const vaporizationFlowRateFixedValueFvPatchVectorField& rhs,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    parent_bctype(rhs, p, iF, mapper)
{
    if (false)
    {
        printMessage("Construct vaporizationFlowRate : patch/DimensionedField/mapper");
    }
}


Foam::
vaporizationFlowRateFixedValueFvPatchVectorField::
vaporizationFlowRateFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    parent_bctype(p, iF, dict)
{
    if (false)
    {
        printMessage("Construct vaporizationFlowRate : patch/dictionary");
    }
}


Foam::
vaporizationFlowRateFixedValueFvPatchVectorField::
vaporizationFlowRateFixedValueFvPatchVectorField
(
    const vaporizationFlowRateFixedValueFvPatchVectorField& rhs
)
:
    parent_bctype(rhs),
    dictionaryContent(rhs)
{
    if (false)
    {
        printMessage("Copy construct vaporizationFlowRate");
    }
}


Foam::
vaporizationFlowRateFixedValueFvPatchVectorField::
vaporizationFlowRateFixedValueFvPatchVectorField
(
    const vaporizationFlowRateFixedValueFvPatchVectorField& rhs,
    const DimensionedField<vector, volMesh>& iF
)
:
    parent_bctype(rhs, iF)
{
    if (false)
    {
        printMessage("Construct vaporizationFlowRate : copy/DimensionedField");
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::
vaporizationFlowRateFixedValueFvPatchVectorField::
~vaporizationFlowRateFixedValueFvPatchVectorField()
{
    if (false)
    {
        printMessage("Destroy vaporizationFlowRate");
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void
Foam::
vaporizationFlowRateFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        printMessage("updateCoeffs vaporizationFlowRate");
    }

//{{{ begin code
    #line 40 "/home/johan/Documents/PhD/openfoam-tank-mesh/tests/0/U.boundaryField.bottom"
const scalar Q_L =  0.019760904210887362;
            const scalar dHvap = 448916.79;
            const scalar mdot_L = Q_L/dHvap;

            const fvPatch& fvP = this->patch();
            const fvMesh& mesh = fvP.boundaryMesh().mesh();
            const label patchi = fvP.index();
            const vectorField& sf = fvP.Sf();
            const scalarField& Af = fvP.magSf();

            const volScalarField q = mesh.lookupObject<volScalarField>("wallHeatFlux");
            const auto& qpatch = q.boundaryField()[patchi];
            const auto rho = mesh.lookupObject<volScalarField>("rho");
            const auto& rhopatch = rho.boundaryField()[patchi];
            scalarField mflux_LG = -qpatch/dHvap;

            /* scalar Q_patch = gSum(qpatch*Af); */
            /* Info << "Q_patch = " << Q_patch << endl; */

            /* scalar Q_total = Q_L + Q_patch; */
            /* mflux_LG = mflux_LG * (Q_total/Q_patch); */
            /* scalar m_total1 = Q_total/dHvap; */
            /* scalar m_total2 = gSum(mflux_LG*Af); */
            /* Info << "m_total1 = " << m_total1 << endl; */
            /* Info << "m_total2 = " << m_total2 << endl; */
            /* vectorField U = - (sf/Af)*mflux_LG/(rhopatch); */

            vectorField U = -(sf/Af)*mdot_L/(rhopatch*gSum(Af)) - (sf/Af)*mflux_LG/(rhopatch);
            // Info << gSum(Af) << endl;
            // Info << -(sf/Af) << endl;
            // Info << "mdot_L = " << mdot_L << endl;
            // Info << "rhopatch = " << rhopatch*(sf/Af) << endl;
            // Info << "mflux_LG" << mflux_LG << endl;
            // Info << Af << endl;
            // Info << "qpatch = " << qpatch*(sf/Af) << endl;
            // Info << U << endl;

            operator==(U);
//}}} end code

    this->parent_bctype::updateCoeffs();
}


// ************************************************************************* //
