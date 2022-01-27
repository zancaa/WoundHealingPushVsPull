/*

Copyright (c) 2005-2021, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "ContactInhibitionDistanceToWoundCellCycleModel.hpp"
#include "AbstractCellBasedTestSuite.hpp"
#include "CellLabel.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "TransitCellProliferativeType.hpp"
#include "QuiescentCellProperty.hpp"
#include "SmartPointers.hpp"
#include "Warnings.hpp"
#include "Debug.hpp"

ContactInhibitionDistanceToWoundCellCycleModel::ContactInhibitionDistanceToWoundCellCycleModel()
    : UniformCellCycleModel(),
      mMeanQuiescentVolumeFraction(0.9),
      mIntervalWidthQuiescentVolumeFraction(0.05),
      mDistanceToWoundMin(4.0),
      mDistanceToWoundMax(10.0)
{
}

ContactInhibitionDistanceToWoundCellCycleModel::ContactInhibitionDistanceToWoundCellCycleModel(const ContactInhibitionDistanceToWoundCellCycleModel& rModel)
    : UniformCellCycleModel(rModel),
      mMeanQuiescentVolumeFraction(rModel.mMeanQuiescentVolumeFraction),
      mIntervalWidthQuiescentVolumeFraction(rModel.mIntervalWidthQuiescentVolumeFraction),
      mDistanceToWoundMin(rModel.mDistanceToWoundMin),
      mDistanceToWoundMax(rModel.mDistanceToWoundMax)
{
    /*
     * Initialize only those member variables defined in this class.
     *
     * The member variables mBirthTime, mReadyToDivide and mDimension
     * are initialized in the AbstractCellCycleModel constructor.
     */
}

void ContactInhibitionDistanceToWoundCellCycleModel::GenerateStochasticQuiescentVolume()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double mean = mMeanQuiescentVolumeFraction;
    double half_width = 0.5*mIntervalWidthQuiescentVolumeFraction;

    double qvf = (mean - half_width) + (2.0*half_width) * p_gen->ranf(); // U[MinCA,MaxCA]

    mpCell->GetCellData()->SetItem("quiescent volume fraction", qvf);
}

void ContactInhibitionDistanceToWoundCellCycleModel::GenerateStochasticTargetArea()
{
    RandomNumberGenerator* p_gen = RandomNumberGenerator::Instance();

    double min_cell_area = 0.9;
    double max_cell_area = 1.1;

    double cell_target_area = min_cell_area + (max_cell_area - min_cell_area) * p_gen->ranf(); // U[MinCA,MaxCA]
    mpCell->GetCellData()->SetItem("target area", cell_target_area);
}

void ContactInhibitionDistanceToWoundCellCycleModel::InitialiseDaughterCell()
{
    UniformCellCycleModel::InitialiseDaughterCell();
    GenerateStochasticQuiescentVolume();
    GenerateStochasticTargetArea();
}

void ContactInhibitionDistanceToWoundCellCycleModel::Initialise()
{
    UniformCellCycleModel::Initialise();
    GenerateStochasticQuiescentVolume();
    GenerateStochasticTargetArea();
}

bool ContactInhibitionDistanceToWoundCellCycleModel::ReadyToDivide()
{
    assert(mpCell != nullptr);

    MAKE_PTR(QuiescentCellProperty, p_quiescent);
    unsigned cell_id = mpCell->GetCellId();

    int wound_dist = mpCell->GetCellData()->GetItem("boundary distance");
    if ( (wound_dist < mDistanceToWoundMin) || (wound_dist > mDistanceToWoundMax))
    {
        // Cells outside of the proliferative hub are not proliferative
        boost::shared_ptr<AbstractCellProperty> p_differentiated_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<DifferentiatedCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_differentiated_type);
        mpCell->RemoveCellProperty<CellLabel>();
        mpCell->RemoveCellProperty<QuiescentCellProperty>();
        mReadyToDivide = false;
    }
    else
    {
        boost::shared_ptr<AbstractCellProperty> p_transit_type =
            mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<TransitCellProliferativeType>();
        mpCell->SetCellProliferativeType(p_transit_type);
    }

    if (!(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
    {
        // Check that cells are large enough to divide. 
        // Get cell volume
        double cell_volume = mpCell->GetCellData()->GetItem("volume");
        double qvf = mpCell->GetCellData()->GetItem("quiescent volume fraction");
        double target_area = mpCell->GetCellData()->GetItem("target area");
        
        // Update cell cycle duration based on cell volume
        double quiescent_volume = target_area * qvf;

        if (cell_volume < quiescent_volume)
        {
            // Update labels
            mpCell->RemoveCellProperty<CellLabel>();
            mpCell->RemoveCellProperty<QuiescentCellProperty>();
            boost::shared_ptr<AbstractCellProperty> p_label =
                mpCell->rGetCellPropertyCollection().GetCellPropertyRegistry()->Get<CellLabel>();
            mpCell->AddCellProperty(p_label);
            mpCell->AddCellProperty(p_quiescent);

            double dt = SimulationTime::Instance()->GetTimeStep();
            mCellCycleDuration += dt;

        }
        else
        {
            // Remove cell properties
            mpCell->RemoveCellProperty<CellLabel>();
            mpCell->RemoveCellProperty<QuiescentCellProperty>();
        }
    }

    if (!mReadyToDivide)
    {
        if ((mpCell->GetAge() >= mCellCycleDuration) && !(mpCell->HasCellProperty<QuiescentCellProperty>()) && !(mpCell->GetCellProliferativeType()->IsType<DifferentiatedCellProliferativeType>()))
        {
            mReadyToDivide = true;
        }
    }

    return mReadyToDivide;
}

AbstractCellCycleModel* ContactInhibitionDistanceToWoundCellCycleModel::CreateCellCycleModel()
{
    return new ContactInhibitionDistanceToWoundCellCycleModel(*this);
}

void ContactInhibitionDistanceToWoundCellCycleModel::SetMeanQuiescentVolumeFraction(double meanQuiescentVolumeFraction)
{
    mMeanQuiescentVolumeFraction = meanQuiescentVolumeFraction;
}

double ContactInhibitionDistanceToWoundCellCycleModel::GetMeanQuiescentVolumeFraction() const
{
    return mMeanQuiescentVolumeFraction;
}

void ContactInhibitionDistanceToWoundCellCycleModel::SetIntervalWidthQuiescentVolumeFraction(double intervalWidthQuiescentVolumeFraction)
{
    mIntervalWidthQuiescentVolumeFraction = intervalWidthQuiescentVolumeFraction;
}

double ContactInhibitionDistanceToWoundCellCycleModel::GetIntervalWidthQuiescentVolumeFraction() const
{
    return mIntervalWidthQuiescentVolumeFraction;
}

void ContactInhibitionDistanceToWoundCellCycleModel::SetDistanceToWoundMin(unsigned distanceToWoundMin)
{
    mDistanceToWoundMin = distanceToWoundMin;
}

unsigned ContactInhibitionDistanceToWoundCellCycleModel::GetDistanceToWoundMin() const
{
    return mDistanceToWoundMin;
}

void ContactInhibitionDistanceToWoundCellCycleModel::SetDistanceToWoundMax(unsigned distanceToWoundMax)
{
    mDistanceToWoundMax = distanceToWoundMax;
}

unsigned ContactInhibitionDistanceToWoundCellCycleModel::GetDistanceToWoundMax() const
{
    return mDistanceToWoundMax;
}

void ContactInhibitionDistanceToWoundCellCycleModel::OutputCellCycleModelParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MeanQuiescentVolumeFraction>" << mMeanQuiescentVolumeFraction << "</MeanQuiescentVolumeFraction>\n";
     *rParamsFile << "\t\t\t<IntervalWidthQuiescentVolumeFraction>" << mIntervalWidthQuiescentVolumeFraction << "</IntervalWidthQuiescentVolumeFraction>\n";

    // Call method on direct parent class
    UniformCellCycleModel::OutputCellCycleModelParameters(rParamsFile);
}

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
CHASTE_CLASS_EXPORT(ContactInhibitionDistanceToWoundCellCycleModel)
