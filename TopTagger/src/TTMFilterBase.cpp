#include "TopTagger/TopTagger/interface/TTMFilterBase.h"

#include "TopTagger/TopTagger/interface/Constituent.h"
#include "TopTagger/CfgParser/include/Context.hh"
#include "TopTagger/CfgParser/include/CfgDocument.hh"

#include "Math/VectorUtil.h"

bool TTMFilterBase::constituentsAreUsed(const std::vector<const Constituent*>& constituents, const std::set<const Constituent*>& usedConsts, const double dRMax, const double dRMaxAK8) const
{
    // loop over set of used constituents for testing
    printf("List of used constituents:\n");
    for(const auto& usedConstituent : usedConsts)
    {
        printf(" - used constituent: (pt=%.5lf, eta=%.5lf, phi=%.5lf, mass=%.5lf), type=%d\n", usedConstituent->p().Pt(), usedConstituent->p().Eta(), usedConstituent->p().Phi(), usedConstituent->p().M(), usedConstituent->getType());
    }
    for(const auto& constituent : constituents)
    {
        if(usedConsts.count(constituent) > 0)
        {
            //First return true if constituent is found (this covers all AK4 and most AK8 jets)
            printf("In %s: constituent already used: (pt=%.5lf, eta=%.5lf, phi=%.5lf, mass=%.5lf), type=%d\n", __func__, constituent->p().Pt(), constituent->p().Eta(), constituent->p().Phi(), constituent->p().M(), constituent->getType());
            printf("In %s: return true\n", __func__);
            return true;
        }
        //else if(constituent->getType() == Constituent::AK8JET)
        //{
        //    //If the constituent is AK8 we will also check its subjets are not used
        //    if(constituent->getSubjets().size() <= 1)
        //    {
        //        // If this jet has only one subjet, use matching to the overall AK8 jet instead 
        //        for(const auto& usedConstituent : usedConsts)
        //        {
        //            if(ROOT::Math::VectorUtil::DeltaR(constituent->p(), usedConstituent->p()) < dRMaxAK8)
        //            {
        //                //we found a match
        //                //printf("In %s: return true 2; dR=%f, dRMaxAK8=%f\n", __func__, ROOT::Math::VectorUtil::DeltaR(constituent->p(), usedConstituent->p()), dRMaxAK8);
        //                return true;
        //            }
        //        }
        //    }
        //    else
        //    {
        //        for(const auto& subjet : constituent->getSubjets())
        //        {
        //            for(const auto& usedConstituent : usedConsts)
        //            {
        //                if(ROOT::Math::VectorUtil::DeltaR(subjet.p(), usedConstituent->p()) < dRMax)
        //                {
        //                    //we found a match
        //                    //printf("In %s: return true 3; dR=%f, dRMax=%f\n", __func__, ROOT::Math::VectorUtil::DeltaR(subjet.p(), usedConstituent->p()), dRMax);
        //                    return true;
        //                }
        //            }
        //        }
        //    }
        //}
    }

    //if nothing is found then we have an unused jet
    printf("In %s: return false\n", __func__);
    return false;
}

void TTMFilterBase::markConstituentsUsed(const std::vector<const Constituent *>& constituents, const std::vector<Constituent>& allConstituents, std::set<const Constituent*>& usedConstituents, const double dRMax, const double dRMaxAK8) const
{
    for(const auto& constituent : constituents)
    {
        //No matter what, add the main constituent
        usedConstituents.insert(constituent);

        //If the constituent is an AK8JET, then add AK4JETs matching its subjets as well 
        if(constituent->getType() == Constituent::AK8JET)
        {
            if(constituent->getSubjets().size() <= 1)
            {
                //If there is one or fewer subjets, instead match to the overall AK8 jet
                for(const auto& matchConst : allConstituents) 
                {
                    if(constituent->getType() == Constituent::AK4JET && ROOT::Math::VectorUtil::DeltaR(constituent->p(), matchConst.p()) < dRMaxAK8)
                    {
                        usedConstituents.insert(&matchConst);
                    }
                }
            }
            else
            {
                for(const auto& subjet : constituent->getSubjets())
                {
                    for(const auto& matchConst : allConstituents) 
                    {
                        if(constituent->getType() == Constituent::AK4JET && ROOT::Math::VectorUtil::DeltaR(subjet.p(), matchConst.p()) < dRMax)
                        {
                            usedConstituents.insert(&matchConst);
                        }
                    }
                }
            }
        }
    }
}


