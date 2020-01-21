#include "TopTagger/TopTagger/interface/TTMConstituentReqs.h"

#include "TopTagger/TopTagger/interface/Constituent.h"
#include "TopTagger/CfgParser/include/Context.hh"
#include "TopTagger/CfgParser/include/CfgDocument.hh"

#include <iostream>

void TTMConstituentReqs::getParameters(const cfg::CfgDocument* cfgDoc, const std::string& localContextName)
{
    //Construct contexts
    cfg::Context commonCxt("Common");
    cfg::Context localCxt(localContextName);

    //common parameter
    dRMatch_      = cfgDoc->get("dRMatch",      commonCxt, -999.9);
    dRMatchAK8_   = cfgDoc->get("dRMatchAK8",   commonCxt, -999.9);
    
    //monojet parameters
    minAK8TopMass_    = cfgDoc->get("minAK8TopMass",    localCxt, -999.9);
    maxAK8TopMass_    = cfgDoc->get("maxAK8TopMass",    localCxt, -999.9);
    maxTopTau32_      = cfgDoc->get("maxTopTau32",      localCxt, -999.9);
    minAK8TopPt_      = cfgDoc->get("minAK8TopPt",      localCxt, -999.9);
    deepAK8TopDisc_   = cfgDoc->get("deepAK8TopDisc",   localCxt, -999.9);

    //mono-W parameters
    deepAK8WDisc_     = cfgDoc->get("deepAK8WDisc",     localCxt, -999.9);

    //dijet parameters
    minAK8WMass_      = cfgDoc->get("minAK8WMass",      localCxt, -999.9);
    maxAK8WMass_      = cfgDoc->get("maxAK8WMass",      localCxt, -999.9);
    maxWTau21_        = cfgDoc->get("maxWTau21",        localCxt, -999.9);
    minAK8WPt_        = cfgDoc->get("minAK8WPt",        localCxt, -999.9);
    minAK4WPt_        = cfgDoc->get("minAK4WPt",        localCxt, -999.9);
    
    // testing
    std::cout << "In " << __func__ << ": minAK8WMass_ = " << minAK8WMass_ << std::endl;
    std::cout << "In " << __func__ << ": maxAK8WMass_ = " << maxAK8WMass_ << std::endl;
    std::cout << "In " << __func__ << ": minAK8TopMass_ = " << minAK8TopMass_ << std::endl;
    std::cout << "In " << __func__ << ": maxAK8TopMass_ = " << maxAK8TopMass_ << std::endl;
    //trijet parameters
}

bool TTMConstituentReqs::passAK8WReqs(const Constituent& constituent) const
{
    //check that it is an AK8 jet
    if(constituent.getType() != Constituent::AK8JET) return false;

    //check that tau1 and 2 are valid
    if(constituent.getTau1() <= 0 || constituent.getTau2() <= 0) return false;

    double tau21 = constituent.getTau2()/constituent.getTau1();

    return float(constituent.p().Pt()) > minAK8WPt_ &&
           float(constituent.getSoftDropMass()) * float(constituent.getWMassCorr()) > minAK8WMass_  && 
           float(constituent.getSoftDropMass()) * float(constituent.getWMassCorr()) < maxAK8WMass_ &&
           float(tau21) < maxWTau21_;
}

bool TTMConstituentReqs::passAK4WReqs(const Constituent& constituent, const Constituent& constituentAK8) const
{
    //basic AK4 jet requirements 
    bool basicReqs = constituent.getType() == Constituent::AK4JET && constituent.p().Pt() > minAK4WPt_;

    //check that the AK4 jet does not overlap with the selected AK8 subjets
    //to keep the algorithm from needing to do unnecessary calculations on all matches
    //we have to have the AK8 jet in the first position, not the second arguement
    bool overlapCheck = constituentsAreUsed({&constituentAK8}, {&constituent}, dRMatch_, dRMatchAK8_);

    return basicReqs && !overlapCheck;
}

bool TTMConstituentReqs::passAK8TopReqs(const Constituent& constituent) const
{
    //check that it is an AK8 jet
    if(constituent.getType() != Constituent::AK8JET) return false;

    //check that tau2 and 3 are valid
    if(constituent.getTau2() <= 0 || constituent.getTau3() <= 0) return false;

    double tau32 = constituent.getTau3()/constituent.getTau2();

    return float(constituent.p().Pt()) > minAK8TopPt_ &&
           float(constituent.getSoftDropMass()) > minAK8TopMass_  && 
           float(constituent.getSoftDropMass()) < maxAK8TopMass_ &&
           float(tau32) < maxTopTau32_;
}

bool TTMConstituentReqs::passDeepAK8WReqs(const Constituent& constituent) const
{
    //check that it is an AK8 jet
    if(constituent.getType() != Constituent::AK8JET) return false;

    return float(constituent.p().Pt()) > minAK8WPt_ &&
           float(constituent.getSoftDropMass()) > minAK8WMass_  && 
           float(constituent.getSoftDropMass()) < maxAK8WMass_ &&
           float(constituent.getWDisc()) > deepAK8WDisc_;
}

bool TTMConstituentReqs::passDeepAK8TopReqs(const Constituent& constituent) const
{
    //check that it is an AK8 jet
    if(constituent.getType() != Constituent::AK8JET) return false;

    return float(constituent.p().Pt()) > minAK8TopPt_ &&
           float(constituent.getSoftDropMass()) > minAK8TopMass_  && 
           float(constituent.getSoftDropMass()) < maxAK8TopMass_ &&
           float(constituent.getTopDisc()) > deepAK8TopDisc_;
}

bool TTMConstituentReqs::passAK4ResolvedReqs(const Constituent& constituent, const double minPt) const
{
    return constituent.getType() == Constituent::AK4JET && constituent.p().Pt() > minPt;
}
