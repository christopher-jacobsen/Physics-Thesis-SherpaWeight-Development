////////////////////////////////////////////////////////////////////////////////////////////////////
// SherpaWeight.cpp
////////////////////////////////////////////////////////////////////////////////////////////////////

#include "SherpaWeight.h"
#include "SherpaMECalculator.h"

#include <SHERPA/Main/Sherpa.H>
#include <ATOOLS/Math/Vector.H>
//#include <ATOOLS/Phys/Cluster_Amplitude.H>

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeight::SherpaWeight()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
SherpaWeight::~SherpaWeight() throw()
{
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::Initialize( const std::vector<const char *> & argv )
{
    m_argv = argv;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
void SherpaWeight::ProcessEvent( int32_t eventId, size_t nInParticles, const std::vector<int> & particleCodes, const ATOOLS::Vec4D_Vector & particleMomenta )
{
    if (!nInParticles || (nInParticles > particleCodes.size()))
    {
        throw int(-2);  // TODO
    }
    
    if (particleCodes.size() != particleMomenta.size())
    {
        throw int(-2);  // TODO
    }

    SHERPA::Sherpa sherpa;
    
    if (!sherpa.InitializeTheRun( static_cast<int>(m_argv.size()), const_cast<char **>(m_argv.data()) ))
    {
        LogMsgError( "Failed to initialize Sherpa framework. Check Run.dat file." );
        throw int(-2);  // TODO
    }
    
    SherpaMECalculator meCalc( &sherpa );

    for (size_t i = 0; i < particleCodes.size(); ++i)
    {
        if (i < nInParticles)
            meCalc.AddInFlav( particleCodes[i] );
        else
            meCalc.AddOutFlav( particleCodes[i] );
    }

    try
    {
        meCalc.Initialize();
    }
    catch (const ATOOLS::Exception & error)
    {
        LogMsgError( "Sherpa exception caught: \"%hs\" in %hs::%hs",
            FMT_HS(error.Info().c_str()), FMT_HS(error.Class().c_str()), FMT_HS(error.Method().c_str()) );
        
        LogMsgInfo( "Sherpa process name: \"%hs\"", FMT_HS(meCalc.Name().c_str()) );
        
        const ATOOLS::Cluster_Amplitude * pAmp = meCalc.GetAmp();
        if (!pAmp)
            LogMsgInfo("Sherpa process cluster amplitude: null");
        else
        {
            //const ATOOLS::ClusterLeg_Vector & legs = pAmp->Legs();
            //int bob = 5;
        }

        throw;
    }

    meCalc.SetMomenta( particleMomenta );
    
    double me = meCalc.MatrixElement();
    LogMsgInfo( "Event %i: ME=%E", FMT_I(eventId), FMT_F(me) );

#if 0
    {
        MODEL::ScalarConstantsMap * pConstants = meCalc.GetModelScalarConstants();
        if (pConstants)
        {
            /**/
            for (const auto & entry : *pConstants)
            {
                std::cout << entry.first << ": " << entry.second << std::endl;
            }
            /**/
            
            double aEWM1 = 227.9;
            
            MODEL::ScalarConstantsMap::iterator itrFind = pConstants->find( "aEWM1" );
            if (itrFind != pConstants->end())
            {
                MODEL::ScalarConstantsMap::value_type & constantPair = *itrFind;
                constantPair.second = aEWM1;
            
                double sw = pConstants->at("sw");
                double cw = pConstants->at("cw");

                double aEW = pow(aEWM1,-1.);
                double ee = 2.*sqrt(aEW)*sqrt(M_PI);
                double gw = ee*pow(sw,-1.);
                double g1 = ee*pow(cw,-1.);

                pConstants->at("aEW") = aEW;
                pConstants->at("ee")  = ee;
                pConstants->at("gw")  = gw;
                pConstants->at("g1")  = g1;
                
            }

            /**/
            std::cout << "------------" << std::endl;
            for (const auto & entry : *pConstants)
            {
                std::cout << entry.first << ": " << entry.second << std::endl;
            }
            /**/
        }
    }
    
    double me2 = meCalc.MatrixElement();
    LogMsgInfo( "Event %i: ME2=%E", FMT_I(eventId), FMT_F(me2) );
#endif
}
