////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeight.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPAWEIGHT_H
#define SHERPAWEIGHT_H

#include "common.h"

#include <TMatrixDfwd.h> // forward declaration only

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

namespace ATOOLS
{

template<typename>
class Vec4;

typedef Vec4<double>            Vec4D;
typedef std::vector<Vec4D>      Vec4D_Vector;

}

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaWeight
{
public:
    SherpaWeight();
    ~SherpaWeight() throw();

    void Initialize( const std::vector<const char *> & argv );

    void ProcessEvent( int32_t eventId, size_t nInParticles, const std::vector<int> & particleCodes, const ATOOLS::Vec4D_Vector & particleMomenta );

private:
    struct ReweightParameter
    {
        std::string     name;
        double          scale   = 1.0;
        double          offset  = 0.0;
    };

private:
    static void GetBilinearMatrices( const std::vector<ReweightParameter> & parameters, TMatrixD & evalMatrix, TMatrixD & invCoefMatrix );

    static void CreateParamCard( const std::string & srcFilePath, const std::string & dstFilePath,
                                 const std::vector<ReweightParameter> & parameters,
                                 const std::vector<double> & paramValues );

private:
    std::vector<const char *>           m_argv;
    std::vector<ReweightParameter>      m_parameters;
    size_t                              m_nEvaluations      = 0;    // number of evaluations
    std::vector<std::string>            m_paramCards;               // size = m_nEvaluations

private:
    SherpaWeight(const SherpaWeight &) = delete;
    SherpaWeight & operator=(const SherpaWeight &) = delete;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPAWEIGHT_H
