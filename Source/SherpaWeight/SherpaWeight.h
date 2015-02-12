////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaWeight.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPAWEIGHT_H
#define SHERPAWEIGHT_H

#include "common.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

namespace SHERPA
{
class Sherpa;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaWeight
{
public:
    struct ReweightParameter
    {
        std::string     name;
        double          scale   = 1.0;
        double          offset  = 0.0;
    };
    
    typedef std::vector<ReweightParameter>  ParameterVector;
    typedef std::vector<std::string>        StringVector;
    typedef std::vector<double>             DoubleVector;
    typedef std::vector<DoubleVector>       DoubleMatrix;
    
public:
    SherpaWeight();
    ~SherpaWeight() throw();

    void Initialize( const std::string & eventFileName, const std::vector<const char *> & argv );
    
    const ParameterVector &  Parameters()       const           { return m_parameters; }
    const StringVector &     CoefficientNames() const           { return m_coefNames;  }

    size_t NParameters()   const throw()                        { return m_parameters.size(); }
    size_t NCoefficients() const throw()                        { return m_coefNames.size();  }
    
    const DoubleMatrix & EvaluationMatrix() const               { return m_evalMatrix;    }
    const DoubleMatrix & InverseCoefficientMatrix() const       { return m_invCoefMatrix; }

    void EvaluateEvents();

    const DoubleVector & MatrixElements(    int32_t eventId ) const;
    const DoubleVector & CoefficientValues( int32_t eventId ) const;

    
    static void GetBilinearMatrices( const ParameterVector & parameters, DoubleMatrix & evalMatrix,
                                                                         DoubleMatrix & invCoefMatrix,
                                                                         StringVector & coefNames );

private:
    static void CreateFeynRulesParamCard( const std::string & srcFilePath,    const std::string & dstFilePath,
                                          const ParameterVector & parameters, const DoubleVector & paramValues );

private:
    std::unique_ptr<SHERPA::Sherpa>     m_upSherpa;
    std::string                         m_eventFileName;
    ParameterVector                     m_parameters;
    DoubleMatrix                        m_evalMatrix;
    DoubleMatrix                        m_invCoefMatrix;
    StringVector                        m_coefNames;                // size = m_nCoefficients
    StringVector                        m_paramCards;               // size = m_nCoefficients
    std::map< int32_t, DoubleVector >   m_matrixElements;

private:
    SherpaWeight(const SherpaWeight &) = delete;
    SherpaWeight & operator=(const SherpaWeight &) = delete;
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPAWEIGHT_H
