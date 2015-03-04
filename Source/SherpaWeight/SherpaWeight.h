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

namespace ATOOLS
{
class Data_Reader;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight
////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaWeight
{
public:     ////// public types //////

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
    
    struct ModelInterface
    {
        virtual ~ModelInterface() = default;

        virtual void Initialize( SherpaWeight & parent, ATOOLS::Data_Reader & modelFileSectionReader ) = 0;

        virtual void ValidateParameters( const ParameterVector & params ) = 0;

        virtual std::string CommandLineArgs( const ParameterVector & params, const DoubleVector & paramValues )  = 0;
    };

    class FeynRulesModel;

public:     ////// public methods //////

    SherpaWeight();
    ~SherpaWeight() throw();

    void Initialize( const std::string & eventFileName, const std::vector<const char *> & argv,
                     bool bReadParametersFromFile = true );

    const std::string & ApplicationRunPath() const throw()  { return m_appRunPath;    }
    const std::string & SherpaRunPath()      const throw()  { return m_sherpaRunPath; }
    const std::string & TemporaryPath()      const throw()  { return m_tmpPath;       }
    
    void ReadParametersFromFile( const char * filePath = nullptr );  // filePath can contain section definition
    void SetParameters( const ParameterVector & params );
    
    size_t NParameters()   const throw()                        { return m_parameters.size();    }
    size_t NEvaluations()  const throw()                        { return m_evalMatrix.size();    }  // equal to NCoefficients()
    size_t NCoefficients() const throw()                        { return m_invCoefMatrix.size(); }  // equal to NEvaluations()
    
    const ParameterVector &  Parameters()               const   { return m_parameters;    }
    const DoubleMatrix &     EvaluationMatrix()         const   { return m_evalMatrix;    }
    const DoubleMatrix &     InverseCoefficientMatrix() const   { return m_invCoefMatrix; }
    const StringVector &     CoefficientNames()         const   { return m_coefNames;     }

    void EvaluateEvents();

    const DoubleVector & MatrixElements(    int32_t eventId ) const;
    const DoubleVector   CoefficientValues( int32_t eventId ) const;

    
    static void GetBilinearMatrices( const ParameterVector & parameters, DoubleMatrix & evalMatrix,
                                                                         DoubleMatrix & invCoefMatrix,
                                                                         StringVector & coefNames );

private:    ////// private types //////

    typedef std::map< int32_t, DoubleVector > EventMatrixElementMap;

private:    ////// private methods //////

    void AddMatrixElementsFromFile( const char * filePath );
    void AddMatrixElement( int32_t eventId, double me );
    
private:    ////// private data //////

    std::unique_ptr<SHERPA::Sherpa>     m_upSherpa;
    std::vector<const char *>           m_argv;
    std::string                         m_appRunPath;
    std::string                         m_sherpaRunPath;
    std::string                         m_tmpPath;
    std::string                         m_sherpaWeightFileSection;

    ModelInterface *                    m_pModel    = nullptr;
    bool                                m_bOwnModel = true;

    ParameterVector                     m_parameters;
    DoubleMatrix                        m_evalMatrix;
    DoubleMatrix                        m_invCoefMatrix;
    StringVector                        m_coefNames;

    std::string                         m_eventFileName;
    EventMatrixElementMap               m_matrixElements;

private:
    SherpaWeight(const SherpaWeight &)              = delete;   // disable copy constructor
    SherpaWeight & operator=(const SherpaWeight &)  = delete;   // disable assignment operator
};

////////////////////////////////////////////////////////////////////////////////////////////////////
// class SherpaWeight::FeynRulesModel
////////////////////////////////////////////////////////////////////////////////////////////////////

class SherpaWeight::FeynRulesModel : public SherpaWeight::ModelInterface
{
public:
    FeynRulesModel();
    virtual ~FeynRulesModel() throw();

public:  // ModelInterface overrides

    virtual void Initialize( SherpaWeight & parent, ATOOLS::Data_Reader & modelFileSectionReader );

    virtual void ValidateParameters( const ParameterVector & params );

    virtual std::string CommandLineArgs( const ParameterVector & params, const DoubleVector & paramValues );

protected:
    static void CreateFeynRulesParamCard( const std::string & srcFilePath,    const std::string & dstFilePath,
                                          const ParameterVector & parameters, const DoubleVector & paramValues );

private:
    SherpaWeight *  m_pParent               = nullptr;
    std::string     m_sourceParamCardFile;

private:
    FeynRulesModel(const FeynRulesModel &)              = delete;   // disable copy constructor
    FeynRulesModel & operator=(const FeynRulesModel &)  = delete;   // disable assignment operator
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPAWEIGHT_H
