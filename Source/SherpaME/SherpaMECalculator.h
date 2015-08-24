////////////////////////////////////////////////////////////////////////////////////////////////////
//  SherpaMECalculator.h
////////////////////////////////////////////////////////////////////////////////////////////////////

#ifndef SHERPAMECALCULATOR_H
#define SHERPAMECALCULATOR_H

#include <vector>
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////
// Sherpa include files

#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Phys/Flavour.H"
#include "MODEL/Main/Model_Base.H"

////////////////////////////////////////////////////////////////////////////////////////////////////
// forward declarations

namespace SHERPA
{
    class Sherpa;
}

namespace ATOOLS
{
    class Cluster_Amplitude;
    class ColorID;
}

namespace PHASIC
{
    class Process_Base;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
class SherpaMECalculator
{
private:
    std::string                     m_name;
    ATOOLS::Cluster_Amplitude *     p_amp;
    SHERPA::Sherpa *                p_gen;
    PHASIC::Process_Base *          p_proc;

    size_t                          m_ncolinds;
    std::vector<std::vector<int> >  m_colcombinations;
    std::vector<int>                m_gluinds, m_quainds, m_quabarinds;
    std::vector<int>                m_inpdgs, m_outpdgs;
    std::vector<size_t>             m_mom_inds;

    size_t                          m_npsp;
    size_t                          m_nin;
    size_t                          m_nout;

    void SetMomentumIndices(const std::vector<int> &pdgs);

    PHASIC::Process_Base * FindProcess();

public:
    SherpaMECalculator(SHERPA::Sherpa * Generator);
    ~SherpaMECalculator() throw();

    void AddInFlav(  const int & id);
    void AddOutFlav( const int & id);
    void AddInFlav(  const int & id, const int & col1, const int & col2);
    void AddOutFlav( const int & id, const int & col1, const int & col2);

    double GenerateColorPoint();
    bool HasColorIntegrator();
    void SetColors();

    void Initialize();

    size_t NumberOfPoints();

    void SetMomenta(size_t n);
    void SetMomenta(const std::vector<double *> & p);
    void SetMomenta(const ATOOLS::Vec4D_Vector & p);
    void SetMomentum(const size_t & id, const double & E , const double & px,
                                        const double & py, const double & pz);
    void SetMomentum(const size_t & id, const ATOOLS::Vec4D & p);

  //MODEL::ScalarConstantsMap * GetModelScalarConstants();

    double MatrixElement();
    double CSMatrixElement();

    double GetFlux();
  
    std::string GeneratorName();

    ATOOLS::Flavour GetFlav(size_t i);
    ATOOLS::Flavour GetInFlav(size_t i)             { return GetFlav(i); }
    ATOOLS::Flavour GetOutFlav(size_t i)            { return GetFlav(i+m_nin); }

    ATOOLS::Cluster_Amplitude * GetAmp()            { return p_amp; }

    std::string Name()                              { return m_name; }
};

////////////////////////////////////////////////////////////////////////////////////////////////////

#endif // SHERPAMECALCULATOR_H
