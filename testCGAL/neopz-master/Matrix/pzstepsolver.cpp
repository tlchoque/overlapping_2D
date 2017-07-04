/**
 * @file
 * @brief Contains the implementation of the TPZStepSolver methods.
 */

#include "pzstepsolver.h"
#include "pzmatrixid.h"
#include <stdlib.h>
using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.converge"));
#endif

template <class TVar>
TPZStepSolver<TVar>::TPZStepSolver(TPZAutoPointer<TPZMatrix<TVar> > refmat) : TPZMatrixSolver<TVar>(refmat), fNumIterations(-1) {
	fPrecond = 0;
	ResetSolver();
}

template <class TVar>
TPZStepSolver<TVar>::TPZStepSolver(const TPZStepSolver<TVar> & copy) : TPZMatrixSolver<TVar>(copy), fNumIterations(copy.fNumIterations) , fSingular(copy.fSingular){
    fSolver = copy.fSolver;
    fDecompose = copy.fDecompose;
    fMaxIterations = copy.fMaxIterations;
    fTol = copy.fTol;
    fOverRelax = copy.fOverRelax;
    fPrecond = 0;
    if(copy.fPrecond) fPrecond = copy.fPrecond->Clone();
    fFromCurrent = copy.fFromCurrent;
    fNumVectors = copy.fNumVectors;//Cedric: 24/04/2003 - 12:39
}

template <class TVar>
TPZStepSolver<TVar>::~TPZStepSolver() {
	if(fPrecond) delete fPrecond;
}

/**
 * This method will reset the matrix associated with the solver \n
 * This is useful when the matrix needs to be recomputed in a non linear problem
 */
template <class TVar>
void TPZStepSolver<TVar>::ResetMatrix()
{
	TPZMatrixSolver<TVar>::ResetMatrix();
}

/**
 * @brief Decompose the system of equations if a direct solver is used
 */
template <class TVar>
void TPZStepSolver<TVar>::Decompose()
{
    if (fSolver == this->EDirect) {
        this->Matrix()->Decompose(fDecompose,fSingular);
    }
}

template <class TVar>
void TPZStepSolver<TVar>::Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual){
	if(!this->Matrix()) {
		cout << "TPZMatrixSolver::Solve called without a matrix pointer\n";
		DebugStop();
	}
	
	TPZAutoPointer<TPZMatrix<TVar> > mat = this->Matrix();
    // update the matrix to which the preconditioner refers
    if(fPrecond)
    {
        
        fPrecond->UpdateFrom(this->Matrix());
    }
    
	if(result.Rows() != mat->Rows() || result.Cols() != F.Cols()) {
		result.Redim(mat->Rows(),F.Cols());
	}
	
	if(this->fScratch.Rows() != result.Rows() || this->fScratch.Cols() != result.Cols()) {
		this->fScratch.Redim(result.Rows(),result.Cols());
	}
	
	REAL tol = fTol;
	long numiterations = fMaxIterations;
	switch(fSolver) {
		case TPZStepSolver::ENoSolver:
		default:
			cout << "TPZMatrixSolver::Solve called without initialized solver, Jacobi used\n";
			SetJacobi(1,0.,0);
		case TPZStepSolver::EJacobi:
			//    cout << "fScratch dimension " << fScratch.Rows() << ' ' << fScratch.Cols() << endl;
			mat->SolveJacobi(numiterations,F,result,residual,this->fScratch,tol,fFromCurrent);
            fNumIterations = numiterations;
			break;
		case TPZStepSolver::ESOR:
			mat->SolveSOR(numiterations,F,result,residual,this->fScratch,fOverRelax,tol,fFromCurrent);
            fNumIterations = numiterations;
			break;
		case TPZStepSolver::ESSOR:
			mat->SolveSSOR(numiterations,F,result,residual,this->fScratch,fOverRelax,tol,fFromCurrent);
            fNumIterations = numiterations;
			break;
		case TPZStepSolver::ECG:
			mat->SolveCG(numiterations,*fPrecond,F,result,residual,tol,fFromCurrent);
			cout << "Number of equations " << mat->Rows() << std::endl;
			cout << "Number of CG iterations " << numiterations << " tol = " << tol << endl;
            fNumIterations = numiterations;
            fTol = tol;
#ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Number of equations " << mat->Rows() << std::endl;
                sout << "Number of CG iterations " << numiterations << " tol = " << tol;
                LOGPZ_DEBUG(logger,sout.str().c_str());
            }
#endif
			break;
		case TPZStepSolver::EGMRES: {
			TPZFMatrix<TVar> H(fNumVectors+1,fNumVectors+1,0.);
			mat->SolveGMRES(numiterations,*fPrecond,H,fNumVectors,F,result,residual,tol,fFromCurrent);
            fNumIterations = numiterations;
            cout << "Number of GMRES iterations " << numiterations << " tol = " << tol;
			if(numiterations == fMaxIterations || tol >= fTol)
			{
				std::cout << "GMRes tolerance was not achieved : numiter " << numiterations <<
				" tol " << tol << endl;
			}
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Number of GMRES iterations " << numiterations << " tol = " << tol;
				if(logger->isDebugEnabled()) LOGPZ_DEBUG(logger,sout.str().c_str());
			}
#endif
		}
			break;
		case TPZStepSolver::EBICGSTAB: 
			mat->SolveBICGStab(numiterations, *fPrecond, F, result,residual,tol,fFromCurrent);
            fNumIterations = numiterations;
			
			if(numiterations == fMaxIterations || tol >= fTol)
			{
				std::cout << "BiCGStab tolerance was not achieved : numiter " << numiterations <<
				" tol " << tol << endl;
			}
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Number of BiCGStab iterations " << numiterations << " tol = " << tol;
			LOGPZ_DEBUG(logger,sout.str().c_str());
		}
#endif
			break;
		case TPZStepSolver::EDirect:
			result = F;
			mat->SolveDirect(result,fDecompose,fSingular);
			if(residual) residual->Redim(F.Rows(),F.Cols());
			break;
		case TPZStepSolver::EMultiply:
			mat->Multiply(F,result);
			if(residual) mat->Residual(result,F,*residual);
			
	}
}
template<class TVar>
void TPZStepSolver<TVar>::ResetSolver() {
	fSolver = this->ENoSolver;
	fDecompose  = ENoDecompose;
	fMaxIterations = 0;
    fNumIterations = -1;
	fTol = 0.;
	fNumVectors = 0;
	fOverRelax = 0.;
	if(fPrecond) delete fPrecond;
	fPrecond = 0;
	fFromCurrent = 0;
}
template <class TVar>
void TPZStepSolver<TVar>::SetDirect (const DecomposeType decomp){
	ResetSolver();
	fSolver = this->EDirect;
	fDecompose = decomp;
}
template <class TVar>
void TPZStepSolver<TVar>::SetCG(const long numiterations, const TPZMatrixSolver<TVar> &pre, const REAL tol, const long FromCurrent){
	ResetSolver();
	fSolver = this->ECG;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	//	fPrecond = &pre;
	if(fPrecond) delete fPrecond;
	fPrecond = pre.Clone();
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetGMRES(const long numiterations, const int numvectors, const TPZMatrixSolver<TVar> &pre, const REAL tol, const long FromCurrent){
	ResetSolver();
	fSolver = this->EGMRES;
	fNumVectors = numvectors;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	//	fPrecond = &pre;
	if(fPrecond) delete fPrecond;
	fPrecond = pre.Clone();
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetBiCGStab(const long numiterations, const TPZMatrixSolver<TVar>&pre,const REAL tol,const long FromCurrent){
	ResetSolver();
	fSolver = this->EBICGSTAB;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	//	fPrecond = &pre;
	if(fPrecond) delete fPrecond;
	fPrecond = pre.Clone();
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetJacobi(const long numiterations, const REAL tol, const long FromCurrent) {
	ResetSolver();
	fSolver = this->EJacobi;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	fFromCurrent = FromCurrent;
}
template <class TVar>void TPZStepSolver<TVar>::SetSSOR(const long numiterations,const REAL overrelax,const REAL tol,const long FromCurrent) {
	ResetSolver();
	fSolver = this->ESSOR;
	fOverRelax = overrelax;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fTol = tol;
	fFromCurrent = FromCurrent;
}
template <class TVar>
void TPZStepSolver<TVar>::SetSOR(const long numiterations,const REAL overrelax,const REAL tol,const long FromCurrent){
	ResetSolver();
	fSolver = this->ESOR;
	fMaxIterations = numiterations;
    fNumIterations = -1;
	fOverRelax = overrelax;
	fTol = tol;
	fFromCurrent = FromCurrent;
}
template<class TVar>
void TPZStepSolver<TVar>::SetMultiply() {
	ResetSolver();
	fSolver = this->EMultiply;
}


/*!
 \fn TPZStepSolver::SetPreconditioner(TPZSolver &solve);
 */
template <class TVar>
void TPZStepSolver<TVar>::SetPreconditioner(TPZSolver<TVar> &solve)
{
    if (fSolver == this->EDirect) {
        DebugStop();
    }
	if(fPrecond) delete fPrecond;
	fPrecond = solve.Clone();
}

template <class TVar>
void TPZStepSolver<TVar>::Write(TPZStream &buf, int withclassid)
{
	TPZMatrixSolver<TVar>::Write(buf, withclassid);
    if (fPrecond) {
        fPrecond->Write(buf, 1);
    }
    else {
        int zero = -1;
        buf.Write(&zero );
    }
	int lfSolver = fSolver;
	buf.Write(&lfSolver, 1);
	int lfDT = fDecompose;
	buf.Write(&lfDT, 1);
	buf.Write(&fMaxIterations, 1);
	buf.Write(&fNumVectors, 1);
	buf.Write(&fTol, 1);
	buf.Write(&fOverRelax, 1);
	buf.Write(&fFromCurrent, 1);
	long size = fSingular.size();
	buf.Write(&size, 1);
	std::list<long>::iterator it = fSingular.begin();
	for(;it != fSingular.end(); it++)
	{
		buf.Write(&*it, 1);
	}
}

template <class TVar>
void TPZStepSolver<TVar>::Read(TPZStream &buf, void *context)
{
	TPZMatrixSolver<TVar>::Read(buf, context);
	fPrecond = dynamic_cast<TPZSolver<TVar> *>(TPZSaveable::Restore(buf, context));
	
	int lfSolver = 0;
	buf.Read(&lfSolver, 1);
	fSolver = (typename TPZMatrixSolver<TVar>::MSolver)lfSolver;
	int lfDT = 0;
	buf.Read(&lfDT, 1);
	fDecompose = (DecomposeType)lfDT;
	buf.Read(&fMaxIterations, 1);
	buf.Read(&fNumVectors, 1);
	buf.Read(&fTol, 1);
	buf.Read(&fOverRelax, 1);
	buf.Read(&fFromCurrent, 1);
	long size = 0;
	buf.Read(&size, 1);
	fSingular.resize(size);
	std::list<long>::iterator it = fSingular.begin();
	for(;it != fSingular.end(); it++)
	{
		buf.Read(&*it, 1);
	}
}

/** @brief Serialization methods */
template <>
int TPZStepSolver<float>::ClassId() const
{
    return TPZSTEPSOLVERFLOAT_ID;
}
template <>
int TPZStepSolver<double>::ClassId() const
{
    return TPZSTEPSOLVERDOUBLE_ID;
}
template <>
int TPZStepSolver<long double>::ClassId() const
{
    return TPZSTEPSOLVERLONGDOUBLE_ID;
}

template <>
int TPZStepSolver<std::complex<float> >::ClassId() const
{
    return TPZSTEPSOLVERCOMPLEXFLOAT_ID;
}
template <>
int TPZStepSolver<std::complex<double> >::ClassId() const
{
    return TPZSTEPSOLVERCOMPLEXDOUBLE_ID;
}
template <>
int TPZStepSolver<std::complex<long double> >::ClassId() const
{
    return TPZSTEPSOLVERCOMPLEXLONGDOUBLE_ID;
}


template class TPZStepSolver<float>;
template class TPZStepSolver<double>;
template class TPZStepSolver<long double>;

template class TPZStepSolver<std::complex<float> >;
template class TPZStepSolver<std::complex<double> >;
template class TPZStepSolver<std::complex<long double> >;

#ifndef BORLAND
template class TPZRestoreClass< TPZStepSolver<float>, TPZSTEPSOLVERFLOAT_ID>;
template class TPZRestoreClass< TPZStepSolver<double>, TPZSTEPSOLVERDOUBLE_ID>;
template class TPZRestoreClass< TPZStepSolver<double>, TPZSTEPSOLVERCOMPLEXDOUBLE_ID>;
#endif
