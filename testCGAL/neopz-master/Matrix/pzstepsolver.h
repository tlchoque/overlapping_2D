/**
 * @file
 * @brief Contains TPZStepSolver class which defines step solvers class.
 */

#ifndef TPZSTEPSOLVER_H
#define TPZSTEPSOLVER_H
#include "pzsolve.h"

#include "pzfmatrix.h"

#include "pzstream.h"

#include <list>

/**
 * @brief Defines step solvers class. \ref solver "Solver"
 * @ingroup solver
 */
template<class TVar>
class TPZStepSolver: public TPZMatrixSolver<TVar>
{
public:
	TPZStepSolver(TPZAutoPointer<TPZMatrix<TVar> > refmat = 0);
	
	TPZStepSolver(const TPZStepSolver<TVar> & copy);
	
	virtual ~TPZStepSolver();
	
	void SetSOR(const long numiterations, const REAL overrelax, const REAL tol, const long FromCurrent);
	
	void SetSSOR(const long numiterations, const REAL overrelax, const REAL tol, const long FromCurrent);
	
	void SetJacobi(const long numiterations, const REAL tol, const long FromCurrent);
	
	void SetCG(const long numiterations, const TPZMatrixSolver<TVar> &pre, const REAL tol, const long FromCurrent);
	
	void SetGMRES(const long numiterations, const int numvectors, const TPZMatrixSolver<TVar> &pre, const REAL tol, const long FromCurrent);
	
	void SetBiCGStab(const long numiterations, const TPZMatrixSolver<TVar> &pre, const REAL tol, const long FromCurrent);
	
	void SetDirect(const DecomposeType decomp);
	
	void SetMultiply();
	
	virtual TPZSolver<TVar> *Clone() const
	{
		return new TPZStepSolver<TVar>(*this);
	}
	
	void SetTolerance(REAL tol)
	{
		fTol = tol;
	}
    
    /** @brief return the value of tolerance from the solver */
    REAL GetTolerance() const
    {
        return fTol;
    }
	
    /** @brief reset the data structure of the solver object */
	void ResetSolver();
    
    virtual typename TPZMatrixSolver<TVar>::MSolver Solver()
    {
        return fSolver;
    }
	
	/** @brief returns the equations for which the equations had zero pivot */
	std::list<long> &Singular()
	{
		return fSingular;
	}
	
	/** @brief This method will reset the matrix associated with the solver */
	/** This is useful when the matrix needs to be recomputed in a non linear problem */
	virtual void ResetMatrix();

	/** @brief Updates the values of the current matrix based on the values of the matrix */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > matrix)
	{
		if (fPrecond)
			fPrecond->UpdateFrom(matrix);
		TPZMatrixSolver<TVar>::UpdateFrom(matrix);
	}
	
    
	void Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual = 0);
    
    /** @brief Decompose the system of equations if a direct solver is used */
    virtual void Decompose();
    
    /** @brief Define the preconditioner as a solver object */
	void SetPreconditioner(TPZSolver<TVar> &solve);
    
    /** @brief Number of iterations of last solve */
    int NumIterations()
    {
        return fNumIterations;
    }
    
    /** @brief access method to the preconditioner */
    TPZSolver<TVar> *PreConditioner()
    {
        return fPrecond;
    }
	
	/** @brief Serialization methods */
	virtual int ClassId() const;
	virtual void Write(TPZStream &buf, int withclassid);
	virtual void Read(TPZStream &buf, void *context);
	
	
private:
	typename TPZMatrixSolver<TVar>::MSolver fSolver;
	DecomposeType fDecompose;
    
    /// Maximum number of iterations
    long fMaxIterations;
    
    /// Number of iterations of last solve
	long fNumIterations;
	int fNumVectors;
	REAL fTol;
	REAL fOverRelax;
	
	/** @brief Solver using preconditioner matrix */
	TPZSolver<TVar> *fPrecond;
	long fFromCurrent;
	
	std::list<long> fSingular;
};

#endif
