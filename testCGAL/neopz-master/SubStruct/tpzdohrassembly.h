/**
 * @file
 * @brief Contains the TPZDohrAssembly class which implements assembling using Dohrmann algorithm.
 */

#ifndef TPZDOHRASSEMBLYH
#define TPZDOHRASSEMBLYH

#include "pzvec.h"
#include "pzsave.h"
#include "pzfmatrix.h"

/**
 * @ingroup substructure
 * @brief Assembling using Dohrmann algorithm. \ref substructure "Sub structure"
 * @author Philippe Devloo
 * @since 04/03/2009
 */
template<class TVar>
class TPZDohrAssembly
// @TODO Implement the methods to make the class actually saveable
{
public:
	/** @brief For each substructure the equation numbering of the substructures */
	/** The order of the equations follows the ordering of the connects */
	TPZVec< TPZVec< int > > fFineEqs;
	
	/** @brief For each substructure the equation numbering of the coarse equations */
	TPZVec< TPZVec< int > > fCoarseEqs;
	
	/** @brief Sum the values in the local matrix into the global matrix */
	void Assemble(int isub, const TPZFMatrix<TVar> &local, TPZFMatrix<TVar> &global);
	
	/** @brief Extract the values from the global matrix into the local matrix */
	void Extract(int isub, const TPZFMatrix<TVar> &global, TPZFMatrix<TVar> &local);
	
	/** @brief Sum the values in the local matrix into the global matrix */
	void AssembleCoarse(int isub, const TPZFMatrix<TVar> &local, TPZFMatrix<TVar> &global);
	
	/** @brief Extract the values from the global matrix into the local matrix */
	void ExtractCoarse(int isub, const TPZFMatrix<TVar> &global, TPZFMatrix<TVar> &local);
    
    /** @brief method for streaming the object to a stream */
    void Write(TPZStream &out);
    
    /** @brief method for reading the object for a stream */
    void Read(TPZStream &input);
};

#endif
