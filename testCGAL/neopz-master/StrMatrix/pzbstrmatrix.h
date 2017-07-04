/**
 * @file
 * @brief Contains the TPZBandStructMatrix class which implements Banded Structural Matrices.
 */

#ifndef TPZBANDSTRUCTMATRIX_H
#define TPZBANDSTRUCTMATRIX_H

#include "pzmatrix.h"
#include "pzstrmatrix.h"
#include "pzfmatrix.h"
#include "pzcmesh.h"

/**
 * @brief Implements Banded Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZBandStructMatrix : public TPZStructMatrix {
public:    
	
    TPZBandStructMatrix(TPZCompMesh *);
    
    ~TPZBandStructMatrix();
	
    TPZBandStructMatrix(const TPZBandStructMatrix &copy) : TPZStructMatrix(copy)
    {
    }
	
    virtual TPZMatrix<STATE> * Create();
	
    virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone();
	
public:
	
};

#endif //TPZBANDSTRUCTMATRIX_H
