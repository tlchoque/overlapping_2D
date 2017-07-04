/**
 * @file
 * @brief Contains the TPZSymetricSpStructMatrix class which implements sparse structural matrices.
 */

#ifndef TPZSymetricSpStructMatrix_H
#define TPZSymetricSpStructMatrix_H

#include "pzstrmatrix.h"
#include "pzysmp.h"

#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzelmat.h"

#include "pzmatrix.h"
#include "pzfmatrix.h"

/**
 * @brief Implements Sparse Structural Matrices. \ref structural "Structural Matrix"
 * @ingroup structural
 */
class TPZSymetricSpStructMatrix : public TPZStructMatrix {
public:    
	
    TPZSymetricSpStructMatrix(TPZCompMesh *);
	
    virtual TPZMatrix<STATE> * Create();
	
	virtual TPZMatrix<STATE> * CreateAssemble(TPZFMatrix<STATE> &rhs, TPZAutoPointer<TPZGuiInterface> guiInterface);
	
    virtual TPZStructMatrix * Clone(); 	
	
    /** Used only for testing */
	static int main();
	
};

#endif //TPZSymetricSpStructMatrix_H
