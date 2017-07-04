/**
 * @file
 * @brief Contains the TPZAbstractFrontMatrix class which implements a matrix stored in a frontal decomposition scheme.
 */

#ifndef TPZFRONTMATRIX_H
#define TPZFRONTMATRIX_H


#include "TPZFront.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontSym.h"


#include "TPZStackEqnStorage.h"
#include "pzvec.h"
#include "TPZFileEqnStorage.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzlog.h"

/**
 * \addtogroup frontal
 * @{
 */
/**
 @brief Implements a matrix stored in a frontal decomposition scheme. \ref frontal "Frontal"
 */
template<class TVar>
class TPZAbstractFrontMatrix : public TPZMatrix<TVar>
{
public:
	
	TPZAbstractFrontMatrix() : TPZMatrix<TVar>()
	{
	}
	
	TPZAbstractFrontMatrix(long ieq, long jeq) : TPZMatrix<TVar>(ieq,jeq)
	{
	}
	
	virtual TPZFront<TVar> & GetFront() = 0;
	
};

/**
 * @brief Responsible for the frontal method as a whole. \ref frontal "Frontal"
 */
/** 
 * Manages the remaining classes connecting them
 */
template <class TVar, class store, class front>
class TPZFrontMatrix : public TPZAbstractFrontMatrix<TVar> {
protected:
	
    /** @brief Indicates storage schema. Assumes values TPZFileEqnStorage for binary file and TPZStackEqnStorage for a TPZStack storage */
	store fStorage;
	
    /** @brief Indicates Front matrix type. Assumes values TPZFrontSym for symmetric front and TPZFrontNonSym for non symmetric front matrix */
	front fFront;
public:
 	int Work();
	/** @brief Finishes writing of a binary file and closes it*/
	void FinishWriting();
	/** @brief Reopen the binary file*/
	void ReOpen();
	
	/** @brief Reinitialize the structure of the frontal matrix */
	virtual int Zero();
    /** @brief Allocates data for the FrontMatrix */
	void AllocData();
	
	/** @brief returns a pointer to the front matrix */
	TPZFront<TVar> &GetFront() { return fFront;}
    /** @brief Checks if FrontMatrix needs a compression,
	 if so calls Compress method */
	void CheckCompress();
    /** Static main for testing */
	static void main();
    /** @brief Prints a FrontMatrix object */
	void Print(const char * name, std::ostream & out ,const MatrixOutputFormat form = EFormatted) const;
    /** @brief Simple Destructor */
    ~TPZFrontMatrix();
    /** @brief Simple Constructor */
    TPZFrontMatrix();
    /** 
	 * @brief Constructor with a globalsize parameter 
	 * @param globalsize Indicates initial global size
	 */
	TPZFrontMatrix(long globalsize);
	
	TPZFrontMatrix(const TPZFrontMatrix &cp) : TPZAbstractFrontMatrix<TVar>(cp), fStorage(cp.fStorage),
	fFront(cp.fFront),fNumEq(cp.fNumEq),fLastDecomposed(cp.fLastDecomposed), fNumElConnected(cp.fNumElConnected),fNumElConnectedBackup(cp.fNumElConnectedBackup)
    {
    }
    
  //  CLONEDEF(TPZFrontMatrix)
	virtual TPZMatrix<TVar>*Clone() const { return new TPZFrontMatrix(*this); }
    /** 
	 * @brief Sends a message to decompose equations from lower_eq to upper_eq, according to destination index
	 * @param destinationindex Contains destination indexes of equations
	 * @param lower_eq Starting index
	 * @param upper_eq Finishing index
	 */
    void EquationsToDecompose(TPZVec<long> &destinationindex, long &lower_eq, long &upper_eq);
	
	
    /** Add a matrix to the frontal matrix */
	//    void AssembleMatrix(TPZVec < int > & eqnumbers, TPZFMatrix<REAL> & ek, TPZFMatrix<REAL> & ef);
	
    /** 
	 * @brief Initializes the number of elements connected to each equation 
	 * @param numelconnected Indicates number of elements connected to that equation
	 */
    void SetNumElConnected(TPZVec < int > &numelconnected);
	
    /** 
	 * @brief Add a contribution of a stiffness matrix putting it on destination indexes position 
	 * @param elmat Indicates number of elements connected to that equation
	 * @param destinationindex Positioning of such members on global stiffness matrix
	 */
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec < long > & destinationindex);
	
    /** 
	 * @brief Add a contribution of a stiffness matrix using the indexes to compute the frontwidth. It does it symbolicaly
	 * @param destinationindex Array containing destination indexes.
	 */
    void SymbolicAddKel(TPZVec < long > & destinationindex);
	
    /** 
	 * @brief Add a contribution of a stiffness matrix 
	 * @param elmat Member stiffness matrix beeing added
	 * @param sourceindex Source position of values on member stiffness matrix
	 * @param destinationindex Positioning of such members on global stiffness matrix
	 */
    virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec < long > & sourceindex, TPZVec < long > & destinationindex);
	
    virtual int SolveDirect( TPZFMatrix<TVar> &B , DecomposeType dt, std::list<long> &singular);
    /**
	 * @brief Forward substitution and result is on b
	 * @param b Result of the substitution
	 */
	int Subst_Forward(TPZFMatrix<TVar> *b) const;
	
    /** @brief Backward substitution and result is on b*/
	int Subst_Backward(TPZFMatrix<TVar> *b) const;
    /** @brief Executes a substitution on a TPZFMatrix<REAL> object
	 applies both forward and backward substitution automaticaly */
	int Substitution(TPZFMatrix<TVar> *) const;
    /*
	 void SetFileName(
	 const char *name = SetTempFileName() //! Name of the file
	 );
	 */
	
    /** 
	 * @brief Sets a file name for generic input or output 
	 * @param option It can be either 'w' or 'r' for writeing or reading respectively
	 * @param name Name of the file
	 */
    void SetFileName(char option, const char *name);
	
protected:
    /** @brief Indicates number of equations */
	long fNumEq;
    /** @brief Indicates last decomposed equation */
	long fLastDecomposed;
	
    /** \ link aggregationByValue */
	//    TPZFront fFront;
	
    /** \ link aggregationByValue */
    //TPZStackEqnStorage fStorage;
	
    /** @brief Contains the number of elements which still need to contribute to a given equation */
    TPZVec <int> fNumElConnected;
	
    /** @brief Contains the number of elements which still need to contribute to a given equation. */
	/** 
	 * This backup copy is needed to be able to reinitialize the matrix through the Zero() method
	 */
    TPZVec <int> fNumElConnectedBackup;
	
    /**
     * @ link aggregationByValue
     */
    //TPZFileEqnStorage fFileEqnStorage;
};

/** @} */

#endif //TPZFRONTMATRIX_H
