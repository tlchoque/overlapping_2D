/**
 * @file
 * @brief Contains the declaration of the TPZCondensedCompEl class, which implements an computational element which condenses the internal connects.
 */

#ifndef TPZCONDENSEDCOMPELH
#define TPZCONDENSEDCOMPELH

#include "pzcompel.h"
#include "pzmatred.h"
#include "pzmanvector.h"
#include "pzelmat.h"

#ifdef USING_BLAS
#define USING_DGER
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#elif USING_MKL
#include <mkl.h>
#else
#include "cblas.h"
//#define USING_DGER
#endif
#endif

#ifdef USING_MKL
#include <mkl.h>
//#include <omp.h>
//#define USING_DGER
#endif


/**
 * @brief Class which implements an element which condenses the internal connects
 * @author Philippe Devloo
 * @ingroup CompElement
 */
class TPZCondensedCompEl : public TPZCompEl
{

    //TPZMatRed<REAL, TPZFMatrix<REAL> > fCondensed;
	TPZMatRed<STATE, TPZFMatrix<STATE> > fCondensed;
    TPZCompEl *fReferenceCompEl;
    TPZManVector<long,27> fIndexes; 
    
    void Resequence();

public:
    
    TPZCondensedCompEl(TPZCompEl *ref);
    
    /** @brief create a copy of the condensed computational element in the other mesh */
    TPZCondensedCompEl(const TPZCondensedCompEl &copy, TPZCompMesh &mesh);
    
    virtual ~TPZCondensedCompEl();

    /**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const;

    /** @brief unwrap the condensed element from the computational element and delete the condensed element */
    void Unwrap();
    
    /**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, long index);
    
    /** @brief Returns the number of nodes of the element */
	virtual int NConnects() const 
    {
        return fReferenceCompEl->NConnects();
    }
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual long ConnectIndex(int i) const 
    {
        return fReferenceCompEl->ConnectIndex(fIndexes[i]);
    }

    TPZCompEl * ReferenceCompEl(){
        return fReferenceCompEl;
    }
    
    virtual void LoadElementReference()
    {
        if(fReferenceCompEl)
        {
            fReferenceCompEl->LoadElementReference();
        }
    }

    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<long> &connectindexes) const;
	
	/** @brief Dimension of the element */
	virtual int Dimension() const 
    {
        return fReferenceCompEl->Dimension();
    }
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const 
    {
        return new TPZCondensedCompEl(*this,mesh);
    }

    /** @brief Loads the solution within the internal data structure of the element */ 
	/** 
	 * Is used to initialize the solution of connect objects with dependency \n
	 * Is also used to load the solution within SuperElements
	 */
	virtual void LoadSolution();

    virtual void TransferMultiphysicsElementSolution()
    {
        if(fReferenceCompEl)
        {
            fReferenceCompEl->TransferMultiphysicsElementSolution();
        }
    }

	/**
	 * @brief Method for creating a copy of the element in a patch mesh
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the computational elements
     */
	/**
	 * Otherwise of the previous clone function, this method don't
	 * copy entire mesh. Therefore it needs to map the connect index
	 * from the both meshes - original and patch
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<long,long> & gl2lcConMap,
									std::map<long,long> & gl2lcElMap) const;

private:
    /**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi master element coordinate
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 * @param axes axes associated with the derivative of the solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
    
public:
    /**
	 * @brief Computes solution and its derivatives in the local coordinate qsi.
	 * @param qsi [in] master element coordinate
     * @param data [in] stores all input data
	 */

    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
    
	/**
	 * @brief Computes solution and its derivatives in the local coordinate qsi. \n
	 * This method will function for both volumetric and interface elements
	 * @param qsi master element coordinate of the interface element
	 * @param normal vector
	 * @param leftsol finite element solution
	 * @param dleftsol solution derivatives
	 * @param leftaxes axes associated with the left solution
	 * @param rightsol finite element solution
	 * @param drightsol solution derivatives
	 * @param rightaxes axes associated with the right solution
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi,
								 TPZVec<REAL> &normal,
								 TPZSolVec &leftsol, TPZGradSolVec &dleftsol,TPZFMatrix<REAL> &leftaxes,
								 TPZSolVec &rightsol, TPZGradSolVec &drightsol,TPZFMatrix<REAL> &rightaxes);
	
	/**
	 * @brief Computes solution and its derivatives in local coordinate qsi
	 * @param qsi master element coordinate
	 * @param phi matrix containing shape functions compute in qsi point
	 * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
	 * @param axes [in] axes indicating the direction of the derivatives
	 * @param sol finite element solution
	 * @param dsol solution derivatives
	 */
	virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
								 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);

	/**
	 * @brief Computes the element stifness matrix and right hand side
	 * @param ek element stiffness matrix
	 * @param ef element load vector
	 */
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
	
	
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	virtual void CalcResidual(TPZElementMatrix &ef);
    
    void EvaluateError(void (* fp)(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv),
                                  TPZVec<REAL> &errors,TPZBlock<REAL> * flux) {
        fReferenceCompEl->EvaluateError(fp, errors, flux);
    }
    
    /**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param graphmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension)
    {
        fReferenceCompEl->CreateGraphicalElement(graphmesh, dimension);
    }



};

#endif
