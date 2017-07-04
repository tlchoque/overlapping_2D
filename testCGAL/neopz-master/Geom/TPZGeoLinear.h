/**
 * @file
 * @brief Contains the TPZGeoLinear class which implements the geometry of a one dimensional linear element.
 */
#ifndef TPZGEOLINEARH
#define TPZGEOLINEARH

#include "pznoderep.h"

#include "pzvec.h"
#include "pzeltype.h"
#include "pzfmatrix.h"
#include "tpzline.h"

#include <string>

class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of a one dimensional linear element. \ref geometry "Geometry"
	 */
	class TPZGeoLinear : public TPZNodeRep<2, pztopology::TPZLine> {
		
	public:
		/** @brief Number of corner nodes */
		enum {NNodes = 2};
		
		/** @brief Constructor with list of nodes */
		TPZGeoLinear(TPZVec<long> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZLine>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoLinear() : TPZNodeRep<NNodes, pztopology::TPZLine>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoLinear(const TPZGeoLinear &cp, std::map<long,long> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoLinear(const TPZGeoLinear &cp) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoLinear(const TPZGeoLinear &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
		{
		}
        
        /** @brief answer if the element side is a linear map */
        static bool IsLinearMapping(int side)
        {
            return true;
        }
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Linear";}
		
		/* @brief Computes the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
        
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            for (int i=0; i<3; i++) {
                gradx(i,0) = (coord(i,1)-coord(i,0))*0.5;
            }
        }
		
        /* @brief Computes the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
		static void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);
        
        template<class T>
        static void GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc,TPZVec<T> &result);
		
		static void Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
							 TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
		
		
	};
	
	// VERSAO ORIGINAL
	//   TPZFNMatrix<9> phi(NNodes, pztopology::TPZLine,1);
	//   TPZFNMatrix<18> dphi(1,NNodes, pztopology::TPZLine);
	//   Shape(param,phi,dphi);
	
	//   int ic;
	//   TPZManVector<REAL,3> v1(3,0.);
	//   REAL mod1 = 0.;
	
	//   for(int i=0; i < NNodes, pztopology::TPZLine; i++) {
	//     for(ic = 0; ic < 3; ic++) {
	//       v1[ic] += coord(ic,i)*dphi(0,i);
	//     }
	//   }
	
	//   for(ic=0; ic<3; ic++) {
	//     mod1 += v1[ic]*v1[ic];
	//   }
	//   mod1 = sqrt(mod1);
	//   jacobian(0,0) = mod1;
	//   detjac = mod1;
	//   jacinv(0,0) = 1./detjac;
	
	
	//  axes.Zero();
	//   for(ic=0; ic<3; ic++) {
	//     axes(0,ic) = v1[ic]/mod1;
	//   }
	
	inline void TPZGeoLinear::X(const TPZFMatrix<REAL> &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result){
		
		int ic;
		REAL xi = loc[0];
		int nrow = coord.Rows();
		for(ic=0; ic<nrow; ic++)
        {
            result[ic] = coord.GetVal(ic,0)*(1.-xi)*0.5+coord.GetVal(ic,1)*(1.+xi)*0.5;
        }
		
		
	}
	
};

#endif
