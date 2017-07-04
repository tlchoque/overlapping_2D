/**
 * @file
 * @brief Contains the TPZGeoTriangle class which implements the geometry of a triangle element.
 */

#ifndef TPZGEOTRIANGLEH
#define TPZGEOTRIANGLEH

#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"
#include "tpztriangle.h"
#include "pzfmatrix.h"

#include <string>
#include <map>

class TPZGeoEl;
class TPZGeoMesh;


namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of a triangle element. \ref geometry "Geometry"
	 */
	class TPZGeoTriangle : public TPZNodeRep<3, pztopology::TPZTriangle> 
	{
	public:

		/** @brief Number of corner nodes */
		enum {NNodes = 3};
		
		/** @brief Constructor with list of nodes */
		TPZGeoTriangle(TPZVec<long> &nodeindexes) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoTriangle() : TPZNodeRep<NNodes,pztopology::TPZTriangle>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoTriangle(const TPZGeoTriangle &cp,
					   std::map<long,long> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTriangle(const TPZGeoTriangle &cp) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTriangle(const TPZGeoTriangle &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp)
		{
		}
        
        static bool IsLinearMapping(int side)
        {
            return true;
        }
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Triangle";}
		
		/** @brief Implementation of two-dimensional bilinear interpolation*/
		static  void Shape(TPZVec<REAL> &x,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		/* brief compute the vectors for defining an HDiv approximation space */
        void VecHdiv(const TPZGeoEl &gel,TPZFMatrix<REAL> &NormalVec,TPZVec<int> & VectorSide) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            VecHdiv(coord,NormalVec,VectorSide);
        }
		/* brief compute the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
        
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            DebugStop();
        }
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
		/** @brief Computes the jacobian*/
		static  void Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		/** @brief Computes the geometric location*/
		static  void X(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);
        
        template<class T>
        static void GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc,TPZVec<T> &result);
        
		/** @brief Implementation of Hdiv space*/
		static	void ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result);
		
		static	void VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result);
		
		static void VecHdiv(TPZFMatrix<REAL> & coord, TPZFMatrix<REAL> & fNormalVec,TPZVec<int> &sidevector);
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, \n
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
	protected:
		/**
		 * @brief This method apply an infinitesimal displacement in some points
		 * to fix singularity problems when using MapToSide() method! \n
		 * This points are CornerNodes, when projected in the opposing side
		 */
		static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes, int matid, long& index);

	};
	
};

#endif 
