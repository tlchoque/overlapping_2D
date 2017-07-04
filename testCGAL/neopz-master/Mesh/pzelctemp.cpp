/**
 * @file
 * @brief Contains the implementation of the TPZIntelGen methods.
 */

#include "pzelctemp.h"
#include "pzquad.h"
#include "pzgeoel.h"
#include "pzmaterial.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzintelgen"));
#endif

template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, long &index) :
TPZInterpolatedElement(mesh,gel,index), fConnectIndexes(TSHAPE::NSides,-1) {

	for(int i=0; i<TSHAPE::NSides; i++) fConnectIndexes[i]=-1;
	//  RemoveSideRestraintsII(EInsert);
	gel->SetReference(this);
	for(int i=0;i<TSHAPE::NSides;i++) {
		fConnectIndexes[i] = CreateMidSideConnect(i);
		mesh.ConnectVec()[fConnectIndexes[i]].IncrementElConnected();
//		IdentifySideOrder(i);
	}
	
    AdjustIntegrationRule();

	
}

template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, long &index, int nocreate) :
TPZInterpolatedElement(mesh,gel,index),fConnectIndexes(TSHAPE::NSides,-1)
{
	//int ic;
	//for(ic=0; ic<TSHAPE::NSides; ic++)
//	{
//		fConnectIndexes[ic] = -1;
//	}
	fPreferredOrder = -1;
}

template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh, const TPZIntelGen<TSHAPE> &copy) :
TPZInterpolatedElement(mesh,copy), fConnectIndexes(copy.fConnectIndexes),fIntRule(copy.fIntRule) {
	fPreferredOrder = copy.fPreferredOrder;
}


template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen(TPZCompMesh &mesh,
								 const TPZIntelGen<TSHAPE> &copy,
								 std::map<long,long> & gl2lcConMap,
								 std::map<long,long> & gl2lcElMap) :
TPZInterpolatedElement(mesh,copy,gl2lcElMap), fConnectIndexes(TSHAPE::NSides,-1), fIntRule(copy.fIntRule)
{
	
	fPreferredOrder = copy.fPreferredOrder;
	int i;
	for(i=0;i<TSHAPE::NSides;i++)
	{
		long lcIdx = -1;
		long glIdx = copy.fConnectIndexes[i];
		if (gl2lcConMap.find(glIdx) != gl2lcConMap.end()) lcIdx = gl2lcConMap[glIdx];
		else
		{
			std::stringstream sout;
			sout << "ERROR in : " << __PRETTY_FUNCTION__
			<< " trying to clone the connect index: " << glIdx
			<< " wich is not in mapped connect indexes!";
			LOGPZ_ERROR(logger, sout.str().c_str());
			fConnectIndexes[i] = -1;
			return;
		}
		fConnectIndexes[i] = lcIdx;
	}
}


template<class TSHAPE>
TPZIntelGen<TSHAPE>::TPZIntelGen() :
TPZInterpolatedElement(), fConnectIndexes(TSHAPE::NSides,-1), fIntRule() {
	fPreferredOrder = -1;
	int i;
	for(i=0;i<TSHAPE::NSides;i++) {
		fConnectIndexes[i] = -1;
	}
}

template<class TSHAPE>
TPZIntelGen<TSHAPE>::~TPZIntelGen(){
    TPZGeoEl *gel = Reference();
	TPZCompEl *cel = gel->Reference();
	if(gel) {
		if(cel == this) {
			RemoveSideRestraintsII(EDelete);
		}
		Reference()->ResetReference();
	}
    TPZStack<long > connectlist;
    BuildConnectList(connectlist);
    long nconnects = connectlist.size();
    for (int ic=0; ic<nconnects ; ic++) {
        fMesh->ConnectVec()[connectlist[ic]].DecrementElConnected();
    }
}

template<class TSHAPE>
MElementType TPZIntelGen<TSHAPE>::Type() {
	return TSHAPE::Type();
}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetConnectIndex(int i, long connectindex){
#ifndef NODEBUG
	if(i<0 || i>= TSHAPE::NSides) {
		std::cout << " TPZIntelGen<TSHAPE>::SetConnectIndex index " << i <<
		" out of range\n";
		return;
	}
#endif
	fConnectIndexes[i] = connectindex;
}

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::NConnectShapeF(int connect, int order) const{
    
    if(connect < TSHAPE::NCornerNodes) return TSHAPE::NConnectShapeF(connect,0);
	if(order < 0) return 0;
    int nshape = TSHAPE::NConnectShapeF(connect, order);
#ifdef PZDEBUG
    if(nshape < 0 )
    {
        nshape = TSHAPE::NConnectShapeF(connect, order);
        DebugStop();
    }
#endif
	return nshape;
}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetIntegrationRule(int ord) {
	TPZManVector<int,3> order(TSHAPE::Dimension,ord);
	fIntRule.SetOrder(order);
}

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::NSideConnects(int side) const {
	return TSHAPE::NContainedSides(side);
}

template<class TSHAPE>
int TPZIntelGen<TSHAPE>::SideConnectLocId(int node, int side) const {
	return TSHAPE::ContainedSideLocId(side,node);
}

/**Sets the interpolation order for the interior of the element*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetInterpolationOrder(int order) {
	fPreferredOrder = order;
}

/**Identifies the interpolation order on the interior of the element*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::GetInterpolationOrder(TPZVec<int> &ord) {
	ord.Resize(TSHAPE::NSides-TSHAPE::NCornerNodes);
	int i;
	for(i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
		ord[i] = Connect(i+TSHAPE::NCornerNodes).Order();
	}
}

/**return the preferred order of the polynomial along side iside*/
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::PreferredSideOrder(int side) {
	if(side < TSHAPE::NCornerNodes) return 0;
	if(side<TSHAPE::NSides) {
		int order =fPreferredOrder;
		return AdjustPreferredSideOrder(side,order);
	}
	PZError << "TPZIntelgen::PreferredSideOrder called for side = " << side << "\n";
	return 0;
	
}

template<class TSHAPE>
long TPZIntelGen<TSHAPE>::ConnectIndex(int con) const{
	
#ifndef NODEBUG
	if(con<0 || con>= NConnects()) {
		std::cout << "TPZIntelgen::ConnectIndex wrong parameter con " << con <<
		" NSides " << TSHAPE::NSides << " NConnects " << NConnects() << std::endl;
		DebugStop();
	}
	
#endif
	return fConnectIndexes[con];
}



/**Sets the preferred interpolation order along a side
 This method only updates the datastructure of the element
 In order to change the interpolation order of an element, use the method PRefine
 */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetPreferredOrder(int order)
{
	fPreferredOrder = order;
}

/**sets the interpolation order of side to order*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SetSideOrder(int side, int order) {
	if(side<0 || side >= TSHAPE::NSides || (side >= TSHAPE::NCornerNodes && order <1)) {
		PZError << "TPZIntelGen::SetSideOrder. Bad paramenter side " << side << " order " << order << std::endl;
		DebugStop();
#ifdef LOG4CXX
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Bad side or order " << side << " order " << order;
		LOGPZ_ERROR(logger,sout.str())
#endif
		return;
	}
	if(side>= TSHAPE::NCornerNodes) {
		if(fConnectIndexes[side] == -1) return;
        int prevorder = fPreferredOrder;
        if (ConnectIndex(TSHAPE::NSides-1) != -1) {
            prevorder = EffectiveSideOrder(TSHAPE::NSides-1);
        }
		TPZConnect &c = Connect(side);
        int previousconnectorder = c.Order();
        if (order != previousconnectorder)
        {
            c.SetOrder(order,ConnectIndex(side));
            long seqnum = c.SequenceNumber();
            int nvar = 1;
            TPZMaterial * mat = Material();
            if(mat) nvar = mat->NStateVariables();
            int nshape = TSHAPE::NConnectShapeF(side, order);
            c.SetNShape(nshape);
            c.SetNState(nvar);
            Mesh()->Block().Set(seqnum,nshape*nvar);
            TPZGeoElSide gelside(Reference(),side);
            TPZStack<TPZCompElSide> equal;
            gelside.EqualLevelCompElementList(equal, 1, 0);
            long neq = equal.size();
            for (long eq=0; eq<neq; eq++) {
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(equal[eq].Element());
                intel->AdjustIntegrationRule();
            }
        }
        
        int neworder = prevorder;
        if (ConnectIndex(TSHAPE::NSides-1) != -1) {
            neworder = EffectiveSideOrder(TSHAPE::NSides-1);
        }
		if(neworder != prevorder) {
            AdjustIntegrationRule();
		}
	}
}

/**returns the actual interpolation order of the polynomial along the side*/
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::EffectiveSideOrder(int side) const {
	if(side < TSHAPE::NCornerNodes || side >= TSHAPE::NSides) return 0;
	if(fConnectIndexes[side] == -1)
	{
		std::stringstream sout ;
		sout << __PRETTY_FUNCTION__ << " side " << side << std::endl;
		//Print(sout);
#ifdef LOG4CXX
		LOGPZ_ERROR(logger,sout.str());
#else
		std::cout << sout.str() << std::endl;
#endif
		DebugStop();
		return -1;
	}
    TPZStack<int> lowdim;
    Reference()->LowerDimensionSides(side,lowdim);
	TPZConnect &c = Connect(side);
	int order = c.Order();
    for (int is=0; is<lowdim.size(); is++) {
        TPZConnect &c = MidSideConnect(lowdim[is]);
        if(c.Order() > order) order = c.Order();
    }
    return order;
}
/**returns the actual interpolation order of the polynomial for a connect*/
template<class TSHAPE>
int TPZIntelGen<TSHAPE>::ConnectOrder(int connect) const {
	return Connect(connect).Order();
}

/**compute the values of the shape function of the side*/
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
	
	int nc = TSHAPE::NContainedSides(side);
	int nn = TSHAPE::NSideNodes(side);
	TPZManVector<long,27> id(nn);
	TPZManVector<int,27> order(nc-nn);
	int n,c;
	TPZGeoEl *ref = Reference();
	for (n=0;n<nn;n++){
		int nodloc = TSHAPE::SideNodeLocId(side,n);
		id [n] = ref->NodePtr(nodloc)->Id();
	}
	for (c=nn;c<nc;c++){
		int conloc = TSHAPE::ContainedSideLocId(side,c);
        order[c-nn] = Connect(conloc).Order();
	}
	TSHAPE::SideShape(side, point, id, order, phi, dphi);
}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
	TPZManVector<long,TSHAPE::NCornerNodes> id(TSHAPE::NCornerNodes,0);
	TPZManVector<int, TSHAPE::NSides-TSHAPE::NCornerNodes+1> ord(TSHAPE::NSides-TSHAPE::NCornerNodes,0);
	int i;
	TPZGeoEl *ref = Reference();
	for(i=0; i<TSHAPE::NCornerNodes; i++) {
		id[i] = ref->NodePtr(i)->Id();
	}
	for(i=0; i<TSHAPE::NSides-TSHAPE::NCornerNodes; i++) {
		ord[i] = Connect(i+TSHAPE::NCornerNodes).Order();
	}
	TSHAPE::Shape(pt,id,ord,phi,dphi);
}

/** Returns the transformation which transform a point from the side to the interior of the element */
template<class TSHAPE>
TPZTransform TPZIntelGen<TSHAPE>::TransformSideToElement(int side) {
	return TSHAPE::TransformSideToElement(side);
}

/** Save the element data to a stream */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Write(TPZStream &buf, int withclassid)
{
	TPZInterpolatedElement::Write(buf,withclassid);
	TPZManVector<int,3> order(3,0);
	fIntRule.GetOrder(order);
	WriteObjects(buf,order);
	buf.Write(fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Write(&fPreferredOrder,1);
	int classid = this->ClassId();
	buf.Write ( &classid, 1 );
}

/** Read the element data from a stream */
template<class TSHAPE>
void TPZIntelGen<TSHAPE>::Read(TPZStream &buf, void *context)
{
	TPZInterpolatedElement::Read(buf,context);
	TPZManVector<int,3> order;
	ReadObjects(buf,order);
	fIntRule.SetOrder(order);
	buf.Read(fConnectIndexes.begin(),TSHAPE::NSides);
	buf.Read(&fPreferredOrder,1);
	int classid = -1;
	buf.Read( &classid, 1 );
	if ( classid != this->ClassId() )
	{
		std::stringstream sout;
		sout << "ERROR - " << __PRETTY_FUNCTION__
        << " trying to restore an object id " << this->ClassId() << " for an package of id = " << classid;
		LOGPZ_ERROR ( logger, sout.str().c_str() );
	}
}

using namespace pzshape;

#include "pzshapepoint.h"
template<>
void TPZIntelGen<pzshape::TPZShapePoint>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == 0) std::cout << "A point element has no graphical representation\n";
}

template<class TSHAPE>
void TPZIntelGen<TSHAPE>::CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) {
	if(dimension == TSHAPE::Dimension /* && Material()->Id() > 0 */) {
		new typename TSHAPE::GraphElType(this,&grafgrid);
	}
}

template<>
int TPZIntelGen<pzshape::TPZShapePoint>::ClassId() const
{
	return TPZINTELPOINTID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapePoint>, TPZINTELPOINTID>;
#endif

#include "pzshapelinear.h"
template<>
int TPZIntelGen<pzshape::TPZShapeLinear>::ClassId() const
{
	return TPZINTELLINEARID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapeLinear>, TPZINTELLINEARID>;
#endif

#include "pzshapetriang.h"
template<>
int TPZIntelGen<pzshape::TPZShapeTriang>::ClassId() const
{
	return TPZINTELTRIANGLEID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapeTriang>, TPZINTELTRIANGLEID>;
#endif

#include "pzshapequad.h"
template<>
int TPZIntelGen<pzshape::TPZShapeQuad>::ClassId() const
{
	return TPZINTELQUADID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapeQuad>, TPZINTELQUADID>;
#endif

#include "pzshapecube.h"
template<>
int TPZIntelGen<pzshape::TPZShapeCube>::ClassId() const
{
	return TPZINTELCUBEID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapeCube>, TPZINTELCUBEID>;
#endif

#include "pzshapetetra.h"
template<>
int TPZIntelGen<pzshape::TPZShapeTetra>::ClassId() const
{
	return TPZINTELTETRAID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapeTetra>, TPZINTELTETRAID>;
#endif

#include "pzshapeprism.h"
template<>
int TPZIntelGen<TPZShapePrism>::ClassId() const
{
	return TPZINTELPRISMID;
}

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapePrism>, TPZINTELPRISMID>;
#endif

#include "pzshapepiram.h"
template<>
int TPZIntelGen<TPZShapePiram>::ClassId() const
{
	return TPZINTELPYRAMID;
}

#include "pzshapepiramHdiv.h"
template<>
int TPZIntelGen<TPZShapePiramHdiv>::ClassId() const
{
    return TPZINTELPYRAMIDHDIV;
}


#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapePiram>, TPZINTELPYRAMID>;
#endif

#ifndef BORLAND
template class
TPZRestoreClass< TPZIntelGen<TPZShapePiramHdiv>, TPZINTELPYRAMIDHDIV>;
#endif

#include "TPZRefCube.h"

#include "TPZRefLinear.h"
#include "pzrefquad.h"

#include "pzgeoquad.h"

#include "pzreftriangle.h"
#include "pzgeotriangle.h"

#include "pzrefprism.h"
#include "pzgeoprism.h"

#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"

#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzgraphelq2dd.h"
#include "tpzgraphelt3d.h"
#include "pzgraphel1dd.h"
#include "pztrigraphd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"
#include "tpzgraphelpyramidmapped.h"
#include "tpzgraphelt2dmapped.h"



#include "tpzpoint.h"
#include "tpzline.h"
#include "tpzquadrilateral.h"
#include "tpztriangle.h"
#include "tpzcube.h"
#include "tpztetrahedron.h"
#include "tpzprism.h"
#include "tpzpyramid.h"

template class TPZIntelGen<TPZShapeTriang>;
template class TPZIntelGen<TPZShapePoint>;
template class TPZIntelGen<TPZShapeLinear>;
template class TPZIntelGen<TPZShapeQuad>;
template class TPZIntelGen<TPZShapeTetra>;
template class TPZIntelGen<TPZShapePrism>;
template class TPZIntelGen<TPZShapePiram>;
template class TPZIntelGen<TPZShapePiramHdiv>;
template class TPZIntelGen<TPZShapeCube>;

