/**
 * @file
 * @brief Implements a tutorial example using shape NeoPZ module
 */
#include "pzreal.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgengrid.h"
#include "TPZExtendGridDimension.h"
#include "pzintel.h"

#include "pzl2projection.h"
#include "pzpoisson3d.h"

#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzgeoelbc.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzshapeprism.h"
#include "pzshapecube.h"

#include "pzvtkmesh.h"
#include "pzanalysis.h"

#include "pzeltype.h"

#include "pzlog.h"

// integration rule
#include <pzquad.h>

#include <iostream>
#include <math.h>

using namespace pzshape;
using namespace std;


TPZCompMesh *InitialMesh(int order,int nsubdiv,int dim,MElementType elem_type);
TPZAutoPointer<TPZCompMesh> CompMesh();

void TestShape(TPZInterpolatedElement *el,int order,ostream &out);
void TestShapeWithPrint(TPZInterpolatedElement *el,int order,ostream &saida);

// Para verifica a linearidade, 'e preciso especificar o lado e a direcao da linearidade
// todas as funcoes associadas ao lado devem ser lineares com respeito a variacao de pontos de integracao na direcao indicada
// uma alternativa seria verificar se grad(phi).direction 'e constante
void TestShapeIsLinear(TPZAutoPointer<TPZCompMesh> mesh, MElementType type,int side, REAL *direction,ostream &saida);

int nquadsides = 4;
int quadsides[] =
{
    4,5,6,7
};
REAL quaddir[][3] =
{
    {0.,1.},
    {-1.,0},
    {0,-1.},
    {1.,0.}
};

int ncubesides = 6;
int cubesides[] =
{
    20,21,22,23,24,25
};
REAL cubedir[][3] =
{
    {0.,0.,1.},
    {0.,1.,0.},
    {-1.,0.,0.},
    {0.,-1.,0.},
    {1.,0.,0.},
    {0.,0.,-1.}
};
int nprismsides = 2;
int prismsides[] =
{
    15,19
};
REAL prismdir[][3] =
{
    {0.,0.,1.},
    {0.,0.,-1.}
};

/**
 * Example of basis functions on non-conformal mesh
 */

int MaterialId = 3;
int idbc0 = -3;

void MakingQuadrilateral(TPZGeoMesh *gmesh,REAL radio);
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction);
void ComputeShape(TPZInterpolatedElement* intel,TPZVec<REAL> &qsi, TPZFNMatrix<220> &phi,TPZFNMatrix<660> &dphi,TPZFNMatrix<9> &axes);

int main() {
    
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    /** To view a basis function on quadrilateral */
    int dim = 2;
    REAL radio = 2.;
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    MakingQuadrilateral(gmesh,radio);
    
    // Refine some elements until two level
    TPZVec<TPZGeoEl *> sub(4);
    TPZVec<TPZGeoEl *> subsub(4);
    gmesh->ElementVec()[1]->Divide(sub);
    sub[1]->Divide(subsub);
    
    // Creating computational mesh (approximation space and materials)
    int p = 2;
    TPZCompEl::SetgOrder(p);
    TPZCompMesh *compmesh;
    compmesh = CreateMesh(gmesh,dim,0);
    TPZAnalysis an(compmesh);
    
    long neq = compmesh->NEquations();
    TPZManVector<long> indices(neq);
    for (int i=0; i<neq; i++) {
        indices[i] = i;
    }
    an.ShowShape("Shapes.vtk", indices);
    
    
    /** Another test for shape functions */
    long elem;
    //int side, dim = 3; // se dim igual a 3 sera construida malha com elementos tridimensionais, se diferente elementos bidimensionais
    int max_order, elem_type;
    TPZInterpolatedElement *el;
    TPZVec<REAL> normal(3,0.);
    
    ofstream saida("malha.txt");
    elem_type = 0; // se igual a 1 triangulo, se diferente quadrilatero
    
    // Insere malha geometrica e computacional
    max_order = 15;
    // generate a mesh that contains all element topologies
    TPZAutoPointer<TPZCompMesh> cmesh = CompMesh();
    
    // Calcula os valores das funcoes shape no lado e no elemento
    for(elem=0;elem<cmesh->NElements();elem++) {
        el = (TPZInterpolatedElement *)cmesh->ElementVec()[elem];
        // Testing side shape function values with shape function values
        TestShapeWithPrint(el,max_order-1,saida);
    }
    
    // verify if the evolution of the prismatic shape functions is linear in the direction oposite of the side
    // it is a nice exercise in understanding the concept of sides, side shape functions and element topology
    for(elem=0; elem<nquadsides; elem++)
    {
        TestShapeIsLinear(cmesh, EQuadrilateral,quadsides[elem], quaddir[elem],saida);
    }
    for(elem=0; elem<ncubesides; elem++)
    {
        TestShapeIsLinear(cmesh, ECube,cubesides[elem], cubedir[elem],saida);
    }
    for(elem=0; elem<nprismsides; elem++)
    {
        TestShapeIsLinear(cmesh, EPrisma,prismsides[elem], prismdir[elem],saida);
    }
    
    saida.close();
    return 0;
}

void ComputeShape(TPZInterpolatedElement* intel,TPZVec<REAL> &qsi, TPZFNMatrix<220> &phi,TPZFNMatrix<660> &dphi,TPZFNMatrix<9> &axes) {
    int nshape = intel->NShapeF();
    TPZGeoEl* ref = intel->Reference();
    const int dim = ref->Dimension();
    
    TPZFNMatrix<660> dphix(dim,nshape);
    TPZFNMatrix<9> jacobian(dim,dim);
    TPZFNMatrix<9> jacinv(dim,dim);
    REAL detjac;
    
    ref->Jacobian( qsi, jacobian, axes, detjac , jacinv);
    
    intel->Shape(qsi,phi,dphi);
}

//**********   Creating computational mesh with materials    *************
TPZCompMesh *CreateMesh(TPZGeoMesh *gmesh,int dim,int hasforcingfunction) {
    
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(TPZCompEl::GetgOrder());
    cmesh->SetAllCreateFunctionsDiscontinuous();
    
    // Creating Poisson material
    TPZMaterial *mat = new TPZMatPoisson3d(MaterialId,dim);
    cmesh->InsertMaterialObject(mat);
    // Make compatible dimension of the model and the computational mesh
    cmesh->SetDimModel(mat->Dimension());
    
    // Creating four boundary condition
    TPZFMatrix<STATE> val1(dim,dim,0.),val2(dim,1,0.);
    TPZBndCond *bc = 0;
    bc = mat->CreateBC(mat,idbc0,0,val1,val2);
    if(bc) cmesh->InsertMaterialObject((TPZMaterial *)bc);
    
    cmesh->AutoBuild();
    cmesh->ApproxSpace().CreateInterfaces(*cmesh);
    cmesh->ExpandSolution();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    return cmesh;
}

void MakingQuadrilateral(TPZGeoMesh *gmesh,REAL radio) {
    REAL co[6][2] = {{1.,1.},{(1.+radio),1.},{(1.+2.*radio),1.},{1.,(1.+radio)},{(1.+radio),(1.+radio)},{(1.+2*radio),(1.+radio)}};
    int indices[2][4] = {{0,1,4,3},{1,2,5,4}};
    TPZGeoEl *elvec[2];
    int nnode = 6;
    int nod;
    for(nod=0; nod<nnode; nod++) {
        int nodind = gmesh->NodeVec().AllocateNewElement();
        TPZVec<REAL> coord(3,0.0);
        coord[0] = co[nod][0];
        coord[1] = co[nod][1];
        gmesh->NodeVec()[nodind] = TPZGeoNode(nod,coord,*gmesh);
    }
    
    int el;
    int nelem = 2;
    for(el=0; el<nelem; el++) {
        TPZVec<long> nodind(4);
        for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
        long index;
        elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,MaterialId,index);
    }
    
    gmesh->BuildConnectivity();
    // bc -1 -> Dirichlet
    TPZGeoElBC gbc1(elvec[0],4,idbc0);
    TPZGeoElBC gbc3(elvec[0],6,idbc0);
    TPZGeoElBC gbc4(elvec[0],7,idbc0);
    TPZGeoElBC gbc5(elvec[1],4,idbc0);
    TPZGeoElBC gbc6(elvec[1],5,idbc0);
    TPZGeoElBC gbc7(elvec[1],6,idbc0);
}

TPZCompMesh *InitialMesh(int order,int nsubdiv,int dim,MElementType type) {
    if(dim != 3)
        dim = 2;
    
    TPZAutoPointer<TPZGeoMesh> g = new TPZGeoMesh;
    TPZGeoMesh *g_extended = 0;
    TPZVec<int> nx(3);
    TPZVec<REAL> X0(3);
    TPZVec<REAL> X1(3);
    for(int i=0;i<3;i++) {
        nx[i] = nsubdiv;
        X0[i] = -3.;
        X1[i] = 3.;
    }
    TPZGenGrid gen(nx,X0,X1,1);
    if(type==ETriangle) gen.SetElementType(type);  // type = 1 para elementos triangulares
    else gen.SetElementType(EQuadrilateral);  // type = 0 para elementos quadrilateros
    gen.Read(g);
    if(dim == 3) {
        TPZExtendGridDimension genext(g,1.);
        g_extended = genext.ExtendedMesh();
    }
    else
        g_extended = g.operator->();
    TPZCompMesh *c = new TPZCompMesh(g_extended);
    TPZVec<STATE> sol(1,0.);
    TPZMaterial * material = new TPZL2Projection(1,2,1,sol);
    c->InsertMaterialObject(material);
    c->SetDefaultOrder(order);
    c->AutoBuild();
    return c;
}

void TestShapeWithPrint(TPZInterpolatedElement *el,int order,ostream &saida) {
    TPZGeoEl *gel = el->Reference();
    int nsides = gel->NSides();
    int npoints;
    TPZVec<REAL> p(3,0);
    TPZVec<REAL> p_side(3,0);
    int nsideconnects, nsideshapef;
    int nconnects, nshapef;
    int side, con, point, shape, conside;
    TPZIntPoints *pointIntRule = 0;
    TPZVec<REAL> peso(1000,0);
    TPZTransform transform;
    
    nconnects = el->NConnects();
    nshapef = el->NShapeF();
    TPZFMatrix<REAL> phiint(nshapef,1,0.), dphiint(3,nshapef,0.);
    TPZVec<int> firstindexinternal(nconnects+1,0);
    for(con=0;con<nconnects;con++)
        firstindexinternal[con+1] = firstindexinternal[con]+el->Connect(con).NShape();
    // Imprime dados do elemento, numero de connects e funçoes shape dependendo da ordem escolhida
    saida << "Elemento: " << el->Index() << "  Order = " << order << endl;
    saida << " " << "connects: " << nconnects << "  n shape: " << nshapef << endl;
    for(side=0;side<nsides;side++) {
        nsideconnects = el->NSideConnects(side);
        nsideshapef = el->NSideShapeF(side);
        TPZFMatrix<REAL> phi(nsideshapef,1,0.), dphi(3,nsideshapef,0.);
        TPZVec<int> firstindexside(nsideconnects+1,0);
        for(conside=0;conside<nsideconnects;conside++) {
            firstindexside[conside+1] = firstindexside[conside] + el->Connect(el->SideConnectLocId(conside,side)).NShape();
        }
        saida << "  Side " << side << ": " << "  side connects: " << nsideconnects << "  n side shape: " << nsideshapef << endl;
        transform = gel->SideToSideTransform(side,(nsides-1));
        pointIntRule = el->Reference()->CreateSideIntegrationRule(side,1);
        int maxorder = pointIntRule->GetMaxOrder();
        TPZManVector<int,3> sideorder(el->Reference()->SideDimension(side),0);
        if (maxorder > order) {
            sideorder.Fill(order);
        }
        else {
            sideorder.Fill(maxorder);
        }
        pointIntRule->SetOrder(sideorder);
        
        npoints = pointIntRule->NPoints();
        for(point=0;point<npoints;point++) {
            pointIntRule->Point(point,p_side,peso[point]);
            transform.Apply(p_side,p);
            el->SideShapeFunction(side,p_side,phi,dphi);
            el->Shape(p,phiint,dphiint);
            saida << " " << " Point " << point << " (" << p[0] << "," << p[1] << "," << p[2] << ")" << "  peso = " << peso[point] << endl;
            saida << " " << " " << " MATRIX Side Shape: ";
            for(conside=0;conside<nsideconnects;conside++) {
                nshapef = el->Connect(el->SideConnectLocId(conside,side)).NShape();
                for(shape=0;shape<nshapef;shape++)
                    saida << phi((firstindexside[conside]+shape),0) << "  ";
            }
            saida << endl << " " << " " << " MATRIX Shape:      ";
            for(conside=0;conside<nsideconnects;conside++) {
                nshapef = el->Connect(el->SideConnectLocId(conside,side)).NShape();
                for(shape=0;shape<nshapef;shape++)
                    saida << phiint((firstindexinternal[el->SideConnectLocId(conside,side)]+shape),0) << "  ";
            }
            saida << endl;
            // Repassando e verificando se as diferencas sao nao nulas
            for(conside=0;conside<nsideconnects;conside++) {
                nshapef = el->Connect(el->SideConnectLocId(conside,side)).NShape();
                for(shape=0;shape<nshapef;shape++) {
                    if(!IsZero(phi((firstindexside[conside]+shape),0)-phiint((firstindexinternal[el->SideConnectLocId(conside,side)]+shape),0)))
                        saida << "  Ops -> " << side << " " << conside << " " << shape << " side,connectside,shapefconnectside." << endl;
                }
            }
        }
        saida << endl;
    }
    delete pointIntRule;
}

void TestShapeIsLinear(TPZAutoPointer<TPZCompMesh> cmesh, MElementType type ,int side,REAL *normal,ostream &saida) {
    int nel = cmesh->NElements();
    int el;
    for(el=0; el<nel; el++)
    {
        TPZInterpolatedElement *cel = dynamic_cast<TPZInterpolatedElement *>(cmesh->ElementVec()[el]);
        if(!cel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(gel->Type() != type) continue;
        int npoints, i;
        int order = 5;
        TPZVec<REAL> p1(3,0.), p2(3,0.), p3(3,0.);
        REAL prod_int1(0.), prod_int2(0.), prod_int3(0.);
        TPZVec<REAL> p_pointIntRule(3,0.);
        int nshapef = cel->NShapeF();
        TPZFMatrix<REAL> phi1(nshapef,1,0.), dphi1(3,nshapef,0.);
        TPZFMatrix<REAL> phi2(nshapef,1,0.), dphi2(3,nshapef,0.);
        TPZFMatrix<REAL> phi3(nshapef,1,0.), dphi3(3,nshapef,0.);
        int nsideconnects, nconnectshapef;
        int conn, sideconnect, sideconnectshapef, point;
        TPZIntPoints *pointIntRule = 0;
        TPZVec<REAL> peso(1000,0.);
        TPZTransform transform;
        TPZVec<REAL> n_unit(3,0.);
        REAL n_modulo;
        
        TPZVec<int> firstindexinternal(cel->NConnects()+1,0);
        for(conn=0;conn<cel->NConnects();conn++)
            firstindexinternal[conn+1] = firstindexinternal[conn]+cel->Connect(conn).NShape();
        
        // Normalizing the normal vetor
        n_modulo = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
        if(IsZero(n_modulo)) return;
        for(i=0;i<3;i++) {
            n_unit[i] = normal[i]/n_modulo;
        }
        // On the side: connects, shape functions and first indexes
        nsideconnects = cel->NSideConnects(side);
        // Print element index, number of connects and shape functions on the side
        saida << "Elemento: " << cel->Index() << "  Side = " << side;
        saida << " " << "side_connects: " << nsideconnects << "  N_side_shape: " << cel->NSideShapeF(side) << endl;
        // Creating integration rule on the side
        transform = gel->SideToSideTransform(side,((cel->Reference()->NSides())-1));
        pointIntRule = cel->Reference()->CreateSideIntegrationRule(side,order);
        npoints = pointIntRule->NPoints();
        // Over all side points (based on the integration rule) compute the gradient values to shape functions associated to connects on the side
        for(point=0;point<npoints;point++) {
            pointIntRule->Point(point,p_pointIntRule,peso[point]);
            transform.Apply(p_pointIntRule,p1);
            for(i=0;i<p2.NElements();i++) {
                p2[i] = p1[i]+0.5*n_unit[i];
                p3[i] = p1[i]+n_unit[i];
            }
            cel->Shape(p1,phi1,dphi1);
            cel->Shape(p2,phi2,dphi2);
            cel->Shape(p3,phi3,dphi3);
            saida << " Point integration rule: " << point << endl;
            saida << " - Point1 " << " (" << p1[0] << "," << p1[1] << "," << p1[2] << ")";
            saida << " - Point2 " << " (" << p2[0] << "," << p2[1] << "," << p2[2] << ")";
            saida << " - Point3 " << " (" << p3[0] << "," << p3[1] << "," << p3[2] << ")" << endl;
            // Computing dphi.normal to all shape functions associated to connects on this side
            for(sideconnect=0;sideconnect<nsideconnects;sideconnect++) {
                nconnectshapef = cel->Connect(cel->SideConnectLocId(sideconnect,side)).NShape();
                for(sideconnectshapef=0;sideconnectshapef<nconnectshapef;sideconnectshapef++) {
                    int index = firstindexinternal[cel->SideConnectLocId(sideconnect,side)] + sideconnectshapef;
                    prod_int1 = dphi1(0,index)*n_unit[0] + dphi1(1,index)*n_unit[1] + dphi1(2,index)*n_unit[2];
                    prod_int2 = dphi2(0,index)*n_unit[0] + dphi2(1,index)*n_unit[1] + dphi2(2,index)*n_unit[2];
                    prod_int3 = dphi3(0,index)*n_unit[0] + dphi3(1,index)*n_unit[1] + dphi3(2,index)*n_unit[2];
                    REAL diff = phi2(index,0)-0.5*(phi1(index,0)+phi3(index,0));
                    // All the scalar products must to be equals
                    if(!IsZero(prod_int1-prod_int2) || !IsZero(prod_int2-prod_int3) || ! IsZero(diff))
                    {
                        saida << "Shape " << index << " on the side_connect " << sideconnect << " on the side " << side << " is not linear." << endl;
                    }
                }
            }
            saida << endl;
        }
        saida << endl;
        delete pointIntRule;
    }
}

// Philippe
void TestShape(TPZInterpolatedElement *el,int order,ostream &saida) {
    TPZGeoEl *gel = el->Reference();
    int nsides = gel->NSides();
    int npoints;
    TPZVec<REAL> p(3,0);
    TPZVec<REAL> p_int(3,0);
    TPZFMatrix<REAL> phi(100,1,0.), dphi(3,100,0.);
    TPZFMatrix<REAL> phiint(100,1,0.), dphiint(3,100,0.);
    int nsideconnects, nsideshapef;
    int nconnects, nshapef;
    int i, j;
    TPZIntPoints *pointIntRule = 0;
    TPZVec<REAL> peso(1000,0);
    TPZTransform transform;
    bool notequal = false;
    
    nconnects = el->NConnects();
    nshapef = el->NShapeF();
    TPZVec<int> firstindexinternal(el->NConnects()+1,0);
    for(i=0; i<nsides; i++)
    {
        firstindexinternal[i+1] = firstindexinternal[i]+el->Connect(i).NShape();
    }
    saida << endl;
    for(i=0;i<nsides;i++) {
        nsideconnects = el->NSideConnects(i);
        nsideshapef = el->NSideShapeF(i);
        TPZVec<int> firstindexside(el->NSideConnects(i)+1,0);
        int isc;
        for(isc=0; isc<nsideconnects; isc++)
            firstindexside[isc+1] = firstindexside[isc]+el->Connect(el->SideConnectLocId(isc,i)).NShape();
        transform = gel->SideToSideTransform(i,nsides-1);
        pointIntRule = el->Reference()->CreateSideIntegrationRule(i,order);
        npoints = pointIntRule->NPoints();
        for(j=0;j<npoints;j++) {
            pointIntRule->Point(j,p_int,peso[j]);
            transform.Apply(p_int,p);
            el->SideShapeFunction(i,p_int,phi,dphi);
            el->Shape(p,phiint,dphiint);
            // Comparando os valores das funções shape nos correspondentes connects
            int isc;
            for(isc=0; isc<nsideconnects; isc++)
            {
                int sideconnectindex = el->SideConnectLocId(isc,i);
                for(int il=firstindexside[isc], ilint = firstindexinternal[sideconnectindex];il<firstindexside[isc+1]; il++,ilint++)
                {
                    if(fabs(phi(il,0)-phiint(ilint,0)) > 1.e-6)
                    {
                        notequal = true;
                    }
                }
            }
        }
        if(notequal) { 
            saida << "\tValor diferente. Side " << i << endl;
            notequal = false;
        }
    }
    delete pointIntRule;
}

TPZAutoPointer<TPZCompMesh> CompMesh()
{
    REAL corners[][3] = {
        {0.,0.,0.},
        {1.,0.,0.},
        {1.,1.,0.},
        {0.,1.,0.},
        {0.,0.,1.},
        {1.,0.,1.},
        {1.,1.,1.},
        {0.,1.,1.}
    };
    MElementType types[] = 
    {
        EOned, ETriangle, EQuadrilateral, ETetraedro, EPiramide,
        //		       6        7         8           9             10
        EPrisma,  ECube
    };
    int elnodes[][8] =
    {
        {0,1},
        {0,1,2},
        {0,1,2,3},
        {0,1,2,4},
        {0,1,2,3,4},
        {0,1,2,4,5,6},
        {0,1,2,3,4,5,6,7}
    };
    int nnodes[] = {
        2,
        3,
        4,
        4,
        5,
        6,
        8
    };
    TPZGeoMesh *gm = new TPZGeoMesh;
    int in;
    for(in = 0; in<8; in++)
    {
        TPZManVector<REAL,3> coord(3);
        int c;
        for(c=0; c<3; c++) coord[c] = corners[in][c];
        TPZGeoNode node(in,coord,*gm);
        gm->NodeVec()[gm->NodeVec().AllocateNewElement()] = node;
        
    }
    int el;
    for(el = 0; el<7; el++)
    {
        TPZManVector<long> cornernodes(nnodes[el]);
        int in;
        for(in=0; in<nnodes[el]; in++)
        {
            cornernodes[in] = elnodes[el][in];
        }
        long index;
        gm->CreateGeoElement(types[el],cornernodes,el+1,index);
    }
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gm);
    for(el=0; el<7; el++)
    {
        TPZMaterial *mat = new TPZMatPoisson3d(el+1,gm->ElementVec()[el]->Dimension());
        cmesh->InsertMaterialObject(mat);
    }
    TPZCompEl::SetgOrder(4);
    cmesh->AutoBuild();
    return cmesh;
}
