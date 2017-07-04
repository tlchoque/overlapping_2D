/**File : pzgclonemesh.c
 
 Method definition for class TPZGeoCloneMesh.*/

#include "TPZRefPattern.h"

#include "pzgclonemesh.h"
#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgnode.h"
#include "pzmaterial.h"
#include "pzerror.h"
#include "pzgeoel.h"
//#include "pzcosys.h"
#include "pzmatrix.h"
//#include "pzelg1d.h"
//#include "pzelgc3d.h"
//#include "pzelgpi3d.h"
//#include "pzelgpoint.h"
//#include "pzelgpr3d.h"
//#include "pzelgq2d.h"
//#include "pzelgt2d.h"
//#include "pzelgt3d.h"
#include "pzelasmat.h"

//static int zero=0;
//static TPZGeoEl *zeropoint = 0;

using namespace std;

TPZGeoCloneMesh::TPZGeoCloneMesh(TPZGeoMesh *ref) : TPZGeoMesh(), fMapNodes(),fMapElements() {
    if (!ref) {
        cout << "TPZGeoCloneMesh::Error\n Reference mesh and reference element must not be NULL!\n";
    }
    fGeoReference = ref;
    fGeoRoot = 0; //will be setted in SetElements method!
}

void TPZGeoCloneMesh::SetElements(TPZStack <TPZGeoEl *> &patch, TPZGeoEl *ref){
    
    int i;
    if (!ref){
        cout << "TPZGeoCloneMesh::Error\n Reference element must not be null\n";
        DebugStop();
    }
    //fGeoRoot = ref;
    CloneElement(ref);
    fGeoRoot = fMapElements[ref];
    
    //  cout << "\n\n\nTeste\n\n\n";
    //  Print(cout);
    
    int nel = patch.NElements();
    //  fReferenceElement.Resize(nel);
    for (i=0; i<nel; i++){
        if(patch[i]) {
            TPZGeoEl *gel = patch[i];
            TPZGeoEl *father = gel->Father();
            while(father) {
                gel = father;
                father = gel->Father();
            }
            //      cout << "\nElemento a ser clonado:\n";
            //      gel->Print(cout);
            CloneElement(gel);
            // verificar se neighbour.Element ja esta no map --->>>> j� � feito no CloneElement
            TPZGeoEl *localpatch = fMapElements[patch[i]];
#ifdef DEBUG 
			if (localpatch == 0) {
				DebugStop();
			}
#endif
            fPatchElements.insert(localpatch);
            AddBoundaryConditionElements(patch[i]);

//			cout << "Printing fPatchElements("<< i << "_NEl_" << fPatchElements.size() << ") : " << endl;
//			std::set<TPZGeoEl *>::iterator it;
//			for(it=fPatchElements.begin();it!=fPatchElements.end();it++) {
//				(*it)->Print(cout);
//			}
		}
    }
}

// IT IS BAD !!!! IMPROVE IT !!!
void TPZGeoCloneMesh::AddBoundaryConditionElements(TPZGeoEl *eltoadd) {
    int nsides = eltoadd->NSides();
    int is;
    for(is=0; is<nsides; is++) {
        TPZGeoElSide elside(eltoadd,is);
        TPZGeoElSide neighbour = elside.Neighbour();
#ifdef DEBUG
        if (!neighbour.Element()) DebugStop();
#endif
        while(neighbour != elside) {
            if(neighbour.Element()->Dimension() < eltoadd->Dimension() &&
               neighbour.Side() == neighbour.Element()->NSides() - 1
               // && neighbour.Element()->Reference()
			   ) {
			   
                TPZGeoEl *gel = neighbour.Element();
                if (HasElement(gel)) {
                    neighbour = neighbour.Neighbour();
                    continue;
                }
                TPZGeoEl *father = gel->Father();
                while(father) {
                    gel = father;
                    father = gel->Father();
                }
                CloneElement(gel);
                // verificar se neighbour.Element ja esta no map
          //      TPZGeoEl *localpatch = fMapElements[neighbour.Element()];
				TPZGeoEl *localpatch = fMapElements[gel];							// Jorge 2013/03/28
                fPatchElements.insert(localpatch);
            }
            neighbour = neighbour.Neighbour();
#ifdef DEBUG
            if (!neighbour.Exists()) {
                DebugStop();
            }
#endif
        }
    }
}

void TPZGeoCloneMesh::AddElement(TPZGeoEl *eltoadd){
    if (!fGeoRoot){
        cout << "TPZGeoCloneMesh::Error\n Reference element must not be null\n";
    }
    CloneElement(eltoadd);
}

int TPZGeoCloneMesh::Index(TPZGeoEl *gel) {
    return gel->Index();
/*    int nel = ElementVec().NElements();
    int iel = 0;
    while(iel < nel && ElementVec()[iel] != gel) iel++;
    if(iel == nel) return -1;
    return iel;
 */
}

int  TPZGeoCloneMesh::CloneElement(TPZGeoEl *orgel){
    int i,j;
//    int nnod = orgel->NNodes();
    //  cout << "Original element nodes = " << nnod << endl;
    if(HasElement(orgel)) return Index(fMapElements[orgel]);
    
    
    // Create a cloned element
    TPZGeoEl *el = InitializeClone(orgel);
    int elindex = Index(el);
    
    //  cout << "\nClonned element\n";
    //  el->Print(cout);
    
    //fill the map
    fMapElements[orgel] = el;
    if(elindex >= fReferenceElement.NElements()) {
        fReferenceElement.Resize(elindex+1,0);
    }	
    fReferenceElement[elindex] = orgel;
    //fill the neighbours
    for (i=0;i<orgel->NSides();i++){
        el->SetSideDefined(i);
        TPZGeoElSide neig = orgel->Neighbour(i);
        if(!neig.Element()) continue;
        // insert all neighbours which have been cloned as neighbours
        // THIS IS OVERKILL it would be suficient to insert a single neighbour
        while(neig.Element() != orgel) {
            // verify if neig.Element has been cloned
            if (HasElement((neig.Element()))){
                TPZGeoElSide sid(el,i);
                //          sid.SetConnectivity(sid);
                //SetNeighbour(i,neig);
                TPZGeoElSide localneig(fMapElements[neig.Element()],neig.Side());
                if(!sid.NeighbourExists(localneig)) {
                    sid.SetConnectivity(localneig);
                }
            }
            neig = neig.Neighbour();
        }
    }
    //loop over the sons
    if (orgel->HasSubElement()){
        int subel = orgel->NSubElements();
        for (j=0;j<subel;j++){
            TPZGeoEl *gelson = orgel->SubElement(j);
            CloneElement(gelson);
            fMapElements[gelson]->SetFather(el);
            fMapElements[gelson]->SetFather(el->Index());
            el->SetSubElement(j,fMapElements[gelson]);
        }
    }
    
    //  el->Print(cout);
    return elindex;
}

int TPZGeoCloneMesh::HasNode(int nodeindex){
    return (fMapNodes.find(nodeindex) != fMapNodes.end());
}

int TPZGeoCloneMesh::HasElement(TPZGeoEl *el){
    if(!el) return 0;
/*	std::map<TPZGeoEl *,TPZGeoEl *>::iterator it;
	it = fMapElements.find(el);
	if(it!=fMapElements.end()) {
		(*it).first->Print();
		(*it).second->Print();
		for(it=fMapElements.begin();it!=fMapElements.end();it++) {
			(*it).first->Print();
			(*it).second->Print();
		}
	}
*/	
    return (fMapElements.find(el) != fMapElements.end());
}

/*
int TPZGeoCloneMesh::IsPatchElement(TPZGeoEl *gel){
    int iel, nel = fPatchElements.NElements();
    for (iel = 0; iel < nel; iel++) {
        if (fPatchElements[iel] == gel) return 1;
    }
    return 0;
}
*/


int TPZGeoCloneMesh:: IsPatchSon(TPZGeoEl *gel) const {
    if (!gel) return 0;
    while(gel) {
        if (fPatchElements.find(gel) != fPatchElements.end()) return 1;
        gel = gel->Father();
    }
    return 0;
}
int TPZGeoCloneMesh:: IsNeighBCForPatchSon(TPZGeoEl *gel) const {
    if (!gel) return 0;
	TPZGeoElSide gelside(gel,gel->NSides()-1);
	std::set<TPZGeoEl *>::iterator it;
    for(it = fPatchElements.begin();it!=fPatchElements.end();it++) {
		int nsides = (*it)->NSides();
		for(int c=0;c<nsides-1;c++)
		{
			if((*it)->NeighbourExists(c,gelside))
				return 1;
		}
    }
    return 0;
}



int TPZGeoCloneMesh::CloneNode(int orgnodindex){
    if (fMapNodes.find(orgnodindex) != fMapNodes.end()) {
        DebugStop();
    }
    int i;
    TPZVec<REAL> coord(3,0.);
    for (i=0; i<3; i++) coord[i] = fGeoReference->NodeVec()[orgnodindex].Coord(i);
    int nodind = NodeVec().AllocateNewElement();
    int nodid =  fGeoReference->NodeVec()[orgnodindex].Id();
    NodeVec()[nodind].Initialize(nodid,coord,*this);
    //fill the map
    fMapNodes[orgnodindex] = nodind;
    return nodind;
}

TPZGeoCloneMesh::~TPZGeoCloneMesh() {
    CleanUp();
}


TPZGeoEl* TPZGeoCloneMesh::InitializeClone(TPZGeoEl *orgel){
    
    TPZManVector<int,27> clonindex(orgel->NNodes(),-1);
    int nnod = orgel->NNodes();
    for (int i=0; i<nnod; i++){
        int nodindex = orgel->NodeIndex(i);	
        if (!HasNode(nodindex)){
            //  cout << "Cloning node " << nodindex << "  resulting node ";
            clonindex[i] = CloneNode(nodindex);
            //      cout << clonid << endl;
        }
        else{
            clonindex[i] = fMapNodes[nodindex];
        }
//        el->SetNodeIndex(i,clonid);
    }

    // the nodeindexes of the new element are equal to the node indexes of the original mesh?
    int index;
    TPZGeoEl *gel = CreateGeoElement((MElementType)orgel->Type(),clonindex,orgel->MaterialId(),index);
    gel->SetRefPattern(orgel->GetRefPattern());
    return gel;
    //   TPZGeoEl1d *el1d = dynamic_cast<TPZGeoEl1d *>(org);
    //   if (el1d){
    //     TPZVec<int> ni(2,0);
    //     for (i=0;i<2;i++){
    //       ni[i]=el1d->NodeIndex(i);
    //     }
    //     TPZGeoEl1d *gel = new TPZGeoEl1d (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   TPZGeoElC3d *elc3d = dynamic_cast<TPZGeoElC3d *>(org);
    //   if (elc3d) {
    //     TPZVec<int> ni(8,0);
    //     for (i=0;i<8;i++){
    //       ni[i]=elc3d->NodeIndex(i);
    //     }
    //     TPZGeoElC3d *gel = new TPZGeoElC3d (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   TPZGeoElPi3d *elpi3d = dynamic_cast<TPZGeoElPi3d *>(org);
    //   if (elpi3d) {
    //     TPZVec<int> ni(5,0);
    //     for (i=0;i<5;i++){
    //       ni[i]=elpi3d->NodeIndex(i);
    //     }
    //     TPZGeoElPi3d *gel = new TPZGeoElPi3d (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   TPZGeoElPoint *elpoint = dynamic_cast<TPZGeoElPoint *>(org);
    //   if (elpoint) {
    //     TPZVec<int> ni(1,0);
    //     for (i=0;i<1;i++){
    //       ni[i]=elpoint->NodeIndex(i);
    //     }
    //     TPZGeoElPoint *gel = new TPZGeoElPoint (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   TPZGeoElPr3d *elpr3d = dynamic_cast<TPZGeoElPr3d *>(org);
    //   if (elpr3d) {
    //     TPZVec<int> ni(6,0);
    //     for (i=0;i<6;i++){
    //       ni[i]=elpr3d->NodeIndex(i);
    //     }
    //     TPZGeoElPr3d *gel = new TPZGeoElPr3d (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   TPZGeoElQ2d *elq2d = dynamic_cast<TPZGeoElQ2d *>(org);
    //   if (elq2d) {
    //     TPZVec<int> ni(4,0);
    //     for (i=0;i<4;i++){
    //       ni[i]=elq2d->NodeIndex(i);
    //     }
    //     TPZGeoElQ2d *gel = new TPZGeoElQ2d (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   TPZGeoElT2d *elt2d = dynamic_cast<TPZGeoElT2d *>(org);
    //   if (elt2d) {
    //     TPZVec<int> ni(3,0);
    //     for (i=0;i<3;i++){
    //       ni[i]=elt2d->NodeIndex(i);
    //     }
    //     TPZGeoElT2d *gel = new TPZGeoElT2d (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   TPZGeoElT3d *elt3d = dynamic_cast<TPZGeoElT3d *>(org);
    //   if (elt3d) {
    //     TPZVec<int> ni(4,0);
    //     for (i=0;i<4;i++){
    //       ni[i]=elt3d->NodeIndex(i);
    //     }
    //     TPZGeoElT3d *gel = new TPZGeoElT3d (ni,org->MaterialId(),*this);
    //     return gel;
    //   }
    //   return 0;
}

TPZGeoEl* TPZGeoCloneMesh::ReferenceElement(int i) {
    int nref = fReferenceElement.NElements();
    if (i < 0 || i >= nref) return 0;
    else  return fReferenceElement[i];
}

void TPZGeoCloneMesh::Print (ostream & out) {
    out << "\n\t\t GEOMETRIC CLONE TPZGeoCloneMesh INFORMATIONS:\n\n";
    out << "Reference Mesh:\t" << fGeoReference->Name() << endl;
    //TPZGeoMesh::Print(out);
    //fGeoReference->Print(out);
    out << "TITLE-> " << Name() << "\n\n";
    out << "number of nodes               = " << NodeVec().NElements() << "\n";
    out << "number of elements            = " << ElementVec().NElements() << "\n";
    out << "\n\nGeometric Reference Element:\n";
    if (fGeoRoot)
        fGeoRoot->Print(out);
    else {
        cout << "Not defined yet\n";
        return;
    }
    out << "\n\tGeometric Node Information:\n\n";
    int i;
    int nnodes = NodeVec().NElements();
    for(i=0; i<nnodes; i++) {
        NodeVec()[i].Print(out);
        out << "\n";
    }
    
    out << "\n\tReference Mesh Geometric Node Information:\n\n";
    nnodes = fGeoReference->NodeVec().NElements();
    for(i=0; i<nnodes; i++) {
        fGeoReference->NodeVec()[i].Print(out);
        out << "\n";
    }
    
    
    out << "\n\tGeometric Element Information:\n\n";
    int nelem = ElementVec().NElements();
    for(i=0; i<nelem; i++) {
        out << "elelelelelelelelelelelelelelelelelelelelelelelelelel\n";
        if(ElementVec()[i]) ElementVec()[i]->Print(out);
        out << "\n";
        out << "Reference Element:";
        int refind = Index(ElementVec()[i]);
        if (refind < 0 || refind >= fReferenceElement.NElements()){continue;}
        if(fReferenceElement[refind]) out<< fReferenceElement[refind]->Id();//->Print(out);
        out <<"\n-----------------------------------------------------\n\n";
    }
}

int TPZGeoCloneMesh::main(){
	cout << "**************************************" << endl;
  	cout << "******Obten��o de Patches!************" << endl;
	cout << "**************************************" << endl;
    
    /*******************************************************
     * Constru��o da malha
     * *****************************************************/
  	//malha quadrada de nr x nc
	const int numrel = 3;
  	const int numcel = 3;
//  	int numel = numrel*numcel;
  	TPZVec<REAL> coord(2,0.);
    
  	// criar um objeto tipo malha geometrica
  	TPZGeoMesh geomesh;
    
  	// criar nos
  	int i,j;
  	for(i=0; i<(numrel+1); i++) {
        for (j=0; j<(numcel+1); j++) {
            int nodind = geomesh.NodeVec().AllocateNewElement();
            TPZVec<REAL> coord(2);
            coord[0] = j;//co[nod][0];
            coord[1] = i;//co[nod][1];
            geomesh.NodeVec()[nodind] = TPZGeoNode(i*(numrel+1)+j,coord,geomesh);
        }
  	}
  	// cria��o dos elementos
  	int elc, elr;
  	TPZGeoEl *gel[numrel*numcel];
  	TPZVec<int> indices(4);
  	for(elr=0; elr<numrel; elr++) {  
        for(elc=0; elc<numcel; elc++) {
            indices[0] = (numrel+1)*elr+elc;
            indices[1] = indices[0]+1;
            indices[3] = indices[0]+numrel+1;
            indices[2] = indices[1]+numrel+1;
            // O proprio construtor vai inserir o elemento na malha
			int index;
			gel[elr*numrel+elc] = geomesh.CreateGeoElement(EQuadrilateral,indices,1,index);
            //gel[elr*numrel+elc] = new TPZGeoElQ2d(elr*numrel+elc,indices,1,geomesh);
        }
  	}
	//Divis�o dos elementos
  	TPZVec<TPZGeoEl *> sub;
  	gel[0]->Divide(sub);
    //  	gel[1]->Divide(sub);
    //  	gel[3]->Divide(sub); 
  	ofstream output("patches.dat");
  	geomesh.Print(output);
    //  	TPZGeoElBC t3(gel[0],4,-1,geomesh); 
    //  	TPZGeoElBC t4(gel[numel-1],6,-2,geomesh); 
  	geomesh.Print(output);
	geomesh.BuildConnectivity();
    std::set <TPZGeoEl *> patch;
    
	
  	TPZCompMesh *comp = new TPZCompMesh(&geomesh);
 	// inserir os materiais
  	TPZMaterial *meumat = new TPZElasticityMaterial(1,1.e5,0.2,0,0);
  	comp->InsertMaterialObject(meumat);
  	// inserir a condicao de contorno
    //  	TPZFMatrix val1(3,3,0.),val2(3,1,0.);
    //  	TPZMaterial *bnd = meumat->CreateBC (-1,0,val1,val2);
    //  	comp->InsertMaterialObject(bnd);
    //  	TPZFMatrix val3(3,3,1);
    // 	bnd = meumat->CreateBC (-2,1,val3,val2);
    //  	comp->InsertMaterialObject(bnd);
	comp->AutoBuild();
	comp->Print(output);
  	output.flush();
    
    /**********************************************************************
     * Cria��o de uma malha computacional clone
     * ********************************************************************/
 	comp->GetRefPatches(patch);
	
	geomesh.ResetReference();
	TPZStack <int> patchel;
	TPZStack <TPZGeoEl *> toclonegel;
	TPZStack <int> patchindex;
	TPZVec<int> n2elgraph;
	TPZVec<int> n2elgraphid;
	TPZStack<int> elgraph;
	TPZVec<int> elgraphindex;
	int k;
	TPZCompMesh *clonecmesh = new TPZCompMesh(&geomesh);
	cout << "Check 1: number of reference elements for patch before createcompel: " << patch.size() << endl;
    std::set<TPZGeoEl *>::iterator it;
    for (it=patch.begin(); it!=patch.end(); it++)
    {
		//patch[i]->Print(cout);
        int index;
        TPZGeoEl *gel = *it;
        clonecmesh->CreateCompEl(gel, index);
        //		patch[i]->CreateCompEl(*clonecmesh,i);
	}
    //	cout << "Check 2: number of reference elements for patch after createcompel: " << patch.NElements() << endl;
	clonecmesh->CleanUpUnconnectedNodes();
    //	clonecmesh->Print(cout);
	clonecmesh->GetNodeToElGraph(n2elgraph,n2elgraphid,elgraph,elgraphindex);
	int clnel = clonecmesh->NElements();
    //	cout << "Number of elements in clonemessh: " << clnel << endl;
	//o primeiro patch come�a em zero
	patchindex.Push(0);
	for (i=0; i<clnel; i++){
		//cout << endl << endl << "Evaluating patch for element: " << i << endl;
		clonecmesh->GetElementPatch(n2elgraph,n2elgraphid,elgraph,elgraphindex,i,patchel);
		cout << "Patch elements: " << patchel.NElements() << endl;
		/*for (k=0;k<patchel.NElements();k++){
         clonecmesh->ElementVec()[patchel[k]]->Reference()->Print();
         cout << endl;
         }*/
		for (j=0; j<patchel.NElements(); j++){
			//obten��o do elemento geom�trico do patch
			//cout << "Creating geometric clone elements for computational element :" << j << endl;
			TPZGeoEl *gel = clonecmesh->ElementVec()[patchel[j]]->Reference();
			//gel->Print(cout);
			//inserir todos os pais do elemento geom�trico do patch
			int count = 0;
			//cout << "Inserting father element:" << "\t"; 
			while(gel){	
				TPZGeoEl *father = gel->Father();
				if (father){
					//father->Print(cout);
					gel = father;
					continue;
				}
				else toclonegel.Push(gel);
				gel = father;
				//cout <<  count << "\t";
				count ++;
			}
			//cout << endl;
		}
		int sum = toclonegel.NElements()-1;
		//cout << endl << sum << endl;
		patchindex.Push(sum);
        
		/*for (k=patchindex[i];k<patchindex[i+1];k++){
         toclonegel[k]->Print();
         }*/
	}
    
    
	cout <<endl;
	cout << endl;
	
	TPZGeoCloneMesh geoclone(&geomesh);
	TPZStack<TPZGeoEl*> testpatch;
	for (j=0; j<1/*patchindex.NElements()-1*/;j++){
		cout << "\n\n\nClone do Patch do elemento: " << j <<endl;
		k=0;
		cout << patchindex[j] << "\t" << patchindex[j+1] <<endl;
		for (i=patchindex[j];i<=patchindex[j+1];i++){
			testpatch.Push(toclonegel[i]);
			toclonegel[i]->Print();
			cout << k << endl;
			k++;
		}
		geoclone.SetElements(testpatch,testpatch[patchindex[j]]);
		geoclone.Print(cout);
        
	}
	//geoclone.SetElements(testpatch);
	//geoclone.Print(cout);
	
    
    /**************************************************************************
     * Fim da cria��o do clone
     **************************************************************************/
    
    
	
    /*	output <<"Impress�o dos Pathces\nN�mero total de patches encontrados\t" << patchindex.NElements()-1 << endl;
     cout << "\n\n&&&&&&&&&&&&&&&&&&&&&&&&\n N�mero total de patches: " << patchindex.NElements()-1 << endl
     << "&&&&&&&&&&&&&&&&&&&&&&&&" << endl;
     for (i=0;i<patchindex.NElements()-1;i++){
     cout << "Patch do elemento " << i << "\t" << "N�mero de elementos componentes do patch: " << (patchindex[i+1]-patchindex[i]) << endl;
     for (j = patchindex[i]; j<patchindex[i+1]; j++){
     toclonegel[j]->Print();
     cout << "||||||||||||||||||||||||||||||||" << endl;
     }
     cout << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n" <<">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n\n";
     cout.flush();
     }
     */	
	comp->LoadReferences();
	
	cout.flush();
	cout << endl;
	cout.flush();
	delete comp;
	delete clonecmesh;
  	return (0);
    
}

