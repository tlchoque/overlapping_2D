#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysiserror.h"
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "pzmatrix.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzsubcmesh.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "pzmaterial.h"
#include "pzelmat.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "pzhdivpressure.h"

#include "TPZRefPattern.h"


#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("Steklov.main"));

#endif



void ValFunction(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType);
TPZGeoMesh * MalhaGeo(const int h);
TPZGeoMesh * MalhaGeoT(const int h);
TPZGeoMesh * MalhaGeo2(const int h);
/** Resolver o problema do tipo 
 * -Laplac(u) = 0
 * du/dn = lambda u em todo contorno
 */

using namespace std;


TPZCompMeshReferred *CreateMesh2d(TPZGeoMesh &gmesh,int porder);
void PrintMesh(TPZCompMesh *cmesh);
int SubStructure(TPZCompMesh *cmesh, int materialid);

TPZGeoMesh * GeraGMesh(int nrows, int ncols,int h);
void SaddlePermute(TPZCompMesh * cmesh);



int main1()
{
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		InitializePZLOG();
		std::stringstream sout;
		sout<< "Problema de Steklov"<<endl;
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif
	
	for (int porder= 4; porder<5; porder++) {
		
		for(int h=1;h<5;h++){
			
			
			TPZGeoMesh *gmesh2 = MalhaGeo(h);//malha geometrica
			
			
			TPZCompMeshReferred *cmesh = CreateMesh2d(*gmesh2,porder+1);//malha computacional
			
            /*
            TPZCompEl *cel = cmesh->ElementVec()[0];
            TPZManVector<int,5> subindex(4,-1);
            int index;
            cel->Divide(index , subindex);
            cmesh->AdjustBoundaryElements();
            cmesh->ExpandSolution();
            cmesh->CleanUpUnconnectedNodes();
            
            TPZCreateApproximationSpace::MakeRaviartTomas(*cmesh);
            TPZCompEl *cel = cmesh->ElementVec()[0];
            TPZElementMatrix ek(cmesh,TPZElementMatrix::EK),ef(cmesh,TPZElementMatrix::EF);
            cel->CalcStiff(ek, ef);
            ek.ApplyConstraints();
            */
//#ifdef LOG4CXX
//            if (logger->isDebugEnabled())
//            {
//                std::stringstream sout;
//                cmesh->Print(sout);
//                LOGPZ_DEBUG(logger, sout.str())
//            }
//#endif
			int submeshindex = -1;
			TPZSubCompMesh *submesh = 0;
			// Aq faz a condensacao estatica			
			if(h >=0)
				if(-1)
				{
					submeshindex = SubStructure(cmesh,1);//monto a submalha com os elementos q serao condesados (externos) e retorna o numero de elementos computacionais da malha
					submesh = dynamic_cast<TPZSubCompMesh *> (cmesh->ElementVec()[submeshindex]);//converte os elementos computacionais para um objeto do tipo TPZSubCompMesh
					submesh->SetNumberRigidBodyModes(1);//aq defino o numero de pivos nulos que podera ter o sistema?
					cmesh->ExpandSolution();//ajusta o vetor de solucao
					TPZAutoPointer<TPZGuiInterface> guiInter = new TPZGuiInterface;
					int numThreads=0;
					submesh->SetAnalysisSkyline(numThreads,1, guiInter);
						cmesh->SetName("Malha depois de SubStructure-----");
#ifdef LOG4CXX
						{
								std::stringstream sout;
								cmesh->Print(sout);
								LOGPZ_DEBUG(logger,sout.str())
						}
#endif
				}
			

			
			
			cmesh->LoadReferences();//mapeia para a malha geometrica lo
			
			TPZAdmChunkVector<TPZCompEl *> elvec = cmesh->ElementVec();
			
			TPZAnalysis analysis(cmesh);
			//SaddlePermute(cmesh);
				
		cmesh->SetName("Malha depois de Analysis-----");
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
			{
				std::stringstream sout;
				cmesh->Print(sout);
			//	submesh->Block().Print("Block",sout);
				LOGPZ_DEBUG(logger,sout.str())
			}
#endif
	
			TPZFStructMatrix str(cmesh);//estou usando estrutura de matriz cheia
			
			// zero the submesh index so that it won't be assembled
			
			int iel;
			int nel = elvec.NElements();
			for (iel=0; iel<nel; iel++) {
				TPZCompEl *cel = cmesh->ElementVec()[iel];
				if(!cel) continue;
				TPZGeoEl *gel = cel->Reference();
				if(!gel || gel->MaterialId() > 0)
				{
					cmesh->ElementVec()[iel] = 0;//zera os elementos q nao sao de contorno
				}
			}
			
			
			if(submeshindex >= 0) cmesh->ElementVec()[submeshindex] = 0;

			 
			TPZFMatrix<STATE> rhs;//matriz cheia
			// this will be the matrix with only the boundary condition
			
			TPZAutoPointer<TPZGuiInterface> guiInter = new TPZGuiInterface;
			
			TPZAutoPointer<TPZMatrix<STATE> > A = str.CreateAssemble(rhs,guiInter);
			
			
			
			// keep only the subcmesh
			cmesh->ElementVec() = elvec;
			for (iel=0; iel<nel; iel++) {
				TPZCompEl *cel = cmesh->ElementVec()[iel];
				if(!cel) continue;
				TPZGeoEl *gel = cel->Reference();
				if(gel && gel->MaterialId() < 0)
				{
					cmesh->ElementVec()[iel] = 0;//vai zerar todos os elementos de contorno q ja foram assemblados(ou seja os elementos de contorno)
				}
			}			
			 
			
			TPZAutoPointer<TPZMatrix<STATE> > B = str.CreateAssemble(rhs,guiInter);
			
			{
				
				std::stringstream hstr; hstr << h;
				
				std::stringstream pstr; pstr << porder;
				
				
				
				std::string fileNameEigen = "Eigenvalues";
				
				fileNameEigen += pstr.str();
				
				fileNameEigen += "h";
				
				fileNameEigen += hstr.str();
				
				fileNameEigen += ".nb";
				
				std::ofstream eig3(fileNameEigen.c_str());					
				
				
				
				A->Print("Acondense = ",eig3,EMathematicaInput);
				
				B->Print("Bcondense = ",eig3,EMathematicaInput);
				B->Print("1/Eigenvalues[{Bcondense,Acondense}] ",eig3,EMathematicaInput);
				
			}		
			// restore the original state
		
	
	/*	cmesh->ElementVec() = elvec;
	 int autovecsize = cmesh->Solution().Rows();
	 
	 
	 TPZFMatrix<REAL> autovec(autovecsize,1,0.);
	 ifstream file;
	 file.open("Autovec.txt");
	 float vec;
	 int count=0, i;
	 for (i=0; i < autovecsize && file; i++) {
	 file>>vec;
	 std::cout << vec;
	 autovec(i,0)=vec;
	 
	 count++;
	 }
	 cout << "Arquivo Lido -->"<<autovec;
	 if(count == autovecsize)
	 {
	 cmesh->LoadSolution(autovec);
	 }
	 
	 
	 //3. Saida para vtk
	 TPZVec<std::string> scalnames(1), vecnames(1);
	 scalnames[0] = "Pressure";
	 vecnames[0] = "Flux";
	 
	 std::string plotfile("AutovecAprox.vtk");
	 const int dim = 2;
	 int div = 2;
	 analysis.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
	 analysis.PostProcess(div,dim);*/
		}}
	
	return 0;
}

TPZCompMeshReferred *CreateMesh2d(TPZGeoMesh &gmesh,int porder){
	TPZCompEl::SetgOrder(porder);
	TPZCompMeshReferred *comp = new TPZCompMeshReferred(&gmesh);
	
	
	
	// Criar e inserir os materiais na malha
	TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
	TPZMaterial * automat(mat);
	comp->InsertMaterialObject(automat);
	
	
	// Condicoes de contorno
	TPZFMatrix<STATE> val1(1,1,1.),val2(1,1,0.);
	
	TPZMaterial *bnd = automat->CreateBC (automat,-1,2,val1,val2);//misto tbem
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-2,2,val1,val2);
	TPZBndCond *bndcond = dynamic_cast<TPZBndCond *> (bnd);
	bndcond->SetValFunction(ValFunction);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-3,2,val1,val2);
	comp->InsertMaterialObject(bnd);
	
	// Mixed
	val1(0,0) = 1.;
	val2(0,0)=0.;
	bnd = automat->CreateBC (automat,-4,2,val1,val2);
	comp->InsertMaterialObject(bnd);
	
   // comp->SetAllCreateFunctionsHDiv();
		comp->SetAllCreateFunctionsHDivPressure();
	//comp->SetAllCreateFunctionsContinuous();
	
	// Ajuste da estrutura de dados computacional
	comp->AutoBuild();
	
	
	comp->AdjustBoundaryElements();//ajusta as condicoes de contorno
	comp->CleanUpUnconnectedNodes();//deleta os nos que nao tem elemntos conectados
	comp->SetName("Malha Computacional Original");
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		comp->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	
    return comp;
	
}



int SubStructure(TPZCompMesh *cmesh, int materialid)
{
	long index;
	TPZSubCompMesh *submesh = new TPZSubCompMesh(*cmesh,index);//alocacao de memoria...o constructor do tpzsubcompmesh eh inicializado com o parametro index que sera o numero de elementos computacionais da malha
	
	int nelem = cmesh->NElements();
	int iel;
	for(iel = 0; iel<nelem; iel++)
	{
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if(!cel || cel == submesh) continue;
		TPZMaterial * celmat = cel->Material();
		if(!celmat) continue;
		int matid = celmat->Id();
		if(matid == materialid)
		{
			submesh->TransferElement(cmesh, iel);//cada elemento de matId=1 eh transferido para a submalha..cada subelemento tera todas as caracterisitcas de um elemento computacional(?)
		}
	}
	cmesh->ComputeNodElCon();//verifica os no connectados
	submesh->MakeAllInternal();//faz as conexoes da malha com os nos internos
	cmesh->CleanUpUnconnectedNodes();//deleta os nos q nao tem elementos conectados
	
	//	submesh->SetAnalysisSkyline(numThreads4Assemble, guiInterface);
	// submesh->SetAnalysis();
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Mesh after substructuring\n";
		cmesh->Print(sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	return index;
	
}


void ValFunction(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType)
{
	BCType = 2;
	Val1.Redim(1, 1);
	Val1(0,0) = 1.;
	Val2[0] = loc[0];
}

TPZGeoMesh * MalhaGeo2(const int h){//malha quadrilatera
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	TPZGeoEl *elvec[2];
	//Criar ns
	const int nnode = 6;//AQUI
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.},{0.,-1.},{0.,1.}};
	int nelem=2;
	
	int nodindAll[2][nnode]={{0,4,5,3},{4,1,2,5}};//como serao enumerados os nos
	
	
	/*for(int i = 0; i < nnode; i++)
	{
		indices[0][i] = i;
	}
	*/
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		long nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos

	for ( int el=0; el<nelem; el++ )
	{
		TPZVec<long> nodind(4);
		nodind[0]=nodindAll[el][0];
		nodind[1]=nodindAll[el][1];
		nodind[2]=nodindAll[el][2];
		nodind[3]=nodindAll[el][3];
		long index;
		elvec[el] = gmesh->CreateGeoElement (EQuadrilateral,nodind,1,index );
	}
	
	
	gmesh->BuildConnectivity();
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],4,-1);
	TPZGeoElBC gbc2(elvec[1],4,-2);
	TPZGeoElBC gbc3(elvec[1],5,-3);
	TPZGeoElBC gbc4(elvec[1],6,-4);
	TPZGeoElBC gbc5(elvec[0],6,-5);
	TPZGeoElBC gbc6(elvec[0],7,-6);

	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	//gmesh->RefPatternList(ETriangle);
	
		//Refinamento uniforme
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		long n = gmesh->NElements();
		for(long i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			
			
		}
		
		
	}
	
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif 		 		 
	return gmesh;

}

TPZGeoMesh * MalhaGeo(const int h){//malha quadrilatera
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Criar nos
	const int nnode = 4;//AQUI
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};
	long indices[1][nnode]={{0,1,2,3}};//como serao enumerados os nos
	
	
	
	/*for(int i = 0; i < nnode; i++)
	{
		indices[0][i] = i;
	}
	*/
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		long nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
	
	
	TPZVec<long> nodind(4);
	for(int i=0; i<4; i++){
		nodind[i] = indices[0][i];
	}
	
	long index;
	TPZGeoEl *elvec = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index); //AQUI
	
	gmesh->BuildConnectivity();
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec,4,-1);// condicao de fronteira tipo -1: (x,y=0)
	TPZGeoElBC gbc2(elvec,5,-2);// condicao de fronteira tipo -2: (x=1,y)
	TPZGeoElBC gbc3(elvec,6,-3);// condicao de fronteira tipo -3: (x,y=1)
	TPZGeoElBC gbc4(elvec,7,-4);// condicao de fronteira tipo -4: (x=0,y)
	
	const std::string nameref;
	
 TPZAutoPointer<TPZRefPattern> ref;
	//gmesh->RefPatternList(EQuadrilateral);
	
	//	Refinamento uniforme
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		long n = gmesh->NElements();
		for(long i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			
			
		}
		
		
	}
	
	
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif 		 		 
	return gmesh;
} 

TPZGeoMesh * MalhaGeoT(const int h){//malha triangulo
	
	
	
	TPZGeoMesh *gmesh = new TPZGeoMesh();
	
	//Criar ns
	const int nnode = 4;//AQUI
	const int nelem = 2;
	TPZGeoEl *elvec[nelem];	
	const int dim = 2;//AQUI
	
	REAL co[nnode][dim] = {{-1.,-1},{1.,-1},{1.,1.},{-1.,1.}};
	int indices[2][nnode];//como serao enumerados os nos
	
	//el 1
	indices[0][0] = 0;
	indices[0][1] = 1;
	indices[0][2] = 3;
	//el2
	indices[1][0] = 2;
	indices[1][1] = 3;
	indices[1][2] = 1;
	
	
	int nod;
	TPZVec<REAL> coord(dim);
	for(nod=0; nod<nnode; nod++) {
		long nodind = gmesh->NodeVec().AllocateNewElement();
		
		for(int d = 0; d < dim; d++)
		{
			coord[d] = co[nod][d];
		}
		gmesh->NodeVec()[nod].Initialize(nodind,coord,*gmesh);
	}
	//Criacao de elementos
	
	
	TPZVec<long> nodind1(3);
	TPZVec<long> nodind2(3);
	for(int i=0; i<3; i++){
		nodind1[i] = indices[0][i];
		nodind2[i] = indices[1][i];
	}
	
	long index;
	elvec[0] = gmesh->CreateGeoElement(ETriangle,nodind1,1,index); //AQUI
	elvec[1] = gmesh->CreateGeoElement(ETriangle,nodind2,1,index); //AQUI
	
	
	gmesh->BuildConnectivity();
	
	//	como usar os padroes de refinamento do caju..nao posso rodar ainda pq tenho o Pz antigo..nao tenho este gRefDBase.
	//{
	//				TPZAutoPointer< TPZRefPattern > refExemplo = gRefDBase->FinfRefPattern("a_Quad_Rib_Side_4_5_6");
	//				if(refExemplo)
	//				{
	//						TPZGeoEl * elExemplo = gmesh->Elementvec()[10];
	//						elExemplo->SetRefPattern(refExemplo);
	//						elEmemplo->Divide(filhos);
	//				}
	//		}
	
	//Cria as condicoes de contorno
	TPZGeoElBC gbc1(elvec[0],3,-1);// condicao de fronteira tipo -1: 
	TPZGeoElBC gbc2(elvec[0],5,-2);// condicao de fronteira tipo -2: 
	
	TPZGeoElBC gbc3(elvec[1],3,-3);// condicao de fronteira tipo -3: 
	TPZGeoElBC gbc4(elvec[1],5,-4);// condicao de fronteira tipo -4: 
	
	const std::string nameref;
	
	TPZAutoPointer<TPZRefPattern> ref;
	
	
	//	Refinamento uniforme
	/*for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		int n = gmesh->NElements();
		for(int i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			}
		
	}*/
	for(int ref = 0; ref < h; ref++){// h indica o numero de refinamentos
		TPZVec<TPZGeoEl *> filhos;
		long n = gmesh->NElements();
		for(long i = 0; i < n; i++){
			TPZGeoEl * gel = gmesh->ElementVec()[i];
			if(!gel->HasSubElement())
			{
				gel->Divide(filhos);
			}		
			
			
		}
		
		
	}
#ifdef LOG4CXX
    if (logger->isDebugEnabled())
	{
		std::stringstream sout;
		gmesh->Print(sout);
		LOGPZ_DEBUG(logger, sout.str());
	}
#endif 		 		 
	return gmesh;
}

void SaddlePermute(TPZCompMesh * cmesh){
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout<< "Implementando permutacao para problemas de ponto de sela"<< std::endl;
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
	TPZVec<long> permute;
	long numinternalconnects = cmesh->NIndependentConnects();
  	permute.Resize(numinternalconnects,0);

	TPZSubCompMesh *submesh = dynamic_cast<TPZSubCompMesh *> (cmesh);
	if(submesh)
	{
		long nexternal = submesh->NConnects();
		numinternalconnects -= nexternal;
	}
//	else {
//		DebugStop();
//	}

	int jperm=0;
	long nel=cmesh->ElementVec().NElements();
	for (long jel=0; jel<nel; jel++) {
		
		for (long ip=0; ip<permute.NElements(); ip++) {
			permute[ip]=ip;
		}
		
		TPZCompEl *elvec=cmesh->ElementVec()[jel];
	//	int idtroca=0;
		long eqmax=0;
		if(!elvec)continue;
		long ncon=elvec->NConnects();
	//	if(ncon==1) continue;
		long eqpress=elvec->Connect(ncon-1).SequenceNumber();
		for (long icon=0; icon< ncon-1; icon++) {
			TPZConnect &coel=elvec->Connect(icon);
			long eqflux=coel.SequenceNumber();
			if (eqflux >= numinternalconnects) {
				continue;
			}
			eqmax = max(eqmax,eqflux);
		}
		
		
		if(eqpress<eqmax) {

			permute[eqpress]=eqmax;
			
		}
		
		
		for ( jperm = eqpress+1; jperm<=eqmax; jperm++) {
			permute[jperm]=jperm-1;
			
		}
		/*
#ifdef LOG4CXX
         if (logger->isDebugEnabled())
		{
			std::stringstream sout;
			sout << "vetor SaddlePermute  do elemento - "<<jel<< " - " <<permute;
			LOGPZ_DEBUG(logger, sout.str().c_str());
		}
#endif
	*/
		cmesh->Permute(permute);
		
	}		
	
}