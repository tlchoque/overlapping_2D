

#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgnode.h"
#include "pzsolve.h"
#include "pzelg1d.h"
#include "pzmat1dlin.h"
#include "pzelmat.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include <stdlib.h>
#include <iostream.h>

#include "pzelgq2d.h"//Cedric

#define NOTDEBUG

int main() {
   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->NodeVec().Resize(4);
   TPZVec<REAL> coord(2);
   coord[0] = 0.;
   coord[1] = 0.;
   //nos geometricos
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);
   coord[0] = 1.0;
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);
   coord[1] = 1.0;
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);
   coord[0] = 0.0;
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);
   TPZVec<int> nodeindexes(4);
   nodeindexes[0] = 0;//local[i] = global[i] , i=0,1,2,3
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
   nodeindexes[3] = 3;
   //elementos geometricos
   TPZGeoElQ2d *elq1 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);
   int i,sen;;
   cout<<"Entre i = 0,1,2,3 ";
   cin>>i;
   cout<<"Sentido local direito/inverso : 0/1 ?  ";
   cin>>sen;
   if(sen==0) {//direito
        nodeindexes[0] = (0+i)%4;//local[i] = global[j] , i,j em {0,1,2,3}
        nodeindexes[1] = (1+i)%4;
        nodeindexes[2] = (2+i)%4;
        nodeindexes[3] = (3+i)%4;
	} else {//inverso
        nodeindexes[0] = (3+i)%4;//local[i] = global[j] , i,j em {0,1,2,3}
        nodeindexes[1] = (2+i)%4;
        nodeindexes[2] = (1+i)%4;
        nodeindexes[3] = (0+i)%4;
   }
   TPZGeoElQ2d *elq2 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);
   cout<<"Entre i = 0,1,2,3 ";
   cin>>i;
   cout<<"Sentido local direito/inverso : 0/1 ?  ";
   cin>>sen;
   if(sen==0) {//direito
        nodeindexes[0] = (0+i)%4;//local[i] = global[j] , i,j em {0,1,2,3}
        nodeindexes[1] = (1+i)%4;
        nodeindexes[2] = (2+i)%4;
        nodeindexes[3] = (3+i)%4;
	} else {//inverso
        nodeindexes[0] = (3+i)%4;//local[i] = global[j] , i,j em {0,1,2,3}
        nodeindexes[1] = (2+i)%4;
        nodeindexes[2] = (1+i)%4;
        nodeindexes[3] = (0+i)%4;
   }
   TPZGeoElQ2d *elq3 = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);
   //Arquivos de saida
	ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
	ofstream outcm2("outcm2.dat");
   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
 	firstmesh->Print(outgm1);
   outgm1.flush();
  	//teste de divisao geometrica
   //TPZAdmChunkVector<TPZGeoEl *> &elements = firstmesh->ElementVec();
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   //material
   int matindex = secondmesh->MaterialVec().AllocateNewElement();
   TPZFMatrix b(1,1,1.),c(1,1,0.);
   TPZMat1dLin *mat = new TPZMat1dLin(1);
   mat->SetMaterial(c,c,b,b);
   secondmesh->MaterialVec()[matindex] = mat;
   //CC : condicao de contorno
/* 	TPZGeoElBC bcel1(el1,2,-1);
   int bcindex = firstmesh->BCElementVec().AllocateNewElement();
   firstmesh->BCElementVec()[bcindex] = bcel1;
	TPZFMatrix val1(1,1,0.), val2(1,1,1.);
   TPZBndCond *bc = mat->CreateBC(-1,0,val1,val2);
   bcindex = secondmesh->BndCondVec().AllocateNewElement();
   secondmesh->BndCondVec()[bcindex] = bc;        */
   //constroe a malha computacional
   secondmesh->AutoBuild();
   secondmesh->InitializeBlock();
   secondmesh->ComputeConnectSequence();
   secondmesh->Print(outcm1);
   outcm1.flush();
   firstmesh->Print(outgm1);
	//Resolucao do sistema
   TPZFMatrix Rhs(secondmesh->NEquations(),1),Stiff(secondmesh->NEquations(),secondmesh->NEquations()),U;
   Stiff.Zero();
   Rhs.Zero();
   secondmesh->Assemble(Stiff,Rhs);
   Rhs.Print("Rhs teste",outcm2);
   Stiff.Print("Bloco teste",outcm2);
	Rhs.Print("Computational Mesh -> fBlock",outcm2);
   TPZMatrixSolver solver(&Stiff);
   solver.SetDirect(ELU);
   solver.Solve(Rhs,U);
   U.Print("Resultado",outcm2);
   secondmesh->LoadSolution(U);
   secondmesh->Solution().Print("Mesh solution ",outcm2);
//   TPZElementMatrix ek,ef;
//   secondmesh->ElementVec()[0]->CalcStiff(ek,ef);
//	ek.fMat->Print();
//   ef.fMat->Print();
   delete secondmesh;
   delete firstmesh;
   return 0;
}
//FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN FIM MAIN
   //if(!elq) elq = (TPZGeoElQ2d *) firstmesh->ElementVec()[0];
