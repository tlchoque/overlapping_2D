

#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgnode.h"
#include "pzsolve.h"
#include "pzelg1d.h"
#include "pzmat1dlin.h"
#include "pzelmat.h"
#include "pzelasmat.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include <stdlib.h>
#include <iostream.h>
#include "pzvec.h"

#include "pzelgq2d.h"
#include "pzelgt2d.h"
#include "pzmat2dlin.h"
#include "pzanalysis.h"
#include "pzmetis.h"
#include "pzplaca.h"

#include <stdio.h>
#include <time.h>
#include "pzelct2d.h"
//template<class T>
//class TPZVec;
//#define NOTDEBUG

TPZMaterial *LerMaterial(char *filename);
void PressaoHid(TPZVec<REAL> &x,TPZVec<REAL> &force);



int main() {

   //malha geometrica
   TPZGeoMesh *firstmesh = new TPZGeoMesh;
   firstmesh->SetName("Malha Geometrica : N�s e Elementos");
   firstmesh->NodeVec().Resize(10);
   TPZVec<REAL> coord(2);   //,coordtrans(2);
//   REAL ct,st,PI=3.141592654;
//   cout << "\nEntre rotacao do eixo n1 (graus) -> ";
//   REAL g;
//   cin >> g;
//   g = g*PI/180.;
//   ct = cos(g);
//   st = sin(g);
//ct = 1.;
//st = 0.;

   //nos geometricos

   //no 0
   coord[0] = 0;
   coord[1] = 0;
// coordtrans[0] =  ct*coord[0]-st*coord[1];
// coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[0].Initialize(coord,*firstmesh);

   //no 1
   coord[0] = 1.;
   coord[1] = 0;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[1].Initialize(coord,*firstmesh);

   //no 2
   coord[0] = 2.;
   coord[1] = 0.;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[2].Initialize(coord,*firstmesh);

   //no 3
   coord[0] = 3.;
   coord[1] = 0.;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[3].Initialize(coord,*firstmesh);

   //no 4
   coord[0] = 4;
   coord[1] = 0.;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[4].Initialize(coord,*firstmesh);

   //no 5
   coord[0] = 0;
   coord[1] = 1;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[5].Initialize(coord,*firstmesh);

   //no 6
   coord[0] = 1.;
   coord[1] = 1.;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[6].Initialize(coord,*firstmesh);

   //no 7
   coord[0] = 2.;
   coord[1] = 1.;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[7].Initialize(coord,*firstmesh);

   //no 8
   coord[0] = 3.;
   coord[1] = 1;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//  coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[8].Initialize(coord,*firstmesh);

   //no 9
   coord[0] = 4;
   coord[1] = 1;
//   coordtrans[0] =  ct*coord[0]-st*coord[1];
//   coordtrans[1] =  st*coord[0]+ct*coord[1];
   firstmesh->NodeVec()[9].Initialize(coord,*firstmesh);

   /*
   TPZVec<int> nodeindexes(3);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
    //elementos geometricos
   TPZGeoElT2d *elg0 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 0;
   nodeindexes[1] = 2;
   nodeindexes[2] = 3;
   TPZGeoElT2d *elg1 = new TPZGeoElT2d(nodeindexes,1,*firstmesh);
   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 2;
   TPZGeoElT2d *elg2 = new TPZGeoElT2d(nodeindexes,2,*firstmesh);
   nodeindexes[0] = 0;
   nodeindexes[1] = 2;
   nodeindexes[2] = 3;
   TPZGeoElT2d *elg3 = new TPZGeoElT2d(nodeindexes,2,*firstmesh);
*/
   TPZVec<int> nodeindexes(4);
    //elementos geometricos
   TPZGeoEl *elg[4];

   nodeindexes[0] = 0;
   nodeindexes[1] = 1;
   nodeindexes[2] = 6;
   nodeindexes[3] = 5;
   elg[0] = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);

   nodeindexes[0] = 1;
   nodeindexes[1] = 2;
   nodeindexes[2] = 7;
   nodeindexes[3] = 6;
   elg[1] = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);

   nodeindexes[0] = 2;
   nodeindexes[1] = 3;
   nodeindexes[2] = 8;
   nodeindexes[3] = 7;
   elg[2] = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);

   nodeindexes[0] = 3;
   nodeindexes[1] = 4;
   nodeindexes[2] = 9;
   nodeindexes[3] = 8;
   elg[3] = new TPZGeoElQ2d(nodeindexes,1,*firstmesh);

   //Arquivos de saida
   ofstream outgm1("outgm1.dat");
   ofstream outcm1("outcm1.dat");
   ofstream outcm2("outcm2.dat");

   //montagem de conectividades entre elementos
   firstmesh->BuildConnectivity();
   //malha computacional
   TPZCompMesh *secondmesh = new TPZCompMesh(firstmesh);
   secondmesh->SetName("Malha Computacional : Conectividades e Elementos");


   //material
   TPZMaterial *pl = LerMaterial("flavio.dat");
   secondmesh->InsertMaterialObject(pl);
//   pl = LerMaterial("placa2.dat");
//   secondmesh->InsertMaterialObject(pl);
//   pl = LerMaterial("placa3.dat");
//   secondmesh->InsertMaterialObject(pl);

   // carregamento hidrostatico no plano vertica xz
   pl->SetForcingFunction(PressaoHid);

   //CC : condic�es de contorno
   TPZBndCond *bc;
   REAL big = 1.e12;
   TPZFMatrix val1(6,6,0.),val2(6,1,0.);

   // engastes nos lados 4 e 7 do elemento 0
   TPZGeoElBC(elg[0],4,-2,*firstmesh);
   TPZGeoElBC(elg[0],7,-2,*firstmesh);

   // engaste no lado 4 do elemento 1
   TPZGeoElBC(elg[1],4,-2,*firstmesh);

   // engaste no lado 4 do elemento 2
   TPZGeoElBC(elg[2],4,-2,*firstmesh);

   // engaste no lado 4 do elemento 3
   TPZGeoElBC(elg[3],4,-2,*firstmesh);

   // imposicao do valor zero associado a condicao -2 (engaste)
   bc = pl->CreateBC(-2,0,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   // imposicao da condicao de simetria no lado 5 do elemento 4
   val1(0,0)=big;
   val1(1,1)=0.;
   val1(2,2)=0.;
   val1(3,3)=0.;
   val1(4,4)=big;
   val1(5,5)=big;
   TPZGeoElBC(elg[3],5,-3,*firstmesh);
   bc = pl->CreateBC(-3,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   //ordem de interpolacao
   int ord;
   cout << "Entre ordem 1,2,3,4,5 : ";
   cin >> ord;
//   TPZCompEl::gOrder = ord;
   firstmesh.SetDefaultOrder(order);
   //constru��o malha computacional
   TPZVec<int> csub(0);
   TPZManVector<TPZGeoEl *> pv(4);
   int n1=1,level=0;
   cout << "\nDividir ate nivel ? ";
   int resp;
   cin >> resp;
   int nelc = firstmesh->ElementVec().NElements();
   int el;
   TPZGeoEl *cpel;
   for(el=0;el<firstmesh->ElementVec().NElements();el++) {
     cpel = firstmesh->ElementVec()[el];
     if(cpel && cpel->Level() < resp)
		cpel->Divide(pv);

   }
   cout << "\nDividir o elemento esquerdo superior quantas vezes? ";
   cin >> resp;
   cpel = firstmesh->ElementVec()[0];
   for(el=0; el<resp; el++) {
		cpel->Divide(pv);
		cpel = pv[3];
   }
   //analysis
   secondmesh->AutoBuild();
   firstmesh->Print(outgm1);
   outgm1.flush();
   secondmesh->AdjustBoundaryElements();
   secondmesh->InitializeBlock();
   secondmesh->Print(outcm1);
   TPZAnalysis an(secondmesh,outcm1);
   int numeq = secondmesh->NEquations();
   secondmesh->Print(outcm1);
   outcm1.flush();
   TPZVec<int> skyline;
   secondmesh->Skyline(skyline);
   TPZSkylMatrix *stiff = new TPZSkylMatrix(numeq,skyline);
   an.SetMatrix(stiff);
   an.Solver().SetDirect(ECholesky);
   secondmesh->SetName("Malha Computacional :  Connects e Elementos");
   // Posprocessamento
   an.Run(outcm2);
   TPZVec<char *> scalnames(5);
   scalnames[0] = "Mn1";
   scalnames[1] = "Mn2";
   scalnames[2] = "Vn1";
   scalnames[3] = "Vn2";
   scalnames[4] = "Deslocz";
   TPZVec<char *> vecnames(0);
   char plotfile[] =  "placaPos.pos";
   char pltfile[] =  "placaView.plt";
   an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
   an.Print("FEM SOLUTION ",outcm1);
   an.PostProcess(3);
   an.DefineGraphMesh(2, scalnames, vecnames, pltfile);
   an.PostProcess(2);
   firstmesh->Print(outgm1);
   outgm1.flush();
   delete secondmesh;
   delete firstmesh;
   return 0;
}

// ****************** fim do programa principal ****************

// rotina para calculo de pressao hidrostatica parede
// calcula pressao hidrostatica aplicada no plano vertical xz

void PressaoHid(TPZVec<REAL> &x,TPZVec<REAL> &force){
	force[2] = -1.+x[1];
}

// ****************** fim da rotina PressaoHid *****************


// ******************  rotina LerMaterial **********************

TPZMaterial *LerMaterial(char *filename) {
	ifstream input(filename);
   TPZFMatrix naxes(3,3);
   REAL ni1,ni2,h,E1,E2,G12,G13,G23,f;
   REAL n00,n01,n02,n10,n11,n12,n20,n21,n22;
   TPZVec<REAL> xf(6);
   int matindex;
   input >> matindex;
   input >> f  >>  h  >>
   		  E1  >>  E2 >>
           G12 >> G13 >> G23 >>
           ni1 >> ni2;
   input >> n00 >> n01 >> n02 >> n10 >> n11 >> n12 >> n20 >> n21 >> n22;
	input >> xf[0] >> xf[1] >> xf[2] >> xf[3] >> xf[4] >> xf[5];
   naxes(0,0) =  n00;    naxes(0,1) =  n01;    naxes(0,2) =  n02;
   naxes(1,0) =  n10;    naxes(1,1) =  n11;    naxes(1,2) =  n12;
   naxes(2,0) =  n20;    naxes(2,1) =  n21;    naxes(2,2) =  n22;
   return new TPZPlaca(matindex,h,f,E1,E2,ni1,ni2,G12,G13,G23,naxes,xf);


}


/*
TESTE 1 TESTE 1 TESTE 1 TESTE 1 TESTE 1 TESTE 1 TESTE 1

   val1(0,0)=big;
   val1(1,1)=big;
   val1(2,2)=big;
   val1(3,3)=big;
   val1(4,4)=big;
   val1(5,5)=big;
 	TPZGeoElBC(elg0,7,-2,*firstmesh);
   bc = pl->CreateBC(-2,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

// 	TPZGeoElBC(elg0,0,-3,*firstmesh);//equivale a val1(1,1)=big;
//   bc = pl->CreateBC(-3,0,val1,val2);
//   secondmesh->InsertMaterialObject(bc);

   val1.Zero();
   val2(0,0) = 100.;
	TPZGeoElBC(elg0,5,-1,*firstmesh);
   bc = pl->CreateBC(-1,1,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1(4,4)=0;
   val1(2,2)=big;
   val1(3,3)=big;
   val1(5,5)=big;
	TPZGeoElBC(elg0,5,-2,*firstmesh);
   bc = pl->CreateBC(-2,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1(3,3)=0;
   val1(2,2)=big;
   val1(4,4)=big;
   val1(5,5)=big;
 	TPZGeoElBC(elg0,6,-3,*firstmesh);
   bc = pl->CreateBC(-3,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);

   val1(0,0)=big;
   val1(1,1)=big;
   val1(2,2)=big;
   val1(3,3)=big;
   val1(4,4)=big;
   val1(5,5)=big;
	TPZGeoElBC(elg0,7,-4,*firstmesh);
   bc = pl->CreateBC(-4,2,val1,val2);
   secondmesh->InsertMaterialObject(bc);
*/
/*
   //redistribuicao de ordem aos lados do elemento
	int nel = secondmesh->ElementVec().NElements();
   cout << "\nEntre Ordem por Elemento\n";
   for(int index = 0;index<nel;index++) {
      TPZInterpolatedElement *el;
		el = ((TPZInterpolatedElement *) secondmesh->ElementVec()[index]);
      cout << "\nElemento Geometrico de Id " << el->Reference()->Id() << " : ordem -> ";
      cin >> ord;
      if(el) el->PRefine(ord);
   }
 */
/*   int iel = 1;
   while(iel) {
      qsi[0] = 1.;
      qsi[1] = 1.;
      cout << "\nEntre id do elemento\n";
      cin >> iel;
	  if(iel<0) break;
      TPZInterpolatedElement *el = (TPZInterpolatedElement *) secondmesh->ElementVec()[iel];
      if(el) {
		  el->Reference()->X(qsi,x);
		  el->Solution(qsi,4,sol);
	  }
	  cout << "coordenadas x = "<< x[0] << ' ' << x[1] << ' ' << x[2] << endl;
      outcm1 << "\n\n\nFLECHA NO CENTRO DA PLACA\n";
	  outcm1 << "coordenadas x = "<< x[0] << ' ' << x[1] << ' ' << x[2] << endl;
      outcm1 << "PZ        : " << sol[0];
      //REAL x = .0334*(-1./(210000.*.05*.05*.05));

      REAL xs = -sol[0]*210000.*.001*100.0 ;
      outcm1<< "\nTABLA 1.2 : " << xs;
   }
*/