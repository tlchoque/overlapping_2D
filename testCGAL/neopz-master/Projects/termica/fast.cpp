//
// C++ Implementation: fast
//
// Description: 
//
//
// Author: Philippe Remy Bernard Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "fast.h"
#include <fstream>
#include <vector>
#include <map>
#include <string.h>
#include <stdlib.h>
#ifdef VC
#include <time.h>
#else
#include <sys/time.h>
#endif
#include "pzfmatrix.h"
#include "pzysmp.h"
#include "pzstack.h"
#include "pzstepsolver.h"
#include "TPZCopySolve.h"

using namespace std;

void RunFast()
{
  
//  string filename("CursoAltoDes.in");
//  string filename("matriz1089.in");
//  string filename("matriz12.in");
//  string filename("matriz16000.in");
  string filename("matriz260000.in");
  TPZFMatrix<STATE> rhs,sol;
  
  TPZMultiTimer timer(5);
  cout << "Leitura\n";
  timer.processName(0) = "Leitura da matriz";
  timer.start(0);
  TPZFYsmpMatrix<STATE> *faster = ReadMatrix(filename,rhs);
  timer.stop(0);
  
  cout << "Multiplicacao Numero de termos " << faster->NumTerms() << "\n";
  TimeMultiply(faster,timer);
  
  cout << "Jacobi\n";
  SolveJacobi(faster, rhs, 1.e-8,timer);
  
  cout << "SSor\n";
  SolveSSOR(faster,rhs,1.e-8,timer);
  
  cout << "CG\n";
  SolveCG(faster,rhs,1.e-8,timer);
  
  cout << timer;

}
// this function will read the matrix
TPZFYsmpMatrix<STATE> *ReadMatrix(const std::string &filename, TPZFMatrix<STATE> &rhs)
{
  ifstream input(filename.c_str());
    long neq,nelem;
    int issymetric;
  input >> neq >> nelem >> issymetric;
  rhs.Redim(neq,1);
  TPZStack<int,10000> elgraph;
  TPZManVector<int> elgraphindex(nelem+1,-1);
  elgraphindex[0] = 0;
  std::vector<int> numelcon(neq,0);
// leitura do grafo dos elementos
  long iel;
  for(iel=0; iel<nelem; iel++)
  {
    int nnodes;
    input >> nnodes;
    elgraphindex[iel+1] = elgraphindex[iel]+nnodes;
    int in,nod;
    for(in=0; in<nnodes; in++) 
    {
      input >> nod;
      elgraph.Push(nod);
      numelcon[nod]++;
    }    
  }
  TPZStack<int,100000> nodegraph,nodegraphindex;
  TPZStack<REAL,100000> avec;
  nodegraph.Expand(1000000);
  avec.Expand(1000000);
  nodegraphindex.Resize(neq+1,-1);
  nodegraphindex[0] = 0;
  
  
  std::vector< std::map<int,REAL> > matrixvalues(neq);

  int nextassemble = 0;
  for(iel = 0; iel<nelem; iel++) 
  {  
    int nnodes = elgraphindex[iel+1]-elgraphindex[iel];
    int in,jn;
    int start = elgraphindex[iel];
    for(in=0; in<nnodes; in++)
    {
      int ing = elgraph[start+in];
      numelcon[ing]--;
      for(jn=0; jn<nnodes; jn++)
      {
        int jng = elgraph[start+jn];
        REAL value;
        input >> value;
        matrixvalues[ing][jng] += value;
      }
    }
    for(in=0; in<nnodes; in++)
    {
      int ing = elgraph[start+in];
      REAL value;
      input >> value;
      rhs(ing,0) += value;
    }

    while(nextassemble < neq && !numelcon[nextassemble])
    {
      std::map<int,REAL>::iterator it;
      nodegraphindex[nextassemble+1] = nodegraphindex[nextassemble]+matrixvalues[nextassemble].size();  
      for(it=matrixvalues[nextassemble].begin();it!=matrixvalues[nextassemble].end(); it++)
      {
        nodegraph.Push((*it).first);
        avec.Push((*it).second);
      }
      matrixvalues[nextassemble].clear();
      nextassemble++;
    }
  }  
  long *IA, *JA;
  STATE *A;
  IA = new long[neq+1];
  JA = new long[nodegraph.NElements()];
  A = new STATE[nodegraph.NElements()];
  memcpy(IA,&nodegraphindex[0],(neq+1)*sizeof(long));
  memcpy(JA,&nodegraph[0],nodegraph.NElements()*sizeof(long));
  memcpy(A,&avec[0],avec.NElements()*sizeof(REAL));
  
  TPZFYsmpMatrix<STATE> *result = new TPZFYsmpMatrix<STATE>(neq,neq);
  result->SetData(IA,JA,A);
  return result;
}

void TimeMultiply(TPZMatrix<STATE> *mat, TPZMultiTimer &timer)
{
  TPZFMatrix<STATE> sol(mat->Rows(),1),rhs(mat->Rows(),1),y;
  long ieq,neq = mat->Rows();
  for(ieq=0; ieq<neq; ieq++) sol(ieq,0) = rand()/2147483647.;
  
  int nummult = 500, imult;
  
  timer.processName(1) = "Produto matriz vetor";
  timer.start(1);
  //start = clock();
  for(imult = 0; imult<nummult; imult++)
  {
    mat->MultAdd(sol,y,rhs,1.,0.);
  }
  timer.stop(1);
  
  return;
}

void SolveJacobi(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> &rhs, REAL tol,TPZMultiTimer &timer)
{
  TPZFMatrix<STATE> sol(rhs.Rows(),1);
  TPZStepSolver<STATE> prec(mat);
  prec.SetJacobi(1,0.,0);
  TPZStepSolver<STATE> step(0);
  step.ShareMatrix(prec);
  step.SetCG(3000,prec,tol,0);
  timer.processName(3) = "Solucao Jacobi";
  timer.start(3);
  step.Solve(rhs,sol,0);  
  timer.stop(3);
  step.ResetMatrix();
  return;
}

void SolveCG(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> &rhs, REAL tol, TPZMultiTimer &timer)
{
  TPZFMatrix<STATE> sol(rhs.Rows(),1);
  TPZCopySolve<STATE> prec(mat);
  TPZStepSolver<STATE> step(0);
  step.ShareMatrix(prec);
  step.SetCG(3000,prec,tol,0);
  timer.processName(2) = "Solucao CG";
  timer.start(2);
  step.Solve(rhs,sol,0);  
  timer.stop(2);
  step.ResetMatrix();
  return;
}

void SolveSSOR(TPZMatrix<STATE> *mat, TPZFMatrix<STATE> &rhs, REAL tol,TPZMultiTimer &timer)
{
  TPZFMatrix<STATE> sol(rhs.Rows(),1);
  TPZStepSolver<STATE> prec(mat);
  prec.SetSSOR(1,1.,0.,0);
  TPZStepSolver<STATE> step(0);
  step.ShareMatrix(prec);
  step.SetCG(3000,prec,tol,0);
  timer.processName(4) = "Solucao SSOR";
  timer.start(4);
  step.Solve(rhs,sol,0);  
  timer.stop(4);
  step.ResetMatrix();
  return;
}

void Compare(TPZMatrix<STATE> *first, TPZMatrix<STATE> *second)
{
  TPZFMatrix<STATE> sol(first->Rows(),1),rhs1(first->Rows(),1),rhs2(first->Rows(),1),y;
  TPZFMatrix<STATE> sol1(rhs1),sol2(rhs1);
  long ieq,neq = first->Rows();
  for(ieq=0; ieq<neq; ieq++) sol(ieq,0) = rand()/2147483647.;
  
  first->MultAdd(sol,y,rhs1,1.,0.);
  second->MultAdd(sol,y,rhs2,1.,0.);
  rhs1.Print("First");
  rhs2.Print("Second");
  
  rhs1 -= rhs2;
  
  rhs1.Print("Difference");
  
  TPZStepSolver<STATE> prec1(first);
  prec1.SetSSOR(1,1.,0.,0);
  TPZStepSolver<STATE> prec2(second);
  prec2.SetSSOR(1,1.,0.,0);
  
  prec1.Solve(rhs2,sol1);
  
  prec2.Solve(rhs2,sol2);

  sol1.Print("First solution");
  sol2.Print("Second solution");
  
  sol1 -= sol2;
  
  sol1.Print("Difference solution");
  
  prec1.ResetMatrix();
  prec2.ResetMatrix();
}

