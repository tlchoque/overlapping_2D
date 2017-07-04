
#include "meshgen.h"

#include "pzgmesh.h"
#include "pzgengrid.h"
#include "pzgeoel.h"
#include "TPZRefPatternTools.h"
#include "pzcheckgeom.h"
#include "TPZVTKGeoMesh.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZSSpStructMatrix.h"

#include "pzlog.h"

#ifdef _AUTODIFF
#include "fadType.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif

template<class TVar>
void TElasticityExample1::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp)
{
    disp[0] = (TVar)0.01*x[0]*x[0];
    TVar a = TVar(M_PI*10.)*x[1];
    disp[1] = (TVar)0.05*x[0]+((TVar)0.03)*x[1]*x[1]*sin(a);
    disp[0] = x[0]*x[0];
    disp[1] = x[1]*x[1];
}

template<>
void TElasticityExample1::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp)
{
    FADFADREAL tmp(2,Fad<REAL>(2,0.));
    disp[0] = tmp*x[0]*x[0];
    FADFADREAL a = FADFADREAL(M_PI*10.)*x[1];
    FADREAL_ sinaval = sin(a.val());
    FADREAL_ cosaval = cos(a.val());
    FADFADREAL sina(2,sinaval);
    for (int i=0; i<2; i++) {
        sina.fastAccessDx(i) = cosaval*a.dx(i);
    }
    disp[1] = (FADFADREAL)0.05*x[0]+((FADFADREAL)0.03)*x[1]*x[1]*sina;
    disp[0] = x[0]*x[0];
    disp[1] = x[1]*x[1];
}

template<class TVar>
void TElasticityExample1::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu)
{
    Elast.val() = 1000.;
    nu.val() = 0.3;
}

template<>
void TElasticityExample1::Elastic(const TPZVec<Fad<double> > &x, Fad<double>  &Elast, Fad<double>  &nu)
{
    Elast = Fad<double> (2,1000.);
    nu = Fad<double> (2,0.3);
}

template<class TVar>
void TElasticityExample1::graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<TVar> temp = Fad<TVar>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(i,j) = result[i].d(j);
        }
    }
    for(int i=0; i<2; i++)
    {
        std::cout << "result " << result[i] << " dx " << result[i].dx() << std::endl;
        for (int j=0; j<2; j++) {
            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " dx " << grad(i,j).dx() << std::endl;
        }
    }
}

void TElasticityExample1::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu)
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>,3> result(2);
    uxy(xfad,result);
    gradu.Resize(2,2);
    u[0] = result[0].val();
    u[1] = result[1].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            gradu(i,j) = result[i].d(j);
        }
    }
    
}


template<>
void TElasticityExample1::graduxy(const TPZVec<Fad<double> > &x, TPZFMatrix<Fad<double> > &grad)
{
    TPZManVector<Fad<Fad<double> >,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<Fad<double> > temp = Fad<Fad<double> >(2,Fad<double>(2,0.));
        temp.val()= x[i];
        Fad<double> temp2(2,1.);
        temp.fastAccessDx(i) = temp2;
        xfad[i] = temp;
    }
    TPZManVector<Fad<Fad<double> >,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(i,j) = result[i].d(j);
        }
    }
    for(int i=0; i<2; i++)
    {
        std::cout << "result " << result[i] << " dx " << result[i].dx() << std::endl;
        for (int j=0; j<2; j++) {
            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " dx " << grad(i,j).dx() << std::endl;
        }
    }
}

template<class TVar>
void TElasticityExample1::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma)
{
    TPZFNMatrix<4,TVar> grad;
    TVar E, nu;
    Elastic(x, E, nu);
    TVar Fac = E/((TVar) 1.+nu)/(((TVar) 1.-2.*nu));
    graduxy(x,grad);
    sigma.Resize(2,2);
    sigma(0,0) = Fac*((1.-nu)*grad(0,0)+nu*grad(1,1));
    sigma(1,1) = Fac*((1.-nu)*grad(1,1)+nu*grad(0,0));
    sigma(0,1) = E/(2.*(1.+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
    for(int i=0; i<2; i++)
    {
        for (int j=0; j<2; j++) {
            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " sigma " << sigma(i,j) << std::endl;
        }
    }
}

template<>
void TElasticityExample1::Sigma(const TPZVec<Fad<double> > &x, TPZFMatrix<Fad<double> > &sigma)
{
    TPZFNMatrix<4,Fad<double> > grad;
    Fad<double>  E, nu;
    Elastic(x, E, nu);
    Fad<double>  Fac = E/((Fad<double> ) 1.+nu)/(((Fad<double> ) 1.-2.*nu));
    graduxy(x,grad);
    sigma.Resize(2,2);
    sigma(0,0) = Fac*((1.-nu)*grad(0,0)+nu*grad(1,1));
    sigma(1,1) = Fac*((1.-nu)*grad(1,1)+nu*grad(0,0));
    sigma(0,1) = E/(2.*(1.+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
    for(int i=0; i<2; i++)
    {
        for (int j=0; j<2; j++) {
            std::cout << "i = " << i << " j = " << j << " grad " << grad(i,j) << " sigma " << sigma(i,j) << std::endl;
        }
    }
}

template<class TVar>
void TElasticityExample1::DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma)
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZFNMatrix<4, Fad<TVar> > sigma(2,2);
    Sigma(xfad,sigma);
    for(int i=0; i<2; i++)
    {
        for (int j=0; j<2; j++) {
            std::cout << "i = " << i << " j = " << j <<  " sigma " << sigma(i,j) <<  " " << sigma(i,j).dx() << std::endl;
        }
    }

    divsigma[0] = sigma(0,0).dx(0)+sigma(0,1).dx(1);
    divsigma[1] = sigma(1,0).dx(0)+sigma(1,1).dx(1);
}

TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::ForcingFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Force);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
}

TPZAutoPointer<TPZFunction<STATE> > TElasticityExample1::DirichletFunction()
{
    TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(Dirichlet);
    dummy->SetPolynomialOrder(5);
    TPZAutoPointer<TPZFunction<STATE> > result(dummy);
    return result;
    
}


template
void TElasticityExample1::DivSigma<REAL>(const TPZVec<REAL> &x, TPZVec<REAL> &divsigma);
template
void TElasticityExample1::Sigma<Fad<REAL> >(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma);


#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.meshgen"));
#endif

TPZGeoMesh *MalhaGeomFredQuadrada(int nelx, int nely, TPZVec<REAL> &x0, TPZVec<REAL> &x1, TPZVec<long> &coarseindices, int ndiv)
{
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    int dimension = 2;
    gmesh->SetDimension(dimension);
    TPZManVector<int,2> nx(2,3);
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid gengrid(nx, x0, x1);
    gengrid.SetRefpatternElements(true);
    gengrid.Read(gmesh, 1);
    gengrid.SetBC(gmesh, 4, -1);
    gengrid.SetBC(gmesh, 5, -2);
    gengrid.SetBC(gmesh, 6, -3);
    gengrid.SetBC(gmesh, 7, -4);
    
    long nel = gmesh->NElements();
    
    coarseindices.resize(nel);
    long elcount = 0;
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (gel->HasSubElement() ||  gel->Dimension() != dimension) {
            continue;
        }
        coarseindices[elcount] = el;
        elcount++;
    }
    coarseindices.resize(elcount);
    
    if(0)
    {
        
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == EQuadrilateral) {
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
        nel = gmesh->NElements();
        
        for (long el=0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if (!gel->HasSubElement() &&  gel->Type() == EOned) {
                TPZAutoPointer<TPZRefPattern> refpat = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if (!refpat) {
                    DebugStop();
                }
                gel->SetRefPattern(refpat);
                TPZManVector<TPZGeoEl *,12> subs;
                gel->Divide(subs);
            }
        }
    }
    TPZCheckGeom geom(gmesh);
    geom.UniformRefine(ndiv);
    //    InsertInterfaceElements(gmesh,1,2);
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
#ifdef PZDEBUG
    {
        std::ofstream file("GMeshFred.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, file);
    }
#endif
    
    return gmesh;
}

void SolveProblem(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<TPZAutoPointer<TPZCompMesh> > compmeshes, std::string prefix, std::string configuration)
{
    //calculo solution
    bool shouldrenumber = false;
    TPZAnalysis an(cmesh,shouldrenumber);
#ifdef USING_MKL
    TPZSymetricSpStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
    
#else
    TPZSkylineStructMatrix strmat(cmesh.operator->());
    strmat.SetNumThreads(0);
#endif
    
    
    an.SetStructuralMatrix(strmat);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    std::cout << "Assembling\n";
    an.Assemble();
    if(1)
    {
        std::string filename = prefix;
        filename += "_Global.nb";
        std::ofstream global(filename.c_str());
        TPZAutoPointer<TPZStructMatrix> strmat = an.StructMatrix();
        an.Solver().Matrix()->Print("Glob = ",global,EMathematicaInput);
        an.Rhs().Print("Rhs = ",global,EMathematicaInput);
    }
    std::cout << "Solving\n";
    an.Solve();
    std::cout << "Finished\n";
    an.LoadSolution(); // compute internal dofs
                       //    an.Solution().Print("sol = ");
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(compmeshes, cmesh);
    
#ifdef PZDEBUG
    {
        std::ofstream out(prefix+"_MeshWithSol.txt");
        cmesh->Print(out);
    }
#endif
    
    //    TPZBuildMultiphysicsMesh::TransferFromMeshes(cmeshes, an.Mesh());
    //    for (int i=0; i<cmeshes.size(); i++) {
    //        cmeshes[i]->Solution().Print("sol = ");
    //    }
    //    cmeshes[0]->Solution().Print("solq = ");
    //    cmeshes[1]->Solution().Print("solp = ");
    std::stringstream sout;
    sout << prefix << "Approx-";
    sout << configuration << ".vtk";
    std::string plotfile = sout.str();
    std::cout << "plotfile " << plotfile.c_str() << std::endl;
    TPZStack<std::string> scalnames,vecnames;
    TPZMaterial *mat = cmesh->FindMaterial(1);
    if (!mat) {
        DebugStop();
    }
    if (mat->NStateVariables() == 2)
    {
        scalnames.Push("SigmaX");
        scalnames.Push("SigmaY");
        scalnames.Push("TauXY");
        vecnames.Push("Displacement");
    }
    else if(mat->NStateVariables() == 1)
    {
        scalnames.Push("Pressure");
        vecnames.Push("Flux");
        vecnames.Push("Derivative");
    }
    an.DefineGraphMesh(cmesh->Dimension(), scalnames, vecnames, plotfile);
    int resolution = 0;
    an.PostProcess(resolution,cmesh->Dimension());
}

