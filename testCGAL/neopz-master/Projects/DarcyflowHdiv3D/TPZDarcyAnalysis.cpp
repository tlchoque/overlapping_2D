
/*
 *  pznondarcyanalysis.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZDarcyAnalysis.h"
#include "pzcheckgeom.h"
#include "pzlog.h"

#include "TPZReadGIDGrid.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "tpzhierarquicalgrid.h"

#include "TPZVTKGeoMesh.h"
#include "TPZDarcyFlow3D.h"
#include "pzpoisson3d.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzfstrmatrix.h"
#include "math.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"

#ifdef PZDEBUG
    #ifdef LOG4CXX
    static LoggerPtr logger(Logger::getLogger("pz.DarcyFlow3D"));
    #endif
#endif

TPZDarcyAnalysis::TPZDarcyAnalysis(TPZAutoPointer<SimulationData> DataSimulation, TPZVec<TPZAutoPointer<ReservoirData> > Layers, TPZVec<TPZAutoPointer<PetroPhysicData> > PetroPhysic, TPZAutoPointer<ReducedPVT> FluidModel)
{
    
    fmeshvec.Resize(2);
    
    fgmesh=NULL;
    fcmeshdarcy=NULL;
    
    // Vector which will store tha residuum in the last state (n)
    fResidualAtn.Resize(0, 0);
    
    // Vector which will store tha residuum in the last state (n+1)
    fResidualAtnplusOne.Resize(0, 0);
    
    fSimulationData     = DataSimulation;
    fLayers             = Layers;
    fRockPetroPhysic    = PetroPhysic;
    fFluidData          = FluidModel;
    
    /** @brief unknowns for n time step */
    falphaAtn.Resize(0, 0);
    
    /** @brief unknowns for n+1 time step */
    falphaAtnplusOne.Resize(0, 0);
    
}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
    
}

void TPZDarcyAnalysis::SetLastState()
{
    fSimulationData->SetnStep(true);
}

void TPZDarcyAnalysis::SetNextState()
{
    fSimulationData->SetnStep(false);
}

void TPZDarcyAnalysis::Assemble()
{
    
}

void TPZDarcyAnalysis::AssembleLastStep(TPZAnalysis *an)
{
    fcmeshdarcy->LoadSolution(falphaAtn);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    SetLastState();
    an->AssembleResidual();
    fResidualAtn = an->Rhs();
    
    // #ifdef PZDEBUG
    //     #ifdef LOG4CXX
    //         if(logger->isDebugEnabled())
    //         {
    //             std::stringstream sout;
    //             fResidualAtn.Print("fResidualAtn = ", sout,EMathematicaInput);
    //             LOGPZ_DEBUG(logger,sout.str())
    //         }
    //     #endif
    // #endif
    
}

void TPZDarcyAnalysis::AssembleNextStep(TPZAnalysis *an)
{
    fcmeshdarcy->LoadSolution(falphaAtnplusOne);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
    SetNextState();
    an->Assemble();
    fResidualAtnplusOne = an->Rhs();
    
    // #ifdef PZDEBUG
    //     #ifdef LOG4CXX
    //         if(logger->isDebugEnabled())
    //         {
    //             std::stringstream sout;
    //             fResidualAtnplusOne.Print("fResidualAtnplusOne = ", sout,EMathematicaInput);
    //             LOGPZ_DEBUG(logger,sout.str())
    //         }
    //     #endif
    // #endif
    
}

void TPZDarcyAnalysis::UpDateAlphaVec(TPZFMatrix<REAL> &alpha)
{
    falphaAtn = alpha;
    falphaAtnplusOne = alpha;
    
    // #ifdef PZDEBUG
    //     #ifdef LOG4CXX
    //         if(logger->isDebugEnabled())
    //         {
    //             std::stringstream sout;
    //             falphaAtn.Print("falphaAtn = ", sout,EMathematicaInput);
    //             falphaAtnplusOne.Print("falphaAtnplusOne = ", sout,EMathematicaInput);
    //             LOGPZ_DEBUG(logger,sout.str())
    //         }
    //     #endif
    // #endif
    
}

void TPZDarcyAnalysis::Residual(TPZFMatrix<STATE> &residual, int icase)
{
    //    TPZNonLinearAnalysis::Residual(residual, icase);
    //    residual = fResidualLastState + residual;
}

void TPZDarcyAnalysis::ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase)
{
    this->SetNextState();
    TPZDarcyAnalysis::ComputeTangent(tangent, coefs, icase);
}

void TPZDarcyAnalysis::InitializeSolution(TPZAnalysis *an)
{
    
    // Compute the intial saturation distribution
    int nalpha = fcmeshdarcy->Solution().Rows();
    falphaAtn.Resize(nalpha, 1);
    falphaAtnplusOne.Resize(nalpha, 1);
    falphaAtn.Zero();
    falphaAtnplusOne.Zero();
//    an->LoadSolution(falphaAtn);
//
//    int nsoil = fmeshvec[3]->Solution().Rows();
//    TPZFMatrix<REAL> SOil(nsoil,1,1.0);
//    fmeshvec[3]->LoadSolution(SOil);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshdarcy);
    
    falphaAtn = fcmeshdarcy->Solution();
    falphaAtnplusOne = fcmeshdarcy->Solution();
    
    // #ifdef PZDEBUG
    //     #ifdef LOG4CXX
    //         if(logger->isDebugEnabled())
    //         {
    //             std::stringstream sout;
    //             falphaAtn.Print("falphaAtn = ", sout,EMathematicaInput);
    //             falphaAtnplusOne.Print("falphaAtnplusOne = ", sout,EMathematicaInput);
    //             LOGPZ_DEBUG(logger,sout.str())
    //         }
    //     #endif
    // #endif
    
}

void TPZDarcyAnalysis::Run()
{
    
    std::string dirname = PZSOURCEDIR;
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//    gRefDBase.InitializeUniformRefPattern(ETriangle);
//    gRefDBase.InitializeUniformRefPattern(ECube);
    
#ifdef PZDEBUG
    #ifdef LOG4CXX
        
        std::string FileName = dirname;
        FileName = dirname + "/Projects/DarcyflowHdiv3D/";
        FileName += "DarcyFlowLog3D.cfg";
        InitializePZLOG(FileName);
        
    #endif
#endif
    
    //  Reading mesh
    std::string GridFileName;
    GridFileName = dirname + "/Projects/DarcyflowHdiv3D/";
    GridFileName += "SingleLayer.dump";
    //GridFileName += "MixLayer.dump";
    //GridFileName += "BatatacoarseQ.dump";
    //GridFileName += "QUAD4.dump";
    
    if(fLayers[0]->GetIsGIDGeometry())
    {
        ReadGeoMesh(GridFileName);
    }
    else
    {
        int nx = 2;
        int ny = 2;
        int nz = 2;
        Geometry(nx,ny,nz);
        //        CreatedGeoMesh();
    }
    
    
    REAL deg = 0.0;
    int hcont = 0;
    RotateGeomesh(deg * M_PI/180.0);
    this->UniformRefinement(fSimulationData->GetHrefinement());
    
    std::set<int> matidstoRef;
    //    matidstoRef.insert(2);
    //    matidstoRef.insert(3);
    //    matidstoRef.insert(4);
    matidstoRef.insert(3);
    matidstoRef.insert(5);
    this->UniformRefinement(hcont, matidstoRef);
    this->PrintGeoMesh();
    
    
    
    int q = fSimulationData->Getqorder();
    int p = fSimulationData->Getporder();
    int s = 0;
    
    //    if (fSimulationData->GetIsH1approx())
    if (false)
    {
        //        CmeshH1(p);
    }
    else
    {
        CreateMultiphysicsMesh(q,p,s);
//        CreateInterfaces();
    }
    
    
    // Analysis
    bool mustOptimizeBandwidth = true;
    TPZAnalysis *an = new TPZAnalysis(fcmeshdarcy,mustOptimizeBandwidth);
    int numofThreads = 0;
    
    bool IsDirecSolver = fSimulationData->GetIsDirect();
    
    if (IsDirecSolver) {
        
        if (fSimulationData->GetIsBroyden()) {
            TPZFStructMatrix fullMatrix(fcmeshdarcy);
            TPZStepSolver<STATE> step;
            fullMatrix.SetNumThreads(numofThreads);
            step.SetDirect(ELU);
            an->SetStructuralMatrix(fullMatrix);
            an->SetSolver(step);
        }
        else{
            
            TPZSkylineNSymStructMatrix skylnsym(fcmeshdarcy);
            TPZStepSolver<STATE> step;
            skylnsym.SetNumThreads(numofThreads);
            step.SetDirect(ELU);
            an->SetStructuralMatrix(skylnsym);
            an->SetSolver(step);
        }
        
    }
    else
    {
        if (fSimulationData->GetIsBroyden()) {
            TPZFStructMatrix fullMatrix(fcmeshdarcy);
            fullMatrix.SetNumThreads(numofThreads);
            
            TPZAutoPointer<TPZMatrix<STATE> > fullMatrixa = fullMatrix.Create();
            TPZAutoPointer<TPZMatrix<STATE> > fullMatrixaClone = fullMatrixa->Clone();
            
            TPZStepSolver<STATE> *stepre = new TPZStepSolver<STATE>(fullMatrixaClone);
            TPZStepSolver<STATE> *stepGMRES = new TPZStepSolver<STATE>(fullMatrixa);
            TPZStepSolver<STATE> *stepGC = new TPZStepSolver<STATE>(fullMatrixa);
            stepre->SetDirect(ELU);
            stepre->SetReferenceMatrix(fullMatrixa);
            stepGMRES->SetGMRES(10, 20, *stepre, 1.0e-10, 0);
            stepGC->SetCG(10, *stepre, 1.0e-10, 0);
            if (fSimulationData->GetIsCG()) {
                an->SetSolver(*stepGC);
            }
            else{
                an->SetSolver(*stepGMRES);
            }
            
        }
        else{
            
            TPZSkylineNSymStructMatrix skylnsym(fcmeshdarcy);
            skylnsym.SetNumThreads(numofThreads);
            
            TPZAutoPointer<TPZMatrix<STATE> > skylnsyma = skylnsym.Create();
            TPZAutoPointer<TPZMatrix<STATE> > skylnsymaClone = skylnsyma->Clone();
            
            TPZStepSolver<STATE> *stepre = new TPZStepSolver<STATE>(skylnsymaClone);
            TPZStepSolver<STATE> *stepGMRES = new TPZStepSolver<STATE>(skylnsyma);
            TPZStepSolver<STATE> *stepGC = new TPZStepSolver<STATE>(skylnsyma);
            
            stepre->SetDirect(ELU);
            stepre->SetReferenceMatrix(skylnsyma);
            stepGMRES->SetGMRES(10, 20, *stepre, 1.0e-10, 0);
            stepGC->SetCG(10, *stepre, 1.0e-10, 0);
            if (fSimulationData->GetIsCG()) {
                an->SetSolver(*stepGC);
            }
            else{
                an->SetSolver(*stepGMRES);
            }
        }
        
    }
    
    this->InitializeSolution(an);
    this->TimeForward(an);
    
}


void TPZDarcyAnalysis::CreateInterfaces()
{
    fgmesh->ResetReference();
    fcmeshdarcy->LoadReferences();
    
    // Creation of interface elements
    int nel = fcmeshdarcy->ElementVec().NElements();
    for(int el = 0; el < nel; el++)
    {
        TPZCompEl * compEl = fcmeshdarcy->ElementVec()[el];
        if(!compEl) continue;
        TPZGeoEl * gel = compEl->Reference();
        if(!gel) {continue;}
        if(gel->HasSubElement()) {continue;}
        int index = compEl ->Index();
        if(compEl->Dimension() == fcmeshdarcy->Dimension())
        {
            TPZMultiphysicsElement * InterpEl = dynamic_cast<TPZMultiphysicsElement *>(fcmeshdarcy->ElementVec()[index]);
            if(!InterpEl) continue;
            InterpEl->CreateInterfaces();
        }
    }
}


void TPZDarcyAnalysis::PrintLS(TPZAnalysis *an)
{
    TPZAutoPointer< TPZMatrix<REAL> > KGlobal;
    TPZFMatrix<STATE> FGlobal;
    KGlobal =   an->Solver().Matrix();
    FGlobal =   an->Rhs();
    
#ifdef PZDEBUG
    #ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            KGlobal->Print("KGlobal = ", sout,EMathematicaInput);
            FGlobal.Print("FGlobal = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
    #endif
#endif
    
}

void TPZDarcyAnalysis::CreateMultiphysicsMesh(int q, int p, int s)
{
    fmeshvec[0] = CmeshFlux(q);
    fmeshvec[1] = CmeshPressure(p);
//    fmeshvec[2] = CmeshSw(s);
//    fmeshvec[3] = CmeshSo(s);
    
    fcmeshdarcy = CmeshMixed();
    

    
#ifdef PZDEBUG
    std::ofstream dumpfile("ComputationaMeshMultiphysics.txt");
    fcmeshdarcy->Print(dumpfile);
#endif
    
}

void TPZDarcyAnalysis::TimeForward(TPZAnalysis *an)
{
    int nsteps = fSimulationData->GetMaxTime() / fSimulationData->GetDeltaT();
    
    REAL tk = 0;
    
    this->fSimulationData->SetTime(tk);
    this->PostProcessVTK(an);
    
    for (int istep = 1 ; istep <=  nsteps; istep++) {
        tk = istep*this->fSimulationData->GetDeltaT();
        this->fSimulationData->SetTime(tk);
        
        std::cout << "Begin of time (days): " << tk/86400.0 << std::endl;
        std::cout<<  "Time step: " << istep << std::endl;
        
        
        this->AssembleLastStep(an);
        this->AssembleNextStep(an);
        
        
        if (fSimulationData->GetIsBroyden())
        {
            
            const clock_t tinia = clock();
            BroydenIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Broyden: " << timea << std::endl;
            std::cout << "Number of DOF = " << fcmeshdarcy->Solution().Rows() << std::endl;
            this->PostProcessVTK(an);
            
        }
        else
        {
            const clock_t tinia = clock();
            NewtonIterations(an);
            const clock_t tenda = clock();
            const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
            std::cout << "Time for Newton: " << timea << std::endl;
            std::cout << "Number of DOF = " << fcmeshdarcy->Solution().Rows() << std::endl;
            this->PostProcessVTK(an);
        }
        
        
        
    }
    
    
    
}


void TPZDarcyAnalysis::NewtonIterations(TPZAnalysis *an)
{
    
    TPZFMatrix<STATE> Residual(an->Rhs().Rows(),1,0.0);
    
    
    //    TPZManVector<long> Actives(0),NoActives(0);
    //
    //    this->FilterSaturations(Actives, NoActives);
    //    an->StructMatrix()->EquationFilter().Reset();
    //    an->StructMatrix()->EquationFilter().SetActiveEquations(Actives);
    
    Residual = fResidualAtn + fResidualAtnplusOne;
    //    an->Rhs() = Residual;
//    this->PrintLS(an);
    
    TPZFMatrix<STATE> X = falphaAtn;
    TPZFMatrix<STATE> DeltaX = falphaAtn;
    
    STATE error     =   1;
    STATE normdx    =   1;
    int iterations  =   0;
    int centinel    =   0;
    int fixed       =   fSimulationData->GetFixediterations();
    
    while (error >= fSimulationData->GetToleranceRes() && iterations <= fSimulationData->GetMaxiterations()) {
        
        an->Rhs() = Residual;
        an->Rhs() *= -1.0;
        
        //        const clock_t tini = clock();
        
        an->Solve();
        //        const clock_t tend = clock();
        //        const REAL time = REAL(REAL(tend - tini)/CLOCKS_PER_SEC);
        //        std::cout << "Time for solving: " << time << std::endl;
        
        DeltaX = an->Solution();
        normdx = Norm(DeltaX);
        X += DeltaX;
        
        fcmeshdarcy->LoadSolution(X);
        if (true)
        {
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
        }
        
        
        //        const clock_t tinia = clock();
        
        if (((fixed+1) * (centinel) == iterations)) {
            an->Assemble();
            centinel++;
        }
        else{
            an->AssembleResidual();
            
        }
        fResidualAtnplusOne = an->Rhs();
        
        //        const clock_t tenda = clock();
        //        const REAL timea = REAL(REAL(tenda - tinia)/CLOCKS_PER_SEC);
        //        std::cout << "Time for assemble: " << timea << std::endl;
        
        Residual = fResidualAtn + fResidualAtnplusOne;
        error = Norm(Residual);
        iterations++;
        
        
#ifdef PZDEBUG
    #ifdef LOG4CXX
            if(logger->isDebugEnabled())
            {
                std::stringstream sout;
                fResidualAtn.Print("fResidualAtn = ", sout,EMathematicaInput);
                fResidualAtnplusOne.Print("fResidualAtnplusOne = ", sout,EMathematicaInput);
                DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
                X.Print("X = ", sout,EMathematicaInput);
                LOGPZ_DEBUG(logger,sout.str())
            }
    #endif
#endif
        
        if(error < fSimulationData->GetToleranceRes() || normdx < fSimulationData->GetToleranceDX())
        {
            std::cout << "Converged; iterations:  " << iterations << std::endl;
            std::cout << "error norm: " << error << std::endl;
            std::cout << "error of dx: " << normdx << std::endl;
            this->UpDateAlphaVec(X);
            break;
        }
        
        if (iterations == fSimulationData->GetMaxiterations()) {
            std::cout << "Out max iterations " << iterations << std::endl;
            std::cout << "error norm " << error << std::endl;
            this->UpDateAlphaVec(X);
            break;
        }
        
    }
    
    
}

void TPZDarcyAnalysis::BroydenIterations(TPZAnalysis *an)
{
    
    
    int m = an->Solution().Rows();
    
    TPZFMatrix<STATE> Residual(m,1,0.0);
    TPZFMatrix<STATE> Rank(m,m,0.0);
    TPZFMatrix<STATE> X(m,1,0.0);
    TPZFMatrix<STATE> DeltaX(m,1,0.0);
    
    TPZAutoPointer<TPZMatrix<STATE> >  D;
    TPZAutoPointer<TPZMatrix<STATE> >  Dk;
    
    TPZFMatrix<STATE> * DInverse =  new TPZFMatrix<STATE> (m,m,0.0);
    
    bool IsShermanMorrison = false;
    STATE ck = 0.0;
    //    TPZFMatrix<STATE> DInverse(m,m,0.0);
    TPZFMatrix<STATE> Identity(m,m,0.0);
    TPZFMatrix<STATE> DInverseT(m,m,0.0);
    TPZFMatrix<STATE> u(m,m,0.0);
    TPZFMatrix<STATE> du(m,m,0.0);
    TPZFMatrix<STATE> duT(m,m,0.0);
    TPZFMatrix<STATE> Ck(1,1,0.0);
    
    STATE error=1;
    STATE normdx=1;
    int iterations=0;
    int centinel    =0;
    int fixed       =fSimulationData->GetFixediterations();
    
    // Computing the first Newton iteration
    
    
    an->Rhs() = fResidualAtn + fResidualAtnplusOne;       // g(X0)
    
    
    an->Rhs() *= -1.0;
    
    if (IsShermanMorrison) {
        DInverse = ComputeInverse();
        DInverse->Multiply(Residual, DeltaX);
        D.operator->()->Multiply(*DInverse, Identity);
        
#ifdef PZDEBUG
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
            Residual.Print("Residual = ", sout,EMathematicaInput);
            Identity.Print("Identity = ", sout,EMathematicaInput);
            DInverse->Print("DInverse = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
#endif
        
    }
    else{
        D = an->Solver().Matrix();  // J(X0)
        an->Solve();
        DeltaX = an->Solution();    // d(X0)
    }
    normdx = Norm(DeltaX);      // d(X0)*d(X0)
    X += DeltaX;                // X1
    
    // End of newton iteration
    
    // Procedure without Inverse computation
    
    iterations++;
    
    while (error >= fSimulationData->GetToleranceRes() && iterations <= fSimulationData->GetMaxiterations()) {
        
        
        fcmeshdarcy->LoadSolution(X);
        if (!fSimulationData->GetIsH1approx())
        {
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshdarcy);
        }
        
        this->AssembleResidual();
        Residual = fResidualAtn + fResidualAtnplusOne;       // g(Xk)
        error = Norm(Residual);     // g(Xk)*g(Xk)
        
        if(error <= fSimulationData->GetToleranceRes() || normdx <= fSimulationData->GetToleranceDX())
        {
            std::cout << "Converged with iterations:  " << iterations << std::endl;
            std::cout << "error norm: " << error << std::endl;
            std::cout << "error of dx: " << normdx << std::endl;
            break;
        }
#ifdef PZDEBUG
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
            Residual.Print("Residual = ", sout,EMathematicaInput);
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
#endif
        
        if (((fixed+1) * (centinel) == iterations)) {
            
            if (IsShermanMorrison) {
                DInverse->Multiply(Residual, u);
                u.Add(DeltaX, du);
                u.Transpose(&duT);
                DeltaX.Multiply(duT,Ck);
                ck = Ck(0,0);
                Rank = TensorProduct(u, DeltaX);
                Rank *= -(1/ck);
                Rank.Multiply(*DInverse, DInverseT);
                DInverse->Add(DInverseT, *DInverse);
                
            }else
            {
                // Application of the Secant condition
                Rank = TensorProduct(Residual, DeltaX);
                Rank *= 1.0/(normdx*normdx);
                D.operator->()->Add(Rank, *D.operator->());
                an->Solver().Matrix() = D;
            }
            
            
            
            centinel++;
            
        }
        
        
        
        //#ifdef LOG4CXX
        //        if(logger->isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            an->Solver().Matrix().operator->()->Print("*an->Solver().Matrix().operator->() = ", sout,EMathematicaInput);
        //            D->Print("D = ", sout,EMathematicaInput);
        //            Rank.Print("Rank = ", sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger,sout.str())
        //        }
        //#endif
        an->Rhs() = fResidualAtn + fResidualAtnplusOne;
        an->Rhs() *= -1.0;
        
        if (IsShermanMorrison) {
            DInverse->Multiply(Residual, DeltaX);
        }
        else{
            an->Solve();
            DeltaX = an->Solution();    // d(Xk)
        }
        
        
        normdx = Norm(DeltaX);      // d(Xk)*d(Xk)
        X += DeltaX;                // Xk+1
        
        
        //#ifdef LOG4CXX
        //        if(logger->isDebugEnabled())
        //        {
        //            std::stringstream sout;
        //            DeltaX.Print("DeltaX = ", sout,EMathematicaInput);
        //            X.Print("X = ", sout,EMathematicaInput);
        //            LOGPZ_DEBUG(logger,sout.str())
        //        }
        //#endif
        iterations++;
        
        if (iterations == fSimulationData->GetMaxiterations()) {
            std::cout << "Out max iterations " << iterations << std::endl;
            std::cout << "error norm " << error << std::endl;
            break;
        }
        
    }
    
    
}


TPZFMatrix<STATE>  TPZDarcyAnalysis::TensorProduct(TPZFMatrix<STATE> &g, TPZFMatrix<STATE> &d)
{
    TPZFMatrix<STATE> dT=d;
    d.Transpose(&dT);
    TPZFMatrix<STATE> RankOne;
    g.Multiply(dT, RankOne);
    
    //#ifdef LOG4CXX
    //    if(logger->isDebugEnabled())
    //    {
    //        std::stringstream sout;
    //        g.Print("g = ", sout,EMathematicaInput);
    //        dT.Print("dT = ", sout,EMathematicaInput);
    //        d.Print("d = ", sout,EMathematicaInput);
    //        RankOne.Print("RankOne = ", sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logger,sout.str())
    //    }
    //#endif
    
    return RankOne;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshMixed()
{
    int dim = 3;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int NId = fLayers[ilayer]->GetMatIDs()[3];
    int SId = fLayers[ilayer]->GetMatIDs()[1];
    int EId = fLayers[ilayer]->GetMatIDs()[2];
    int WId = fLayers[ilayer]->GetMatIDs()[4];
    int BId = fLayers[ilayer]->GetMatIDs()[5];
    int TId = fLayers[ilayer]->GetMatIDs()[6];
    
    int typeFluxin = 1, typePressurein = 0;
//    int typeFluxout = 3, typePressureout = 2;
    TPZFMatrix<STATE> val1(1,2,0.), val2(1,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    // Material medio poroso
    TPZDarcyFlow3D * mat = new TPZDarcyFlow3D(RockId);
    mat->SetSimulationData(fSimulationData);
    mat->SetReservoirData(fLayers[ilayer]);
    mat->SetPetroPhysicsData(fRockPetroPhysic[ilayer]);
    mat->SetFluidModelData(fFluidData);
    cmesh->InsertMaterialObject(mat);
    
    
    // Rigth hand side function
    TPZDummyFunction<STATE> *dum = new TPZDummyFunction<STATE>(Ffunction);
    TPZAutoPointer<TPZFunction<STATE> > forcef;
    dum->SetPolynomialOrder(0);
    forcef = dum;
    mat->SetForcingFunction(forcef);
    
    // Setting up linear tracer solution
    TPZDummyFunction<STATE> *Analytic = new TPZDummyFunction<STATE>(AnalyticSolution);
    TPZAutoPointer<TPZFunction<STATE> > fAnalytic = Analytic;
    mat->SetTimeDependentFunctionExact(fAnalytic);
    
    
    // Bc N
    val2(0,0) = -1.0;
    TPZBndCond * bcN = mat->CreateBC(mat, NId, typeFluxin, val1, val2);
    
    // Bc S
    val2(0,0) = 1.0;
    TPZBndCond * bcS = mat->CreateBC(mat, SId, typePressurein, val1, val2);
    
    // Bc E
    val2(0,0) = 0.0;
    TPZBndCond * bcE = mat->CreateBC(mat, EId, typeFluxin, val1, val2);
    
    // Bc W
    val2(0,0) = 0.0;
    TPZBndCond * bcW = mat->CreateBC(mat, WId, typeFluxin, val1, val2);
    
    // Bc B
    val2(0,0) = 0.0;
    TPZBndCond * bcB = mat->CreateBC(mat, BId, typeFluxin, val1, val2);
    
    // Bc T
    val2(0,0) = 0.0;
    TPZBndCond * bcT = mat->CreateBC(mat, TId, typeFluxin, val1, val2);
    
    
    cmesh->InsertMaterialObject(bcN);
    cmesh->InsertMaterialObject(bcS);
    cmesh->InsertMaterialObject(bcE);
    cmesh->InsertMaterialObject(bcW);
    cmesh->InsertMaterialObject(bcB);
    cmesh->InsertMaterialObject(bcT);
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    bool StaticCondensation = false;
    if (StaticCondensation)
    {
//        //Creating multiphysic elements containing skeletal elements.
//        TPZBuildMultiphysicsMesh::AddElements(fmeshvec,cmesh);
//        cmesh->Reference()->ResetReference();
//        cmesh->LoadReferences();
//        
//        long nel = cmesh->ElementVec().NElements();
//        
//        std::map<long, long> bctoel, eltowrap;
//        for (long el=0; el<nel; el++) {
//            TPZCompEl *cel = cmesh->Element(el);
//            TPZGeoEl *gel = cel->Reference();
//            int matid = gel->MaterialId();
//            if (matid < 0) {
//                TPZGeoElSide gelside(gel,gel->NSides()-1);
//                TPZGeoElSide neighbour = gelside.Neighbour();
//                while (neighbour != gelside) {
//                    if (neighbour.Element()->Dimension() == dim && neighbour.Element()->Reference()) {
//                        // got you!!
//                        bctoel[el] = neighbour.Element()->Reference()->Index();
//                        break;
//                    }
//                    neighbour = neighbour.Neighbour();
//                }
//                if (neighbour == gelside) {
//                    DebugStop();
//                }
//            }
//        }
//        
//        TPZStack< TPZStack< TPZMultiphysicsElement *,7> > wrapEl;
//        for(long el = 0; el < nel; el++)
//        {
//            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(cmesh->Element(el));
//            if(mfcel->Dimension()==dim) TPZBuildMultiphysicsMesh::AddWrap(mfcel, RockId, wrapEl);//criei elementos com o mesmo matId interno, portanto nao preciso criar elemento de contorno ou outro material do tipo TPZLagrangeMultiplier
//        }
//        
//        for (long el =0; el < wrapEl.size(); el++) {
//            TPZCompEl *cel = wrapEl[el][0];
//            long index = cel->Index();
//            eltowrap[index] = el;
//        }
//        
//        fmeshvec[0]->CleanUpUnconnectedNodes();
//        TPZBuildMultiphysicsMesh::AddConnects(fmeshvec,cmesh);
//        TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, cmesh);
//        
//        std::map<long, long>::iterator it;
//        for (it = bctoel.begin(); it != bctoel.end(); it++) {
//            long bcindex = it->first;
//            long elindex = it->second;
//            if (eltowrap.find(elindex) == eltowrap.end()) {
//                DebugStop();
//            }
//            long wrapindex = eltowrap[elindex];
//            TPZCompEl *bcel = fcmeshdarcy->Element(bcindex);
//            TPZMultiphysicsElement *bcmf = dynamic_cast<TPZMultiphysicsElement *>(bcel);
//            if (!bcmf) {
//                DebugStop();
//            }
//            wrapEl[wrapindex].Push(bcmf);
//            
//        }
//        
//        //------- Create and add group elements -------
//        long index, nenvel;
//        nenvel = wrapEl.NElements();
//        TPZStack<TPZElementGroup *> elgroups;
//        for(long ienv=0; ienv<nenvel; ienv++){
//            TPZElementGroup *elgr = new TPZElementGroup(*wrapEl[ienv][0]->Mesh(),index);
//            elgroups.Push(elgr);
//            nel = wrapEl[ienv].NElements();
//            for(int jel=0; jel<nel; jel++){
//                elgr->AddElement(wrapEl[ienv][jel]);
//            }
//        }
//        
//        cmesh->ComputeNodElCon();
//        // create condensed elements
//        // increase the NumElConnected of one pressure connects in order to prevent condensation
//        for (long ienv=0; ienv<nenvel; ienv++) {
//            TPZElementGroup *elgr = elgroups[ienv];
//            int nc = elgr->NConnects();
//            for (int ic=0; ic<nc; ic++) {
//                TPZConnect &c = elgr->Connect(ic);
//                if (c.LagrangeMultiplier() > 0) {
//                    c.IncrementElConnected();
//                    break;
//                }
//            }
//            
////            TPZCondensedCompEl *condense = new TPZCondensedCompEl(elgr);
//        }
//        
//        cmesh->CleanUpUnconnectedNodes();
//        cmesh->ExpandSolution();
    }
    else{
     
        // Transferindo para a multifisica
        TPZBuildMultiphysicsMesh::AddElements(fmeshvec, cmesh);
        TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, cmesh);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, cmesh);
    }
    
    return cmesh;
}


TPZCompMesh * TPZDarcyAnalysis::CmeshFlux(int qorder)
{
    int dim = 3;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int NId = fLayers[ilayer]->GetMatIDs()[3];
    int SId = fLayers[ilayer]->GetMatIDs()[1];
    int EId = fLayers[ilayer]->GetMatIDs()[2];
    int WId = fLayers[ilayer]->GetMatIDs()[4];
    int BId = fLayers[ilayer]->GetMatIDs()[5];
    int TId = fLayers[ilayer]->GetMatIDs()[6];
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZDarcyFlow3D * mat = new TPZDarcyFlow3D(RockId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, NId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, SId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcS);
    
    // Bc E
    TPZBndCond * bcE = mat->CreateBC(mat, EId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcE);
    
    // Bc W
    TPZBndCond * bcW = mat->CreateBC(mat, WId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcW);
    
    // Bc B
    TPZBndCond * bcB = mat->CreateBC(mat, BId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcB);
    
    // Bc T
    TPZBndCond * bcT = mat->CreateBC(mat, TId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcT);
    
    // Setando Hdiv
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(qorder);
    cmesh->SetAllCreateFunctionsHDiv();
    
    
    cmesh->AutoBuild();
    
    
#ifdef PZDEBUG
    std::ofstream out("cmeshFlux.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * TPZDarcyAnalysis::CmeshPressure(int porder)
{
    
    int dim = 3;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int NId = fLayers[ilayer]->GetMatIDs()[3];
    int SId = fLayers[ilayer]->GetMatIDs()[1];
    int EId = fLayers[ilayer]->GetMatIDs()[2];
    int WId = fLayers[ilayer]->GetMatIDs()[4];
    int BId = fLayers[ilayer]->GetMatIDs()[5];
    int TId = fLayers[ilayer]->GetMatIDs()[6];
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
    cmesh->InsertMaterialObject(mat);
    
    // Bc N
    TPZBndCond * bcN = mat->CreateBC(mat, NId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcN);
    
    // Bc S
    TPZBndCond * bcS = mat->CreateBC(mat, SId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcS);
    
    // Bc E
    TPZBndCond * bcE = mat->CreateBC(mat, EId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcE);
    
    // Bc W
    TPZBndCond * bcW = mat->CreateBC(mat, WId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcW);
    
    // Bc B
    TPZBndCond * bcB = mat->CreateBC(mat, BId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcB);
    
    // Bc T
    TPZBndCond * bcT = mat->CreateBC(mat, TId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcT);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    
#ifdef PZDEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshSw(int Sworder)
{
    
    int dim = 2;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(Sworder);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
#ifdef PZDEBUG
    std::ofstream out("cmeshSw.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}

TPZCompMesh * TPZDarcyAnalysis::CmeshSo(int Soorder)
{
    
    int dim = 2;
    int ilayer = 0;
    int RockId = fLayers[ilayer]->GetMatIDs()[0];
    int bottomId = fLayers[ilayer]->GetMatIDs()[1];
    int rigthId = fLayers[ilayer]->GetMatIDs()[2];
    int topId = fLayers[ilayer]->GetMatIDs()[3];
    int leftId = fLayers[ilayer]->GetMatIDs()[4];
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
    
    TPZMatPoisson3d * mat = new TPZMatPoisson3d(RockId,dim);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Right
    TPZBndCond * bcRight = mat->CreateBC(mat, rigthId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcRight);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Bc Left
    TPZBndCond * bcLeft = mat->CreateBC(mat, leftId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcLeft);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(Soorder);
    
    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    
#ifdef PZDEBUG
    std::ofstream out("cmeshSo.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}


void TPZDarcyAnalysis::ReadGeoMesh(std::string GridFileName)
{
    TPZReadGIDGrid GeometryInfo;
    GeometryInfo.SetfDimensionlessL(1.0);
    fgmesh = GeometryInfo.GeometricGIDMesh(GridFileName);
    fgmesh->SetDimension(3);
}

void TPZDarcyAnalysis::ParametricfunctionX(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = par[0];
    X[1] = 0.0;
    X[2] = 0.0;
}

void TPZDarcyAnalysis::ParametricfunctionY(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = par[0];
    X[2] = 0.0;
}

void TPZDarcyAnalysis::ParametricfunctionZ(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}

void TPZDarcyAnalysis::Geometry(int nx, int ny, int nz)
{
    REAL t=0.0;
    REAL dt;
    int n;
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh0D = new TPZGeoMesh;
    GeoMesh0D->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh0D->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom0D(GeoMesh0D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncX = new TPZDummyFunction<STATE>(ParametricfunctionX);
    CreateGridFrom0D.SetParametricFunction(ParFuncX);
    CreateGridFrom0D.SetFrontBackMatId(5,3);
    
    dt = 1.0;
    n = nx;
    // Computing Mesh extruded along the parametric curve Parametricfunction
    TPZGeoMesh * GeoMesh1D = CreateGridFrom0D.ComputeExtrusion(t, dt, n);
    
    TPZHierarquicalGrid CreateGridFrom1D(GeoMesh1D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncY = new TPZDummyFunction<STATE>(ParametricfunctionY);
    CreateGridFrom1D.SetParametricFunction(ParFuncY);
    CreateGridFrom1D.SetFrontBackMatId(2,4);
    CreateGridFrom1D.SetTriangleExtrusion();
    
    dt = 1.0;
    n = ny;
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    TPZGeoMesh * GeoMesh2D = CreateGridFrom1D.ComputeExtrusion(t, dt, n);
    
    
    TPZHierarquicalGrid CreateGridFrom2D(GeoMesh2D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncZ = new TPZDummyFunction<STATE>(ParametricfunctionZ);
    CreateGridFrom2D.SetParametricFunction(ParFuncZ);
    CreateGridFrom2D.SetFrontBackMatId(6,7);
    CreateGridFrom2D.SetTriangleExtrusion();
    
    dt = 1.0;
    n = nz;
    // Computing Mesh extruded along the parametric curve Parametricfunction2
    fgmesh = CreateGridFrom2D.ComputeExtrusion(t, dt, n);
    
}

void TPZDarcyAnalysis::PrintGeoMesh()
{
    
    
#ifdef PZDEBUG
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    fgmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fgmesh,Dummyfile, true);
    
#endif
}

void TPZDarcyAnalysis::RotateGeomesh(REAL CounterClockwiseAngle)
{
    REAL theta = CounterClockwiseAngle;
    // It represents a 3D rotation around the z axis.
    TPZFMatrix<STATE> RotationMatrix(3,3,0.0);
    RotationMatrix(0,0) =   +cos(theta);
    RotationMatrix(0,1) =   -sin(theta);
    RotationMatrix(1,0) =   +sin(theta);
    RotationMatrix(1,1) =   +cos(theta);
    RotationMatrix(2,2) = 1.0;
    TPZVec<STATE> iCoords(3,0.0);
    TPZVec<STATE> iCoordsRotated(3,0.0);
    
    RotationMatrix.Print("Rotation = ");
    
    int NumberofGeoNodes = fgmesh->NNodes();
    for (int inode = 0; inode < NumberofGeoNodes; inode++)
    {
        TPZGeoNode GeoNode = fgmesh->NodeVec()[inode];
        GeoNode.GetCoordinates(iCoords);
        // Apply rotation
        iCoordsRotated[0] = RotationMatrix(0,0)*iCoords[0]+RotationMatrix(0,1)*iCoords[1]+RotationMatrix(0,2)*iCoords[2];
        iCoordsRotated[1] = RotationMatrix(1,0)*iCoords[0]+RotationMatrix(1,1)*iCoords[1]+RotationMatrix(1,2)*iCoords[2];
        iCoordsRotated[2] = RotationMatrix(2,0)*iCoords[0]+RotationMatrix(2,1)*iCoords[1]+RotationMatrix(2,2)*iCoords[2];
        GeoNode.SetCoord(iCoordsRotated);
        fgmesh->NodeVec()[inode] = GeoNode;
    }
}

void TPZDarcyAnalysis::UniformRefinement(int nh)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = fgmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fgmesh->ElementVec() [i];
            if (gel->Dimension() == 2 || gel->Dimension() == 1) gel->Divide (filhos);
        }//for i
    }//ref
    fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::UniformRefinement(int nh, std::set<int> &MatToRef)
{
    for ( int ref = 0; ref < nh; ref++ ){
        TPZVec<TPZGeoEl *> filhos;
        long n = fgmesh->NElements();
        for ( long i = 0; i < n; i++ ){
            TPZGeoEl * gel = fgmesh->ElementVec() [i];
            if(!gel){continue;}
            //            int reflevel = gel->Level();
            //            if (reflevel == ref + 1) {
            //                continue;
            //            }
            TPZRefPatternTools::RefineDirectional(gel,MatToRef);
        }//for i
    }//ref
    fgmesh->BuildConnectivity();
}

void TPZDarcyAnalysis::UniformRefinement(int nh, int MatId)
{
    //    for ( int ref = 0; ref < nh; ref++ ){
    //        TPZVec<TPZGeoEl *> filhos;
    //        long n = fgmesh->NElements();
    //        for ( long i = 0; i < n; i++ ){
    //            TPZGeoEl * gel = fgmesh->ElementVec() [i];
    //            if(!gel){continue;}
    //            if (gel->Dimension() == 1){
    //                if (gel->MaterialId() == MatId) {
    //                    gel->Divide(filhos);
    //                }
    //
    //            }
    //        }//for i
    //    }//ref
    
    ///Refinamento
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    
    for (int idivide = 0; idivide < nh; idivide++){
        const int nels = fgmesh->NElements();
        TPZVec< TPZGeoEl * > allEls(nels);
        for(int iel = 0; iel < nels; iel++){
            allEls[iel] = fgmesh->ElementVec()[iel];
        }
        
        for(int iel = 0; iel < nels; iel++){
            TPZGeoEl * gel = allEls[iel];
            if(!gel) continue;
            if(gel->HasSubElement()) continue;
            int nnodes = gel->NNodes();
            int found = -1;
            for(int in = 0; in < nnodes; in++){
                if(gel->NodePtr(in)->Id() == MatId){
                    found = in;
                    break;
                }
            }///for in
            if(found == -1) continue;
            
            MElementType gelT = gel->Type();
            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
            if(!uniform){
                DebugStop();
            }
            gel->SetRefPattern(uniform);
            TPZVec<TPZGeoEl*> filhos;
            gel->Divide(filhos);
            
        }///for iel
    }//idivide
    
    fgmesh->BuildConnectivity();
    
}

////refinamento uniforme em direcao ao no
//void DirectionalRef(int nh, int MatId){
//
//    ///Refinamento
//    gRefDBase.InitializeUniformRefPattern(EOned);
//    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
//
//
//    for (int idivide = 0; idivide < nh; idivide++){
//        const int nels = fgmesh->NElements();
//        TPZVec< TPZGeoEl * > allEls(nels);
//        for(int iel = 0; iel < nels; iel++){
//            allEls[iel] = gmesh->ElementVec()[iel];
//        }
//
//        for(int iel = 0; iel < nels; iel++){
//            TPZGeoEl * gel = allEls[iel];
//            if(!gel) continue;
//            if(gel->HasSubElement()) continue;
//            int nnodes = gel->NNodes();
//            int found = -1;
//            for(int in = 0; in < nnodes; in++){
//                if(gel->NodePtr(in)->Id() == nodeAtOriginId){
//                    found = in;
//                    break;
//                }
//            }///for in
//            if(found == -1) continue;
//
//            MElementType gelT = gel->Type();
//            TPZAutoPointer<TPZRefPattern> uniform = gRefDBase.GetUniformRefPattern(gelT);
//            if(!uniform){
//                DebugStop();
//            }
//            gel->SetRefPattern(uniform);
//            TPZVec<TPZGeoEl*> filhos;
//            gel->Divide(filhos);
//
//        }///for iel
//    }//idivide
//
//    gmesh->BuildConnectivity();
//
//#ifdef LOG4CXX
//    if (logger->isDebugEnabled())
//    {
//        std::stringstream sout;
//        sout<<"gmesh depois de refinar direcionalmente\n";
//        gmesh->Print(sout);
//        LOGPZ_DEBUG(logger, sout.str());
//    }
//#endif
//
//}///void



void TPZDarcyAnalysis::PostProcessVTK(TPZAnalysis *an)
{
    const int dim = 3;
    int div = fSimulationData->GetHPostrefinement();
    TPZStack<std::string> scalnames, vecnames;
    std::string plotfile;
    if (fSimulationData->GetIsH1approx()) {
        plotfile = "3DH1Darcy.vtk";
    }
    else{
        plotfile = "3DMixedDarcy.vtk";
    }
    
    scalnames.Push("WeightedPressure");
    scalnames.Push("WaterSaturation");
    scalnames.Push("OilSaturation");
    scalnames.Push("WaterDensity");
    scalnames.Push("OilDensity");
    scalnames.Push("Porosity");
    scalnames.Push("DivOfBulkVeclocity");
//    scalnames.Push("ExactSaturation");
    vecnames.Push("BulkVelocity");
    an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
    an->PostProcess(div);
}

void TPZDarcyAnalysis::Ffunction(const TPZVec<REAL> &pt, TPZVec<STATE> &ff)
{
    
    ff[0] = 0.0*0.0000001;
}

TPZFMatrix<STATE> * TPZDarcyAnalysis::ComputeInverse()
{
    int neq = fcmeshdarcy->NEquations();
    TPZFMatrix<STATE> * PreInverse =  new TPZFMatrix<STATE> (neq,neq,0.0);
    TPZFStructMatrix skyl(fcmeshdarcy);
    std::set<int> matids; // to be computed
    matids.insert(1);
    matids.insert(2);
    matids.insert(3);
    matids.insert(4);
    matids.insert(5);
    skyl.SetMaterialIds(matids);
    TPZFMatrix<STATE> rhsfrac;
    TPZFMatrix<STATE> Identity;
    TPZAutoPointer<TPZGuiInterface> gui = new TPZGuiInterface;
    TPZAutoPointer<TPZMatrix<STATE> > MatG = skyl.CreateAssemble(rhsfrac, gui);
    TPZFMatrix<STATE> oldmat = *MatG.operator->();
    oldmat.Inverse( * PreInverse);
    oldmat.Multiply(*PreInverse, Identity);
    
#ifdef PZDEBUG
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        sout << "Is decomposed=  " << MatG->IsDecomposed() << std::endl;
        oldmat.Print("oldmat = ", sout,EMathematicaInput);
        PreInverse->Print("PreInverse = ", sout,EMathematicaInput);
        Identity.Print("Identity = ", sout,EMathematicaInput);
        LOGPZ_DEBUG(logger,sout.str())
    }
#endif
#endif
    
    return PreInverse;
    
}

void TPZDarcyAnalysis::FilterSaturations(TPZManVector<long> &active, TPZManVector<long> &nonactive)
{
    int ncon_flux       = fmeshvec[0]->NConnects();
    int ncon_pressure   = fmeshvec[1]->NConnects();
    int ncon_sw = fmeshvec[2]->NConnects();
    int ncon_so    = fmeshvec[3]->NConnects();
    int ncon = fcmeshdarcy->NConnects();
    
    for(int i = 0; i < ncon-ncon_sw-ncon_so; i++)
    {
        TPZConnect &con = fcmeshdarcy->ConnectVec()[i];
        int seqnum = con.SequenceNumber();
        int pos = fcmeshdarcy->Block().Position(seqnum);
        int blocksize = fcmeshdarcy->Block().Size(seqnum);
        
        int vs = active.size();
        active.Resize(vs+blocksize);
        for(int ieq = 0; ieq<blocksize; ieq++)
        {
            active[vs+ieq] = pos+ieq;
        }
    }
    
    //    for(int i = ncon-ncon_gravity; i < ncon; i++)
    //    {
    //        TPZConnect &con = mphysics->ConnectVec()[i];
    //        int seqnum = con.SequenceNumber();
    //        if (seqnum == -1) {
    //            continue;
    //        }
    //        int pos = mphysics->Block().Position(seqnum);
    //        int blocksize = mphysics->Block().Size(seqnum);
    //        int vs = active.size();
    //        active.Resize(vs+blocksize);
    //        for(int ieq = 0; ieq<blocksize; ieq++){
    //            active[vs+ieq] = pos+ieq;
    //        }
    //
    //    }
    //
    //    for(int i = ncon-ncon_saturation-ncon_gravity; i<ncon-ncon_gravity; i++)
    //    {
    //        TPZConnect &con = mphysics->ConnectVec()[i];
    //        int seqnum = con.SequenceNumber();
    //        int pos = mphysics->Block().Position(seqnum);
    //        int blocksize = mphysics->Block().Size(seqnum);
    //        int vs = active.size();
    //        active.Resize(vs+1);
    //
    //        int ieq = blocksize-1;
    //        active[vs] = pos+ieq;
    //    }
    //
    //    for(int i = ncon-ncon_saturation-ncon_gravity; i<ncon-ncon_gravity; i++)
    //    {
    //        TPZConnect &con = mphysics->ConnectVec()[i];
    //        int seqnum = con.SequenceNumber();
    //        int pos = mphysics->Block().Position(seqnum);
    //        int blocksize = mphysics->Block().Size(seqnum);
    //        int vs = nonactive.size();
    //        nonactive.Resize(vs+blocksize-1);
    //        for(int ieq = 0; ieq<blocksize-1; ieq++)
    //        {
    //            nonactive[vs+ieq] = pos+ieq;
    //        }
    //
    //    }
    
}

void TPZDarcyAnalysis::AnalyticSolution(const TPZVec< REAL >& pt, REAL time, TPZVec< STATE >& Saturation, TPZFMatrix< STATE >& Grad)
{
    REAL x = pt[0];
    REAL v = 0.00001;
    REAL Porosity = 0.1;
    REAL xshock = v*time/Porosity;
    
    if(x <= xshock)
    {
        Saturation[0] = 1.0;
    }
    else
    {
        Saturation[0] = 0.0;
    }
    
    return;
}

