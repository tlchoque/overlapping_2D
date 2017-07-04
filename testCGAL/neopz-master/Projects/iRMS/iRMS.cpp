#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZIntQuadQuarterPoint.h"
#include <time.h>
#include "pzgmesh.h"
#include "TPZRefPatternTools.h"

#include "TRMRawData.h"
#include "TRMSimworxMeshGenerator.h"
#include "TRMOrchestra.h"
#include "TRMSpaceOdissey.h"

#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"


// iRMS (i rise for / innovatory Reservoir Muli-scale Simulator)

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.iRMS"));
#endif

void LinearTracerPrimal();
void LinearTracerDual();
void BoxLinearTracerDual();
void CheckQuarterPoint();
void CreateExampleRawData(TRMRawData &data);

int main()
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    // This code use normalized piola contravariant mapping for nonlinear mappings
    HDivPiola = 1;
    
    int ele_id=0;
    int mat_id=0;

    TPZGeoMesh *  geometry = new TPZGeoMesh;
    geometry->NodeVec().Resize(8);
    geometry->ElementVec().Resize(1);
    
    TPZManVector<long, 2> topology_l(2,0);
    
    TPZGeoElRefPattern< pzgeom::TPZGeoLinear > * line = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear >(ele_id,topology_l,mat_id,*geometry);
    
    

    // Running primal problem
//    LinearTracerPrimal();

    // Running dual problem on box shape
//    BoxLinearTracerDual();
    
//    // Running dual problem on Reservoir
    LinearTracerDual();
    
    
    std::cout << "Process complete normally." << std::endl;
    return 0;
}

void LinearTracerPrimal()
{

    TRMOrchestra  * SymphonyX = new TRMOrchestra;
    SymphonyX->CreateAnalysisPrimal();
    std::cout << "Primal complete normally." << std::endl;
    
}

void LinearTracerDual()
{
    
    TRMOrchestra  * SymphonyX = new TRMOrchestra;
    SymphonyX->CreateAnalysisDual();
    std::cout << "Dual complete normally." << std::endl;       
    
}

void BoxLinearTracerDual()
{
    
    TRMOrchestra  * SymphonyX = new TRMOrchestra;
    SymphonyX->CreateAnalysisDualonBox();
    std::cout << "Dual complete normally." << std::endl;
    
}

void CheckQuarterPoint()
{
    TPZIntQuadQuarterPoint qt(10);
    qt.SetCorner(1);
//    TPZIntQuad qt(60);
    int np = qt.NPoints();
    TPZManVector<REAL,2> pt(2,0.);
    std::cout << "Numpoints = " << np << std::endl;
    REAL weight;
    REAL integral = 0.;
    for (int ip = 0; ip<np; ip++) {
        qt.Point(ip, pt, weight);
        std::cout << "ip " << ip << " pt " << pt << " weight " << weight << std::endl;
        REAL r = sqrt((pt[0]-1)*(pt[0]-1)+(pt[1]+1)*(pt[1]+1));
        integral += weight/r;
    }
    std::cout << "Integral " << integral << std::endl;
    
}

