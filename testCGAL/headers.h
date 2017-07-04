
#define PI 3.14159265
const int nroLabels = 8;

double total_area=0;

#define HOLE 7

#include <iostream>

#include <functional>  
using namespace std;

#include <cstring>
#include <string>

#include <math.h>
#include <fstream> 
#include <ctime>
#include "imesh2D.h"
#include "IOimesh.h"
#include "globalVTK.h"
#include "outVTK.h"//  read and write files
#include "generation.h"
#include "segmentation.h"

#include "verify.h"
#include "operator.h"
#include "movement.h"
#include "singular_vertices.h"

#include "Hausdorff.h"

#include "winding.h"

#include "criteria.h"
#include "join.h"
#include "point_insertion.h"
#include "repairing.h"
#include "repairing_mesh.h"
#include "simulated.h"
#include "peak_removal.h"

#include "histogram.h"

#include "relabeling.h"
#include "smoothing.h"
#include "optimal.h"