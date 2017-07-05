#define vtkRenderingCore_AUTOINIT 4(vtkRenderingOpenGL2, vtkInteractionStyle,vtkRenderingVolumeOpenGL2,vtkRenderingFreeType)

#include <Windows.h>// for threads

#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkCellType.h>
#include <vtkDataSetMapper.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLine.h>
#include <vtkTetra.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSphereSource.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkProperty.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkElevationFilter.h>
#include <vtkVectorText.h>
#include <vtkCommand.h>
#include <vtkAutoInit.h>

#include "vtkShrinkFilter.h"
#include "vtkStructuredGridReader.h"
#include "vtkExtractEdges.h"
#include "vtkGenericCell.h"
#include "vtkPolyhedron.h"
#include "vtkCubeSource.h"
#include "vtkIdList.h"
#include "vtkDataArray.h"
#include "vtkPointLocator.h"

#include <vtkRegularPolygonSource.h>
#include <vtkFloatArray.h>
#include <vtkGlyph3D.h>
#include <vtkGlyph2D.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkLookupTable.h>

#include <vtkPolyDataWriter.h>


void* handler = 0;
void* mtx;

#define VTK_CRT(type, name) \
    vtkSmartPointer<type> name = vtkSmartPointer<type>::New()

vtkSmartPointer<vtkRenderer> renderer =  vtkSmartPointer<vtkRenderer>::New();
vtkSmartPointer<vtkRenderWindow> renderWindow =  vtkSmartPointer<vtkRenderWindow>::New();	
vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
vtkInteractorStyleTrackballCamera *styleDisplay =         vtkInteractorStyleTrackballCamera::New();
vtkSmartPointer<vtkRenderWindowInteractor> irenDisplay = vtkSmartPointer<vtkRenderWindowInteractor>::New();

vtkSmartPointer<vtkPolyDataMapper> mapper1 =  vtkSmartPointer<vtkPolyDataMapper>::New();//for tetrahedrons
vtkSmartPointer<vtkPolyDataMapper> mapper2 =  vtkSmartPointer<vtkPolyDataMapper>::New();//for triangles
vtkSmartPointer<vtkPolyDataMapper> mapper3 =  vtkSmartPointer<vtkPolyDataMapper>::New();//for points
vtkSmartPointer<vtkDataSetMapper> mapper4 =  vtkSmartPointer<vtkDataSetMapper>::New();//for set cells

vtkSmartPointer<vtkPolyDataMapper> mapper5 =  vtkSmartPointer<vtkPolyDataMapper>::New();//for circles

vtkSmartPointer<vtkActor> actor1 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor3 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor4 = vtkSmartPointer<vtkActor>::New();
vtkSmartPointer<vtkActor> actor5 = vtkSmartPointer<vtkActor>::New();//for circles

unsigned char grey[3] = {192,192,192};
unsigned char red[3] = {255, 0, 0};	
unsigned char blue[3] = {0, 0, 255};
unsigned char green[3] = {0, 255, 0};
unsigned char purple[3] = {153,0,153};
unsigned char bluewhite[3] = {0, 255, 255};
unsigned char brown[3] = {153, 76, 0};	
unsigned char black[3] = {0, 0, 0};
unsigned char pink[3] = {255, 102, 255};
unsigned char orange[3] = {255, 102, 255};
unsigned char *colorsArray[10] = {grey,red,blue,green,purple,bluewhite,brown,black,pink,orange};

class CommandSubclass2 : public vtkCommand{
public:
    vtkTypeMacro(CommandSubclass2, vtkCommand);
    static CommandSubclass2 *New()    {
        return new CommandSubclass2;
    }
    void Execute(vtkObject *caller, unsigned long vtkNotUsed(eventId), 
        void *vtkNotUsed(callData))    {
        vtkRenderWindowInteractor *iren = static_cast<vtkRenderWindowInteractor*>(caller);
        iren->Render();
    }
};

unsigned long __stdcall displayVTK(void* param){
    renderWindow->SetSize(800, 600);
    renderWindow->AddRenderer(renderer);
    renderWindow->Render();
    irenDisplay->SetRenderWindow(renderWindow);
    irenDisplay->SetInteractorStyle(styleDisplay);// move slowly
    irenDisplay->Initialize();
    vtkSmartPointer<CommandSubclass2> timerCallback = vtkSmartPointer<CommandSubclass2>::New();
    irenDisplay->AddObserver ( vtkCommand::TimerEvent, timerCallback );
    irenDisplay->CreateRepeatingTimer(100);
    irenDisplay->Start();
    return 0;
}

void begin_VTK(Delaunay & dt){
	unsigned long id_thread;
    VTK_CRT(vtkVectorText, text);
    text->SetText("...");
    VTK_CRT(vtkElevationFilter, elevation);
    elevation->SetInputConnection(text->GetOutputPort());
	mapper1->SetInputConnection(elevation->GetOutputPort());
    mapper1->Update();
	
    actor1->SetMapper(mapper1);
	//actor1->GetProperty()->SetOpacity(0.7);
	actor2->SetMapper(mapper2);
	actor3->SetMapper(mapper3);
	actor3->GetProperty()->SetPointSize(15);
	actor4->SetMapper(mapper4);
	actor5->SetMapper(mapper5);
	actor5->GetProperty()->SetLineWidth(1.5);

    renderer->AddActor(actor1);
	renderer->AddActor(actor2);
	renderer->AddActor(actor3);
	renderer->AddActor(actor4);
	renderer->AddActor(actor5);
    renderer->SetBackground(1, 1, 1);	
    handler = CreateThread(0, 0, displayVTK, 0, 0, &id_thread);
    if(!handler){
        printf("Cannot create thread. Error code = %d\n", GetLastError());
        getchar();
        return;
    }	
}

void drawMesh(Delaunay & dt, vector<Point> &vertices,bool sa);

void drawMesh(Delaunay & dt);

void drawMesh(Delaunay & dt, vector<Point> &vertices);

void drawMesh(Delaunay & dt, bool sa);

void drawMesh(Delaunay & dt){
	bool sa = false;
	vector<Point> vertices;
	drawMesh(dt,vertices,sa);	
}

void drawMesh(Delaunay & dt, vector<Point> &vertices){
	bool sa = false;
	drawMesh(dt,vertices,sa);	
}

void drawMesh(Delaunay & dt, bool sa){
	vector<Point> vertices;
	drawMesh(dt,vertices,sa);
}

void drawMesh(Delaunay & dt, vector<Point> &vertices,bool sa){
	actor1->GetProperty()->SetOpacity(0.5);
	Vertex v;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k=0;
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) {
		points->InsertNextPoint(vi->point().x(), vi->point().y(), 0.0);
		vi->info().m_index=k;
		++k;
	} 
	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");	

	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		if( sa == false ){
			if(it->info().m_label > 0){
				colors->InsertNextTypedTuple(*(colorsArray + it->info().m_label));
				triangle->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
				triangle->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
				triangle->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
				triangles->InsertNextCell(triangle);
			}
		}
		else{
			if(it->info().m_solution > 0){
				colors->InsertNextTypedTuple(*(colorsArray + it->info().m_solution));
				triangle->GetPointIds()->SetId(0, it->vertex(0)->info().m_index);
				triangle->GetPointIds()->SetId(1, it->vertex(1)->info().m_index);
				triangle->GetPointIds()->SetId(2, it->vertex(2)->info().m_index);
				triangles->InsertNextCell(triangle);
			}
		}		
	}	
	//draw lines
	vtkSmartPointer<vtkUnsignedCharArray> colorLines = vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorLines->SetNumberOfComponents(3);
	colorLines->SetName("ColorsLines");

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();

	for(Edge_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		Face f = (ei->first);
		int i = ei->second;
		Vertex vs = f->vertex(f->cw(i));
		Vertex vt = f->vertex(f->ccw(i));
		line->GetPointIds()->SetId(0,vs->info().m_index); 
		line->GetPointIds()->SetId(1,vt->info().m_index);
		lines->InsertNextCell(line);
		colorLines->InsertNextTypedTuple(*(colorsArray + 7));
	}

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
	polydata->GetCellData()->SetScalars(colors);
	
	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);
	linesPolyData->GetCellData()->SetScalars(colorLines);

	vtkSmartPointer<vtkPoints> pointsVertex =   vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");

	for(unsigned int i=0;i<vertices.size();++i){
		pointsVertex->InsertNextPoint ( vertices[i].x(), vertices[i].y(), 0);
		//colorsp->InsertNextTypedTuple( *(colorsArray + (i+5)%10 ) );
		colorsp->InsertNextTypedTuple( *(colorsArray + i%10) );
	}
	vtkSmartPointer<vtkPolyData> pointsPolydata =   vtkSmartPointer<vtkPolyData>::New(); 
	pointsPolydata->SetPoints(pointsVertex);

	vtkSmartPointer<vtkVertexGlyphFilter> vertexFilter =    vtkSmartPointer<vtkVertexGlyphFilter>::New();
	vertexFilter->SetInputData(pointsPolydata);	
	vertexFilter->Update();

	vtkSmartPointer<vtkPolyData> polydatapc =    vtkSmartPointer<vtkPolyData>::New();
	polydatapc->ShallowCopy(vertexFilter->GetOutput());
	polydatapc->GetPointData()->SetScalars(colorsp);

	mapper1->SetInputData(polydata);
	mapper2->SetInputData(linesPolyData);
	mapper3->SetInputData(polydatapc);

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) 
		vi->info().m_index=-1;
}

void draw_relocating(vector<pair<Point,double>> & circles){	
	actor1->GetProperty()->SetOpacity(1);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	//cout<<"drawing circles size: "<<circles.size()<<endl;
	for(unsigned int i = 0; i < circles.size(); ++i){	
		points->InsertNextPoint( circles[i].first.x() , circles[i].first.y() , 0 ); // sphere in circle
        scales->InsertNextValue( 2*circles[i].second );
		//scales->InsertNextValue( circles[i].second );
		colorsp->InsertNextTypedTuple( *(colorsArray + 2) );
	}

	vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	grid->SetPoints(points);
    grid->GetPointData()->AddArray(scales);    
    grid->GetPointData()->SetActiveScalars("scales"); // !!!to set radius first
    grid->GetPointData()->AddArray(colorsp);	

	vtkSmartPointer<vtkRegularPolygonSource> polygonSource = vtkSmartPointer<vtkRegularPolygonSource>::New();
	polygonSource->GeneratePolygonOff();
	polygonSource->SetNumberOfSides(100);	

	vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
    glyph3D->SetInputData(grid);
    glyph3D->SetSourceConnection(polygonSource->GetOutputPort());
	mapper5->SetInputConnection(glyph3D->GetOutputPort());
}