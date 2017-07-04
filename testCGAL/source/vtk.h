#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkLine.h>


void drawMesh(Delaunay & dt){
	unsigned char grey[3] = {192,192,192};
	unsigned char red[3] = {255, 0, 0};
	unsigned char blue[3] = {0, 0, 255};
	unsigned char green[3] = {0, 255, 0};
	unsigned char purple[3] = {153,0,153};
	unsigned char bluewhite[3] = {0, 255, 255};
	unsigned char brown[3] = {153, 76, 0};	
	unsigned char black[3] = {0, 0, 0};
	unsigned char pink[3] = {255, 102, 255};
	unsigned char *colorsArray[9] = {grey,red,blue,green,purple,bluewhite,brown,black,pink};

	Vertex v;
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	int k=0;
	for (Vertex_iterator vi = dt.vertices_begin(); vi != dt.vertices_end(); vi++) { 
		points->InsertNextPoint(vi->point().x(), vi->point().y(), 0.0);
		vi->info().index=k;
		++k;
	} 

	vtkSmartPointer<vtkUnsignedCharArray> colors =  vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetNumberOfComponents(3);
	colors->SetName("Colors");

	vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
	vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		colors->InsertNextTupleValue(*(colorsArray + it->info().m_label));
		triangle->GetPointIds()->SetId(0, it->vertex(0)->info().index);
		triangle->GetPointIds()->SetId(1, it->vertex(1)->info().index);
		triangle->GetPointIds()->SetId(2, it->vertex(2)->info().index);
		triangles->InsertNextCell(triangle);
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
		line->GetPointIds()->SetId(0,vs->info().index); 
		line->GetPointIds()->SetId(1,vt->info().index);
		lines->InsertNextCell(line);
		colorLines->InsertNextTupleValue(*(colorsArray + 7));
	}

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->SetPolys(triangles);
	polydata->GetCellData()->SetScalars(colors);
	
	vtkSmartPointer<vtkPolyData> linesPolyData = vtkSmartPointer<vtkPolyData>::New();
	linesPolyData->SetPoints(points);
	linesPolyData->SetLines(lines);
	linesPolyData->GetCellData()->SetScalars(colorLines);

	vtkSmartPointer<vtkPolyDataMapper> mapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
	vtkSmartPointer<vtkPolyDataMapper> mapper2 =  vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION <= 5
	mapper->SetInputConnection(polydata->GetProducerPort());
	mapper2->SetInputConnection(linesPolyData->GetProducerPort());
#else
	mapper->SetInputData(polydata);
	mapper2->SetInputData(linesPolyData);
#endif
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkActor> actor2 = vtkSmartPointer<vtkActor>::New();
	actor2->SetMapper(mapper2);
 
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->AddActor(actor2);
	/*renderer->SetBackground(.3, .6, .3);*/
	renderer->SetBackground(1, 1, 1);
	renderer->ResetCamera();
 
	renderWindow->Render();
	renderWindowInteractor->Start();
}