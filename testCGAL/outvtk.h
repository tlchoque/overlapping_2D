
#include <fstream>
#include <sstream>
#include <string>

void saveMesh(Delaunay &dt, const char* filename,string file ){
	vector<vector<Face>> faceByLabel(nroLabels);
	vector<vector<Vertex>> vertexByLabel(nroLabels);

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;

	//cout<<"save stop 0"<<endl;

	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		if(it->info().m_label  < 0) {
			//cout<<"label "<<it->info().m_label<<endl;

			vector<Point> pts;
			pts.push_back(it->vertex(0)->point() );
			pts.push_back(it->vertex(1)->point() );
			pts.push_back(it->vertex(2)->point() );
			drawMesh(dt,pts);
			getchar();
			continue;
		}
		faceByLabel[it->info().m_label].push_back(it);		
	}

	//cout<<"save stop 1"<<endl;
	for(unsigned int i = 0; i < faceByLabel.size(); ++i){
		for(unsigned int j = 0; j < faceByLabel[i].size(); ++j){
			Vertex v0 = faceByLabel[i][j]->vertex(0);
			Vertex v1 = faceByLabel[i][j]->vertex(1);
			Vertex v2 = faceByLabel[i][j]->vertex(2);
			if(v0->info().m_index != i){
				vertexByLabel[i].push_back(v0);
				v0->info().m_index=i;
			}
			if(v1->info().m_index != i){
				vertexByLabel[i].push_back(v1);
				v1->info().m_index=i;
			}
			if(v2->info().m_index != i){
				vertexByLabel[i].push_back(v2);
				v2->info().m_index=i;
			}
		}
	}

	//cout<<"save stop 2"<<endl;
	for(unsigned int i = 0; i < nroLabels; ++i){
		std::ofstream os;
		std::stringstream ss;
		ss<<filename<<"/"<<file<<i<<".vtk";
		//cout<<ss.str()<<endl;
		os.open( ss.str() );
		os << "# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID" << std::endl;
		size_t n,m;
		n = vertexByLabel[i].size();
		m = faceByLabel[i].size();
		os << "POINTS" << ' ' << n << ' ' << "float\n";
		int k = 0;	
		for(unsigned int j = 0; j < vertexByLabel[i].size(); ++j){
			vertexByLabel[i][j]->info().m_index=k;
			os << vertexByLabel[i][j]->point()<<" "<<0<< "\n";
			++k;
		}		

		os << "CELLS" << ' ' << m << ' ' << m*4<< '\n';
		for(unsigned int j = 0; j < faceByLabel[i].size(); ++j){
			os <<"3";
            for(unsigned int n=0; n < 3; ++n)
                os << ' ' << faceByLabel[i][j]->vertex(n)->info().m_index;
            os << '\n';
		}

		os << "CELL_TYPES" << ' ' << m << '\n';
		for(unsigned int j = 0; j < faceByLabel[i].size(); ++j)
			os << "5 ";
		os << "\n\n";

		os.close();
	}	
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;
	cout<<"save finished"<<endl;
}

void savePoints(vector<Point> &points, const char* filename){
	ofstream os(filename);
	os << "# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID" << std::endl;
	os << "POINTS" << ' ' << points.size() << ' ' << "float\n";

	for(unsigned int i = 0; i < points.size(); ++i)
		os << points[i]<<" "<<0<< "\n";

	os << "CELLS" << ' ' << points.size() << ' ' << points.size()*2<< '\n';
	for(unsigned int i = 0; i < points.size(); ++i){
		os <<"1"<< ' ' << i<< '\n';
	}
	os << "CELL_TYPES" << ' ' << points.size() << '\n';
	for(unsigned int i = 0; i < points.size(); ++i){
		os << "1 ";
	}
	os << "\n\n";
	os.close();

}

void savePoints2(vector<Point> &circles, const char* filename,double radio){
	std::stringstream ss;
	ss<<filename<<radio<<".vtk";

	actor1->GetProperty()->SetOpacity(1);
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkFloatArray> scales = vtkSmartPointer<vtkFloatArray>::New();
    scales->SetName("scales");
    vtkSmartPointer<vtkUnsignedCharArray> colorsp =   vtkSmartPointer<vtkUnsignedCharArray>::New();
	colorsp->SetNumberOfComponents(3);
	colorsp->SetName ("Colors");
	for(unsigned int i = 0; i < circles.size(); ++i){	
		points->InsertNextPoint( circles[i].x() , circles[i].y() , 0 ); // sphere in circle
        scales->InsertNextValue( radio );
		colorsp->InsertNextTypedTuple( *(colorsArray + 4) );
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

	vtkSmartPointer<vtkPolyDataWriter> writer =    vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetInputConnection(glyph3D->GetOutputPort());

	

	writer->SetFileName((ss.str() ).c_str());
	writer->Write();
}
