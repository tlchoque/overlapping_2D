

void begin_verification(Delaunay &dt){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
		it->info().m_original_label = it->info().m_label;
}

void end_verification(Delaunay &dt, const char* filename){
	vector<Face> faces;
	double partial=0,total=0;
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		total += CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point() );
		if( it->info().m_original_label  != it->info().m_label && it->info().m_original_label !=-1){
			faces.push_back(it);
			partial += CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point() );
		}
	}	
	cout<<"total area: "<<total<<endl;
	cout<<"partial area: "<<partial<<endl;
	cout<<"ratio: "<<partial/total<<endl;

	// get a mesh of modified cells

	/*cout<<"faces size "<<faces.size()<<endl;
	getchar();*/

	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++)
		vi->info().m_index=-1;

	vector<Vertex> vertices;
	for(unsigned int i = 0; i < faces.size(); ++i){
		Face f = faces[i];
		for( unsigned int j = 0; j < 3; j++){
			Vertex v = f->vertex(j);
			if(v->info().m_index == -1){
				v->info().m_index=1;
				vertices.push_back(v);				
			}
		}
	}

	std::ofstream os;
	std::stringstream ss;
	ss<<filename<<"/"<<"relabeled_cells"<<".vtk";
	cout<<ss.str()<<endl;
	os.open( ss.str() );
	os << "# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID" << std::endl;
	size_t n,m;
	n = vertices.size();
	m = faces.size();
	os << "POINTS" << ' ' << n << ' ' << "float\n";
	int k = 0;	
	for(unsigned int j = 0; j < vertices.size(); ++j){
		vertices[j]->info().m_index=k;
		os << vertices[j]->point()<<" "<<0<< "\n";
		++k;
	}		
	os << "CELLS" << ' ' << m << ' ' << m*4<< '\n';
	for(unsigned int j = 0; j < faces.size(); ++j){
		os <<"3";
        for(unsigned int n=0; n < 3; ++n)
            os << ' ' << faces[j]->vertex(n)->info().m_index;
        os << '\n';
	}
	os << "CELL_TYPES" << ' ' << m << '\n';
	for(unsigned int j = 0; j < faces.size(); ++j)
		os << "5 ";
	os << "\n\n";

	os.close();
}

void get_inserted_points(Delaunay &dt,int old, const char* filename,int size){
	vector<Point> nwpoints;
	int count = 0;
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++){
		++count;
		if( count > old){
			nwpoints.push_back(vi->point() );
		}
	}
	savePoints2(nwpoints,filename,size); 
}