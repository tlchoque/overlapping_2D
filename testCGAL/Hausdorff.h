
double point_edge(Point x, Edge e){
	Point y;
	Point p = e.first->vertex( e.first->ccw( e.second ) )->point();
	Point q = e.first->vertex( e.first->cw( e.second ) )->point();
	Vector w = unit_vector(q - p);
	double n = norm(q - p);
	Vector v = x - p;
	double par = (v*w)/n;

	if( par >= 0 && par <= 1){// is in the edge
		y = p + (q - p)*par;
	}
	else {
		if( par < 0 )	{y = p;}
		else y = q;
	}

	/*vector<Point> pts;
	pts.push_back(p);
	pts.push_back(q);
	pts.push_back(x);
	pts.push_back(y);
	drawMesh(dt,pts);
	getchar();*/
	
	return euclidian(x,y);
}

double min_distance(Point p, vector<Edge> edges){
	double min = 100000000000;
	for( unsigned int i = 0; i < edges.size() ; ++i){
		Edge e = edges[i];
		double d = point_edge(p,e);
		if( d < min ){
			min = d;
		}
	}
	return min;
}

vector<Edge> boundary_edges(Delaunay &dt){
	vector<Edge> edges;
	for(Finite_edges_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		if( ei->first->info().m_label == ei->first->neighbor(ei->second)->info().m_label )		continue;
		if( dt.is_infinite(ei->first) ||  dt.is_infinite(ei->first->neighbor(ei->second)) )		continue;
		Face c = ei->first;
		int index_facet = ei->second;
		if( ei->first->neighbor(ei->second)->info().m_label <  c->info().m_label ) {
			c = ei->first->neighbor(ei->second);
			index_facet = dt.mirror_index(ei->first, ei->second );
		}
		edges.push_back(Edge(c,index_facet));
	}
	return edges;
}

void simliarty(Delaunay &a, Delaunay &b,double sampling, const char* filename,string file, bool save){
	vector<Edge> aedges = boundary_edges(a);
	vector<Edge> bedges = boundary_edges(b);
	int index = 0;
	vector<pair<Point,double>> vtk_points;
	vector<pair<int,int>> sampled_edges;

	
	Vertex initial = a.finite_vertices_begin();
	double left = initial->point().x(), right = initial->point().x() , top  = initial->point().y() ,  bottom= initial->point().y();

	double total_lenght= 0;
	double max_H =0;
	for(unsigned int i =  0; i < aedges.size(); ++i){//edges in a
		Edge e = aedges[i];
		Vertex v = e.first->vertex( e.first->ccw( e.second ) );
		Vertex w = e.first->vertex( e.first->cw( e.second ) );

		Point p = v->point();
		Point q = w->point();

		if( p.x() < left ) left = p.x();
		if( p.x() > right ) right = p.x();
		if( p.y() > top ) top = p.y();
		if( p.y() < bottom ) bottom = p.y();

		if( q.x() < left ) left = q.x();
		if( q.x() > right ) right = q.x();
		if( q.y() > top ) top = q.y();
		if( q.y() < bottom ) bottom = q.y();

		total_lenght += euclidian(v->point(), w->point() );

		if(v->info().m_index == -1){
			double H = min_distance(v->point(),bedges);
			vtk_points.push_back( make_pair( v->point() ,H) );
			v->info().m_index = index; ++index;
			if( H > max_H) max_H = H;

		}
		if(w->info().m_index == -1){
			double H = min_distance(w->point(),bedges);
			vtk_points.push_back( make_pair( w->point() ,H) );
			w->info().m_index = index; ++index;
			if( H > max_H) max_H = H;
		}

		int s_index = v->info().m_index;
		int t_index;


		for(double j =  1; j < sampling + 1; ++j){//point created
			if( j == sampling) {
				t_index = w->info().m_index;
			}
			else{
				Point target = CGAL::ORIGIN +( (p - CGAL::ORIGIN) +  (j/sampling)*( q - p ) ) ;
				//cout<<"target: "<<target<<" j and sampling "<<j/sampling<<endl;
				double H = min_distance(target,bedges);
				vtk_points.push_back( make_pair(target,H) );
				t_index = index;
				++index;
				if( H > max_H) max_H = H;
			}

			//cout<<vtk_points[s_index].first<<" *** "<<vtk_points[t_index].first<<endl;
			sampled_edges.push_back( make_pair(s_index,t_index) );
			s_index = t_index;			
		}
	}
	//cout<<" l r t b : "<<left<<" * "<<right<<" * "<<top<<" * "<<bottom<<endl;
	double diag = euclidian(Point(left,bottom),Point(right,top));
	
	double RMS_p = 0;
	for(unsigned int i = 0; i < sampled_edges.size(); ++i){
		pair<int,int> t = sampled_edges[i];
		double mean = 0.5*( vtk_points[ t.first ].second + vtk_points[ t.second ].second );
		RMS_p += mean * euclidian( vtk_points[ t.first ].first, vtk_points[ t.second ].first); 
	}
	
	double RMS = RMS_p/total_lenght;
	/*cout<<"total_lenght: "<<total_lenght<<endl;
	cout<<"RMS_p: "<<RMS_p<<endl;*/
	cout<<"max: "<<max_H<<endl;
	cout<<"RMS: "<<RMS<<endl;
	cout<<"diag max: "<<max_H/diag<<endl;
	cout<<"diag RMS: "<<RMS/diag<<endl;

	if(save){
		std::ofstream os;
		std::stringstream ss;
		ss<<filename<<"/"<<file<<".vtk";
		//cout<<ss.str()<<endl;
		os.open( ss.str() );
		os << "# vtk DataFile Version 1.0\nUnstructured Grid Example\nASCII\n\nDATASET UNSTRUCTURED_GRID" << std::endl;
		size_t n,m;
		n = vtk_points.size();
		m = sampled_edges.size();
		os << "POINTS" << ' ' << n << ' ' << "float\n";
	
		for(unsigned int i = 0; i < vtk_points.size(); ++i){
			os << vtk_points[i].first<<" "<<0<< "\n";
		}		

		os << "CELLS" << ' ' << m << ' ' << m*3<< '\n';
		for(unsigned int i = 0; i < sampled_edges.size(); ++i){
			os <<"2"<< ' ' << sampled_edges[i].first<< ' ' << sampled_edges[i].second<< '\n';
		}

		os << "CELL_TYPES" << ' ' << m << '\n';
		for(unsigned int j = 0; j < sampled_edges.size(); ++j)
			os << "3 ";
		os << "\n";

		os << "POINT_DATA" << ' ' << n << "\n";
		os << "SCALARS scalars float 1"<< "\n";
		os << "LOOKUP_TABLE default"<< "\n";

		for(unsigned int i = 0; i < vtk_points.size(); ++i){
			os<< vtk_points[i].second<<' ';
		}
		os << "\n";
		os.close();
	}
	cout<<"hasudorff fnished"<<endl;
	for (Finite_vertices_iterator vi = a.finite_vertices_begin(); vi != a.finite_vertices_end(); vi++) 
		vi->info().m_index = -1;
}