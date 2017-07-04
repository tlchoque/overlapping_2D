
double sum_norm_square(Delaunay &dt,Edge &e){
	Vertex v0 = e.first->vertex( (e.second) );
	double sum=0;
	for(int i = 1; i < 3 ; ++i){
		Vector vec = e.first->vertex( (e.second+i)%3 )->point() - v0->point();
		sum += vec*vec;
	}
	return sum;
}

double area_around_vertex(Delaunay &dt, Vertex &v){
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	double sum=0;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		sum += abs( dt.triangle(fc).area() );
		++fc;
	}while (fc!= fc_start);
	return sum;
}

Vector unit_vector_perpendicular(Delaunay &dt,Edge &e){
	Point p1 = e.first->vertex((e.second + 1) % 3)->point();
	Point p2 = e.first->vertex((e.second + 2) % 3)->point();
	Vector v = p2 - p1;
	v = Vector(- v.y() , v.x() );
	return unit_vector(v);
}

bool is_embedded_after_displacement(Delaunay &dt, Vertex v, const Point &p){
	if(v->point() == p) return true;
	Point ant = v->point();
	v->set_point(p);
	// are incident faces well-oriented
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		if( dt.orientation(fc->vertex(0)->point(), fc->vertex(1)->point(), fc->vertex(2)->point() ) != CGAL::POSITIVE ){
			v->set_point(ant);
			//cout<<"oriented "<<endl;
			return false;
		}
		++fc;
	}while (fc!= fc_start);
	v->set_point(ant);
	return true;
}

Point inner_smoothed_point(Delaunay &dt, Vertex &v){
	Vector total = Vector(0,0);
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		Edge e = Edge(fc,fc->index(v) );
		double sumNS= sum_norm_square(dt,e);
		double S = sqrt( dt.segment(e).squared_length() );
		//cout<<"S: "<<S<<endl;
		Vector n = unit_vector_perpendicular(dt,e);
		//cout<<"n: "<<n<<endl;
		total = total + sumNS*S*n;
		//cout<<"total: "<<total<<endl;
		++fc;
	}while (fc!= fc_start);
	double ohmio = area_around_vertex(dt,v);
	//cout<<"ohmio: "<<ohmio<<endl;
	
	Point x = v->point() - total/(4*ohmio);
	/*cout<<"new point: "<<x<<endl;
	cout<<"old point: "<<v->point()<<endl;*/
	//getchar();
	return x;
}

bool inner_smooth(Delaunay &dt, Vertex &v){
	Point p = inner_smoothed_point(dt,v);
	/*cout<<"new point: "<<p<<endl;
	cout<<"old point: "<<v->point()<<endl;

	vector<Point> points;
	points.push_back(p);
	points.push_back(v->point());
	drawMesh(dt,points);
	getchar();*/

	if( is_embedded_after_displacement(dt,v,p) ){
		v->set_point(p);
		//cout<<"setted point"<<endl;
		return true;
	}
	//cout<<"not setted point"<<endl;
	return false;
}

double det(Vector v, Vector w){
	return v.x()*w.y() - v.y()*w.x();
}

Point displacement(Delaunay &dt, Vertex &v){
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	double t = 1, t_aux;
	Vector p = v->point() - CGAL::ORIGIN;
	Vector r = v->info().m_smoothed_point - v->point();

	/*vector<Point> points;
	points.push_back(v->point());
	points.push_back(v->info().m_smoothed_point);
	drawMesh(dt,points);
	getchar();*/

	Vector q,s;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		int i = fc->index(v);
		q = fc->vertex(fc->ccw(i))->point() - CGAL::ORIGIN;
		s = fc->vertex(fc->cw(i))->point() - fc->vertex(fc->ccw(i))->point();
		t_aux = det( ( q - p) , s ) / det(r,s);
		if( det( r, s) != 0 ){// not parallel, not colinear
			if( t_aux >= 0 && t_aux <= 1){// intersection
				if( t_aux < t)
					t = t_aux;
			}
		}
		++fc;
	}while (fc!= fc_start);
	Vector intersection;
	if(t == 1)
		intersection = p + t*r;
	else
		intersection = p + t*0.98*r;
	/*cout<<"t: "<<t<<endl;
	points.push_back(CGAL::ORIGIN + intersection);
	drawMesh(dt,points);
	getchar();*/
	return CGAL::ORIGIN + intersection;
}

bool boundary_smooth( Delaunay &dt, Vertex &v ){
	if(v->info().m_smoothed_point == v->point() ) {
		//cout<<" problem"<<endl;
		return false;
	}
	//cout<<"not problem "<<endl;
	Point p = displacement(dt,v);
	/*cout<<"new point: "<<p<<endl;
	cout<<"old point: "<<v->point()<<endl;

	vector<Point> points;
	points.push_back(v->point());
	points.push_back(p);	
	drawMesh(dt,points);
	getchar();*/

	if( is_embedded_after_displacement(dt,v,p) ){
		v->set_point(p);
		//cout<<"setted boundary point "<<endl;
		return true;
	}
	//cout<<"not setted boundary point"<<endl;
	return false;
}

void ODT(Delaunay &dt, int num){
	mark_boundary(dt);
	vector<Vertex> boundary;
	boundary_vertices(dt,boundary);
	for( int i = 0 ; i < num ; ++i){
		for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) {
			Vertex v = vi;
			if( v->info().is_corner() ){
			}
			else if( v->info().is_boundary() ){
				//cout<<"stop 1"<<endl;
				boundary_smooth(dt,v);
			}		
			else{//inner
				inner_smooth(dt,v);
			}
		}
	}
}