
double SafeAcos (double x){
	if (x < -1.0) x = -1.0 ;
	else if (x > 1.0) x = 1.0 ;
	return acos (x) ; 
}

double normVector(Vector v){
	return sqrt(v*v);
}

double getAngle(Vertex a, Vertex b,Vertex c){
	Vector x =a->point() - b->point();
	Vector y =c->point() - b->point();
	return acos( ( x*y )/ ( normVector( x ) * normVector( y ) ) )* 180.0 / PI;
	//return acos( ( x*y )/ ( normVector( x ) * normVector( y ) ) );
}

double get_angle(Vector x, Vector y){
	//return acos( ( x*y )/ ( normVector( x ) * normVector( y ) ) )* 180.0 / PI;

	//cout<<"angle: "<<acos( ( x*y )/ ( normVector( x ) * normVector( y ) ) )* 180.0 / PI<<endl;
	return SafeAcos( ( x*y )/ ( normVector( x ) * normVector( y ) ) );
}

double curvature(Vertex & prev_v, Vertex &v, Vertex &next_v){
	return  normVector( ( next_v->point() - v->point() ) - ( v->point() - prev_v->point() ) ) / normVector( (next_v->point() - prev_v->point() ) );
}

bool is_obtuse(Vector a , Vector b){
	return a.x()*b.y() - a.y()*b.x() < 0;
}

double curvature2(Vertex & prev_v, Vertex &v, Vertex &next_v){
	Vector t = next_v->point() - v->point();
	Vector prev_t = v->point() - prev_v->point();
	double angle2 = get_angle( prev_t,t );
	if( is_obtuse(prev_t,t) )
		angle2 = -angle2;

	//double angle = atan2( t.y(),t.x() ) - atan2( prev_t.y() , prev_t.x() ) * 180.0 / PI;
	//cout<<"angle: "<<angle<<" " <<angle2<<endl;
	return  angle2 / normVector( (next_v->point() - prev_v->point() ) );
}

double unsigned_curvature(Vertex & prev_v, Vertex &v, Vertex &next_v){
	Vector t = next_v->point() - v->point();
	Vector prev_t = v->point() - prev_v->point();
	double angle2 = get_angle( prev_t,t );
	return  angle2 / normVector( (next_v->point() - prev_v->point() ) );
}

vector<Face> faces_in_vertex_by_label(Delaunay &dt, Vertex v, int label){
	vector<Face> faces;

	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;

	do{
		if( dt.is_infinite(fc) ){++fc; continue;}
		if( fc->info().m_label == label)
			faces.push_back(fc);
		++fc;
	}while (fc!= fc_start);
	return faces;
}

bool vertex_curvature(Delaunay &dt, Vertex &v, int label, double &curvature){
	vector<Face> faces = faces_in_vertex_by_label(dt,v,label);
	/*Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		if( fc->info().m_label == label )
			faces.push_back( (Face)fc );
		++fc;
	}while (fc!= fc_start);*/

	vector< pair<Face,int> > boundary_faces;
	for(unsigned int i = 0; i < faces.size(); ++i ){
		int indexv = faces[i]->index(v);
		Face f = faces[i];
		Face f1 = f->neighbor( f->ccw(indexv));
		Face f2 = f->neighbor( f->cw(indexv));
		if( f1->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->cw(indexv) ) );
		if( f2->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->ccw(indexv) ) );
	}

	if(boundary_faces.size() == 2){
		Face bf1 = boundary_faces[0].first;
		Face bf2 = boundary_faces[1].first;
		Vertex prev_v,next_v;

		if( bf1->index( v ) ==  bf1->ccw(boundary_faces[0].second) ){// v is after
			prev_v = bf1->vertex(boundary_faces[0].second);
			next_v = bf2->vertex(boundary_faces[1].second);
		}
		else{
			prev_v = bf2->vertex(boundary_faces[1].second);
			next_v = bf1->vertex(boundary_faces[0].second);
		}
		curvature = curvature2( prev_v, v, next_v );

		/*cout<<test<<" "<<cur<<endl;
		vector<Point> vtkpoints;
		vtkpoints.push_back(v->point());
		vtkpoints.push_back(bf1->vertex(boundary_faces[0].second)->point());
		vtkpoints.push_back(bf2->vertex(boundary_faces[1].second)->point());
		drawMesh(dt,vtkpoints);
		getchar();*/

		return 1;
	}
	//cout<<"not possible"<<endl;
	return 0;
}

bool unsigned_vertex_curvature(Delaunay &dt, Vertex &v, int label, double &curvature){
	vector<Face> faces = faces_in_vertex_by_label(dt,v,label);

	vector< pair<Face,int> > boundary_faces;
	for(unsigned int i = 0; i < faces.size(); ++i ){
		int indexv = faces[i]->index(v);
		Face f = faces[i];
		Face f1 = f->neighbor( f->ccw(indexv));
		Face f2 = f->neighbor( f->cw(indexv));
		if( f1->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->cw(indexv) ) );
		if( f2->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->ccw(indexv) ) );
	}

	if(boundary_faces.size() == 2){
		Face bf1 = boundary_faces[0].first;
		Face bf2 = boundary_faces[1].first;
		Vertex prev_v,next_v;

		if( bf1->index( v ) ==  bf1->ccw(boundary_faces[0].second) ){// v is after
			prev_v = bf1->vertex(boundary_faces[0].second);
			next_v = bf2->vertex(boundary_faces[1].second);
		}
		else{
			prev_v = bf2->vertex(boundary_faces[1].second);
			next_v = bf1->vertex(boundary_faces[0].second);
		}
		curvature = unsigned_curvature( prev_v, v, next_v );
		return 1;
	}
	//cout<<"not possible"<<endl;
	return 0;
}

bool curvature_and_neighbors(Delaunay &dt, Vertex &v, int label, double &curvature,vector<Face> &faces, Vertex &prev_v, Vertex &next_v){
	faces.clear();
	faces = faces_in_vertex_by_label(dt,v,label);
	vector< pair<Face,int> > boundary_faces;
	for(unsigned int i = 0; i < faces.size(); ++i ){
		int indexv = faces[i]->index(v);
		Face f = faces[i];
		Face f1 = f->neighbor( f->ccw(indexv));
		Face f2 = f->neighbor( f->cw(indexv));
		if( f1->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->cw(indexv) ) );
		if( f2->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->ccw(indexv) ) );
	}

	if(boundary_faces.size() == 2){
		Face bf1 = boundary_faces[0].first;
		Face bf2 = boundary_faces[1].first;
		//Vertex prev_v,next_v;

		//cout<<"start if"<<endl;
		if( bf1->index( v ) ==  bf1->ccw(boundary_faces[0].second) ){// v is after
			prev_v = bf1->vertex(boundary_faces[0].second);
			next_v = bf2->vertex(boundary_faces[1].second);
			//cout<<"operator prev_v "<<bf1->vertex(boundary_faces[0].second)->point();
		}
		else{
			prev_v = bf2->vertex(boundary_faces[1].second);
			next_v = bf1->vertex(boundary_faces[0].second);
			//cout<<"operator prev_v "<<bf2->vertex(boundary_faces[1].second)->point();
		}
		//cout<<"operator prev_v "<<prev_v->point();
		curvature = curvature2( prev_v, v, next_v );
		return 1;
	}
	return 0;
}

bool get_adyacents(Delaunay &dt, Vertex &v, int label, Vertex &prev_v, Vertex &next_v){
	vector<Face> faces = faces_in_vertex_by_label(dt,v,label);
	vector< pair<Face,int> > boundary_faces;
	for(unsigned int i = 0; i < faces.size(); ++i ){
		int indexv = faces[i]->index(v);
		Face f = faces[i];
		Face f1 = f->neighbor( f->ccw(indexv));
		Face f2 = f->neighbor( f->cw(indexv));
		if( f1->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->cw(indexv) ) );
		if( f2->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->ccw(indexv) ) );
	}

	if(boundary_faces.size() == 2){
		Face bf1 = boundary_faces[0].first;
		Face bf2 = boundary_faces[1].first;
		//Vertex prev_v,next_v;

		//cout<<"start if"<<endl;
		if( bf1->index( v ) ==  bf1->ccw(boundary_faces[0].second) ){// v is after
			prev_v = bf1->vertex(boundary_faces[0].second);
			next_v = bf2->vertex(boundary_faces[1].second);
			//cout<<"operator prev_v "<<bf1->vertex(boundary_faces[0].second)->point();
		}
		else{
			prev_v = bf2->vertex(boundary_faces[1].second);
			next_v = bf1->vertex(boundary_faces[0].second);
			//cout<<"operator prev_v "<<bf2->vertex(boundary_faces[1].second)->point();
		}
		return 1;
	}
	return 0;
}

bool singular_vertex_curvature(Vertex &v, int label, vector<Face> &faces , double &curvature, Vertex &prev_v, Vertex &next_v  ){
	vector< pair<Face,int> > boundary_faces;
	for(unsigned int i = 0; i < faces.size(); ++i ){
		int indexv = faces[i]->index(v);
		Face f = faces[i];
		Face f1 = f->neighbor( f->ccw(indexv));
		Face f2 = f->neighbor( f->cw(indexv));
		if( f1->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->cw(indexv) ) );
		if( f2->info().m_label != label )
			boundary_faces.push_back( make_pair( f, f->ccw(indexv) ) );
	}

	if(boundary_faces.size() == 2){
		Face bf1 = boundary_faces[0].first;
		Face bf2 = boundary_faces[1].first;
		double test = 0;
		if( bf1->index( v ) ==  bf1->ccw(boundary_faces[0].second) ){// v is after
			prev_v = bf1->vertex(boundary_faces[0].second);
			next_v = bf2->vertex(boundary_faces[1].second);
		}
		else{
			prev_v = bf2->vertex(boundary_faces[1].second);
			next_v = bf1->vertex(boundary_faces[0].second);
		}
		curvature = curvature2( prev_v, v, next_v );
		return 1;
	}
	//cout<<"not possible"<<endl;
	return 0;
}

double total_vertex_curvature(Delaunay &dt, Vertex v){
	double sum = 0;
	for( unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		int label = v->info().m_regions_around[i];
		double k;
		if ( unsigned_vertex_curvature(dt,v,label,k) )
			sum += k;
	}
	return sum;
}