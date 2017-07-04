int new_map(map<Point,int> &m, Point& p){
    std::map<Point,int>::const_iterator it = m.find( p );
    if (it == m.end())
        return -1;
    return it->second;
}

void update_regions_around_vertex(Delaunay &dt, Vertex &v){
	vector<int> &regions = v->info().m_singular_regions_around;
	vector<int> &labels = v->info().m_regions_around;
	regions.clear();
	labels.clear();
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;

	//to advance the begining of label
	bool isInternal = true;
	int faceLabel =  fc->info().m_label;
	do{
		if( faceLabel != fc->info().m_label){
			fc_start = fc;
			isInternal = false;
			break;
		}
		++fc;
	}while (fc!= fc_start);
	if( isInternal ) return;
	// end advance

	int labelcirculator = -2;
	do{
		int label =  fc->info().m_label;
		if( label == labelcirculator ){ fc++;continue;}
		labelcirculator = label;
		if( !dt.is_infinite(fc) && label != -1 ){			
			if( std::find(labels.begin(), labels.end(), label ) != labels.end() ){//found in labels
				if( !(std::find(regions.begin(), regions.end(), label ) != regions.end()) )//not foun in regions
					regions.push_back( label );
			}
			else
				labels.push_back( label );
		}		
		++fc;
	}while (fc!= fc_start);
}	

void decrease_decimals(Point &p){
	double dec= 1000;
	double x = floor(p.x()  * dec)/dec;
	double y = floor(p.y()  * dec)/dec;
	p = Point(x,y);
}

bool is_delaunay_after_displacement(Delaunay &dt, Vertex v, const Point &p){
	if(v->point() == p) return true;
	Point ant = v->point();
	v->set_point(p);
	// are incident faces well-oriented

	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		if( dt.orientation(fc->vertex(0)->point(), fc->vertex(1)->point(), fc->vertex(2)->point() ) ){
			v->set_point(ant);
			return false;
		}
		++fc;
	}while (fc!= fc_start);

	// are incident bi-Faces Delaunay?

	Edge_circulator ec_start = dt.incident_edges(v);
	Edge_circulator ec = ec_start;
	do{
		Face f = ec->first;
		int j = ec->second;
		Face fj = f->neighbor(j);
		int mj = dt.mirror_index(f, j);
		Vertex h1 = f->vertex(j);
		if(dt.is_infinite(h1)) {
			if( dt.side_of_oriented_circle(f, fj->vertex(mj)->point(), true) != CGAL::ON_UNBOUNDED_SIDE) {
				v->set_point(ant);
				return false;
			}
		} else {
			if(dt.side_of_oriented_circle(fj, h1->point(), true) != CGAL::ON_UNBOUNDED_SIDE) {
				v->set_point(ant);
				return false;
			}
		}
		++ec;
	}while (ec!= ec_start);
	v->set_point(ant);
	return true;
}

Point average_neighbor_boundary_vertices(Delaunay &dt,vector<Face> &faces, Vertex &v){
	Point p = v->point();
	vector<Vertex> vertices;
	vertices.push_back(v);
	Vertex v1;
	for(unsigned int i = 0; i < faces.size(); ++i){
		int idx = faces[i]->index(v);
		for(unsigned int j = 1; j < 3; ++j){
			v1 = faces[i]->vertex( (idx+j)%3 );
				if( v1->info().m_auxiliar == false ){
				p = Point(p.x() + v1->point().x(), p.y() + v1->point().y() );
				v1->info().m_auxiliar = true;
				vertices.push_back(v1);
			}
		}
	}
	for(unsigned int i = 0; i < vertices.size(); ++i)
		vertices[i]->info().m_auxiliar = false;
	if(vertices.size() != 0 )
		p = Point( p.x() / vertices.size(),p.y() / vertices.size() );
	return p; 
}

void hash_centroids_to_labels(Delaunay &dt, Vertex &v, vector<Face> &faces, map<Point,int> &pointLabel,Point p){// Faces doesn't have to be hashed
	Vertex v1,v2;
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		if ( std::find(faces.begin(), faces.end(), fc) != faces.end() )	{++fc;continue;}
		int idx = fc->index(v);
		v1 = fc->vertex( fc->ccw(idx) );
		v2 = fc->vertex( fc->cw(idx ) );
		Point cen = CGAL::centroid(v1->point(),v2->point(),p);
		decrease_decimals(cen);
		pointLabel[cen] = fc->info().m_label;
		++fc;
	}while (fc!= fc_start);

}

void check_if_Faces_set_is_out(Delaunay &dt, vector<Face> &faces, Vertex v,int label, deque<Vertex> &singularVertices){
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		if ( std::find(faces.begin(), faces.end(), fc) == faces.end() )	{++fc;continue;}// c is not in Faces
		double a2 = 0;int lbl = 0;int idx=0; // there will always be a facet 
		for(unsigned int j = 1; j < 3;++j){
			idx = (fc->index(v)+j)%3;
			if( dt.segment( Edge(fc, idx) ).squared_length() > a2 && fc->neighbor( idx )->info().m_label != label && fc->neighbor( idx )->info().m_label != -1) {
				a2 = dt.segment( Edge(fc, idx) ).squared_length();
				lbl = fc->neighbor( idx )->info().m_label;
			}
		}// at the end, we got the label & index			
		fc->info().m_label = lbl;
		fc->info().m_state = 1;	

		int index = fc->index(v);
		for(unsigned int j = 1; j < 3;++j){
			Vertex vc = fc->vertex( (index+j)%3 );
			//update_regions_around_vertex( dt, vc );
			update_vertex_info(dt,vc);
			if( vc->info().m_state != 2  && vc->info().is_singular() ){
				singularVertices.push_front( vc );
				vc->info().m_state=2;
			}
		}
		++fc;
	}while (fc!= fc_start);
}

int adjacent_label(Delaunay &dt,Face &c, Vertex &v,int label){
	Face c1 = c->neighbor( (c->index(v)+1)%3 );
	Face c2 = c->neighbor( (c->index(v)+2)%3 );

	if( c1->info().m_label != -1 && c1->info().m_label != label)
		return c1->info().m_label;

	else if( c2->info().m_label != -1 && c2->info().m_label != label)
		return c2->info().m_label;
	else return -1;
}

void relabel_with_neighbor(Delaunay &dt, vector<Face> &faces,Vertex &v,Vertex &s, int label){
	deque<Face> facesd;
	for(unsigned int i = 0; i <faces.size();++i)
		facesd.push_back(faces[i]);

	unsigned int sizec = 0;
	size_t faceSize = facesd.size();	
	while(facesd.size()!=0 && sizec < faceSize ){
		Face main = facesd.back();
		facesd.pop_back();
		Face c = main->neighbor(main->index(v));
		if(c->info().m_label != -1 && c->info().m_label != label){//get the oposite label			
			main->info().m_label = c->info().m_label;
		}
		else{// get the adjacent, but no oposite			
			int idx = adjacent_label(dt,main,v,label);
			if(idx != -1)	{main->info().m_label = idx;}
			else	{facesd.push_front(main);++sizec;}
		}
	}

	if( sizec >= faceSize){// there not exist and adjacent cel with different label 
		//cout<<" problem insert label"<<endl;
		int alternateLabel;
		for(unsigned int i = 0; i < s->info().m_regions_around.size(); ++i){
			int reg = s->info().m_regions_around[i];
			if( reg != label)
				alternateLabel = i;
		}	
		for(unsigned int i = 0; i <facesd.size();++i)
			facesd[i]->info().m_label = alternateLabel;
	}
}

void relabel_with_hash(Delaunay &dt,Vertex &v,Vertex &w, int label, map<Point,int> &pointLabel, deque<Vertex> &singularVertices){
	vector<Face> withoutLabel;
	Face_circulator fc_start = dt.incident_faces(w);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		Point cenc = CGAL::centroid(dt.triangle(fc));
		decrease_decimals(cenc);
		int mapLabel = new_map(pointLabel,cenc);
		if( mapLabel != -1){
			fc->info().m_label = mapLabel;
		}
		else
			withoutLabel.push_back(fc);
		++fc;
	}while (fc!= fc_start);
	if(withoutLabel.size()!=0){
		relabel_with_neighbor(dt,withoutLabel,w,v,label);
	}
}

void detect_new_singularities(Delaunay &dt,Vertex w,deque<Vertex> &singularVertices){
	Vertex_circulator vc_start = dt.incident_vertices(w);
	Vertex_circulator vc = vc_start;
	do{
		if(dt.is_infinite(vc)) {++vc;continue;}
		Vertex v1 = vc;
		//update_regions_around_vertex(dt, v1 );
		update_vertex_info(dt, v1 );
		if( v1->info().m_state != 2 && v1->info().is_singular_2() ){
			v1->info().m_state=2;
			singularVertices.push_back(v1);
		}
		++vc;
	}while (vc!= vc_start);

	//update_regions_around_vertex(dt, w );
	update_vertex_info(dt, w );
	if( w->info().is_singular_2() && w->info().m_state!= 2)
			singularVertices.push_back(w);
}

bool insert_point(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Face>> &faceSets,Vertex &v,int label,bool draw ){
	vector<double> ratio(faceSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,faceSets[l],v) , area_face_set(faceSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	for(unsigned int l = 0; l < ratio.size(); ++l)
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
	
	int max=0,min=1;
	if(ratio[min] > ratio[max]){
		max = 1;
		min = 0;
	}

	if(!v->info().m_boundary){
		Point pv1 = average_neighbor_boundary_vertices(dt,faceSets[max],v);		
		vector<pair<Point,double>> circles1;

		int step1 = getDisplacedPoint(dt,v,pv1,circles1);
		if(step1){//exist space to move
			if( is_delaunay_after_displacement(dt,v,pv1 ) ){
				//v->set_point(pv1);
			}
		}
	}		
	vector<Point> todraw;
	Point pv2 = average_neighbor_boundary_vertices(dt,faceSets[min],v);

	/*vector<Point> points;
	points.push_back(v->point());
	points.push_back(pv2);
	drawMesh(dt,points);
	getchar();*/

	vector<pair<Point,double>> circles;
	int step2 = getDisplacedPoint(dt,v,pv2,circles);

	/*Edge_circulator ec = dt.incident_edges(v);
	Face c1 = ec->first;
	Vertex v2 = c1->vertex( (ec->second+1)%3 );
	Vertex v3 = c1->vertex( (ec->second+2)%3 );
	double rec = sqrt( pow(v2->point().x() - v3->point().x() , 2 )  + pow(v2->point().y() - v3->point().y() , 2 )  );*/
	//circles.push_back(make_pair(v->point(),rec) );

	/*points.push_back(pv2);
	drawMesh(dt,points);
	draw_relocating(circles);
	getchar();*/

	if(step2){
		map<Point,int> pointLabel;
		hash_centroids_to_labels(dt,v,faceSets[max],pointLabel,pv2);	
		Vertex w = dt.insert(pv2);	
		check_if_Faces_set_is_out(dt,faceSets[min],v,label,singularVertices);		
		//cout<<"stop 4"<<endl;		
		relabel_with_hash(dt,v,w,label,pointLabel,singularVertices);		
		detect_new_singularities(dt,w,singularVertices);		
		//update_regions_around_vertex(dt, v );
		update_vertex_info(dt, v );
		if(draw){
			todraw.push_back(v->point());
			todraw.push_back(w->point());

			/*circles.clear();
			draw_relocating(circles);*/

			drawMesh(dt,todraw);
			//getchar();
		}
		return true;
	}
	else{
		cout<<"not possible inserted"<<endl;
		return false;
	}
}

// new insertion point

vector<Face> boundary_faces(Delaunay &dt,vector<Face> &erode, Vertex &v){
	int label = erode[0]->info().m_label;
	vector<Face> faces;
	Face vc,c0;
	for( unsigned int i = 0; i < erode.size(); ++i){
		vc = erode[i];
		int index = vc->index(v);
		for(unsigned int i = 1; i < 3;++i){
			c0 =  vc->neighbor( ( index+i )%3 );
			if( c0->info().m_label != label ){
				faces.push_back(vc);
				break;
			}
		}
	}
	return faces;
}

Point intersection_1(Delaunay &dt, vector<Face> &erode, Vertex &v){// one cell
	int index = erode[0]->index(v);
	Face c = erode[0];
	Vector a = c->vertex( (index+1)%3 )->point() - CGAL::ORIGIN;
	Vector b = c->vertex( (index+2)%3 )->point() - CGAL::ORIGIN;
	Vector m = (a + b)/2;

	Vector v_vec  =  v->point()- CGAL::ORIGIN;
	Vector p = m - v_vec ;
	Vector mid = v_vec + 0.5*p;
	//return  CGAL::ORIGIN + m;
	return  CGAL::ORIGIN + mid;
}

Point intersection_2(Delaunay &dt, vector<Face> &erode, Vertex &v){
	Face a = erode[0];
	Face b = erode[1];
	int index = a->index(v);
	int idx1 = (index+1)%3;
	int idx2 = (index+2)%3;
	Vector p;
	if( a->neighbor(idx1) == b )
		p = a->vertex(idx2)->point()- CGAL::ORIGIN;
	else
		p = a->vertex(idx1)->point()- CGAL::ORIGIN;	
	Vector vec = v->point() - CGAL::ORIGIN;
	//Vector s = vec + 0.75*( p - vec );
	Vector s = vec + 0.5*( p - vec );
	return CGAL::ORIGIN + s;
}

Point intersection_3(Delaunay &dt, vector<Face> &erode, Vertex &v){
	vector<Face> boundary = boundary_faces(dt,erode,v);
	Point c,d,m,x;
	Vector cv,cd,cm,ce,vx,test;
	Face a = boundary[0];
	Face b = boundary[1];
	c = dt.dual(a);
	cv = v->point() - c;
	d = dt.dual(b);
	cd = unit_vector(d - c);
	cm = ( cv*cd )*cd;		
	ce = cm-cv;// ce -> vector parallel to mx
	m = c + cm;// m-> middle point between v and x  HERE i can short the computation
	x = m + ce;	
	Point intersection = x;	
	intersection = (CGAL::ORIGIN + ( v->point() - CGAL::ORIGIN ) + 0.8*(intersection - v->point() ) );

	for( unsigned int i = 0; i < erode.size(); ++i){
		Face c = erode[i];
		if( dt.side_of_oriented_circle(c,intersection,true) != CGAL::ON_BOUNDED_SIDE){
			return v->point();
		}
	}
	return intersection;
}

bool circumcircle_intersection(Delaunay &dt,vector<Face> &erode, Vertex v, Point &collision){
	if( erode.size() == 1 )
		collision = intersection_1(dt,erode,v);
	else if( erode.size() == 2  )
		collision = intersection_2(dt,erode,v);
	else
		collision = intersection_3(dt,erode,v);
	if(collision == v->point() ) return false;
	return true;
}

void unvisited(Delaunay &dt, Vertex v, vector<Face> &preserved_faces ){
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		if( !fc->info().m_visited ){
			int index = fc->index(v);
			Face c = fc->neighbor(index);
			if( !dt.is_infinite(c) )
				preserved_faces.push_back(c);
		}
		else
			fc->info().m_visited = false;
		++fc;
	}while (fc!= fc_start);
}

vector<Face> mark_visited_boundary_faces(Delaunay &dt,vector<Face> &erode, Vertex &v){
	int label = erode[0]->info().m_label;
	vector<Face> faces;
	Face vc,c0;
	for( unsigned int i = 0; i < erode.size(); ++i){
		vc = erode[i];
		vc->info().m_visited=true;
		int index = vc->index(v);
		for(unsigned int i = 1; i < 3;++i){
			c0 = vc->neighbor( ( index+i )%3  );
			if( c0->info().m_label != label ){
				faces.push_back(vc);
				break;
			}
		}
	}
	return faces;
}

bool check_circle_collision(Delaunay &dt, vector<Face> &collision_faces, Point &collision,Vertex v){
	double t = 1,ta;
	for( unsigned int i = 0; i < collision_faces.size(); ++i){
		Face c = collision_faces[i];
		Point center = dt.dual(c);
		Vector r = center - c->vertex(0)->point();
		double radio = sqrt( r * r );
		int step = getNearPoint(v->point(),collision,center,radio,ta);
		if(step){
			if(ta < t)  t = ta;
		}
	}
	if(t > 0 && t <= 1){	
		t = 0.98*t;
		collision = v->point() + t*(collision - v->point()); 
		if( collision == v->point() )
			return 0;
		else	return 1;
	}
	else{
		collision = v->point();
		return 0;
	}
}

bool check_preserved_faces(Delaunay &dt, vector<Face> &preserved_faces, Point &collision,Vertex v){
	for( unsigned int i = 0; i < preserved_faces.size(); ++i){
		Face c = preserved_faces[i];
		if( dt.side_of_oriented_circle(c,collision,false) == CGAL::ON_BOUNDED_SIDE ){
			return false;
		}
	}
	return true;
}

bool safe_point(Delaunay &dt, vector<Face> &erode, vector<Face> &keep, Vertex v, Point &collision){	
	//cout<<"safe stop 0 "<<endl;
	if( !circumcircle_intersection(dt,erode,v, collision) ) { /*cout<<"not intersection: "<<endl;*/ return false;}
	//cout<<"safe stop 1 "<<endl;
	vector<Face> preserved_faces = mark_visited_boundary_faces(dt,keep,v);
	vector<Face> collision_faces;
	unvisited(dt,v,collision_faces);
	//cout<<"safe stop 2 "<<endl;
	if( !check_circle_collision(dt,collision_faces,collision,v) ) { /*cout<<"collision: "<<endl;*/ return false;}
	//cout<<"safe stop 3 "<<endl;}
	if( !check_preserved_faces(dt,preserved_faces,collision,v) ) { /*cout<<"preserved: "<<endl;*/ return false;}
	return true;
}

bool insert_point_2(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Face>> &faceSets,Vertex &v,int label,int min,int max ){
	Point p;
	if( safe_point(dt,faceSets[max],faceSets[min],v,p) ){
		//double total = area_face_set(faceSets[max] );

		//cout<<" point inserted 1"<<endl;
		map<Point,int> pointLabel;
		map<Point,int> pointOriginaLabel;
		/*hash_centroids_to_labels(dt,v,cellSets[min],pointLabel,p);	
		Vertex w = dt.insert(p, cellSets[min][0]);	
		relabel_with_hash(dt,v,w,cellSets[min][0]->info().m_label,pointLabel,singularVertices,holes);*/
		hash_centroids_to_labels(dt,v,faceSets[min],pointLabel,p);	
		Vertex w = dt.insert(p, faceSets[max][0]);	
		relabel_with_hash(dt,v,w,faceSets[min][0]->info().m_label,pointLabel,singularVertices);
		detect_new_singularities(dt,w,singularVertices);
		//update_regions_around_vertex(dt, v );
		update_vertex_info(dt, v );
		/*if( v->info().is_singular_2 && v->info().m_state!= 2)
			singularVertices.push_back(v);*/


		/*vector<Face> faces =faces_in_vertex_by_label(dt,w,faceSets[min][0]->info().m_label);
		double partial = area_face_set(faces);
		total_area +=total - partial;*/
		
		return true;
	}	
	else if( safe_point(dt,faceSets[min],faceSets[max],v,p) ){
		//double total = area_face_set(faceSets[min] );

		//cout<<" point inserted 2"<<endl;
		map<Point,int> pointLabel;
		map<Point,int> pointOriginaLabel;
		/*hash_centroids_to_labels(dt,v,cellSets[max],pointLabel,p);	
		Vertex w = dt.insert(p, cellSets[max][0]);	
		relabel_with_hash(dt,v,w,cellSets[max][0]->info().m_label,pointLabel,singularVertices,holes);*/
		hash_centroids_to_labels(dt,v,faceSets[max],pointLabel,p);	
		Vertex w = dt.insert(p, faceSets[min][0]);	
		relabel_with_hash(dt,v,w,faceSets[max][0]->info().m_label,pointLabel,singularVertices);

		detect_new_singularities(dt,w,singularVertices);
		//update_regions_around_vertex( dt, v );
		update_vertex_info(dt, v );


		/*vector<Face> faces =faces_in_vertex_by_label(dt,w,faceSets[max][0]->info().m_label);
		double partial = area_face_set(faces);
		total_area +=total - partial;*/

		return true;
	}
	else{
		//cout<<"not possible inserted"<<endl;
		return false;
	}
}

bool insert_to_erode(Delaunay &dt, deque<Vertex> &singularVertices, vector<vector<Face>> &faceSets,Vertex &v,int label,int min, int max,bool draw ){		
	//return false;
	return insert_point_2(dt,singularVertices,faceSets,v,label,min,max);
}
