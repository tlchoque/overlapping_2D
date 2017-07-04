bool produce_singularities( Delaunay &dt, vector<Vertex> &vertices ){
	for(unsigned int i = 0; i < vertices.size(); ++i ){
		//if( vertices[i]->info().has_more_than_two_regions() || vertices[i]->info().is_corner() ){
			if( vertices[i]->info().m_state != 2   )
				if( is_singular_after_relabeling_2(dt, vertices[i]) )
					return true;
		//}
	}
	return false;
}

bool hole_face(Face &c,deque<Face> &holes){
	c->info().m_label = HOLE;
	holes.push_back(c);
	//if( !c->info().m_hole_state ){// state point f the cell is in the deque		
	//	c->info().m_hole_state = true;
	//	holes.push_back(c);
	//	return true;
	//}		
	return false;
}

void make_hole(Delaunay &dt, vector<Face> &faces, deque<Vertex> &singularities,Vertex &v, deque<Face> &holes){
	for(unsigned int i = 0; i < faces.size(); ++i ){
		Face c = faces[i];
		hole_face(c,holes);
	}
	vector<Vertex> neighbors;
	/*vector<Point> pts;
	pts.push_back(v->point());*/
	neighbor_vertices(dt,faces,neighbors,v);
	for(unsigned int i = 0; i < neighbors.size(); ++i ){
		//pts.push_back(neighbors[i]->point());
		update_vertex_info( dt,neighbors[i] );
		if( neighbors[i]->info().is_singular_2() && neighbors[i]->info().m_state != 2 ){
			singularities.push_back( neighbors[i] );
			neighbors[i]->info().m_state = 2;
			//cout<<"is singular"<<endl;
		}
		/*drawMesh(dt,pts);
		getchar();*/
	}
	update_vertex_info( dt, v );
	/*if( v->info().is_singular_2() && v->info().m_state != 2 ){
		singularities.push_back( v );
		v->info().m_state = 2;
	}*/

	
	/*drawMesh(dt,pts);
	getchar();*/
}

bool relabel_for_simulated(Delaunay & dt, vector<Face> &faces, int other_label,Vertex v, vector<Face> &modified){
	vector<Vertex> neighbors;
	neighbor_vertices(dt,faces,neighbors,v);

	//cout<<"relabel stop 0"<<endl;

	int old_label = faces[0]->info().m_label; 
	//cout<<"relabel stop 1"<<endl;
	for(unsigned int i = 0; i < faces.size(); ++i )
		faces[i]->info().m_label = other_label;

	//cout<<"relabel stop 2"<<endl;
	if( produce_singularities(dt,neighbors) ){
		//cout<<"relabel stop 3 1"<<endl;
		for(unsigned int i = 0; i < faces.size(); ++i )
			faces[i]->info().m_label = old_label;
		return false;
	}
	else{
		//cout<<"relabel stop 3 2"<<endl;
		for(unsigned int i = 0; i < neighbors.size(); ++i ){
			update_vertex_info( dt,neighbors[i] );
		}
		modified.insert(modified.end(), faces.begin(), faces.end());//append the c faces to modified faces
		update_vertex_info( dt, v );
		return true;
	}
}

bool erode_sa( Delaunay &dt,vector<vector<Face>> &facesSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Face> &holes, vector<Face> &modified){	
	int min,max;
	bool change = false;

	//cout<<"stop 0 , size "<<facesSets.size()<<endl;
	sort_by_criteria(dt,facesSets,v,min,max,change);
	int newlabel = new_label_for_relabeling(dt,facesSets[min],v,label); // change because it should not be HOLE
	int index;

	//cout<<"stop 1"<<endl;
	double r = ((double) rand() / (RAND_MAX));
	if( r < 0.4) index = min;
	else index = max;

	//cout<<"stop 2"<<endl;
	if( newlabel == HOLE )	{
		make_hole(dt,facesSets[index], singularities,v,holes);
		return false;
	}
	//cout<<"stop 3"<<endl;
	if( !relabel_for_simulated(dt,facesSets[index],newlabel,v,modified ) ){ // test if we canerode change to MIN
		//if( !insert_point_2(dt,singularities,facesSets,v,label,min,max,holes,draw) ){
			//cout<<"make hole"<<endl;
		//cout<<"stop 4 "<<index<<endl;
			make_hole(dt,facesSets[index], singularities,v,holes); //change to MIN
			return false;
		//}
	}
	return true;
}

bool relabel_path_sa(Delaunay & dt, vector<Face> &faces, int alabel,int blabel,Vertex v, deque<Vertex> &singularities){
	if( !check_faces_state(faces) ){
		return false;
	}

	deque<Face> queue;
	for( unsigned int i = 0; i< faces.size(); ++i)
		queue.push_back(faces[i]);	

	while( queue.size() != 0 ){
		Face c = queue.front();
		queue.pop_front();
		int lbl = adyacent_face_label(dt,c,alabel,blabel);
		if( lbl != -1){
			c->info().m_label = lbl;
			//state of the cell
			c->info().m_state = true;
		}
		else{
			queue.push_back(c);
		}
	}
	vector<Vertex> restartIndex;
	int k = 1;
	for(unsigned int i = 0; i < faces.size(); ++i){
		Face c = faces[i];
		int vindex = c->index(v);
		for( unsigned int j = 1; j < 3; ++j){//			
			int index = ( vindex + j )%3;
			Vertex w = c->vertex(index);
			if(w->info().m_index == -1){
				w->info().m_index = k;
				restartIndex.push_back(w);
				update_vertex_info(dt,w);
				if( w->info().is_singular_2() && w->info().m_state!= 2)
					singularities.push_back(w);
			}
		}
	}
	update_vertex_info(dt,v);
	if( v->info().is_singular_2() && v->info().m_state!= 2)
		singularities.push_back(v);

	for(unsigned int i = 0 ; i < restartIndex.size(); ++i)
		restartIndex[i]->info().m_index=-1;
	return true;
}

bool dilate_sa( Delaunay &dt,vector<vector<Face>> &faceSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Face> &holes, vector<Face> &modified){	
	if(faceSets.size() != 2 ) return false;
	vector<Face> path;
	path_to_join(dt,faceSets,v,path);

	double s1 = criteria_1(dt,faceSets[0],v);
	double s2 = criteria_1(dt,faceSets[1],v);
	double s3 = criteria_1(dt,path,v);

	double s = s1;
	if( s2 > s1) s = s2;
	if( s/s3 > 2 ) return false;	
	return relabel_for_simulated(dt,path,label,v,modified) ; 
}

bool solve_2(Delaunay &dt,vector<vector<Face>> &faceSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Face> &holes, vector<Face> &modified){
	if( ( !is_erodible(dt,faceSets[0],v,label) || !is_erodible(dt,faceSets[1],v,label) ) 
		&& !v->info().is_corner() )
	{	
	//if( 1 == 2){//changed
		
		//cout<<"dilate decision"<<endl;
		if( v->info().has_two_regions() ){
			//cout<<"dilate 2 region"<<endl;
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			vector<vector<Face>> faceSets_2;
			set_faces_around_singular_vertex(dt,v,opposite_label,faceSets_2);
			//cout<<"after set"<<endl;
			if( faceSets_2.size() == 2 && !erode_sa(dt, faceSets_2,singularities,v,opposite_label,draw,holes,modified) ){
				//cout<<"pass erode"<<endl;
				dilate_sa(dt,faceSets,singularities,v,label,draw,holes,modified);
			}
		}
		else{
			//cout<<"dilate + 2 region"<<endl;
			dilate_sa(dt,faceSets,singularities,v,label,draw,holes,modified);
		}
	}
	set_faces_around_singular_vertex(dt,v,label,faceSets);

	if( faceSets.size() == 1 )
		return true;
	else{
		if( !erode_sa(dt, faceSets,singularities,v,label,draw,holes,modified) ){
			return dilate_sa(dt,faceSets,singularities,v,label,draw,holes,modified);
		}
		return false;
	}

	// here we go
	/*if( faceSets.size() == 1 )
		return true;
	else{
		if( !dilate_sa(dt,faceSets,singularities,v,label,draw,holes,modified) ){
			return erode_sa(dt, faceSets,singularities,v,label,draw,holes,modified);
		}
		return false;
	}*/
}

bool solve_3(Delaunay &dt,vector<vector<Face>> &faceSets,deque<Vertex> &singularities,Vertex &v,int label,bool draw, deque<Face> &holes, vector<Face> &modified){
	vector<double> ratio(faceSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	//cout<<"solve 3 stop -1: "<<faceSets.size()<<endl;
	for(unsigned int i = 0; i < ratio.size(); ++i){
		ratio[i] = criteria_1(dt,faceSets[i],v);
	}
	int v1=0,v2=1;
	if(ratio[v2] < ratio[v1]){
		v1 = 1;
		v2 = 0;
	}
	//cout<<"solve 3 stop 0"<<endl;
	for(unsigned int l = 2; l < ratio.size(); ++l){// choose the 2 greater volumes v1 < v2
		if(ratio[l] < ratio[v1]){
			int aux = v1;
			v1 = l;
			v2 = aux;
		}
		else if(ratio[l] < ratio[v2])   v2 = l;
	}
	//cout<<"solve 3 stop 1"<<endl;

	for(unsigned int l = 0; l < faceSets.size(); ++l){
		if( l != v1 && l!=v2){
			int newlabel = new_label_for_relabeling(dt,faceSets[l],v,label);
			if( newlabel == HOLE) {
				make_hole(dt,faceSets[l],singularities,v,holes);
			}
			else if(!relabel_for_simulated(dt,faceSets[l],newlabel,v,modified ) ){
				make_hole(dt,faceSets[l],singularities,v,holes);
			}
			faceSets.erase(faceSets.begin()+l);
			--l;
			--v1;--v2;
		}
	}
	cout<<"solve 3 stop 2"<<endl;
	return solve_2(dt,faceSets,singularities,v,label,draw,holes,modified);
}

void choose_holes(Delaunay &dt,deque<Face> &holes){
	deque<Face> aux = holes;
	for(unsigned int i = 0 ; i < holes.size(); ++i){
		Face c = holes[i];
		if( c->info().m_label != HOLE || !dt.is_face(c->vertex(0), c->vertex(1), c->vertex(2)) || c->info().m_hole_visited ){
			c->info().m_hole_state = false;
			if( c->info().m_hole_visited )
				c->info().m_hole_state = true;
			holes.erase(holes.begin() + i);
			--i;
		}		
		else{
			c->info().m_hole_visited = true;
			c->info().m_hole_state = true;
		}
	}
}

void separate_hole(Delaunay &dt, deque<Face> &holes,vector<Face> &modified){
	for( unsigned int i = 0; i < holes.size(); ++i ){
		Face c = holes[i];
		c->info().m_hole_visited = false; // to return original state

		if( c->info().m_label != HOLE){ // if the cell is not a hole
			holes.erase(holes.begin() + i);
			--i;
			continue;
		}
		vector<Vertex> neighbors;
		for( unsigned int j = 0; j < 3; ++j){
			neighbors.push_back(c->vertex(j) );			
		}

		c->info().m_label = c->info().m_original_label;
		if( c->info().m_original_label == HOLE )		continue;// if the new label is hole

		if( produce_singularities(dt,neighbors) )
			c->info().m_label = HOLE;
		else{
			modified.push_back(c);
			c->info().m_hole_state = false;

			holes.erase( holes.begin() + i );
			--i;
			for( unsigned int j = 0; j < 3; ++j)
				update_vertex_info( dt,neighbors[j] );
		}
	}
}

void solve_holes(Delaunay &dt,deque<Vertex> &singularVertices,deque<Face> &holes,vector<Face> &modified){
	bool draw = 0;
	int count = 1;
	vector<Point> pts;
	while(singularVertices.size() != 0 ){
		Vertex v = singularVertices.front();
		singularVertices.pop_front();
		//if( v->point() == Point(150,	110) ) draw = 1;

		v->info().m_state=1;

		vector<pair<int,int>> &sing_regions = v->info().m_singular_regions_and_subgraphs;
		while( sing_regions.size() != 0){
			int label= sing_regions.front().first;
			vector<vector<Face>> faceSets;
			set_faces_around_singular_vertex(dt,v,label,faceSets);
			if(draw){				
				pts.clear();
				pts.push_back(v->point() );
				drawMesh(dt,pts);
				getchar();
			}
			if( faceSets.size() == 1 ) break;

			if( faceSets.size() == 2 ){
				//cout<<"stop 2"<<endl;
				solve_2(dt,faceSets,singularVertices,v,label,draw,holes,modified);
			}
			else{
				//cout<<"stop 3"<<endl;
				solve_3(dt,faceSets,singularVertices,v,label,draw,holes,modified);
			}
			if(draw){
				pts.clear();
				pts.push_back(v->point() );
				drawMesh(dt,pts);
				getchar();
			}
		}
		//v->info().m_state=1;
	}

	choose_holes(dt,holes);
	separate_hole(dt,holes,modified);
	/*cout<<"after separate"<<endl;
	drawMesh(dt);
	getchar();*/
}

int label_for_cell(Delaunay &dt, Face &c){
	vector<int> neighbor;
	vector<double> areas;
	Face vc;
	for(unsigned int j = 0; j < 3;++j){
		vc = c->neighbor( j );
		if( dt.is_infinite(vc) ) continue;
		if( vc->info().m_label != -1 && vc->info().m_label != HOLE){// if neighbor label is differen of the c original label
			Edge f = Edge(vc,j);
			double squaredLength = dt.segment( f ).squared_length();
			vector<int>::iterator ite = std::find(neighbor.begin(), neighbor.end(), vc->info().m_label ); // locate iterator in neighbor
			if( ite != neighbor.end() ){// found
				int position = std::distance( neighbor.begin(), ite );										
				areas[position] += squaredLength;
			}				
			else{// not found
				neighbor.push_back( vc->info().m_label );
				areas.push_back(squaredLength);
			}
		}
	}
	if( areas.size() != 0){
		int poslabel = std::distance( areas.begin() , max_element( areas.begin(),areas.end() ) );
		return neighbor[poslabel];
	}
	return -1;
}

void neighbor_state(Delaunay &dt,deque<Face> &holes,vector<Face> &modified){
	Face c;
	int new_label;
	while(true){
		unsigned int index = rand()%holes.size();
		c = holes[index];
		if( c->info().m_original_label != HOLE){
			new_label = c->info().m_original_label;
			break;
		}
		else{
			int lbl = label_for_cell(dt,c);
			if(lbl == -1) {continue;}
			new_label = lbl;
			break;
		}
	}
	c->info().m_label = new_label;
	modified.push_back(c);

	deque<Vertex> singularities;
	for( unsigned int i = 0; i < 4; ++i){
		Vertex v = c->vertex(i);
		update_vertex_info(dt,v);
		if( v->info().is_singular_2() ){
			singularities.push_back(v);
			v->info().m_state = 2;
		}
	}
	solve_holes(dt,singularities,holes,modified);	
}

void neighbor_state_2(Delaunay &dt,deque<Face> &holes,vector<Face> &modified){
	deque<Vertex> singularities;
	for(unsigned int i = 0; i < holes.size(); ++i){
		Face c = holes[i];
		c->info().m_label = c->info().m_original_label;
		double r = ((double) rand() / (RAND_MAX));
		if( r < 0.4) {
			int lbl = label_for_cell(dt,c);
			if(lbl != -1) {
				c->info().m_label =lbl;
			}
		}
		modified.push_back(c);
		for( unsigned int j = 0; j < 3; ++j){
			Vertex v = c->vertex(j);
			update_vertex_info(dt,v);
			if( v->info().is_singular_2() && v->info().m_state != 2 ){
				singularities.push_back(v);
				v->info().m_state = 2;
			}
		}
	}
	solve_holes(dt,singularities,holes,modified);
}

double area_set(deque<Face> &cells){
	double vol=0;
	for(unsigned int i = 0; i < cells.size(); ++i)
		vol+= CGAL::area( cells[i]->vertex(0)->point(),cells[i]->vertex(1)->point(),cells[i]->vertex(2)->point()) ;
	return vol;	
}

double cost(Delaunay &dt,deque<Face> &holes){
	double c1 = area_set(holes);
	//double c2 = mesh_curvature(dt);
	//return c1 + c2;
	return c1;
	//return holes.size();
}

void initiate_simulated(Delaunay &dt){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		Face c = it;
		c->info().m_original_label = c->info().m_label;
		c->info().m_solution = c->info().m_label;
	}
}

void accept_state(deque<Face> &holes, vector<Face> &modified){	
	for( unsigned int i = 0; i < modified.size(); ++i){
		modified[i]->info().m_solution = modified[i]->info().m_label;
	}
	for( unsigned int i = 0; i < holes.size(); ++i){
		holes[i]->info().m_solution = holes[i]->info().m_label;
	}
}

void not_accept_state(deque<Face> &holes, vector<Face> &modified){
	for( unsigned int i = 0; i < modified.size(); ++i){
		modified[i]->info().m_label = modified[i]->info().m_solution;
	}
	for( unsigned int i = 0; i < holes.size(); ++i){
		holes[i]->info().m_label = holes[i]->info().m_solution;
	}
}

int check_process(Delaunay &dt){
	vector<Face> incorrect;
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		Face c = it;
		if( c->info().m_label != c->info().m_solution ){
			cout<<"label and solution "<<c->info().m_label<<" "<<c->info().m_solution<<" "<<c->info().m_hole_state <<endl;
			incorrect.push_back(c);
		}
	}
	return (int)incorrect.size();
}

deque<Face> get_holes_solution(Delaunay &dt){
	deque<Face> holes;
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		if( it->info().m_solution == HOLE){
			holes.push_back(it);
		}
	}
	return holes;
}

double acceptance_probability(double old_cost,double new_cost, double T ){
	return exp(old_cost - new_cost)/T;
}

double random_sa(){
	return (double) (rand()/ (double) RAND_MAX);
}

void simulated_annealing_whithout_points(Delaunay &dt){
	double count = 0;
	std::stringstream manifold_simulated;
	cout<<"begin"<<endl;
	vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;
	deque<Face> holes;
	vector<Face> modifiedFaces;	

	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);

	clock_t begin = clock();

	initiate_simulated(dt);// equal soluton to m_label
	solve_holes(dt,singularVertices,holes,modifiedFaces);
	accept_state(holes,modifiedFaces);

	/*manifold_simulated<<count<<"manifold_simulated";
	saveMesh(dt,"tweety/manifold_simulated", ( manifold_simulated.str() ).c_str()  );*/

	//cout<<"stop 1 "<<endl;
	double old_cost = cost(dt,holes);
	double T = 1.0;
	double T_min = 0.00001;
    double alpha = 0.9;

	
	while ( T > T_min && holes.size() != 0){
        int i = 1;
        while (i <= 100 && holes.size() != 0){	
			++count;
			manifold_simulated.str("");
			manifold_simulated<<count<<"manifold_simulated";
			

			modifiedFaces.clear();
			deque<Face> new_holes = holes;
			neighbor_state_2(dt,new_holes,modifiedFaces);// it return less holes 
			//if(new_holes.size() == 0) { T_min = 2; break;}
            double new_cost = cost(dt,new_holes);
            double ap = acceptance_probability(old_cost, new_cost, T);

			/*deque<Face> verify = get_holes(dt);
			cout<<"verified holes before "<<verify.size()<<endl;*/

			cout<<"old and new cost: "<<old_cost<<"  "<<new_cost<<endl;
			cout<<"old and new size holes "<<holes.size()<<"  "<<new_holes.size()<<endl;
			
			/*drawMesh(dt,1);
			getchar();*/

            if( ap > random_sa() ){
				cout<<"accepted"<<endl;
				accept_state(new_holes,modifiedFaces);
				holes = new_holes;
                old_cost = new_cost;
				//saveMesh(dt,"tweety/manifold_simulated", ( manifold_simulated.str() ).c_str()  );
			}
			else{
				not_accept_state(new_holes,modifiedFaces);
				
			}
			//int check = check_process(dt);

			/*verify = get_holes(dt);
			cout<<"verified holes "<<verify.size()<<endl;*/
			
			vector<Point> pts;
            i += 1;			

			cout<<"iterations "<<count<<endl;
		}
        T = T*alpha;
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time simulated: "<<elapsed_secs<<endl;

	restart_vertex_state(dt);
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices simulated "<<singularVertices.size()<<endl;
	cout<<"number of vertices simulated "<<dt.number_of_vertices()<<endl;

	/*vector<Point> vtkPoints;
	for(unsigned int i = 0 ; i < singularVertices.size();++i){
		vtkPoints.push_back(singularVertices[i]->point());
		cout<<singularVertices[i]->point()<<endl;
	}
	drawMesh(dt,vtkPoints);
	getchar();*/

	restart_vertex_state(dt);
}
