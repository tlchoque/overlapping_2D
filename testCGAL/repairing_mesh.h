
vector<Face> faces_of_label(Delaunay &dt, Vertex v, int label){
	vector<Face> faces;
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;

	bool isInternal = true;
	int faceLabel =  fc->info().m_label;
	do{
		if( dt.is_infinite(fc) ){++fc; continue;}
		if( fc->info().m_label == label )
			faces.push_back(fc);
		++fc;
	}while (fc!= fc_start);
	return faces;
}

bool is_isolated(vector<Face> &faces,Vertex v, int label){
	bool cen = true;
	for(unsigned int i=0; i < faces.size(); ++i ){
		Face f =faces[i];
		int vindex = f->index(v);
		Face nf = f->neighbor(vindex);
		if( nf->info().m_label == label){
			cen = false;
			break;
		}
	}
	return cen;
}

bool is_erodible(Delaunay &dt, vector<Face> &faces,Vertex v, int label){
	double c,c_prev,c_next;
	Vertex prev_v,next_v,a,b;
	singular_vertex_curvature(v,label,faces,c,prev_v,next_v);

	bool is_isolated = true;
	for(unsigned int i = 0 ; i < faces.size(); ++i){
		Face f =faces[i];
		int vindex = f->index(v);
		Face nf = f->neighbor(vindex);
		if( nf->info().m_label == label){
			is_isolated = false;
			break;
		}
	}
	if(is_isolated)	return true;

	double sum = 0;
	if( !prev_v->info().is_singular_2() ){
		vector<Face> faces = faces_of_label(dt,prev_v,label);
		singular_vertex_curvature(prev_v,label,faces,c_prev,a,b);	
		sum+= c_prev;
	}
	if( !next_v->info().is_singular_2() ){
		vector<Face> faces = faces_of_label(dt,next_v,label);
		singular_vertex_curvature(next_v,label,faces,c_next,a,b);	
		sum+= c_next;
	}

	//cout<<" sum and c: "<<sum<<" "<<c<<endl;
	if( sum*c < 0 )
		return true;
	else
		return false;
}

bool relabel_if_possible_3(Delaunay & dt, vector<Face> &faces, int other_label,Vertex v, deque<Vertex> &singularities){
	if( !check_faces_state(faces) ){
		cout<<"not checked"<<endl;
		return false;
	}

	int old_label = faces[0]->info().m_label; 
	for(unsigned int i = 0; i < faces.size(); ++i ){
		faces[i]->info().m_label = other_label;
		faces[i]->info().m_state = 1;	
	}

	vector<Vertex> neighbors;
	neighbor_vertices(dt,faces,neighbors,v);
	for(unsigned int i = 0; i < neighbors.size(); ++i ){
		Vertex w = neighbors[i];
		update_vertex_info( dt, w );
		if(  w->info().m_state!=2 && w->info().is_singular_2()  ){
			w->info().m_state = 2;
			if( w->info().has_two_regions() )
				singularities.push_back(w);
			else
				singularities.push_front(w);
		}
	}
	update_vertex_info( dt, v );
	return true;	
}

bool erode( Delaunay &dt,vector<vector<Face>> &faceSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw, deque<Face> &holes ){	
	if(faceSets.size() != 2 ) return false;
	int min,max;
	bool change = false;
	//sort_by_criteria_w(dt,faceSets,v,min,max,change);	
	sort_by_criteria(dt, faceSets, v, min, max, change);

	//if( !insert_to_erode(dt,singulsarities,faceSets,v,label,min,max,draw) ){
	//	int newlabel = new_label_for_relabeling(dt,faceSets[max],v,label); 
	//	if( !relabel_if_possible_3(dt,faceSets[max],newlabel,v,singularities ) ){
	//		if( change){
	//			//cout<<"eroded changed"<<endl;
	//			newlabel = new_label_for_relabeling(dt,faceSets[min],v,label); 
	//			if( !relabel_if_possible_3(dt,faceSets[min],newlabel,v,singularities) ){
	//				return false;
	//			}
	//		}
	//		else return false;
	//	}
	//}

	int newlabel = new_label_for_relabeling(dt,faceSets[max],v,label); 
	if( !relabel_if_possible_3(dt,faceSets[max],newlabel,v,singularities) ){	
		if( change){
			newlabel = new_label_for_relabeling(dt,faceSets[min],v,label); 
			//cout<<"new label 2 "<<newlabel<<endl;
			if( !relabel_if_possible_3(dt,faceSets[min],newlabel,v,singularities ) ){						
				if( !insert_to_erode(dt,singularities,faceSets,v,label,min,max,draw) ){
					return false;
				}
			}
		}
		else{
			if( !insert_to_erode(dt,singularities,faceSets,v,label,min,max,draw) ){
				return false;
			}
		}
	}

	return true;
}

int adyacent_face_label(Delaunay &dt, Face &c, int alabel, int blabel){
	vector<int> neighbor;
	vector<double> areas;
	Face vc;
	for(unsigned int j = 0; j < 3;++j){
		vc = c->neighbor( j );
		if( dt.is_infinite(vc) ) continue;
		if( vc->info().m_label == alabel || vc->info().m_label == blabel){// if neighbor label is equl to alabel or blabel
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

bool relabel_path(Delaunay & dt, vector<Face> &faces, int alabel,int blabel,Vertex v, deque<Vertex> &singularities){
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
		for( unsigned int j = 1; j < 3; ++j){			
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

bool dilate_with_checking( Delaunay &dt,vector<vector<Face>> &faceSets,deque<Vertex> &singularities,Vertex &v,int label,bool &draw){	
	if(faceSets.size() != 2 ) return false;
	vector<Face> path;
	path_to_join(dt,faceSets,v,path);

	double s1 = criteria_1(dt,faceSets[0],v);
	double s2 = criteria_1(dt,faceSets[1],v);
	double s3 = criteria_1(dt,path,v);

	double s = s1;
	if( s2 > s1) s = s2;
	if( s/s3 > 2 ) return false;
	return relabel_path(dt,path,faceSets[0][0]->info().m_label,faceSets[1][0]->info().m_label,v,singularities) ; 
}

bool repairing_2(Delaunay &dt, vector< vector<Face> > &faceSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Face> &holes ) {
	if( ( !is_erodible(dt,faceSets[0],v,label) || !is_erodible(dt,faceSets[1],v,label) )  
		//&& ( !is_isolated(faceSets[0],v,label) && !is_isolated(faceSets[1],v,label) )
		&& !v->info().is_corner() 
		){	
		//cout<<"dilate decision"<<endl;
		if( v->info().has_two_regions() ){
			int opposite_label = v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
			vector<vector<Face>> faceSets_2;
			set_faces_around_singular_vertex(dt,v,opposite_label,faceSets_2);
			if( erode(dt, faceSets_2,singularities,v,opposite_label,draw,holes) ){
				if(draw){				
					drawMesh(dt);
					getchar();
				}
			}	
			else{
				dilate_with_checking(dt,faceSets,singularities,v,label,draw);
			}
		}
		else{
			//cout<<"dilate stop 3"<<endl;
			dilate_with_checking(dt,faceSets,singularities,v,label,draw);
		}
	}
	set_faces_around_singular_vertex(dt,v,label,faceSets);

	if( faceSets.size() == 1 )
		return true;
	else{
		if( erode(dt, faceSets,singularities,v,label,draw,holes) ){
			return true;
		}
		else{
			return dilate_with_checking(dt,faceSets,singularities,v,label,draw);
		}
		//return erode_with_checking(dt, faceSets,singularities,v,label,draw,holes);	
	}
}

bool dilate_reparing_3(Delaunay &dt, vector< vector<Face> > &faceSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Face> &holes ){
	int cen = true ;
	while( faceSets.size() > 1 ){
		cen = true;
		for( unsigned int i = 0; i < faceSets.size() - 1; ++i){
			vector< vector<Face> > sets;
			sets.push_back( faceSets[i] );
			sets.push_back( faceSets[i + 1] );
			if( dilate_with_checking(dt,sets,singularities,v,label,draw) ){// dilated
				cen = false;
				set_faces_around_singular_vertex(dt,v,label,faceSets);
				break;
			}
		}
		if (cen) return false;
	}
	return true;
}

bool set_is_erodible(Delaunay &dt, vector< vector<Face> > &faceSets,Vertex v,int label){
	double count = 0;
	for(unsigned int i = 0; i < faceSets.size(); ++i){
		if( !is_erodible(dt,faceSets[i],v,label) )
			++count;
	}
	if( count > 0) return false;
	return true;
}

bool repairing_3(Delaunay &dt, vector< vector<Face> > &faceSets, deque<Vertex> &singularities, Vertex v, int label, bool draw, deque<Face> &holes ) {
	if(draw){		
		cout<<"size: "<<faceSets.size()<<endl;
		drawMesh(dt);
		getchar();
	}

	bool result = false;
	if( !set_is_erodible(dt,faceSets,v,label) || v->info().has_more_than_two_regions() )
		result = dilate_reparing_3(dt,faceSets,singularities,v,label,draw,holes);

	if ( result ) return true;

	vector<double> ratio(faceSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;

	for(unsigned int i = 0; i < ratio.size(); ++i){
		ratio[i] = criteria_1(dt,faceSets[i],v);
	}

	int v1=0,v2=1;
	if(ratio[v2] < ratio[v1]){
		v1 = 1;
		v2 = 0;
	}
	for(unsigned int l = 2; l < ratio.size(); ++l){// choose the 2 greater volumes v1 < v2
		if(ratio[l] < ratio[v1]){
			int aux = v1;
			v1 = l;
			v2 = aux;
		}
		else if(ratio[l] < ratio[v2])   v2 = l;
	}
	for(unsigned int l = 0; l < faceSets.size(); ++l){		
		if( l != v1 && l!=v2){
			int newlabel = new_label_for_relabeling(dt,faceSets[l],v,label);
			if(! relabel_if_possible_3(dt,faceSets[l],newlabel,v,singularities ) ){
				continue;
			}
			faceSets.erase(faceSets.begin()+l);
			--l;
			--v1;--v2;
		}
	}

	set_faces_around_singular_vertex(dt,v,label,faceSets);
	if(draw){		
		//cout<<"size: "<<faceSets.size()<<endl;
		drawMesh(dt);
		getchar();
		// = 0;
	}
	return repairing_2(dt,faceSets,singularities,v,label,draw,holes);
}

void repairing(Delaunay &dt){
	deque<Face> holes;
	vector<Vertex> boundaries;
	deque<Vertex> singularities;

	clock_t begin0 = clock();

	boundary_vertices(dt,boundaries);
	singular_vertices(dt,boundaries,singularities);

	clock_t end0 = clock();
	double elapsed_secs0 = double(end0 - begin0) / CLOCKS_PER_SEC;
	cout<<"time get singularities: "<<elapsed_secs0<<endl;

	cout<<"number of non manifold vertices before "<<singularities.size()<<endl;
	cout<<"number of vertices "<<dt.number_of_vertices()<<endl;
	
	//getchar();
	vector<Point> pts;
	for(unsigned int i = 0 ; i < singularities.size();++i){
		pts.push_back( singularities[i]->point() );
	}
	/*savePoints2(pts,"titicaca/original/singularities7.vtk",7);
	savePoints2(pts,"titicaca/original/singularities6.vtk",5);
	savePoints2(pts,"titicaca/original/singularities5.vtk",5);
	savePoints2(pts,"titicaca/original/singularities4.vtk",4);
	savePoints2(pts,"titicaca/original/singularities3.vtk",3); */

	/*cout<<"finish save singularitites"<<endl;
	drawMesh(dt,pts);
	getchar();*/
	clock_t begin = clock();


	bool draw = 0;	
	int count = 0;
	while(singularities.size() != 0 ){	
		++count;		
		//if( count == 174 ) draw = 1;// 3570 - 3590    ------> 3587 bad s for singularity
		//if( count == 215 ) draw = 1;// 3570 - 3590    ------> 3587 bad s for singularity

		/*draw = 0;
		if( count == 150 ) draw = 1;
		if( count == 155 ) draw = 1;
		if( count == 160 ) draw = 1;
		if( count == 165 ) draw = 1;
		if( count == 170 ) draw = 1;
		if( count == 175 ) draw = 1;*/

		//cout<<"count: "<<count<<endl;
		Vertex v = singularities.front();

		//if( v->point() == Point(83.225,201.98) ) draw = 1;
		//if( v->point() == Point(64,278) ) draw = 1;

		singularities.pop_front();
		v->info().m_state = 1;				
		vector<pair<int,int>> &sing_regions = v->info().m_singular_regions_and_subgraphs;
		while( sing_regions.size() != 0){
			//cout<<"size sing: "<<sing_regions.size()<<" point: " <<v->point()<<endl;
			int label = sing_regions.front().first;
			//cout<<"label: "<<label<<endl;
			vector<vector<Face>> faceSets;
			set_faces_around_singular_vertex(dt,v,label,faceSets);

			if(draw){		
				pts.clear();
				pts.push_back(v->point() );
				drawMesh(dt,pts);
				getchar();
			}

			if( faceSets.size() == 2 ){
				//cout<<"case 2"<<endl;
				if( !repairing_2(dt,faceSets,singularities,v,label,draw,holes) ){
					//cout<<"stop 2"<<endl;
					v->info().m_state=2;
					break;
				}
			}
			else{
				//cout<<"case 3"<<endl;
				if( !repairing_3(dt,faceSets,singularities,v,label,draw,holes) ){
					//cout<<"stop 3"<<endl;
					v->info().m_state=2;
					break;
				}
			}
			if(draw){
				pts.clear();
				pts.push_back(v->point() );
				drawMesh(dt,pts);
				getchar();
			}
		}
	}

	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time repairing: "<<elapsed_secs<<endl;

	restart_vertex_state(dt);
	boundary_vertices(dt,boundaries);
	singular_vertices(dt,boundaries,singularities);
	cout<<"number of non manifold vertices after "<<singularities.size()<<endl;
	cout<<"number of vertices "<<dt.number_of_vertices()<<endl;

	//vector<Point> pts;
	pts.clear();
	for(unsigned int i = 0 ; i < singularities.size();++i){
		pts.push_back( singularities[i]->point() );
		cout<<"sing "<<singularities[i]->point() <<endl;
	}
	/*drawMesh(dt,pts);
	getchar();*/
	//savePoints2(pts,"tweety/manifold/remaining3.vtk",3); 	

	restart_vertex_state(dt);
}