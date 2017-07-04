
void mark_boundary(Delaunay &dt){
	Vertex_circulator vc_start = dt.incident_vertices(dt.infinite_vertex());
	Vertex_circulator vc = vc_start;
	do{
		if( !dt.is_infinite(vc) )	vc->info().m_type = 0;
		++vc;
	}while (vc!= vc_start);
}

void boundary_vertices(Delaunay &dt,vector<Vertex> & boundaryVertices){
	boundaryVertices.clear();
	for(Finite_edges_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		if( ei->first->info().m_label == ei->first->neighbor(ei->second)->info().m_label )		continue;
		int i = ei->second;
		Vertex v1 = ei->first->vertex(ei->first->cw(i));
		Vertex v2 = ei->first->vertex(ei->first->ccw(i));

		if (v1->info().m_state == 0)	{
			//update_regions_around_vertex( dt, v1 );
			update_vertex_info(dt,v1);
			v1->info().m_state=1;
			boundaryVertices.push_back(v1);
		}

		if (v2->info().m_state == 0)	{
			//update_regions_around_vertex( dt, v2 );
			update_vertex_info(dt,v2);
			v2->info().m_state=1;
			boundaryVertices.push_back(v2);
		}
	}
}

void boundary_vertices2(Delaunay &dt,vector<Vertex> & boundaryVertices){
	boundaryVertices.clear();
	for(Finite_edges_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		if( ei->first->info().m_label == ei->first->neighbor(ei->second)->info().m_label )		continue;
		if( dt.is_infinite(ei->first ) || dt.is_infinite( ei->first->neighbor(ei->second) ) )		continue;

		int i = ei->second;
		Vertex v1 = ei->first->vertex(ei->first->cw(i));
		Vertex v2 = ei->first->vertex(ei->first->ccw(i));

		if (v1->info().m_state == 0)	{
			update_regions_around_vertex( dt, v1 );
			//update_vertex_info( dt, v1 );
			v1->info().m_state=1;
			boundaryVertices.push_back(v1);
		}

		if (v2->info().m_state == 0)	{
			update_regions_around_vertex( dt, v2 );
			//update_vertex_info( dt, v2 );
			v2->info().m_state=1;
			boundaryVertices.push_back(v2);
		}
	}
}

void singular_vertices(Delaunay &dt,vector<Vertex> & boundaryVertices, deque<Vertex> & singularVertices){
	for(unsigned int i = 0; i < boundaryVertices.size(); ++i){
		Vertex v = boundaryVertices[i];
		if( v->info().is_singular_2() ){
		//if( v->info().is_singular() ){
			if( v->info().has_two_regions() )
				singularVertices.push_front(v);
			else
				singularVertices.push_back(v);
			v->info().m_state = 2;// singular vertex
		}
	}
}

//for new relabeling

bool is_label_around(Delaunay &dt, Vertex &v, int label){
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if( !dt.is_infinite(fc) && fc->info().m_label != -1 ){
			if( fc->info().m_label == label )
				return true;
		}		
		++fc;
	}while (fc!= fc_start);
	return false;
}

bool add_triangle(Delaunay &dt, Face &f, int label){
	if( f->info().m_label == label) return false;
	Face neighborFace;
	vector<int> facets;
	for(unsigned int i = 0; i < 3; ++i){
		neighborFace = f->neighbor(i);
		if(neighborFace->info().m_label == label)
			facets.push_back(i);
	}

	if( facets.size() == 0 ){
		for(unsigned int i = 0; i < 3; ++i){
			if( is_label_around(dt,f->vertex(i),label ) )
				return false;
		}
		return true;
	}		
	else if( facets.size() == 1 )
		return !is_label_around(dt,f->vertex(facets[0]),label);
	else	return true;
}

bool is_singular_after_relabeling(Delaunay &dt, Vertex &v){
	vector<int> labels;
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
	if( isInternal ) return false;
	// end advance

	int labelcirculator = -2;
	do{
		int label =  fc->info().m_label;
		if( label == labelcirculator ){ fc++;continue;}
		labelcirculator = label;
		if( !dt.is_infinite(fc) && label != -1 ){			
			if( std::find(labels.begin(), labels.end(), label ) != labels.end() )//found in labels
				return true;
			else
				labels.push_back( label );
		}		
		++fc;
	}while (fc!= fc_start);
	return false;
}

void neighbor_vertices(Delaunay &dt, vector<Face> &faces,vector<Vertex> &vertices,Vertex &v) {
	Vertex v1;
	for(unsigned int i = 0; i < faces.size(); ++i){
		int idx = faces[i]->index(v);
		for(unsigned int j = 1; j < 3; ++j){
			v1 = faces[i]->vertex( (idx+j)%3 );
			if( v1->info().m_auxiliar == false ){
				v1->info().m_auxiliar = true;
				vertices.push_back(v1);
			}
		}
	}
	for(unsigned int i = 0; i < vertices.size(); ++i)
		vertices[i]->info().m_auxiliar = false;
}

bool check_neighbor_vertices( Delaunay &dt, vector<Vertex> &vertices ){
	for(unsigned int i = 0; i < vertices.size(); ++i ){
		if( vertices[i]->info().has_more_than_two_regions() || vertices[i]->info().m_boundary){
			if( is_singular_after_relabeling(dt, vertices[i] ) )
				return false;
		}
	}
	return true;
}

bool relabel_if_possible(Delaunay & dt, vector<Face> &faces, int label,Vertex v){
	deque<Face> cpyFaces(faces.begin(),faces.end());
	size_t cen = cpyFaces.size();
	Face f;
	int old = faces[0]->info().m_label;
	while(cpyFaces.size() != 0 && cen !=0){		
		f = cpyFaces.back();
		cpyFaces.pop_back();
		if( add_triangle(dt,f,label) ){
			f->info().m_label = label;
			f->info().m_state = 1;
			cen = cpyFaces.size();	
		}
		else{
			cpyFaces.push_front(f);
			--cen;
		}
	}

	vector<Vertex> neighbors;
	neighbor_vertices(dt,faces,neighbors,v);

	if( cpyFaces.size() != 0  ){
		cout<<"this case "<<endl;
		for(unsigned int i = 0; i < faces.size(); ++i){
			faces[i]->info().m_label = old;
			faces[i]->info().m_state = 0;
		}
		return false;
	}
	else{		
		if( !check_neighbor_vertices(dt,neighbors) ){
			cout<<"generate singularities "<<endl;
			for(unsigned int i = 0; i < faces.size(); ++i){
				faces[i]->info().m_label = old;
				faces[i]->info().m_state = 0;
			}
			return false;
		}
		else{
			for(unsigned int i = 0; i < neighbors.size(); ++i ){
				update_regions_around_vertex( dt,neighbors[i] );
			}
			update_regions_around_vertex( dt, v );
			return true;
		}		
	}	
}

//for fix singularities

void set_faces_around_singular_vertex(Delaunay &dt, Vertex &v, int label, vector<vector<Face>> &faceSets){
	faceSets.clear();
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
	if( isInternal ) {return;}
	// end advance


	vector<Face> faces;
	int previous = -2;

	
	do{
		faceLabel = fc->info().m_label;
		if( faceLabel == previous ){
			if( faceLabel == label )
				faces.push_back(fc);
		}
		else{
			if( faceLabel == label )
				faces.push_back(fc);
			else{
				if( previous == label ){
					faceSets.push_back(faces);
					faces.clear();
				}
			}
			previous = faceLabel;
		}
		/*cout<<"face label "<<fc->info().m_label<<endl;
		cout<<"size faces "<<faces.size()<<endl;
		cout<<"size faceSets "<<faceSets.size()<<endl;
		getchar();*/
		++fc;
		
	}while (fc!= fc_start);

	if( faces.size() > 0 && previous == label)	faceSets.push_back(faces);	
}

bool check_faces_state(vector<Face> &faces){
	for(unsigned int i = 0; i < faces.size(); ++i){
		if( faces[i]->info().m_state )
			return false;
	}
	return true;
}

void relabel_Faces_1(Delaunay &dt, deque<Vertex> &singularVertices, vector<Face> &faces,Vertex &v, int label){
	int newlabel = v->info().m_regions_around[0];
	if( newlabel == label )	newlabel = v->info().m_regions_around[1];
	if( newlabel == -1 ) return;

	for(unsigned int i = 0; i < faces.size();++i){
		Face c = faces[i];
		c->info().m_label = newlabel;
		c->info().m_state = 1;	
		int index = c->index(v);
		for(unsigned int j = 1; j < 3;++j){
			Vertex vc = c->vertex( (index+j)%3 );
			update_regions_around_vertex( dt, vc );
			if( vc->info().m_state != 2  && vc->info().is_singular() ){
				singularVertices.push_front( vc );
				vc->info().m_state=2;
			}
		}
	}
	update_regions_around_vertex( dt, v );
}

void relabel_Faces_2(Delaunay &dt, deque<Vertex> &singularVertices, vector<Face> &faces,Vertex &v, int label){
	for(unsigned int i = 0; i < faces.size();++i){
		Face c = faces[i];
		double a2 = 0;int lbl = 0;int idx=0; // there will always be an edge 
		for(unsigned int j = 1; j < 3;++j){
			idx = ( c->index(v) + j )%3;
			Edge e = Edge(c, idx);

			double squaredLengthEdge = dt.segment( e ).squared_length();
			if( squaredLengthEdge > a2 && c->neighbor( idx )->info().m_label != label && c->neighbor( idx )->info().m_label != -1) {
				a2 = squaredLengthEdge;
				lbl = c->neighbor( idx )->info().m_label;
			}
		}// at the end, we got the label & index	
		c->info().m_label = lbl;
		c->info().m_state = 1;	

		int index = c->index(v);
		for(unsigned int j = 1; j < 3;++j){
			Vertex vc = c->vertex( ( c->index(v) + j )%3 );
			update_regions_around_vertex( dt, vc );
			if( vc->info().m_state != 2  && vc->info().is_singular() ){
				singularVertices.push_front( vc );
				vc->info().m_state=2;
			}
		}			
	}
	update_regions_around_vertex( dt, v );
}

void relabeling_if_checked( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Face>> &faceSets,Vertex &v,int label ){
	vector<double> ratio(faceSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxMass=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,faceSets[l],v) , area_face_set(faceSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxMass)		maxMass = data[l].second;
	}
	int max=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxMass; 
		if( ratio[max] < ratio[l])  max = l;
	}

	if(v->info().has_more_than_two_regions() ){
		for(unsigned int l = 0; l < faceSets.size(); ++l){
			if( l != max && check_faces_state(faceSets[l]) ){
				relabel_Faces_2(dt,singularVertices,faceSets[l],v,label);//update all affected vertices
				faceSets.erase(faceSets.begin()+l);
				--l;
				--max;
			}
		}
	}
	else{
		for(unsigned int l = 0; l < faceSets.size(); ++l){
			if( l != max && check_faces_state(faceSets[l]) ){
				relabel_Faces_1(dt,singularVertices,faceSets[l],v,label);//update all affected vertices
				faceSets.erase(faceSets.begin()+l);
				--l;
				--max;
			}
		}
	}
}

void relabeling( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Face>> &faceSets,Vertex &v,int label ){
	vector<double> ratio(faceSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,faceSets[l],v) , area_face_set(faceSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	int max=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
		if( ratio[max] < ratio[l])  max = l;
	}

	for(unsigned int l = 0; l < faceSets.size(); ++l){
		//if( l != max ){
			relabel_Faces_2(dt,singularVertices,faceSets[l],v,label);//update all affected vertices
			//faceSets.erase(faceSets.begin()+l);
		//}
	}
}

void relabeling_decrease_set_to_two( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Face>> &faceSets,Vertex &v,int label ){
	vector<double> ratio(faceSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,faceSets[l],v) , area_face_set(faceSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
	}
	int v1=0,v2=1;
	if(ratio[v2] > ratio[v1]){
		v1 = 1;
		v2 = 0;
	}
	for(unsigned int l = 2; l < ratio.size(); ++l){// choose the 2 greater volumes v1 > v2
		if(ratio[l] > ratio[v1]){
			int aux = v1;
			v1 = l;
			v2 = aux;
		}
		else if(ratio[l] > ratio[v2])   v2 = l;
	}					
	for(unsigned int l = 0; l < faceSets.size(); ++l){
		if( l != v1 && l!=v2){
			relabel_Faces_2(dt,singularVertices,faceSets[l],v,label);
			faceSets.erase(faceSets.begin()+l);
		}
	}
}

void restart_vertex_state(Delaunay &dt){
	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		vi->info().m_state = 0;
	} 
}

int new_label_2(Delaunay &dt, vector<Face> &faces,Vertex &v,int label){
	Face vc,c;
	vector<int> neighbor;
	vector<double> distances;
	for(unsigned int i = 0; i < faces.size(); ++i ){
		vc = faces[i];
		int index = vc->index(v);
		for(unsigned int j = 1; j < 3;++j){
			int index2 = (index+j)%3; 
			c = vc->neighbor( index2 );
			if(dt.is_infinite(c)) continue;
			if( c->info().m_label != label ){
				Edge f = Edge(vc,index2);
				double squaredSegment = dt.segment( f ).squared_length();
				vector<int>::iterator ite = std::find(neighbor.begin(), neighbor.end(), c->info().m_label );
				if( ite != neighbor.end() ){// found
					int position = std::distance( neighbor.begin(), ite );										
					distances[position] += squaredSegment;
				}				
				else{// not found
					neighbor.push_back( c->info().m_label );
					distances.push_back(squaredSegment);
				}
			}
		}
	}
	int poslabel = std::distance( distances.begin() , max_element( distances.begin(),distances.end() ) );
	return neighbor[poslabel];
}

int new_label_for_relabeling(Delaunay &dt, vector<Face> &faces,Vertex &v,int label){
	if(v->info().has_more_than_two_regions() )
		return new_label_2(dt,faces,v,label);
	else
		return v->info().m_regions_around[ v->info().m_regions_around[0] == label ];
}

//relabeling singular vertex, 2 set of triangles, 2 regions around v
bool relabeling_with_curvature(Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Face>> &faceSets,Vertex &v,int &label,bool draw ){
	double c1=0,c2=0;
	Vertex prev_v_c1,next_v_c1;
	Vertex prev_v_c2,next_v_c2;
	singular_vertex_curvature(v,label,faceSets[0],c1,prev_v_c1,next_v_c1);
	singular_vertex_curvature(v,label,faceSets[1],c2,prev_v_c2,next_v_c2);
	vector<Face> selected_set;
	double prev_curvature,next_curvature;
	bool cen=false;

	//cout<<"curvatures: " <<c1 << " * "<<c2<<endl;
	if( c1 > c2 ){
		selected_set = faceSets[0];
		if( ! prev_v_c1->info().is_singular() && ! next_v_c1->info().is_singular()){
			vertex_curvature(dt,prev_v_c1,label,prev_curvature);
			vertex_curvature(dt,next_v_c1,label,next_curvature);
		}
		else cen = true;
	}
	else{
		selected_set = faceSets[1];
		if( ! prev_v_c2->info().is_singular() && ! next_v_c2->info().is_singular()){
			vertex_curvature(dt,prev_v_c2,label,prev_curvature);
			vertex_curvature(dt,next_v_c2,label,next_curvature);
		}
		else cen = true;
	}

	bool relabel = true;
	int newlabel = new_label_for_relabeling(dt,selected_set,v,label);
	if( ( cen == 0 && prev_curvature < 0 && next_curvature < 0 ) || cen == 1 ){
		if( !relabel_if_possible(dt,selected_set,newlabel,v ) ){
			relabel = false;
		}
	}
	else relabel = false;

	if( relabel ) return true;

	//if(v->info().is_corner()) return false;
	//it happens that the opposite lable has 1 set


	vector<vector<Face>> aux =faceSets;
	faceSets.clear();
	int aux_label = label;
	label = newlabel;
	set_faces_around_singular_vertex(dt,v,newlabel,faceSets);


	cout<<"face sets size: "<<faceSets.size() <<" "<<faceSets[0].size()<<endl;
	if(faceSets.size() == 1){
		cout<<"label "<<aux_label<<endl;
		if(!relabel_if_possible(dt,faceSets[0],aux_label,v ) ){
			faceSets = aux;
			return false;
		}
		else return true;
	}

	singular_vertex_curvature(v,label,faceSets[0],c1,prev_v_c1,next_v_c1);
	singular_vertex_curvature(v,label,faceSets[1],c2,prev_v_c2,next_v_c2);

	//cout<<"2 curvatures: " <<c1 << " * "<<c2<<endl;
	if( c1 > c2 )
		selected_set = faceSets[0];
	else
		selected_set = faceSets[1];

	/*if( !insert_to_erode(dt,singularVertices,faceSets,v,label,draw )  ){
		if( !relabel_if_possible(dt,selected_set,aux_label,v )  ){
			return false;
		}
	}*/
	return true;
}

bool relabeling_without_making_singularities( Delaunay &dt, deque<Vertex> &singularVertices,vector<vector<Face>> &faceSets,Vertex &v,int &label,bool draw ){
	/*if( v->info().has_two_regions() && !v->info().is_corner() )	
		return  relabeling_with_curvature( dt, singularVertices,faceSets,v,label,draw );*/

	vector<double> ratio(faceSets.size());
	vector<pair<double,double>> data;
	double maxEdge=0,maxVolume=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		data.push_back( make_pair( larger_edge(dt,faceSets[l],v) , area_face_set(faceSets[l]) )  );
		if(data[l].first > maxEdge)		maxEdge = data[l].first ;
		if(data[l].second > maxVolume)		maxVolume = data[l].second;
	}
	int max=0;
	for(unsigned int l = 0; l < ratio.size(); ++l){
		ratio[l] = 0.7*data[l].first/maxEdge + 0.3*data[l].second/maxVolume; 
		if( ratio[max] < ratio[l])  max = l;
	}

	bool cen = true;
	for(unsigned int l = 0; l < faceSets.size(); ++l){
		//if( l != max && check_cells_state(faceSets[l]) ){
		if( l != max ){
			int newlabel = new_label_for_relabeling(dt,faceSets[l],v,label);
			if( !relabel_if_possible(dt,faceSets[l],newlabel,v ) ){
				cen = false;
			}
			else{
				faceSets.erase(faceSets.begin()+l);
				--l;
				--max;
			}
		}
	}
	return cen;
}

void fix_singularities2(Delaunay &dt){
	mark_boundary(dt);
	vector<Vertex> boundaryVertices;
	deque<Vertex> singularVertices;		
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices before "<<singularVertices.size()<<endl;
	//vector<Point> points;
	//getchar();
	vector<Point> vtkPoints;
	vtkPoints.push_back(Point(255,163));
	/*for(unsigned int i = 0 ; i < singularVertices.size();++i)
		vtkPoints.push_back(singularVertices[i]->point());*/

	
	/*drawMesh(dt,vtkPoints);
	getchar();*/
	
	//savePoints(vtkPoints,"tauro/original/singularities.vtk");
	/*std::stringstream ss,ss1;
	ss<<file<<"/original/singularities_circles_9.vtk";
	savePoints2(vtkPoints,ss.str().c_str(), 9);
	ss1<<file<<"/original/singularities_circles_11.vtk";
	savePoints2(vtkPoints,ss1.str().c_str(), 11);
	cout<<"finish "<<endl;*/
	//getchar();

	bool draw = 0;
	int counting = 0;
	while(singularVertices.size() != 0 ){
		++counting;
		//cout<<"counting: "<<counting<<endl;
		Vertex v = singularVertices.front();
		singularVertices.pop_front();
		//cout<<"size singularities: "<<singularVertices.size()<<endl;

		while(v->info().m_singular_regions_around.size() != 0){
			int label = v->info().m_singular_regions_around.back();
			v->info().m_singular_regions_around.pop_back();
			
			vector<vector<Face>> faceSets;
			set_faces_around_singular_vertex(dt,v,label,faceSets);
			if(draw){				
				vtkPoints.clear();
				vtkPoints.push_back(v->point() );
				drawMesh(dt, vtkPoints);
				getchar();
				//aroundVertex(dt,v,v,label);

				//getchar();
			}

			if(faceSets.size() > 2)
				relabeling_decrease_set_to_two(dt,singularVertices,faceSets,v,label );

			if( !relabeling_without_making_singularities(dt,singularVertices,faceSets,v,label,draw) ){
				//cout<<"show image"<<endl;
				if(draw){						
					vtkPoints.clear();
					vtkPoints.push_back(v->point() );
					drawMesh(dt, vtkPoints);
					//getchar();
				}
				cout<<"insertion "<<endl;
				if( !insert_point(dt,singularVertices,faceSets,v,label,draw ) ){
					cout<<"break"<<endl;
					break;
				}	
			}
			if(draw){
				vtkPoints.clear();
				vtkPoints.push_back(v->point() );
				drawMesh(dt, vtkPoints);
				getchar();
			}
		}
		v->info().m_state=1;
	}

	restart_vertex_state(dt);
	boundaryVertices.clear();
	singularVertices.clear();
	boundary_vertices(dt,boundaryVertices);
	singular_vertices(dt,boundaryVertices,singularVertices);
	cout<<"number of non manifold vertices, after "<<singularVertices.size()<<endl;

	vtkPoints.clear();
	for(unsigned int i = 0 ; i < singularVertices.size();++i)
		vtkPoints.push_back(singularVertices[i]->point());

	/*drawMesh(dt,vtkPoints);
	getchar();*/

	//for smoothing
	restart_vertex_state(dt);
}