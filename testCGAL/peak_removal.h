
void curvature_polygonal_curve(Delaunay &dt){
	vector<Vertex> boundaryVertices;
	boundary_vertices(dt,boundaryVertices);
	double totalcurvature=0;
	for(unsigned int i = 0 ; i < boundaryVertices.size(); ++i ){
		double curvature=0;
		if(vertex_curvature(dt, boundaryVertices[i],1 , curvature))//label = 1
			totalcurvature+=curvature;
	}
	cout<<totalcurvature;
	restart_vertex_state(dt);
}

bool remove_pick(Delaunay & dt, Vertex &v, double minAngle){
	int label0 = v->info().m_regions_around[0];
	int label1 = v->info().m_regions_around[1];
	vector<Face> Faces0;
	vector<Face> Faces1;
	list<Face> incident;

	Vertex v1,v2;
	double total = 0;

	vector<Face> faces;
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		Face f = fc;
		if( fc->info().m_label == label0 ){
			Faces0.push_back( f );			
			int indx = f->index(v);	
			v1 = f->vertex( f->ccw(indx) );
			v2 = f->vertex( f->cw(indx) );
			total += getAngle( v2 , v, v1 );
		}else
			Faces1.push_back(f);
		++fc;
	}while (fc!= fc_start);
	
	//cout<<total<<endl;
	if( total > 180){
		total = 360 - total;
		if( total < minAngle)
			return relabel_if_possible(dt,Faces1,label0,v);
	}
	else{
		if( total < minAngle)
			return relabel_if_possible(dt,Faces0,label1,v) ;
	}
	return false;
}

void remove_pick2(Delaunay & dt, Vertex &v, double minAngle){
	map<int,int> label_to_index;

	cout<<"stop 0 "<<endl;
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		label_to_index[ v->info().m_regions_around[i] ] = i;
	}

	cout<<"stop 1 "<<endl;
	vector<vector<Face>> FacesSet( v->info().m_regions_around.size() );	
	vector<double> regionAngle( v->info().m_regions_around.size() ,0 );	

	Vertex v1,v2,v3;
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		cout<<"stop 1 1/2 "<<endl;
		if(dt.is_infinite(fc)) {++fc;continue;}
		Face f = fc;
		FacesSet[ label_to_index[ f->info().m_label ]  ].push_back(f);
		int indx = f->index(v);	
		v1 = f->vertex( f->ccw(indx) );
		v2 = f->vertex( f->cw(indx) );
		regionAngle[ label_to_index[ f->info().m_label ] ] += getAngle( v2 , v, v1 );
	}while (fc!= fc_start);

	cout<<"stop 2 "<<endl;

	for(unsigned int i = 0; i < regionAngle.size(); ++i){
		cout<<"angle: "<<i<<" "<<regionAngle[i]<<endl;
		if( regionAngle[i] < minAngle ){
			int newlabel = new_label_for_relabeling(dt,FacesSet[i],v,v->info().m_regions_around[i]);
			relabel_if_possible(dt,FacesSet[i],v->info().m_regions_around[i],v);
		}
	}
}

void remove_boundary_peaks(Delaunay &dt, double minAngle){
	vector<Vertex> boundaryVertices;
	boundary_vertices(dt,boundaryVertices);
	cout<<boundaryVertices.size()<<endl;
	for(unsigned int i = 0 ; i< boundaryVertices.size(); ++i){
		Vertex v = boundaryVertices[i];
		if( v->info().is_singular() || v->info().is_corner() ) continue;
		if( v->info().has_two_regions() ){
			cout<<"2 regions"<<endl;
			remove_pick( dt,v,minAngle );
		}
		else{
			if( v->info().has_more_than_two_regions() ){
				cout<<"has_more_than_two_regions"<<endl;
				remove_pick2( dt,v,minAngle );
			}
		}
	}
	restart_vertex_state(dt);
	cout<<"ends remove boundary picks"<<endl;
}
