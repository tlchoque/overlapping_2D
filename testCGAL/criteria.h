

double larger_edge(Delaunay &dt,vector<Face> &faces, Vertex v){
	double max=0;
	double squaredLength;
	for(unsigned int i = 0; i < faces.size(); ++i){
		for(unsigned int j = 0 ; j < 3 ; ++j){
			//if( j == faces[i]->index(v) ) continue;//important to change
			squaredLength = dt.segment(faces[i],j).squared_length();
			if( squaredLength > max)
				max = squaredLength;
		}
	}
	return max;
}

double area_face_set(vector<Face> &faces){
	double A=0;
	for(unsigned int i = 0; i < faces.size(); ++i)
		A += CGAL::area( faces[i]->vertex(0)->point(), faces[i]->vertex(1)->point(), faces[i]->vertex(2)->point() );
	return A;	
}

double criteria_1( Delaunay &dt,vector<Face> &faces, Vertex v){
	/*Vertex prev_v,next_v;double k;
	singular_vertex_curvature(v,faces[0]->info().m_label,faces,k,prev_v,next_v);
	return k;*/
	//return 1/area_face_set(faces);
	return 1/larger_edge(dt,faces,v);
}

double criteria_2( Delaunay &dt,vector<Face> &faces, Vertex v){
	//return 1/larger_edge(dt,faces,v);
	return 1/area_face_set(faces);
	//return larger_edge(dt,faces);
}


void sort_by_criteria( Delaunay &dt,vector<vector<Face>> &faceSets, Vertex v,int &min, int &max, bool &changeable){
	vector<double> ratio(faceSets.size());
	for(unsigned int i = 0; i < faceSets.size(); ++i)
		ratio[i] = criteria_2(dt,faceSets[i],v);
	max=0;min=1;
	if(ratio[min] > ratio[max]){
		max = 1;
		min = 0;
	}
	if( ratio[max]/ratio[min] < 1 )//for chest is different 2.65 3
		changeable = true;
	else
		changeable = false;
}