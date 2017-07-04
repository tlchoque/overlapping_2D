
bool IsNumber(double x) {
    return (x == x); 
}

void set_adyacent_vertices(Delaunay & dt){
	vector<Vertex> boundaryVertices;
	boundary_vertices2(dt,boundaryVertices);
	cout<<"boundry vertices "<<boundaryVertices.size()<<endl;
	for(unsigned int i = 0; i < boundaryVertices.size(); ++i){
		Vertex v = boundaryVertices[i];
		vector<int> &regions = v->info().m_regions_around;
		vector<vector<Vertex>> &adyacent = v->info().m_adyacent_vertices;
		adyacent.clear();
		for(unsigned int i = 0; i< regions.size(); ++i){
			int label = regions[i];
			vector<Vertex> pair; Vertex a,b;
			//check what we got
			get_adyacents(dt,v,label,a,b);
			pair.push_back(a);
			pair.push_back(b);
			adyacent.push_back(pair);
		}
	}	
	//cout<<boundaryVertices[0]->point()<<" * "<<boundaryVertices[0]->info().m_adyacent_vertices.size()<<" + "<<boundaryVertices[0]->info().m_regions_around.size()<<endl;
	restart_vertex_state(dt);
}

void mesh_boundary_histogram(Delaunay & dt, const char* filename){
	vector<Vertex> boundaryVertices;
	boundary_vertices2(dt,boundaryVertices);
	Vertex v;

	double avgAngleBoundary=0;
	double r=0.004,index;
	double sum=0;
	double limit = 0.2;
	int sizeVector=ceil(limit/r);
	vector<int> hist(sizeVector, 0);
	int nro=0;
	double maxvalue = 0;
	//cout<<boundaryVertices.size()<<endl;
	for(unsigned int i = 0; i < boundaryVertices.size(); ++i){
		v = boundaryVertices[i];
		if( v->info().is_singular_2() )		continue;
		vector<int> &regions = v->info().m_regions_around;
		for(unsigned int j = 0; j < regions.size(); ++j){
			//cout<<"stop 3 1/2 "<<v->point()<< " * "<<v->info().m_adyacent_vertices.size()<<endl;
			int label = regions[j];
			Vertex prev_v = v->info().m_adyacent_vertices[j][0];
			Vertex next_v = v->info().m_adyacent_vertices[j][1];
			double k = unsigned_curvature( prev_v, v, next_v );
			if( k > maxvalue) maxvalue=k;
			if( k < 0 || k > limit || !IsNumber(k) ) continue;
			if(k == 0)	k = 0.0001;
			

			index=ceil( abs(k)/r )-1;
			hist[index]++;
			sum += k;
			++nro;	
		}
	}
	//cout<<"max value"<<maxvalue<<endl;
	avgAngleBoundary=sum/(nro);
	ofstream o(filename);
	for(int i=0;i<sizeVector;i++){
		o<<hist[i]<<endl;
	}
	o<<avgAngleBoundary<<endl;
	o<<sum<<endl;
	o<<maxvalue<<endl;
	o.close();

	restart_vertex_state(dt);
}