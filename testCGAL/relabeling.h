
#include <queue>

struct Ccomponent{
	Vertex v;
	double k;
	int label;
	Ccomponent( Vertex vi,int labeli,double ki ){
		v = vi;
		label = labeli;
		k = ki;
	}
	bool operator<(const Ccomponent& rhs) const {
        return k < rhs.k;
    }
};

bool comparator2(const Ccomponent &a, const Ccomponent &b){
    return (a.k < b.k);
}

bool is_feature(Delaunay &dt,double k, Vertex prev_v, Vertex next_v,int label){
	double sum = 0,ki;
	if( !prev_v->info().is_singular_2() ){
		if( vertex_curvature(dt,prev_v,label,ki) )
			sum += ki;
	}
	if( !next_v->info().is_singular_2() ){
		if( vertex_curvature(dt,next_v,label,ki) )
			sum += ki;
	}
	//cout<<"sum and c: "<<sum<<" * "<<k<<endl;
	//if( sum*k  > 0 )	return true;
	if( sum*k  > 0  || abs(k/sum) > 8)	return true;
	//if( sum*k  > 0 )	return true;
	else	return false;
}

bool is_isolated(Delaunay &dt,vector<Face> &faces,Vertex v){
	int label = faces[0]->info().m_label;
	for(unsigned int i = 0; i< faces.size(); ++i){
		Face f =faces[i];
		Face n = f->neighbor( f->index(v) );
		//cout<<" n and new_label "<<n->info().m_label <<" * " <<label<<endl;
		if( n->info().m_label == label  )
			return false;
	}
	return true;
}

void mesh_curvature(Delaunay &dt, vector<Vertex> &boundary, priority_queue<Ccomponent> &peaks, double &total_curvature){	
	for(unsigned int i = 0 ; i < boundary.size(); ++i){
		vector<Ccomponent> components;
		Vertex v = boundary[i];
		bool add = true;
		//cout<<"vertex "<<v->point()<<endl;
		for(unsigned int j = 0; j < v->info().m_regions_around.size(); ++j){
			Vertex prev_v,next_v;
			vector<Face> faces;
			int label = v->info().m_regions_around[j];
			double k;
			if ( curvature_and_neighbors(dt,v,label,k,faces,prev_v,next_v) ){
				vector<Point> ptos;
				if( is_feature(dt,k,prev_v,next_v,label) && !is_isolated(dt,faces,v) ){
					add = false;
					v->info().m_feature = true;
					//cout<<"is feature "<<k<<endl;
					//ptos.push_back(v->point());
					//drawMesh(dt,ptos),
					////cout<<label<<" * "<<k<<endl;
					//getchar();
					break;
				}

				//cout<<"not feature "<<k<<endl;
				//ptos.push_back(v->point());
				//drawMesh(dt,ptos),
				////cout<<label<<" * "<<k<<endl;
				//getchar();

				if( k > 0 )
					components.push_back(Ccomponent(v,label,k ) ); 
			}			
		}
		if( add ){
			for(unsigned int l = 0; l < components.size(); ++l){
				Ccomponent c = components[l];
				peaks.push(c); 
			}
		}
	}
}
// il manque checking

bool is_relabeled(vector<Face> &faces){
	for(unsigned int i = 0; i < faces.size(); ++i){
		Face f = faces[i];
		if( f->info().m_relabeled )
			return true;
	}
	return false;
}

void add_component(Delaunay &dt,Vertex v,priority_queue<Ccomponent> &peaks){
	if( v->info().m_feature ) return;	
	bool add = true;
	vector<Ccomponent> components;
	for(unsigned int i = 0; i < v->info().m_regions_around.size(); ++i){
		Vertex prev_v,next_v;
		vector<Face> faces;
		int label = v->info().m_regions_around[i];
		double k;
		//cout<<"add stop 0 1"<<endl;
		if ( curvature_and_neighbors(dt,v,label,k,faces,prev_v,next_v) ){
			//cout<<"add stop 0 2"<<endl;
			vector<Point> ptos;
			ptos.push_back(v->point());
			/*drawMesh(dt,ptos),
			cout<<label<<" * "<<k<<endl;
			getchar();*/

			//if( is_feature(dt,k,prev_v,next_v,label)  ){
			/*if( v->info().m_feature  ){	
				cout<<"feat"<<endl;
				add = false;
				break;
			}*/
			//cout<<"not feat"<<endl;
			if( k > 0 && !is_relabeled(faces) )
				components.push_back(Ccomponent(v,label,k ) ); 
		}		
	}
	//cout<<"add stop 1"<<endl;
	if( add ){
		for(unsigned int l = 0; l < components.size(); ++l){
			Ccomponent c = components[l];
			peaks.push(c); 
		}
	}
}

bool right_neigborhood(Delaunay &dt,Vertex v,vector<Face> &faces,int label){
	for(unsigned int i = 0; i< faces.size(); ++i){
		Face f =faces[i];
		Face n = f->neighbor( f->index(v) );
		//cout<<" n and new_label "<<n->info().m_label <<" * " <<label<<endl;
		if( n->info().m_label == label  )
			return false;
	}
	return true;
}

bool is_possible_to_relabel(Delaunay &dt, Ccomponent c,vector<Vertex> &neighbors){
	//is label info of vertex 
	if( !( std::find(c.v->info().m_regions_around.begin(), c.v->info().m_regions_around.end(), c.label ) != c.v->info().m_regions_around.end() ) ) return false;
	vector<Face> faces = faces_in_vertex_by_label(dt,c.v,c.label);
	Vertex prev_v,next_v;
	if( is_relabeled(faces) ) return false;
	int new_label = new_label_for_relabeling(dt,faces,c.v,c.label); 
	if( !right_neigborhood(dt,c.v,faces,new_label) ) {
		return false;
	}

	//cout<<" possible stop 1 "<<endl;
	vector<Vertex> vertices;
	neighbor_vertices(dt,faces,vertices,c.v);
	vertices.push_back(c.v);
	double sum_before=0;
	for(unsigned int i = 0; i< vertices.size(); ++i){
		Vertex v = vertices[i];
		if( v->info().is_boundary() )
			sum_before+=total_vertex_curvature(dt,v);
	}

	//cout<<" possible stop 2 "<<endl;
	for(unsigned int i = 0; i< faces.size(); ++i)
		faces[i]->info().m_label = new_label;

	double sum_after=0;
	for(unsigned int i = 0; i< vertices.size(); ++i){
		Vertex v = vertices[i];
		update_vertex_info(dt,v);
		if( v->info().is_boundary() && !v->info().is_singular_2() )
			sum_after+=total_vertex_curvature(dt,v);
	}	
	//cout<<" possible stop 3 "<<endl;
	//cout<<sum_before<<" * "<<sum_after<<endl;


	//if( produce_singularities(dt,vertices) || sum_before < sum_after ){
	if( produce_singularities(dt,vertices) || sum_before/sum_after < 2){//1.355 change to 3
		//cout<<" possible stop 3 1"<<endl;
		for(unsigned int i = 0; i< faces.size(); ++i)
			faces[i]->info().m_label = c.label;
		
		for(unsigned int i = 0; i< vertices.size(); ++i)
			update_vertex_info(dt,vertices[i]);

		return false;
	}
	else{
		//cout<<" possible stop 3 2"<<endl;
		for(unsigned int i = 0; i< faces.size(); ++i)
			faces[i]->info().m_relabeled = true;

		c.v->info().m_state=4;
		for(unsigned int i = 0; i< vertices.size(); ++i){
			Vertex v = vertices[i];
			if( v->info().m_state == 1 && v->info().is_boundary() )
				neighbors.push_back(v);
		}

		/*vector<Face> test;
		cout<<"curvature_and_neighbors "<<endl;
		curvature_and_neighbors(dt,c.v,c.label,k,test,prev_v,next_v);
		cout<<"prev_v "<<k<<prev_v->point();
		neighbors.push_back(prev_v);
		neighbors.push_back(next_v);*/
		return true;
	}
}

void relabeling_to_smooth(Delaunay &dt){
	vector<Vertex> boundary;	
	boundary_vertices2(dt,boundary);
	double total_curvature;
	std::priority_queue<Ccomponent> peaks;


	clock_t begin = clock();

	mesh_curvature(dt,boundary,peaks,total_curvature);

	vector<Point> pts;
	for(unsigned int i = 0; i < boundary.size(); ++i){
		Vertex v = boundary[i];
		if( v->info().m_feature ) pts.push_back(v->point());
	}
	/*savePoints2(pts,"thundercats/relabeling/features.vtk",3);
	cout<<"finish features"<<endl;*/
	//getchar();
	
	while( peaks.size() !=0 ){
		Ccomponent c = peaks.top();
		peaks.pop();
		vector<Vertex> neighbors;		
		vector<Point> ptos;
		//cout<<" before "<<c.label<<endl;
		/*ptos.push_back(c.v->point());
		drawMesh(dt,ptos),
		getchar();*/

		if( is_possible_to_relabel(dt,c,neighbors) ){
			//cout<<" after"<<endl;
			/*drawMesh(dt,ptos),
			getchar();*/

			/*for( unsigned int i = 0 ; i < neighbors.size(); ++i ){
				add_component(dt,neighbors[i],peaks);
			}*/
		}
	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time relabeling: "<<elapsed_secs<<endl;

	restart_vertex_state(dt);	
}
