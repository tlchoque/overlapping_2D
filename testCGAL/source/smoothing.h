

#define PI 3.14159265

double normVector(Vector a){
	return sqrt( pow(a.x(),2) + pow(a.y(),2) );
}

double getAngle(Vertex a, Vertex b,Vertex c){
	return acos( ( (a->point() - b->point() )*(c->point() - b->point() ) )/ ( normVector( a->point()- b->point() ) * normVector( c->point()- b->point() ) ))* 180.0 / PI;	
}

void getBoundary(Delaunay & dt, vector<pair<Vertex,Vertex>> &boundary){
	boundary.clear();
	for(Edge_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		Delaunay::Face &f1 = *(ei->first);
		int i = ei->second;
		Delaunay::Face &f2 = *f1.neighbor(i);
		if(f1.info().m_label == f2.info().m_label){
			continue;
		}
		if(f1.info().m_label < 0 || f2.info().m_label < 0){
			continue;
		}
		Vertex vs = f1.vertex(f1.cw(i));
		Vertex vt = f1.vertex(f1.ccw(i));
		boundary.push_back(make_pair(vs,vt));
	}
}

void getArrayBoundaries(Delaunay & dt,vector<pair<Vertex,Vertex>> &boundary, vector<deque<Vertex>> &boundaries ){
	boundaries.clear();
	deque<Vertex> pointArray;
	Vertex t=0,s=0,cen=0;
	int k=0;
	bool state=0;
	pointArray.clear();
	while(boundary.size() != 0){		
		if(pointArray.size() == 0){
			state=0;//state=0 we look for t vertex, state=1 look for s vertex
			pointArray.push_back(boundary[0].first);
			s = boundary[0].first;
			pointArray.push_back(boundary[0].second);
			t = boundary[0].second;
			boundary.erase(boundary.begin());
			cen=t;
		}
		//find the next vertex
		k=-1;
		for(unsigned int j=0;j < boundary.size();++j){
			if(cen == boundary[j].first){
				k=j;
				cen = boundary[j].second;
				break;
			}
			else if(cen == boundary[j].second){
				k=j;
				cen = boundary[j].first;
				break;
			}
		}
		if(k!=-1){
			if(state == 0)
				pointArray.push_back(cen);	
			else
				pointArray.push_front(cen);
			boundary.erase(boundary.begin() + k);
		}
		else {
			if(boundary.size() == 0){
				boundaries.push_back(pointArray);
				pointArray.clear();
			}
			else{
				if(state == 0){
					cen=s;
					state=1;
				}
				else{
					boundaries.push_back(pointArray);
					pointArray.clear();
				}
			}
		}
	}
	if(pointArray.size() != 0){
		boundaries.push_back(pointArray);
	}
}

void getBoundaryAngles(vector<deque<Vertex>> & boundaries){
	Vertex v;
	double t;
	for(unsigned int i=0;i<boundaries.size();++i){
		for(unsigned int j=1;j<boundaries[i].size()-1;++j){
			t=getAngle(boundaries[i][j-1],boundaries[i][j],boundaries[i][j+1]);
			boundaries[i][j]->info().m_boundary_angle=t;
		}
		if(boundaries[i][0] == boundaries[i][boundaries[i].size()-1]){
			t=getAngle(boundaries[i][boundaries[i].size()-2],boundaries[i][0],boundaries[i][1]);
			boundaries[i][0]->info().m_boundary_angle=t;
		}
		else{
			boundaries[i][0]->info().m_boundary_angle=180.0;
		}
	}
}

//relabeling
int getLowerRegion(Delaunay & dt,Vertex v){
	Vertex u;
	Vertex_circulator vc_start = dt.incident_vertices(v);
	Vertex_circulator vc = vc_start;

	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;

	double sum=0;
	int label=fc->info().m_label;
	int opositeLabel=0;;
	do{
		u=vc;
		++vc;
		if(fc->info().m_label == label)
			sum+=getAngle(u,v,vc);
		else
			opositeLabel=fc->info().m_label;
		++fc;
	} while (vc!= vc_start);
	if(sum < v->info().m_boundary_angle+2)
		return opositeLabel;
	else
		return label;
}

void setlabel(Delaunay & dt, Vertex v, int label){
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		fc->info().m_label=label;
		++fc;
	} while (fc!= fc_start);
}

void relabeling(Delaunay & dt, vector<deque<Vertex>> & boundaries, int num, int maxAngle){
	for(int k=0; k<num;++k){
		for(unsigned int i=0;i<boundaries.size();++i){
			for(unsigned int j=0;j<boundaries[i].size()-1;++j){
				//if(state==1)	state=0;
				if(boundaries[i][j]->info().m_boundary_angle < maxAngle ){
					cout<<boundaries[i][j]->point()<<endl;				
					cout<<" region label: "<<getLowerRegion(dt,boundaries[i][j])<<endl;
					setlabel( dt,boundaries[i][j],getLowerRegion(dt,boundaries[i][j]));
					++j;
					//state=1;
				}
			}
		}
	}
}

//TSUR
double getEle(Vertex p1, Vertex p3){
	return normVector(p3->point()-p1->point())/3;
}

void TSUR(vector<deque<Vertex>> & boundaries, int num){
	Vector t,n, x1,x2;
	for(int k=0;k < num;++k){//10
		for(unsigned int i=0;i<boundaries.size();++i){
			for(unsigned int j=0;j<boundaries[i].size()-3;++j){
				double L = getEle(boundaries[i][j],boundaries[i][j+3]);
				double A = CGAL::area(boundaries[i][j]->point(),boundaries[i][j+1]->point(),boundaries[i][j+2]->point() ) + CGAL::area(boundaries[i][j]->point(),boundaries[i][j+2]->point(),boundaries[i][j+3]->point() );
				double h = A/(2*L);
				double nr = normVector(boundaries[i][j+3]->point() - boundaries[i][j]->point());
				t =  (boundaries[i][j+3]->point() - boundaries[i][j]->point() )*(1/nr);
				n = Vector(t.y(),-t.x());
				x1 = (boundaries[i][j]->point() - CGAL::ORIGIN)+ L*t + h*n;
				x2 = (boundaries[i][j]->point() - CGAL::ORIGIN)+ 2*L*t + h*n;
				boundaries[i][j+1]->set_point( Point(x1.x(), x1.y() ));
				boundaries[i][j+2]->set_point( Point(x2.x(), x2.y() ));
			}
		}
	}
}