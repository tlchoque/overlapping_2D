
#define PI 3.14159265

double getRatio(Vertex a, Vertex b,Vertex c){
	Vector x =a->point() - b->point();
	Vector y =c->point() - b->point();
	return ( x*y )/ ( normVector( x ) * normVector( y ) );	
}

void getBoundary(Delaunay & dt, vector<Vertex> &boundary){
	boundary.clear();
	for(Edge_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		Delaunay::Face &f1 = *(ei->first);
		int i = ei->second;
		Delaunay::Face &f2 = *f1.neighbor(i);
		
		//if(f1.info().m_label < 0 || f2.info().m_label < 0)	continue;
		if(f1.info().m_label == f2.info().m_label)	continue; //share the same label
		if(f1.info().m_label < 0 && f2.info().m_label == 0)	continue; // limit and background
		if(f1.info().m_label == 0 && f2.info().m_label < 0)	continue; // limit and background

		Vertex vs = f1.vertex(f1.cw(i));
		Vertex vt = f1.vertex(f1.ccw(i));
		if (vs->info().m_state == 0)	{
			vs->info().m_state=1;
			boundary.push_back(vs);
		}

		if (vt->info().m_state == 0)	{
			vt->info().m_state=1;
			boundary.push_back(vt);
		}
	}
}

int getLabel(Face_circulator f){
	Face_circulator F = f;
	Face_circulator prevF = --f;
	if( F->info().m_label == prevF->info().m_label  )	return -1;
	if(F->info().m_label != 0)	return F->info().m_label;
	else return prevF->info().m_label;
}

int compareLabel(Face_circulator f, int label){
	Face_circulator F = f;
	Face_circulator prevF = --f;
	if(F->info().m_label == prevF->info().m_label) return -1;
	if(F->info().m_label == label) return 1;
	if(prevF->info().m_label == label) return 1;
	return 0;
}

void boundary_sets(Delaunay & dt,vector<Vertex> &boundary, vector<deque<Vertex>> &boundaries){
	boundaries.clear();
	for(unsigned int i = 0; i < boundary.size(); ++i){
		if(boundary[i]->info().m_state != 1)		continue;
		vector<Vertex> search;
		deque<Vertex> contour;
		boundary[i]->info().m_state=2;
		search.push_back(boundary[i]);
		contour.push_back(boundary[i]);
		int prevLabel = -1;
		Vertex flag = 0;
		bool push=1;
		while(search.size() != 0){
			Vertex v = search.back();
			search.pop_back();
			if(v == flag)	push=0;
			int count = 0;			
			Vertex_circulator vc_start = dt.incident_vertices(v);
			Vertex_circulator vc = vc_start;
			Face_circulator fc = dt.incident_faces(v);
			do{
				if(vc->info().m_state == 1){
					if(prevLabel == -1){
						int label = getLabel(fc);
						if(label != -1){
							Vertex c = vc;
							c->info().m_state = 2;
							search.push_back(c);
							contour.push_front(c);							
							prevLabel = getLabel(fc);
							flag = c;
						}
					}
					else{
						int labelCheck = compareLabel(fc,prevLabel);
						if( labelCheck == 1 ){
							Vertex c = vc;
							c->info().m_state = 2;
							search.push_back(c);
							if(push == 1)	contour.push_back(c);
							else	contour.push_front(c);	
							++count; //not count when shares the same label
						}
						if(labelCheck == 0)
							++count;						
					}
				}
				else if( dt.is_infinite(fc) ){
					v->info().m_limit=1;
				}
				++vc;
				++fc;
			} while (vc!= vc_start);
			if(count > 1)	v->info().m_restricted=1;
		}
		boundaries.push_back(contour);
	}
}

void getBoundaryAngles(Delaunay & dt,vector<deque<Vertex>> & boundaries){
	Vertex v; double t,t2;	
	for(unsigned int i=0;i<boundaries.size();++i){
		bool isClosed = dt.is_edge( boundaries[i].front(),boundaries[i].back() ); // to walk trhough all points
		if( isClosed ){
			boundaries[i].push_back(boundaries[i][0]);
			boundaries[i].push_back(boundaries[i][1]);
		}
		for(unsigned int j=1;j<boundaries[i].size()-1;++j){
			if( boundaries[i][j]->info().m_limit || boundaries[i][j]->info().m_restricted )	  continue;
			t=getRatio(boundaries[i][j-1],boundaries[i][j],boundaries[i][j+1]);
			t2=getAngle(boundaries[i][j-1],boundaries[i][j],boundaries[i][j+1]);
			boundaries[i][j]->info().m_ratio=t;
			boundaries[i][j]->info().m_boundary_angle=t2;
		}
		if( isClosed ){
			boundaries[i].pop_back();
			boundaries[i].pop_back();
		}
	}
}

//relabeling

int getLowerRegionRatio(Delaunay & dt,Vertex bv, Vertex v, Vertex av){// return label wich replace the lower
	Vertex u;
	Vertex_circulator vc_start = dt.incident_vertices(v);
	Vertex_circulator vc = vc_start;
	Face_circulator fc = dt.incident_faces(v);
	Vector a,b;
	a = bv->point() - v->point();
	b = av->point() - v->point();
	double cross = a.x()*b.y()-a.y()*b.x();
	if(cross < 0){
		Vertex t = bv;
		bv = av;
		av=t;
	}	

	int label,opositeLabel;
	do{
		u=vc;
		++vc;
		if(u == bv){
			label = fc->info().m_label;
			opositeLabel = fc->neighbor(fc->index(vc))->info().m_label;
			break;
		}		
		++fc;
	} while (vc!= vc_start);

	return opositeLabel;
	/*double a = bv->point().x() - v->point().x();
	double b = bv->point().y() - v->point().y();
	double c = v->point().x()*b - v->point().y()*a;


	cout<<label<<" - "<<opositeLabel<<endl;
	cout<<bv->point()<<" * "<<v->point()<<" * "<<av->point()<<endl;
	cout<<a<<" : "<<b<<" : "<<c<<endl;
	cout<<a*av->point().x() - b*av->point().y() + c<<endl;
	getchar();

	if( a*av->point().x() - b*av->point().y() + c >= 0)
			return opositeLabel;
		else 
			return label;*/
}

void setlabel(Delaunay & dt, Vertex v, int label){
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		fc->info().m_label=label;
		++fc;
	} while (fc!= fc_start);
}

void relabeling(Delaunay & dt, vector<deque<Vertex>> & boundaries, int maxAngle){
	double ratio= cos(maxAngle*PI/180);
	for(unsigned int i=0;i<boundaries.size();++i){
		for(int j=0;j<(int)boundaries[i].size();++j){
			if( boundaries[i][j]->info().m_limit || boundaries[i][j]->info().m_restricted )		continue;
			if( boundaries[i][j]->info().m_ratio > ratio ){
				if( j+2 < (int)boundaries[i].size() && boundaries[i][j]->info().m_ratio < boundaries[i][j+1]->info().m_ratio){
					setlabel( dt,boundaries[i][j+1],getLowerRegionRatio(dt,boundaries[i][j],boundaries[i][j+1],boundaries[i][j+2]));
					j+=2;			
				}
				else{
					if(j-1 < 0) continue;
					if(j+1 < (int) boundaries[i].size() ){
					setlabel( dt,boundaries[i][j],getLowerRegionRatio(dt,boundaries[i][j-1],boundaries[i][j],boundaries[i][j+1]));
					++j;
					}
				}
			}
		}
	}
}

void restartBoundary(vector<Vertex> &boundary){
	for(unsigned int i = 0; i < boundary.size();++i){
		boundary[i]->info().m_state=0;
		boundary[i]->info().m_limit=0;
		boundary[i]->info().m_restricted=0;
		boundary[i]->info().m_boundary_angle=180;
		boundary[i]->info().m_ratio=-1;
	}
}

void relabel2D(Delaunay & dt,int num,int maxAngle){
	vector<Vertex> boundary;
	vector<deque<Vertex>> boundaries;
	for(int i=0;i<num;++i){
		getBoundary(dt,boundary);		
		boundary_sets(dt,boundary,boundaries);
		getBoundaryAngles(dt,boundaries);
		relabeling(dt,boundaries,maxAngle);
		restartBoundary(boundary);
	}
}

//TSUR
double getEle(Point p1, Point p3){
	return normVector(p3 - p1 )/3;
}

void TSUR(Delaunay &dt, int num){
	vector<deque<Vertex>> boundaries;
	vector<Vertex> boundary;
	boundary_vertices2(dt,boundary);
	vector<Point> pts;
	

	boundary_sets(dt,boundary,boundaries);
	Vector t,n, x1,x2;

	clock_t begin = clock();

	for(int k=0;k < num;++k){//10
		for(unsigned int i=0;i<boundaries.size();++i){
			if( boundaries[i].size() < 8)	continue;

			bool isClosed = dt.is_edge( boundaries[i].front(),boundaries[i].back() ); // to walk trhough all points
			if( isClosed){
				boundaries[i].push_back(boundaries[i][0]);
				boundaries[i].push_back(boundaries[i][1]);
				boundaries[i].push_back(boundaries[i][2]);
			}
			for(unsigned int j=0;j<boundaries[i].size()-3;++j){		
				Vertex a,b,c,d;
				a = boundaries[i][j];
				b = boundaries[i][j+1];
				c = boundaries[i][j+2];
				d = boundaries[i][j+3];
				if( b->info().m_limit || b->info().m_restricted || c->info().m_limit || c->info().m_restricted || b->info().is_corner() || c->info().is_corner() )	continue;
				double L = getEle( boundaries[i][j]->point() , boundaries[i][j+3]->point() );
				double A = CGAL::area(boundaries[i][j]->point(),boundaries[i][j+1]->point(),boundaries[i][j+2]->point() ) + CGAL::area(boundaries[i][j]->point(),boundaries[i][j+2]->point(),boundaries[i][j+3]->point() );

				if(2*L==0 ) continue;
				double h = A/(2*L);

				double nr = normVector(boundaries[i][j+3]->point() - boundaries[i][j]->point());
				if(nr == 0) continue;
				t =  (boundaries[i][j+3]->point() - boundaries[i][j]->point() )/nr;
				n = Vector(t.y(),-t.x());
				x1 = (boundaries[i][j]->point() - CGAL::ORIGIN)+ L*t + h*n;
				x2 = (boundaries[i][j]->point() - CGAL::ORIGIN)+ 2*L*t + h*n;
				/*cout<<x1<<" * "<<x2<<endl;
				getchar();*/
				boundaries[i][j+1]->set_point( Point(x1.x(), x1.y() ));
				boundaries[i][j+2]->set_point( Point(x2.x(), x2.y() ));
			}
			if( isClosed ){
				boundaries[i].pop_back();
				boundaries[i].pop_back();
				boundaries[i].pop_back();
			}
		}
		/*drawMesh(dt);
		getchar();*/

	}
	clock_t end = clock();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	cout<<"time smoothing: "<<elapsed_secs<<endl;

	restartBoundary(boundary);
}

void initialize_smoothed_point(vector<Vertex> &boundary){
	for(unsigned int i = 0; i < boundary.size(); ++i){
		boundary[i]->info().m_smoothed_point = boundary[i]->point();
	}
}

void TSUR_points(Delaunay &dt, int num){
	//vector<deque<Vertex>> boundaries;
	//vector<Vertex> boundary;
	//boundary_vertices(dt,boundary);
	//initialize_smoothed_point(boundary);
	//boundary_sets(dt,boundary,boundaries);	

	//Vector t,n, x1,x2;
	//for(int k=0;k < num;++k){//10
	//	for(unsigned int i=0;i<boundaries.size();++i){
	//		if( boundaries[i].size() < 4)	continue;
	//		bool isClosed = dt.is_edge( boundaries[i].front(),boundaries[i].back() ); // to walk trhough all points
	//		if( isClosed){
	//			boundaries[i].push_back(boundaries[i][0]);
	//			boundaries[i].push_back(boundaries[i][1]);
	//			boundaries[i].push_back(boundaries[i][2]);
	//		}
	//		for(unsigned int j=0;j<boundaries[i].size()-3;++j){
	//			Vertex v1 = boundaries[i][j];
	//			Vertex v2 = boundaries[i][j+1];
	//			Vertex v3 = boundaries[i][j+2];
	//			Vertex v4 = boundaries[i][j+3];
	//			if( v2->info().m_limit || v2->info().m_restricted || v3->info().m_limit || v3->info().m_restricted)	continue;
	//			Point p1 = v1->info().m_smoothed_point;
	//			Point p2 = v2->info().m_smoothed_point;
	//			Point p3 = v3->info().m_smoothed_point;
	//			Point p4 = v4->info().m_smoothed_point;
	//			double L = getEle( p1, p4 );
	//			double A = CGAL::area( p1 , p2 , p3 ) + CGAL::area( p1 , p3 , p4 );
	//			double h = A/(2*L);
	//			double nr = normVector( p4 - p1 );
	//			t =  ( p4 - p1 )*( 1/nr );
	//			n = Vector(t.y(),-t.x());
	//			x1 = ( p1 - CGAL::ORIGIN ) + L*t + h*n;
	//			x2 = ( p1 - CGAL::ORIGIN ) + 2*L*t + h*n;
	//			boundaries[i][j+1]->info().m_smoothed_point =  Point(x1.x(), x1.y() );
	//			boundaries[i][j+2]->info().m_smoothed_point =  Point(x2.x(), x2.y() );
	//		}
	//		if( isClosed ){
	//			boundaries[i].pop_back();
	//			boundaries[i].pop_back();
	//			boundaries[i].pop_back();
	//		}
	//	}
	//}
	//restartBoundary(boundary);

	/*vector<deque<Vertex>> boundaries;
	vector<Vertex> boundary;
	boundary_vertices2(dt,boundary);*/
	
	vector<deque<Vertex>> boundaries;
	vector<Vertex> boundary;
	boundary_vertices2(dt,boundary);
	initialize_smoothed_point(boundary);
	boundary_sets(dt,boundary,boundaries);	


	Vector t,n, x1,x2;
	for(int k=0;k < num;++k){//10
		for(unsigned int i=0;i<boundaries.size();++i){
			/*cout<<"here"<<endl;
			getchar();*/
			if( boundaries[i].size() < 10)	continue;

			bool isClosed = dt.is_edge( boundaries[i].front(),boundaries[i].back() ); // to walk trhough all points
			if( isClosed){
				boundaries[i].push_back(boundaries[i][0]);
				boundaries[i].push_back(boundaries[i][1]);
				boundaries[i].push_back(boundaries[i][2]);
			}
			for(unsigned int j=0;j<boundaries[i].size()-3;++j){		
				Vertex a,b,c,d;
				a = boundaries[i][j];
				b = boundaries[i][j+1];
				c = boundaries[i][j+2];
				d = boundaries[i][j+3];
				if( b->info().m_limit || b->info().m_restricted || c->info().m_limit || c->info().m_restricted || b->info().is_corner() || c->info().is_corner() )	continue;
				double L = getEle( boundaries[i][j]->point() , boundaries[i][j+3]->point() );
				double A = CGAL::area(boundaries[i][j]->point(),boundaries[i][j+1]->point(),boundaries[i][j+2]->point() ) + CGAL::area(boundaries[i][j]->point(),boundaries[i][j+2]->point(),boundaries[i][j+3]->point() );

				if(2*L==0 ) continue;
				double h = A/(2*L);

				double nr = normVector(boundaries[i][j+3]->point() - boundaries[i][j]->point());
				if(nr == 0) continue;
				t =  (boundaries[i][j+3]->point() - boundaries[i][j]->point() )/nr;
				n = Vector(t.y(),-t.x());
				x1 = (boundaries[i][j]->point() - CGAL::ORIGIN)+ L*t + h*n;
				x2 = (boundaries[i][j]->point() - CGAL::ORIGIN)+ 2*L*t + h*n;
				/*cout<<x1<<" * "<<x2<<endl;
				getchar();*/
				//cout<<"before "<<boundaries[i][j+1]->info().m_smoothed_point<<endl;
				boundaries[i][j+1]->info().m_smoothed_point =  Point(x1.x(), x1.y() );
				boundaries[i][j+2]->info().m_smoothed_point =  Point(x2.x(), x2.y() );

				/*boundaries[i][j+1]->set_point( Point(x1.x(), x1.y() ));
				boundaries[i][j+2]->set_point( Point(x2.x(), x2.y() ));*/
				//cout<<"before "<<boundaries[i][j+1]->info().m_smoothed_point<<endl;
				//getchar();
			}
			if( isClosed ){
				boundaries[i].pop_back();
				boundaries[i].pop_back();
				boundaries[i].pop_back();
			}
		}
		/*drawMesh(dt);
		getchar();*/

	}
	restartBoundary(boundary);
}

int loorForLabel(vector<pair<int,double>> &labels, int lbl){
	for(unsigned int i=0;i< labels.size();++i){
		if(labels[i].first==lbl)
			return i;
	}
	return -1;
}

void checkManifold(Delaunay &dt){
	vector<Vertex> boundary;
	getBoundary(dt,boundary);
	vector<pair<int,double>> around;
	int changes=0,lbl=0,sum=0,vx1=0,vx2=0,vx3=0,a=0,b=0,theta=0;
	for(unsigned int i = 0; i < boundary.size();++i){
		sum=0;
		around.clear();
		changes=0;
		Face_circulator fc = dt.incident_faces(boundary[i]);
		Face_circulator fc_start=fc;
		lbl = fc->info().m_label;
		do{
			if( dt.is_infinite(fc) ){
				++changes;
				if(changes == 1){// start to iterate
					fc_start = fc;
				}else{
					around.push_back(make_pair(lbl,sum)); // push the last iterations
				}
				sum = 360;
				lbl=-1;
				++fc;
				continue;
			}
			
			vx1 = fc->index(boundary[i]);
			vx2 = (vx1+1)%3;
			vx3 = (vx1+2)%3;

			Vertex v1 = fc->vertex( vx1 );//angle
			Vertex v2 = fc->vertex( vx2 );
			Vertex v3 = fc->vertex( vx3 );

			if(lbl != fc->info().m_label){
				++changes;		
				if(changes == 1){
					fc_start = fc;	
				}
				else{
					around.push_back(make_pair(lbl,sum));
				}
				sum = getAngle(v2,v1,v3);
				lbl=fc->info().m_label;
			}
			else{
				if(changes != 0)
					sum += getAngle(v2,v1,v3);
			}
			++fc;
		} while (fc!= fc_start);

		if(changes==0)
			around.push_back(make_pair(lbl,360));
		else
			around.push_back(make_pair(lbl,sum));

		vector<pair<int,double>> labels;
		int cen=0;
		for(unsigned int i=0;i< around.size();++i){
			if(around[i].first == -1)	continue;
			if(labels.size() ==  0){labels.push_back(make_pair(around[i].first,around[i].second));continue;}
			int look=loorForLabel(labels,around[i].first);
			if( look < 0) 
				labels.push_back(make_pair(around[i].first,around[i].second));
			else{
				labels[look].second+=around[i].second;
				cen=1;
			}
		}

		if(cen==0)	continue;
		/*cout<<"reaches here"<<endl;
		for(unsigned int i=0;i< around.size();++i){
			cout<<around[i].first<<" " <<around[i].second<<endl;
		}*/

		int greater=0;
		for(unsigned int i=1;i< labels.size();++i){//i=1
			//cout<<labels[i].first<<" " <<labels[i].second<<endl;
			if(labels[greater].second<labels[i].second){
				greater = i;
			}
		}
		fc=fc_start;
		do{			
			if( dt.is_infinite(fc) ){++fc;continue;}
			fc->info().m_label=labels[greater].first;
			++fc;
		}while(fc!=fc_start);

	}
	restartBoundary(boundary);
}

void getAngleBoundaryHistogram(vector<deque<Vertex>> &boundaries){
	double avgAngleBoundary=0;
	int r=20,index;
	double sum=0;
	int sizeVector=ceil(180/r);
	vector<int> hist(sizeVector, 0);

	int count=0;
	
	//cout<<"size : "<<boundaries.size()<<endl;
	for(unsigned int i=0;i<boundaries.size();++i){
		//cout<<"size : "<<boundaries[i].size()<<endl;
		for(unsigned int j=0;j<boundaries[i].size();++j){
			if( boundaries[i][j]->info().m_limit || boundaries[i][j]->info().m_restricted )	  continue;
			if(boundaries[i][j]->info().m_boundary_angle < 0 || boundaries[i][j]->info().m_boundary_angle >180) continue;
			//cout<<boundaries[i][j]->info().m_boundary_angle<<endl;
			index=ceil(boundaries[i][j]->info().m_boundary_angle/r)-1;

			if(index <0 || index >(int)hist.size()) continue;
			//cout<<index<<endl;
			hist[index]++;
			sum+=boundaries[i][j]->info().m_boundary_angle;
			++count;
		}
	}	
	
	//avgAngleBoundary=sum/(boundaries[0].size()-2);
	avgAngleBoundary=sum/(count);
	ofstream o("histoC.txt");
	for(int i=0;i<sizeVector;i++){
		o<<hist[i]<<endl;
	}
	o<<avgAngleBoundary<<endl;
}