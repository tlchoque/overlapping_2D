
#include <algorithm> 
#include <limits>

double squaredHalfAnnulusWidth(Point p1, Point p2, Point p3, Point p4){
	double a1,b1,c1,a2,b2,c2;
	double den,g,nx,ny,xannulus,yannulus;
	double r2,R2,w;

	a1=2*( p2.x() - p4.x() );
	b1=2*( p2.y() - p4.y() );
	c1=( pow( p4.x() , 2) + pow( p4.y() , 2) ) - ( pow( p2.x() , 2) + pow( p2.y() , 2) );

	a2=2*(p1.x() - p3.x());
	b2=2*(p1.y() - p3.y());
	c2=( pow( p3.x() , 2) + pow( p3.y() , 2) ) - ( pow( p1.x() , 2) + pow( p1.y() , 2) );

	g = (a1 * b2) - (b1 * a2);
	nx = (b1 * c2) - (b2 * c1);
	ny = (a2 * c1) - (a1 * c2);

	den = 1/g;
	xannulus = nx*den;
	yannulus = ny*den;

	r2 = pow( (p2.x() - xannulus) , 2 ) + pow( (p2.y() - yannulus), 2 );
	R2 = pow( (p3.x() - xannulus), 2 ) + pow( (p3.y() - yannulus) , 2);
	w = (r2 + R2 - 2 * sqrt(r2 * R2))/4;
	return w;
}

double boundarySquaredHalfAnnulusWidth(Vertex p1, Vertex p2, Vertex p3){
	double a2,b2;
	a2 = pow( CGAL::area(p1->point(),p2->point(),p3->point()) , 2) ;
	b2 = pow( p3->point().x() - p1->point().x() , 2 ) +  pow( p3->point().y() - p1->point().y() , 2 );
	return a2/b2;
}

void orderIndex(int &a,int &b,int &c){
	if(a != 0 && b!= 0)		c = 0;
	else if(a != 1 && b!= 1)	c = 1;
	else if(a != 2 && b!= 2)	c = 2;
}

void initializationTolerance(Delaunay &dt){
	for(Edge_iterator ei=dt.finite_edges_begin();ei!=dt.finite_edges_end(); ei++){
		Face f1 = ei->first;
		Face f2 = f1->neighbor(ei->second);
		Vertex p1 = f1->vertex(f1->cw(ei->second));
		Vertex p3 = f1->vertex(f1->ccw(ei->second));

		if(!dt.is_infinite(f2) && !dt.is_infinite(f1)){
			int a,b,c;
			a = f2->index(p1);
			b = f2->index(p3);
			orderIndex(a,b,c);
			Vertex p2 = f1->vertex(ei->second);
			Vertex p4 = f2->vertex(c);
			double w = squaredHalfAnnulusWidth(p1->point(),p2->point(),p3->point(),p4->point());
			//cout<<p1->point()<<" " << p3->point()<<" w: "<<sqrt(w)<<endl;
			p1->info().m_bicellTolerances.push_back(w);
			p2->info().m_bicellTolerances.push_back(w);
			p3->info().m_bicellTolerances.push_back(w);
			p4->info().m_bicellTolerances.push_back(w);
		}
		else{
			Vertex p2;
			if( dt.is_infinite(f1) ){
				int a,b,c;
				a = f2->index(p1);
				b = f2->index(p3);
				orderIndex(a,b,c);
				p2 = f2->vertex(c);
			}
			else{
				p2 = f1->vertex(ei->second);
			}
			double w = boundarySquaredHalfAnnulusWidth(p1,p2,p3);
			//cout<<p1->point()<<" " << p3->point()<<" w: "<<sqrt(w)<<endl;
			p1->info().m_bicellTolerances.push_back(0);
			p2->info().m_bicellTolerances.push_back(w);
			p3->info().m_bicellTolerances.push_back(0);
		}		
	}

	// get the minimum bicell tolerance for vertex
	for (Vertex_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
		vi->info().m_tolerance = *std::min_element(vi->info().m_bicellTolerances.begin(), vi->info().m_bicellTolerances.end());
		cout<<sqrt(vi->info().m_tolerance)<<endl;
		vi->info().m_moving.x = vi->point().x();
		vi->info().m_moving.y = vi->point().y();
		vi->info().m_fixed.x = vi->point().x();
		vi->info().m_fixed.y = vi->point().y();
	} 
}


bool checkDelaunay(Vertex p1, Vertex p2, Vertex p3, Vertex p4){
	double sr = CGAL::squared_radius(p1->point(),p2->point(),p3->point());
	Point center = CGAL::circumcenter(p1->point(),p2->point(),p3->point());
	double dist = CGAL::squared_distance(p4->point(),center);
	if( dist < sr )		return 1;
	else return 0;
}

void relocation(Delaunay &dt, Vertex v,vector<Delaunay::Edge> &edges){
	Edge_circulator ec_start = dt.incident_edges(v);
	Edge_circulator ec = ec_start;
	Face_circulator fc = dt.incident_faces(v);
	
	do{
		Face f1 = ec->first;
		Face f2 = f1->neighbor(ec->second);

		if( dt.is_infinite(f1) || dt.is_infinite(f2) ) {fc++;ec++;continue;}

		// set f2 as the previous face ; edge is p1 an p3, p2 is par of f1 face ;
		Vertex p1,p2,p3,p4;
		int a,b,c;
		if( f2 == fc){
			f2 = f1;
			f1 = fc;
			p1 = f2->vertex(f2->cw(ec->second));
			p3 = f2->vertex(f2->ccw(ec->second));
			a = f1->index(p1);
			b = f1->index(p3);
			orderIndex(a,b,c);
			p2 = f2->vertex(ec->second);
			p4 = f1->vertex(c);
		}
		else{
			p1 = f1->vertex(f1->cw(ec->second));
			p3 = f1->vertex(f1->ccw(ec->second));
			a = f2->index(p1);
			b = f2->index(p3);
			orderIndex(a,b,c);
			p2 = f1->vertex(ec->second);
			p4 = f2->vertex(c);
		}
		
		cout<<dt.segment(ec) <<endl;
		cout<<dt.triangle(fc) <<endl;
		++ec;
		++fc;
	} while (ec!= ec_start);
}

void filteringRelocation( Delaunay &dt, Vertex v, Point p){
	v->info().m_moving.x = p.x();
	v->info().m_moving.y = p.y();
	deque<Vertex> Q;
	Vertex h;
	double inf = std::numeric_limits<double>::infinity();
	double D = sqrt( pow( p.x() - v->info().m_fixed.x , 2) + pow( p.y() - v->info().m_fixed.y , 2) );
	if( D < v->info().m_tolerance ){
		v->info().m_tolerance -= D;
	}
	else{
		Q.push_back(v);
		while( Q.size() != 0 ){
			double Dh;
			h = Q.front();
			Q.pop_front();
			h->info().m_fixed.x = h->info().m_moving.x;
			h->info().m_fixed.y = h->info().m_moving.y;
			h->info().m_tolerance = inf;
			Dh = 0;
			vector<Delaunay::Edge> edges;

		}
	}
}

void circle_properties(Point p1, Point p2, Point p3, Point p4, Point &c, double &r){
	double a1,b1,c1,a2,b2,c2;
	double den,g,nx,ny,xannulus,yannulus;
	double r2,R2,w;

	a1 = 2*(p2.x() - p4.x());
	b1 = 2*(p2.y() - p4.y());
	c1 = ( pow( p4.x() , 2) + pow( p4.y() , 2) ) - ( pow( p2.x() , 2) + pow( p2.y() , 2) );

	a2 = 2*(p1.x() - p3.x());
	b2 = 2*(p1.y() - p3.y());
	c2 = ( pow( p3.x() , 2) + pow( p3.y() , 2) ) - ( pow( p1.x() , 2) + pow( p1.y() , 2) );

	g = (a1 * b2) - (b1 * a2);
	nx = (b1 * c2) - (b2 * c1);
	ny = (a2 * c1) - (a1 * c2);

	den = 1/g;
	xannulus = nx*den;
	yannulus = ny*den;

	c = Point( xannulus,yannulus );

	r2 = pow( (p2.x() - xannulus) , 2 ) + pow( (p2.y() - yannulus), 2 );
	R2 = pow( (p3.x() - xannulus), 2 ) + pow( (p3.y() - yannulus) , 2);
	w = (r2 + R2 - 2 * sqrt(r2 * R2))/4;

	r = sqrt(r2)+sqrt(w);
}

void boundary_circle_properties(Point &v1, Point &v2, Point &v3, Point &c, double &r){
	//cout<<"boundary circle"<<endl;
	Vector va = v2 - v3;
	Vector nr = Vector( -va.y() , va.x() );
	nr = unit_vector(nr);
	Vector vv = v1 - v3;
	c = Point( (v2.x() + v3.x() )/2 , (v2.x() + v3.x() )/2 );
	r = (nr*vv)/2;
}

int getNearPoint(Point &p0,  Point &p1 , Point &pc, double &r , double &t){
	Vector d = p1 - p0;
	Vector f = p0 - pc;
	double a = d*d;
	double b = 2*f*d;
	double c = f*f- r*r;
	double discriminant = b*b-4*a*c;
	if( discriminant < 0 )
		return 0;
	else{
		discriminant = sqrt( discriminant );
		t = (-b - discriminant)/(2*a);
		//cout<<"t: "<<t<<endl;
		if( t >= 0 && t <= 1 ){
			return 1;
		}
		return 0;
	}
}

int getDisplacedPoint(Delaunay &dt, Vertex &v, Point &p1, vector<pair<Point,double>> &circles){
	double t = 1000;
	Edge_circulator ec_start = dt.incident_edges(v);
	Edge_circulator ec = ec_start;
	do{
		Face c1 = ec->first;
		Face c2 = c1->neighbor(ec->second);
		///if( dt.is_infinite(c1) || dt.is_infinite(c2)) continue;
		Vertex v1 = c1->vertex( ec->second );
		Vertex v2 = c1->vertex( (ec->second+1)%3 );
		Vertex v3 = c1->vertex( (ec->second+2)%3 );
		Vertex v4 = dt.mirror_vertex( c1,ec->second );

		double ta;Point center;double radio;
		if( dt.is_infinite(c1) )
			boundary_circle_properties(v4->point(),v2->point(),v3->point(),center,radio);
		else if( dt.is_infinite(c2) )
			boundary_circle_properties(v1->point(),v2->point(),v3->point(),center,radio);
		else
			circle_properties(v1->point(), v3->point(), v4->point(), v2->point(), center, radio);

		int step = getNearPoint(v->point(),p1,center,radio,ta);

		//circles.push_back( make_pair (center,radio) );
		if(step){
			/*cout<<center<<" "<<radio<<endl;
			cout<<ta<<endl;*/
			if(ta < t)  t = ta;
		}
		++ec;
	}while (ec != ec_start);

	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;
	do{
		if(dt.is_infinite(fc)) {++fc;continue;}
		Face c1 = fc;
		Face c2 = c1->neighbor( c1->index(v) );
		//if(dt.is_infinite(c2)) continue;

		int indexv = c1->index(v);
		Vertex v1 = v;
		Vertex v2 = c1->vertex( (indexv+1)%3 );
		Vertex v3 = c1->vertex( (indexv+2)%3 );
		Vertex v4 = dt.mirror_vertex( c1,indexv );
		double ta;Point center;double radio;
		
		if( dt.is_infinite(c2) )
			boundary_circle_properties(v1->point(),v2->point(),v3->point(),center,radio);
		else
			circle_properties(v1->point(), v3->point(), v4->point(), v2->point(), center, radio);
		int step = getNearPoint(v->point(),p1,center,radio,ta);

		circles.push_back( make_pair (center,radio) );
		//cout<<"center and radio: "<<center<<" "<<radio<<endl;
		if(step){
			/*cout<<center<<" "<<radio<<endl;
			cout<<ta<<endl;*/
			if(ta < t)  t = ta;
		}

		++fc;
	}while (fc!= fc_start);

	//cout<<"t: "<<t<<endl;
	if(t > 0 && t != 1000){		
		//t = 0.7*t;
		p1 = Point( v->point().x() + t*(p1.x() - v->point().x()) , v->point().y() + t*(p1.y() - v->point().y()) );
		return 1;
	}
	else{
		p1 = v->point();
		return 0;
	}
}