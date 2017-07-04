
double wind(Vector a, Vector b) {
	double angle2 = get_angle(a, b);
	if (is_obtuse(a, b))
		angle2 = -angle2;
	//return  angle2* 180 / PI;
	//return  angle2/ (2*PI);
	return angle2;
}

double w_point( Point p, vector<Edge> &facets) {
	double total = 0;
	for (unsigned int i = 0; i < facets.size(); ++i) {
		Edge e = facets[i];
		Point c_i = e.first->vertex(e.first->ccw(e.second))->point();
		Point c_j = e.first->vertex(e.first->cw(e.second))->point();

		if (p == c_i || p == c_j) { 
			continue; 
			total += 0.5;
		}

		Vector a = c_i - p;
		Vector b = c_j - p;
		double angle = wind(a, b);
		total += angle;
	}
	//return total / (2 + PI);
	return total;
}


double w_face(Delaunay &dt,Face f, Vertex v, int n , vector<Edge> &facets) {
	int idx = f->index(v);
	Vector b = f->vertex(f->cw(idx))->point() - v->point();
	Vector c = f->vertex(f->ccw(idx))->point() - v->point();
	b = (1 / ( (double)n - 1) )*b;
	c = (1 / ( (double)n - 1) )*c;
	//cout << "v and b and c: " <<v->point()<<" "<< b << " " << c << endl;
	vector<vector<Point> > samples(n, vector<Point>(n) );
	vector<vector<double> > winding( n, vector<double>(n) );

	vector<Point> pts;
	for ( int i = 0; i < n; ++i) {
		for (int j = 0; j < n - i; ++j) {
			samples[i][j] = v->point() + i*b + j*c;
			winding[i][j] = w_point(samples[i][j], facets);

			/*cout << "winding[i][j]: " << winding[i][j] << endl;
			pts.push_back(samples[i][j]);
			cout << samples[i][j] << endl;
			drawMesh(dt, pts);
			getchar();*/
		}
	}
	/*cout << "ola k ase" << endl;
	drawMesh(dt, pts);
	getchar();*/

	double a_w, b_w, c_w, total = 0;
	for ( int i = 0, imax = n - 1; i<imax; i++) {
		for ( int j = 0, jmax = imax - i; j<jmax; j++) {
			a_w = winding[i][j];
			b_w = winding[i][j + 1];
			c_w = winding[i + 1][j];
			total += b_w + b_w + c_w;
		}
	}

	for ( int i = 1; i<n; i++) {
		for ( int j = 1, jmax = n - i; j<jmax; j++) {
			a_w = winding[i - 1][j];
			b_w = winding[i][j - 1];
			c_w = winding[i][j];
			total += b_w + b_w + c_w;
		}
	}	
	//cout << "total: " << total << endl;
	return total/ ((n - 1)*(n - 1) * 3);
}

vector<Edge> extract_facets(vector<vector<Face>> &faceSets,Vertex v) {
	vector<Edge> facets;
	for (unsigned int i = 0; i < faceSets.size(); ++i) {
		for (unsigned int j = 0; j < faceSets[i].size(); ++j) {
			Face f = faceSets[i][j];
			int idx = f->index(v);
			for (unsigned int j = 1; j < 3; ++j) {
				int iedge = (idx + j) % 3;
				if (f->neighbor(iedge)->info().m_label != f->info().m_label) {
					Edge e = Edge(f, iedge);
					facets.push_back(e);
				}
			}
		}
	}
	return facets;
}

double triangle_w(Delaunay &dt,Face f, vector<Edge> &facets){
	Point p = CGAL::centroid(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
	double total = 0;
	for (unsigned int i = 0; i < facets.size(); ++i) {
		Edge e = facets[i];
		Point c_i = e.first->vertex(e.first->ccw(e.second))->point();
		Point c_j = e.first->vertex(e.first->cw(e.second))->point();
		Vector a = c_i - p;
		Vector b = c_j - p;

		double angle = wind(a,b);
		//double angle = atan2( a.y(), a.x() ) - atan2( b.y(), b.x() ) * 180 / PI;

		//double angle = ( atan2( b.y(), b.x() ) - atan2( a.y(), a.x() ) )* 180 / PI;	
		//double area = CGAL::area( f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point() );
		total += angle;

		/*vector<Point> pts;
		pts.push_back(c_i);
		pts.push_back(c_j);
		pts.push_back(p);
		cout << "b y a: " << b.x()<<" "<<b.y()<<" * "<< a.x()<<" "<<a.y() << endl;
		cout << "angle: " << angle << endl;
		drawMesh(dt,pts);
		getchar();*/
	}
	double area = CGAL::area( f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point() );
	//return total*area/(2+PI);
	return total/(2 + PI);

	/*vector<Point> pts;
	cout<<"total: "<<total<<endl;
	pts.push_back(p);
	drawMesh(dt,pts);
	getchar();*/

}

double criteria_w(Delaunay &dt, vector<Face> &faces, vector<Edge> &facets, Vertex v ) {
	double w_average = 0;
	for (unsigned int i = 0; i < faces.size(); ++i) {
		Face f = faces[i];
		//double w = triangle_w(dt, faces[i], facets); 
		double w = w_face( dt, f,v,4, facets );
		double area = CGAL::area(f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
		w_average += w * area;
		//w_average += w;

		/*vector<Point> pts;
		pts.push_back( f->vertex(0)->point() );
		pts.push_back( f->vertex(1)->point());
		pts.push_back( f->vertex(2)->point());
		cout << "w face and w_average: " << w <<" * "<<w_average<< endl;
		drawMesh(dt, pts);
		getchar();*/
	}
	/*vector<Point> pts;
	pts.push_back(faces[0]->vertex(0)->point() );
	pts.push_back(faces[0]->vertex(1)->point());
	pts.push_back(faces[0]->vertex(2)->point());
	cout << "value w averge: " << w_average << endl;
	drawMesh(dt, pts);
	getchar();*/

	//w_average *=  10;

	//return w_average / faces.size();

	//return faces.size() / w_average;
	return 1/w_average;
}

void sort_by_criteria_w(Delaunay &dt, vector<vector<Face>> &faceSets, Vertex v, int &min, int &max, bool &changeable) {
	vector<Edge> facets = extract_facets(faceSets, v);
	vector<double> ratio(faceSets.size());
	for (unsigned int i = 0; i < faceSets.size(); ++i) {
		ratio[i] = criteria_w(dt, faceSets[i], facets, v);
		//cout << "ratio: " << ratio[i] << endl;
	}
	max = 0; min = 1;
	if (ratio[min] > ratio[max]) {
		max = 1;
		min = 0;
	}
	if (ratio[max] / ratio[min] < 1)//for chest is different 2.65 3
		changeable = true;
	else
		changeable = false;
}

