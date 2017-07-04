

Point centroid_edge(Delaunay &dt, Edge e){
	Vector a = dt.segment(e).source() - CGAL::ORIGIN;
	Vector b = dt.segment(e).target() - CGAL::ORIGIN;
	Vector c = (a + b)/2;
	return CGAL::ORIGIN + c;
}

bool path_to_join(Delaunay &dt, vector<vector<Face>> &faceSets,Vertex &v,vector<Face> &path){
	path.clear();
	if( faceSets.size() != 2 ) return false;	
	vector<Point> agroup,bgroup;
	vector<Face> acells,bcells;

	vector<Edge> afacets;
	for( unsigned int i = 0; i < faceSets[0].size(); ++i){// i can take jus the boundary cells
		Face c = faceSets[0][i];
		int label = c->info().m_label;
		int vindex = c->index(v);
		for( unsigned int j = 1; j < 3; ++j){
			int index = ( j + vindex )%3;
			Face nc = c->neighbor(index);
			if( !dt.is_infinite(nc) && nc->info().m_label != label){
				Edge f = Edge(c,index);
				afacets.push_back(f);
				Point cen = centroid_edge(dt,f);
				agroup.push_back(cen);
			}
		} 
	}

	for( unsigned int i = 0; i < faceSets[1].size(); ++i){// i can take just the boundary cells
		Face c = faceSets[1][i];
		bcells.push_back(c);
		int label = c->info().m_label;
		int vindex = c->index(v);
		for( unsigned int j = 1; j < 3; ++j){
			int index = ( j + vindex )%3;//index of c
			Face nc = c->neighbor(index);
			if( !dt.is_infinite(nc) && nc->info().m_label != label){
				Edge f = Edge(c,index);
				Point cen = centroid_edge(dt,f);
				bgroup.push_back(cen);
			}
		} 
	}
	double distance = 10000000000;
	int origin = 0, destiny = 0;
	for( unsigned int i = 0; i < agroup.size(); ++i){
		Vector a = agroup[i] - CGAL::ORIGIN;
		for( unsigned int j = 0; j < bgroup.size(); ++j){
			Vector b = bgroup[j] - CGAL::ORIGIN;
			Vector c = a - b;			
			double current = c*c;	
			if( current < distance ){
				distance = current;
				origin = i;
				destiny = j;
			}
		}
	}	

	Face c = afacets[origin].first;
	Face next = c;
	Point p = bgroup[destiny];

	bool is_path = false;
	bool try_next_cell = true;
	while(try_next_cell){
		try_next_cell = false;
		int cindex = c->index(v);
		const Point* pts[3] = { &(c->vertex(0)->point()),&(c->vertex(1)->point()),&(c->vertex(2)->point()) };

		for( unsigned int i = 1; !try_next_cell && i != 3; ++i){
			int index = ( i + cindex)%3;
			Face next = c->neighbor( index );
			if( next->info().m_path_state || dt.is_infinite(next) ) continue;// visited cell or infinite cell
			//if( std::find( bcells.begin(), bcells.end(), next) != bcells.end() )	{is_path = true; break;} //we arrive to destiny
			const Point* backup = pts[index];
			pts[index] = &p;
			if ( orientation(*pts[0], *pts[1], *pts[2]) != CGAL::NEGATIVE )
				pts[index] = backup;
			else{
				next->info().m_path_state = true;
				c = next;
				path.push_back(next);
				
				for( unsigned int j = 0 ; j < 3; ++j){
					Face nc = next->neighbor(j);
					if( std::find( bcells.begin(), bcells.end(), nc) != bcells.end() )	{is_path = true; } //we arrive to destiny
				}

				if(is_path)
					break;
				else
					try_next_cell = true;;
			}
		}
	}

	for( unsigned int i = 0; i < faceSets[0].size(); ++i)
		faceSets[0][i]->info().m_path_state = false;

	for( unsigned int i = 0; i < path.size(); ++i){
		if( !path[i]->info().m_path_state ){
			path.erase( path.begin() + i);
			--i;
		}
		else	path[i]->info().m_path_state = false;
	}

	if( is_path) return true;
	else	return false;
}