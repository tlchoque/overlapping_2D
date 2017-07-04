#include <CGAL/basic.h>
#include <map>
#include <set>
#include <CGAL/Unique_hash_map.h>

using namespace std;

std::ostream& operator<<(std::ostream& os, Delaunay &dt){
	Tds::size_type n = dt.tds().number_of_vertices();
	Tds::size_type m = dt.tds().number_of_full_dim_faces(); 
	
	if(CGAL::is_ascii(os))  os << n << ' ' << m << ' ' << dt.tds().dimension() << std::endl;
	else     os << n << m << dt.tds().dimension();
	if (n==0) return os;

	CGAL::Unique_hash_map<Vertex,int> V;
	CGAL::Unique_hash_map<Face,int> F;

	// first vertex 	
	int inum = 0;
	Vertex v = dt.infinite_vertex();
	if ( v != Tds::Vertex_handle()) {
		V[v] = inum++;
	}

	for( Tds::Vertex_iterator vit= dt.tds().vertices_begin(); vit != dt.tds().vertices_end() ; ++vit) {
		if ( v != vit ) {
			V[vit] = inum++;
			os << vit->point();
			if(CGAL::is_ascii(os)) os << "\n";
		}
	}	
	if(CGAL::is_ascii(os)) os << "\n";

	inum = 0;
	int dim = (dt.tds().dimension() == -1 ? 1 :  dt.tds().dimension() + 1);
	for( Tds::Face_iterator ib = dt.tds().face_iterator_base_begin();ib != dt.tds().face_iterator_base_end(); ++ib) {
		F[ib] = inum++;
		for(int j = 0; j < dim ; ++j) {
			os << V[ib->vertex(j)];
			if(CGAL::is_ascii(os)) os << " ";
		}
		//os << *ib ;		
		os << (double)(ib->info().m_color[0])/255<<" "<<(double)ib->info().m_color[1]/255<<" "<<(double)ib->info().m_color[2]/255<<" "<<ib->info().m_label<<" ";
		if(CGAL::is_ascii(os)) os << "\n";
	}
	if(CGAL::is_ascii(os)) os << "\n";

	for( Tds::Face_iterator it = dt.tds().face_iterator_base_begin();it != dt.tds().face_iterator_base_end(); ++it) {
		//std::cout<<cen<<" * ";
		for(int j = 0; j < dt.dimension()+1; ++j){
			//std::cout<<F[it->neighbor(j)]<<" ";
			os << F[it->neighbor(j)];
			if(CGAL::is_ascii(os))  os << " ";
		}
		//std::cout<<std::endl;
		if(CGAL::is_ascii(os)) os << "\n";
	}
	return os;
}


std::istream& operator>>(std::istream& is, Delaunay &dt){
	if(dt.tds().number_of_vertices() != 0)    dt.tds().clear();  
	Tds::size_type n, m;
	int d;
	is >> n >> m >> d;
	if (n==0){ return is;}
	dt.tds().set_dimension(d);

	std::vector<Vertex  > V(n);
	std::vector<Face> F(m);

	Tds::size_type i = 0;
	V[0] =dt.tds().create_vertex();
	++i;
	for( ; i < n; ++i) {
		V[i] = dt.tds().create_vertex();
		is >> *(V[i]);
	}
  
	int index;
	int dim = (dt.tds().dimension() == -1 ? 1 :  dt.tds().dimension() + 1);
	
	double r,g,b;int label;
	for(i = 0; i < m; ++i) {
		F[i] = dt.tds().create_face() ;
		for(int j = 0; j < dim ; ++j){
			is >> index;
			F[i]->set_vertex(j, V[index]);
			V[index]->set_face(F[i]);
		}
		is>>r>>g>>b>>label;
		F[i]->info().m_color[0] = (int)(255*r);
		F[i]->info().m_color[1] = (int)(255*g);
		F[i]->info().m_color[2] = (int)(255*b);
		F[i]->info().m_label = label;
		is >> *(F[i]) ;
	}
		
	for(i = 0; i < m; ++i) {
		for(int j = 0; j < dt.tds().dimension()+1; ++j){
			is >> index;
			F[i]->set_neighbor(j, F[index]);
		}
	}
	 
	dt.set_infinite_vertex(V[0]);
	return is;
}