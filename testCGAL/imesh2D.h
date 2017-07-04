#ifndef IMESH2D_H
#define IMESH2D_H

#include <iostream>
#include <vector>


struct Cposition{
	double x,y;
};
//template
struct CCell2D{	//CName 
	//XType
	//mesh generation
	//mesh segmentation
	//ColorPattern 
	unsigned char m_color[3]; 
	//smoothing mesh
	int m_label;//m_segementation_label unsigned long
	double m_minAngle; 

	//new version
	bool m_auxiliar;
	int m_state;//the cell has been relabeled = 1


	bool m_visited;

	int m_original_label;
	bool m_path_state;

	bool m_hole_state;
	bool m_hole_visited;

	bool m_relabeled;

	int m_solution;

	//for overlapping
	int m_previous_label;
	bool m_overlap_state;
	
	CCell2D(){
		m_label=-1;
		m_color[0]=255;
		m_color[1]=255;
		m_color[2]=255;

		//new version
		m_auxiliar = false;
		m_state=0;
		
		m_visited = false;

		m_path_state = false;
		m_hole_state = false; //is not a hole
		m_hole_visited = false;

		m_relabeled = false;

		m_original_label=-1;

		m_previous_label = -1;
		m_overlap_state = false;
	}
};

struct CVertex2D{
	//smoothing mesh
	double m_boundary_angle;//m_boundary_angle
	double m_ratio;//m_boundary_angle
	int m_index;//unsigned long
	unsigned int m_state;//0 -> normal 1->part of boundary 2-> divided and ordered
	bool m_limit;//point in th order
	bool m_restricted; //intersection point between more than 2 labels

	//moving points
	double m_tolerance;
	Cposition m_fixed;
	Cposition m_moving;
	std::vector<double> m_bicellTolerances;
	double angle;

	//new version 
	bool m_boundary;
	bool m_auxiliar;

	int m_type;
	bool is_corner()		const { return m_type == 0; }
	bool is_internal()	const { return m_type == 1; }
	bool is_boundary()	const { return m_type == 2 || m_regions_around.size() > 1; }
	

	vector<int> m_regions_around;
	bool has_two_regions()		const { return m_regions_around.size() == 2; }
	bool has_more_than_two_regions()		const { return m_regions_around.size() > 2; }

	vector<int> m_singular_regions_around;
	bool is_singular()		const { return m_singular_regions_around.size() > 0; }


	//new function
	

	CVertex2D(){
		m_limit = 0;
		m_state = 0;
		m_restricted = 0;
		m_boundary_angle = 180;//  ideal angle
		m_ratio = -1;

		//new version
		m_type = -1;
		m_boundary= false;
		m_auxiliar = false;
	}
};

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/point_generators_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_face_base_with_info_2< CCell2D,K> Fb;

template < class GT, class Vb = CGAL::Triangulation_vertex_base_2<GT> >
class My_vertex_base
  : public Vb
{

public:
  typedef typename Vb::Vertex_handle  Vertex_handle;
  typedef typename Vb::Face_handle    Face_handle;
  typedef typename Vb::Point          Point;

  struct Info{
		double m_boundary_angle;//m_boundary_angle
		double m_ratio;//m_boundary_angle
		int m_index;//unsigned long
		unsigned int m_state;//0 -> normal 1->part of boundary 2-> divided and ordered
		bool m_limit;//point in th order
		bool m_restricted; //intersection point between more than 2 labels

		//moving points
		double m_tolerance;
		Cposition m_fixed;
		Cposition m_moving;
		std::vector<double> m_bicellTolerances;
		double angle;

		//new version 
		bool m_boundary;
		bool m_auxiliar;
		Point m_smoothed_point;

		int m_type;
		bool is_corner()		const { return m_type == 0; }
		bool is_internal()	const { return m_type == 1; }
		bool is_boundary()	const { return m_type == 2 || m_regions_around.size() > 1; }
	

		vector<int> m_regions_around;
		bool has_two_regions()		const { return m_regions_around.size() == 2; }
		bool has_more_than_two_regions()		const { return m_regions_around.size() > 2; }

		vector<int> m_singular_regions_around;
		bool is_singular()		const { return m_singular_regions_around.size() > 0; }

		vector<pair<int,int>> m_singular_regions_and_subgraphs;
		bool is_singular_2()		const { return m_singular_regions_and_subgraphs.size() > 0; }	

		bool m_feature;
		std::vector<vector<Vertex_handle>> m_adyacent_vertices;

		Vertex_handle m_vertex_in_other_mesh;
		Info(){
			m_limit = 0;
			m_state = 0;
			m_restricted = 0;
			m_boundary_angle = 180;//  ideal angle
			m_ratio = -1;

			//new version
			m_index = -1;
			m_type = -1;
			m_boundary= false;
			m_auxiliar = false;
			m_feature = false;
		}
	};

   Info _info;

  template < class TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other  Vb2;
    typedef My_vertex_base<GT, Vb2>                        Other;
  };
  My_vertex_base() {}
  My_vertex_base(const Point& p)
    : Vb(p) {}
  My_vertex_base(const Point& p, Face_handle c)
    : Vb(p, c) {}

  const Info& info() const { return _info; }
  Info&       info()       { return _info; }
};

//typedef CGAL::Triangulation_vertex_base_with_info_2<CVertex2D, K>    Vb;
//typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;

typedef CGAL::Triangulation_data_structure_2<My_vertex_base<K>,Fb>    Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;

typedef K::Point_2 Point;//template  mesh generation (traits)
typedef K::Vector_2 Vector;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::Edge_iterator Edge_iterator;
typedef Delaunay::Vertex_iterator Vertex_iterator;
typedef Delaunay::Finite_vertices_iterator Finite_vertices_iterator;
typedef Delaunay::Vertex_handle Vertex;
typedef Delaunay::Edge Edge;
typedef Delaunay::Face_handle Face;
typedef Delaunay::Vertex_circulator Vertex_circulator;
typedef Delaunay::Face_circulator Face_circulator;
typedef Delaunay::Edge_circulator Edge_circulator;
typedef Delaunay::Finite_edges_iterator Finite_edges_iterator;
#endif