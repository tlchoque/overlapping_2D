#ifndef IMESH2D_H
#define IMESH2D_H

//template
struct  CCell2D{	//CName 
	//XType
	//mesh generation
	 
	//mesh segmentation

	//ColorPattern 

	//smoothing mesh
	int m_label;//m_segementation_label unsigned long
	double m_minAngle; 
};

struct CVertex2D{
	double m_boundary_angle;//m_boundary_angle
	int m_index;//unsigned long
};

#include <vector>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/Color.h>
#include <CGAL/int.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/point_generators_2.h>


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_with_info_2<CVertex2D, K>    Vb;
typedef CGAL::Triangulation_face_base_with_info_2< CCell2D,K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<K,Tds> Delaunay;

typedef K::Point_2 Point;//template  mesh generation (traits)
typedef K::Vector_2 Vector;
typedef Delaunay::Finite_faces_iterator Finite_faces_iterator;
typedef Delaunay::Edge_iterator Edge_iterator;
typedef Delaunay::Vertex_iterator Vertex_iterator;
typedef Delaunay::Vertex_handle Vertex;
typedef Delaunay::Face_handle Face;
typedef Delaunay::Vertex_circulator Vertex_circulator;
typedef Delaunay::Face_circulator Face_circulator;

#endif