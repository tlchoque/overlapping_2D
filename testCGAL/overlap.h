#include <r3d.h>
#include <r2d.h>
#include <v3d.h>
#include <v2d.h>
#include <iostream>
#include <vector>

#include <tests/test_helpers.h>
#include <tests/utest.h>
#include <math.h>

#define TOL_WARN 1.0e-8
#define TOL_FAIL 1.0e-4

// minimum volume allowed for test polyhedra 
#define MIN_AREA 1.0e-8

// number of trials per test
#define NUM_TRIALS 100

// for recursive tests, the stack capacity and 
// maximum recursion depth allowed
#define STACK_SIZE 64
#define MAX_DEPTH 16

// order of polynomial integration for all tests 
#define POLY_ORDER 2

r2d_rvec2 point2vec(Point &p){
	r2d_rvec2 tmp;
	tmp.x = p.x();
	tmp.y = p.y();
	return tmp;
}

double overlap_cells( Face a, Face b){
	r2d_poly poly; 
	r2d_rvec2 verts[3];
	r2d_plane faces[3];
	vector<r2d_real> om(R2D_NUM_MOMENTS(POLY_ORDER));

	for(unsigned int i = 0 ; i < 3; ++i)
		verts[i] = point2vec( a->vertex(i)->point() );
	r2d_init_poly(&poly, verts, 3);

	// generate the second random tet
	for(unsigned int i = 0 ; i < 3; ++i)
		verts[i] = point2vec( b->vertex(i)->point() );

	r2d_poly_faces_from_verts(faces, verts, 3);
	// clip the first tet against the faces of the second
	r2d_clip(&poly, faces, 3);

	// find the moments (up to quadratic order) of the clipped poly
	r2d_reduce(&poly, om, POLY_ORDER);

	double a_area = CGAL::area(a->vertex(0)->point(),a->vertex(1)->point(),a->vertex(2)->point());
	if( om[0]>a_area){
		cout<<"om and real "<<om[0]<<" "<<a_area<<endl;
	}

	return om[0];
}

void begin_overlapping(Delaunay &dt){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++)
		it->info().m_previous_label = it->info().m_label;
}

vector<Face> involved_triangles(Delaunay &dt, Face c, vector<Vertex> &vertices){
	vector<Vertex> involved_vertices;
	for(unsigned int i = 0 ; i < vertices.size(); ++i){
		Vertex v = vertices[i];
		if( dt.side_of_oriented_circle(c,v->point(),false) == CGAL::ON_BOUNDED_SIDE ){
			involved_vertices.push_back(v);
		}
	}
	vector<Face> involved_cells;
	for(unsigned int i = 0 ; i < involved_vertices.size(); ++i){
		Vertex w = involved_vertices[i];

		Face_circulator fc_start = dt.incident_faces(w);
		Face_circulator fc = fc_start;
		do{
			if(dt.is_infinite(fc)) {++fc;continue;}

			if( fc->info().m_label == c->info().m_label && fc->info().m_overlap_state == false ){
				fc->info().m_overlap_state = true;
				involved_cells.push_back(fc);
			}
			++fc;
		}while (fc!= fc_start);
	}

	for(unsigned int i = 0 ; i < involved_cells.size(); ++i)
		involved_cells[i]->info().m_overlap_state = false;

	return involved_cells;
}

//Vertex is_vertex(Delaunay &dt,Point p){
//	for (Finite_vertices_iterator vi = dt.finite_vertices_begin(); vi != dt.finite_vertices_end(); vi++) { 
//		if( vi->point() == p)
//			return vi;
//	}
//}


void end_overlapping(Delaunay &dt_old, Delaunay &dt_new, double num_vertices){
	double overlapping=0,total=0;
	//int ct=0;
	//int dels=0;
	//for(Finite_faces_iterator it = dt_new.finite_faces_begin(); it != dt_new.finite_faces_end(); it++){
	//	total += CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point() );
	//	if( it->info().m_previous_label == it->info().m_label && it->info().m_previous_label != -1){
	//		if( it->info().m_label  == -1) cout<<"waaa"<<endl;
	//		
	//		Vertex v0 = is_vertex(dt_old,it->vertex(0)->point() );
	//		Vertex v1 = is_vertex(dt_old,it->vertex(1)->point() );
	//		Vertex v2 = is_vertex(dt_old,it->vertex(2)->point() );
	//		Face tt;
	//		if( !dt_old.is_face(v0,v1,v2,tt) ) {
	//			cout<<"problem cell: "<<dt_new.triangle(it)<<endl;
	//			cout<<it->info().m_previous_label<<" " <<it->info().m_label <<endl;
	//			vector<Point> pts;
	//			pts.push_back(it->vertex(0)->point() );
	//			pts.push_back(it->vertex(1)->point() );
	//			pts.push_back(it->vertex(2)->point() );
	//			drawMesh(dt_old,pts);
	//			getchar();
	//		}			++ct;
	//		overlapping += CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point()  );
	//	}
	//	else{
	//		++dels;
	//		/*vector<Point> pts;
	//		pts.push_back(it->vertex(0)->point() );
	//		pts.push_back(it->vertex(1)->point() );
	//		pts.push_back(it->vertex(2)->point() );
	//		drawMesh(dt_old,pts);
	//		getchar();*/
	
	//get the new vertices
	vector<Vertex> new_vertices;
	int count = 0;
	for (Finite_vertices_iterator vi = dt_new.finite_vertices_begin(); vi != dt_new.finite_vertices_end(); vi++) { 
		++count;
		if( count > num_vertices){
			new_vertices.push_back(vi);
		}
	}

	//relate the old and new vertices
	for (Finite_vertices_iterator v_old = dt_old.finite_vertices_begin(),  v_new = dt_new.finite_vertices_begin(); 
		v_old != dt_old.finite_vertices_end(); v_old++, v_new ++) { 
		v_old->info().m_vertex_in_other_mesh = v_new;
	}

	//get the deleted cells
	vector<Face> deleted_cells;
	//int tdel=0;
	for(Finite_faces_iterator it = dt_old.finite_faces_begin(); it != dt_old.finite_faces_end(); it++){
		total += CGAL::area(it->vertex(0)->point(), it->vertex(1)->point(), it->vertex(2)->point() );
		Face c_old = it;		
		Face c_new;
		Vertex v0 = c_old->vertex(0)->info().m_vertex_in_other_mesh;
		Vertex v1 = c_old->vertex(1)->info().m_vertex_in_other_mesh;
		Vertex v2 = c_old->vertex(2)->info().m_vertex_in_other_mesh;
		if(!dt_new.is_face(v0,v1,v2,c_new) ){
			deleted_cells.push_back(c_old);
		}
		else{
			//include the remaining cells  with probably different label
			if( c_new->info().m_label == c_old->info().m_label ){
				overlapping += CGAL::area(c_old->vertex(0)->point(), c_old->vertex(1)->point(), c_old->vertex(2)->point() );
			}
			//++tdel;
		}
	}

	/*cout<<"NEW: changed cells: "<<dels<<" ,faces"<<dt_new.number_of_faces()<<" ,keeped cells: "<<ct<<endl;
	cout<<"OLD: keeped cells "<<tdel<<" ,faces"<<dt_old.number_of_faces()<<" ,deleted cells"<<deleted_cells.size()<<endl;*/

	//get the overlapping
	/*double aux = overlapping;
	double area_deleted=0;*/
	for(unsigned int i = 0 ; i < deleted_cells.size(); ++i){
		Face del = deleted_cells[i];
		vector<Face> involved = involved_triangles(dt_new,del,new_vertices);
		//double tot=0;
		for(unsigned int j = 0 ; j < involved.size();++ j ){
			Face inv = involved[j];
			overlapping +=overlap_cells(del,inv);
			//tot +=overlap_cells(del,inv);
		}
		/*double del_area = CGAL::area(del->vertex(0)->point(), del->vertex(1)->point(), del->vertex(2)->point()  );
		area_deleted+=del_area;*/
		/*if(del_area < tot)
			cout<<"area overl: "<<del_area<<" "<<tot<<endl;*/
	}

	double difference = total - overlapping;
	cout<<"total area: "<<total<<endl;
	cout<<"partial area: "<<difference<<endl;
	cout<<"ratio: "<<difference/total<<endl;
}
