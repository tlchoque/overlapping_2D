
#include "headers.h"
#include "overlap.h"

//#include "test_integral.h"


using namespace std;

const int windowWidth = 1600;//original = 1600
const int windowHeight = 1200;//original=1200
const int numPoints = 3000;//3000

int main(){
	Delaunay dt;
	begin_VTK(dt);	
	/*generation(dt,numPoints,windowWidth,windowHeight);
	segmentation(dt,windowWidth,windowHeight);*/

	std::stringstream ss,original,manifold,relabeling,smoothing,manifold_histo,relabeling_histo,manifold_insertion_histo,
		smoothing_histo,smoothing_rel,smoothing_rel_histo,manifold_volume,manifold_edge,manifold_simulated,manifold_insertion;
	string file = "tweety";

	ss<<"D:/project/samples/meshes2D/"<<file<<".txt";
	original<<file<<"/original";
	manifold<<file<<"/manifold";
	relabeling<<file<<"/relabeling";
	smoothing<<file<<"/smoothing";
	smoothing_rel<<file<<"/smoothing_rel";
	manifold_volume<<file<<"/manifold_volume";
	manifold_edge<<file<<"/manifold_edge";
	manifold_simulated<<file<<"/manifold_simulated";
	manifold_insertion<<file<<"/manifold_insertion";

	manifold_histo<<file<<"/manifold/histo.txt";
	relabeling_histo<<file<<"/relabeling/histo.txt";
	smoothing_histo<<file<<"/smoothing/histo.txt";
	smoothing_rel_histo<<file<<"/smoothing_rel/histo.txt";
	manifold_insertion_histo<<file<<"/manifold_insertion/histo.txt";


	ifstream input(ss.str(), std::ifstream::in);
	input>>(input,dt);	

	/*Delaunay dt_pre;
	ifstream input2(ss.str(), std::ifstream::in);
	input2>>(input2,dt_pre);
	Delaunay dt_pre = dt;*/

	/*double result = atan2(5, -10)- atan2(5,-5);
	result *= 180 / PI;
	cout<<result<<endl;
	getchar();*/

	
	double old = dt.number_of_vertices();
	begin_overlapping(dt);
	Delaunay old_dt = dt;
	
	repairing(dt);
	simulated_annealing_whithout_points(dt);

	drawMesh(dt);
	getchar();

	/*get_inserted_points(dt,old,"tauro/manifold/new_points2",2);
	get_inserted_points(dt,old,"tauro/manifold/new_points3",3);
	get_inserted_points(dt,old,"tauro/manifold/new_points4",4);
	get_inserted_points(dt,old,"tauro/manifold/new_points5",5);
	get_inserted_points(dt,old,"tauro/manifold/new_points6",6);
	drawMesh(dt);
	getchar();*/	

	simliarty(old_dt, dt, 4, (manifold.str()).c_str(), "", 0);
	simliarty(dt,old_dt,4, ( manifold.str() ).c_str(),"rep2",0);
	

	end_overlapping(old_dt,dt,old);
	
	drawMesh(dt);
	getchar();


	/*drawMesh(dt);
	getchar();*/

	//saveMesh(dt,( manifold.str() ).c_str(),"rep1");
	/*get_inserted_points(dt,old,"titicaca/manifold/new_points",2);
	get_inserted_points(dt,old,"titicaca/manifold/new_points",3);
	get_inserted_points(dt,old,"titicaca/manifold/new_points",4);
	get_inserted_points(dt,old,"titicaca/manifold/new_points",5);*/

	drawMesh(dt);
	getchar();
	
	/*end_overlapping(old_dt,dt,old_num_vertices);
	drawMesh(dt);
	getchar();*/

	/*saveMesh(dt,( manifold.str() ).c_str(),"manifold");
	end_verification(dt, ( manifold.str() ).c_str() );	
	get_inserted_points(dt,old,"thundercats/manifold/new_points.vtk");*/

	/*set_adyacent_vertices(dt);	
	mesh_boundary_histogram(dt,( manifold_histo.str() ).c_str() );*/
		

	begin_verification(dt);
	relabeling_to_smooth(dt);
	/*end_verification(dt, ( relabeling.str() ).c_str() );
	saveMesh(dt,( relabeling.str() ).c_str(),"relabeling"  );*/
	//set_adyacent_vertices(dt);
	//mesh_boundary_histogram(dt,( relabeling_histo.str() ).c_str() );
		
	/*drawMesh(dt);
	getchar();*/

	//set_adyacent_vertices(dt);
	TSUR(dt,4);
	//mesh_boundary_histogram(dt,( smoothing_histo.str() ).c_str() );
	//saveMesh(dt, ( smoothing.str() ).c_str(),"smoothing" );	
		
	drawMesh(dt);
	getchar();
	return 0;
}

