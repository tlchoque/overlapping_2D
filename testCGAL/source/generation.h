#include <iostream>
#include "additional.h"

using namespace std;

void generation(Delaunay & dt, int  numPoints, int windowWidth, int windowHeight){
	vector<Point> rPoints;
	int max = 10;
	Vector *p = new Vector[max]; 
	Vector *p_Temp = NULL;   
	int nv = 0;
	double x,y;
	bool ok = false;
	randomize();
	nv = 0;
	p = new Vector[max];
	while (nv != numPoints){
		do{
			ok = true;
			x = (double)(rand()%windowWidth);
			y = (double)(rand()%windowHeight);
			for(int n_Cpt = 0; n_Cpt <= nv; ++n_Cpt){
				if((x == p[n_Cpt].x() ) && (y == p[n_Cpt].y() )) {ok = false;break;}
			}
		}while(!ok);
		if (nv >= max){
			max = max * 2;            
			p_Temp = new Vector[max]; 
			for (int i = 0; i < nv; ++i)
				p_Temp[i] = p[i];  

			delete []p;  
			p = p_Temp; 
		}   
		p[nv] = Vector(x,y); 
		nv++;
	}
	
	for (int i = 0; i < nv; ++i)
		rPoints.push_back(Point(p[i].x() ,p[i].y() ));   
	delete []p_Temp;  

	dt.insert(rPoints.begin(),rPoints.end());
}