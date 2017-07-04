#include <iostream>
#include <time.h>

void randomize(){
	//srand((time_t) time(NULL));  
}

double norm(Vector &v){
	return sqrt( v*v );
}

Vector unit_vector(Vector &v){
	double n = norm(v);
	return v/n;
}

double euclidian(Point a, Point b){
	Vector n = a - b ;
	return sqrt(n*n);
}