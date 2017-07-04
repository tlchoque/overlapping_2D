//#include <iostream>
//#include <math.h>
//#include "imesh2D.h"
//#include "additional.h"

int bigger(double x,double y,double s,double b){
	double y2=s*x+b;
	if(y > y2)	return 1;
	else return 0;
}

void segmentedCircle(Delaunay & dt,double xc,double yc, double r){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		double r2 =r*r;
		if( (pow(dt.triangle(it)[0].x()-xc, 2)+pow(dt.triangle(it)[0].y()-yc , 2)  < r2 ) && 
			(pow(dt.triangle(it)[1].x()-xc, 2)+pow(dt.triangle(it)[1].y()-yc , 2)  < r2  ) && 
			(pow(dt.triangle(it)[2].x()-xc, 2)+pow(dt.triangle(it)[2].y()-yc , 2)  < r2  ) ){
			it->info().m_label = 1;
		}	
		else it->info().m_label = 0;
	}
}

void segmentedSquare(Delaunay & dt, double x, double y,double side){
	double left,right,up,down;
	left=x;right=x+side;up=y;down=y-side;
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		if( (dt.triangle(it)[0].x() < right && dt.triangle(it)[0].x() > left && dt.triangle(it)[0].y() < up && dt.triangle(it)[0].y() > down) &&
			(dt.triangle(it)[1].x() < right && dt.triangle(it)[1].x() > left && dt.triangle(it)[1].y() < up && dt.triangle(it)[1].y() > down) &&
			(dt.triangle(it)[2].x() < right && dt.triangle(it)[2].x() > left && dt.triangle(it)[2].y() < up && dt.triangle(it)[2].y() > down) ){
			it->info().m_label = 2;
		}
	}
}

void segmentedLine(Delaunay & dt, double s, double b){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		if( bigger(dt.triangle(it)[0].x(),dt.triangle(it)[0].y(),s,b) && 
			bigger(dt.triangle(it)[1].x(),dt.triangle(it)[1].y(),s,b) && 
			bigger(dt.triangle(it)[2].x(),dt.triangle(it)[2].y(),s,b) ){
			it->info().m_label = 3;
		}
	}
}

void segmentation(Delaunay & dt, int windowWidth, int windowHeight){
	double xc,yc,r; // center (xc,yc) and radio r
	double x,y,side; // left upper point of square (x,y) and area of square = side^2

	r=(double)( rand()%(int)floor(windowHeight/4)   ) + (windowHeight/6);
	xc=(double)( rand()%windowWidth );
	yc=(double)( rand()%windowHeight );
	segmentedCircle(dt,xc,yc,r);
	
	side=(double)( rand()%(int)floor(windowHeight/4)   ) + (windowHeight/3);
	x=(double)(rand()%windowWidth - (int)(side/2) );
	y=(double)(rand()%windowHeight + (int)(side/2) ) ;	
	segmentedSquare(dt,x,y,side);
}