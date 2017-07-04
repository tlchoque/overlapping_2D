

void segmentedCircle(Delaunay & dt,double xc,double yc, double r){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		double r2 =r*r;
		if( (pow(dt.triangle(it)[0].x()-xc, 2)+pow(dt.triangle(it)[0].y()-yc , 2)  < r2 ) && 
			(pow(dt.triangle(it)[1].x()-xc, 2)+pow(dt.triangle(it)[1].y()-yc , 2)  < r2  ) && 
			(pow(dt.triangle(it)[2].x()-xc, 2)+pow(dt.triangle(it)[2].y()-yc , 2)  < r2  ) ){
			it->info().m_label = 1;
			it->info().m_color[0]=255;
			it->info().m_color[1]=0;
			it->info().m_color[2]=0;
		}	
		else {
			it->info().m_label = 0;
			it->info().m_color[0]=192;
			it->info().m_color[1]=192;
			it->info().m_color[2]=192;
		}
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
			it->info().m_color[0]=0;
			it->info().m_color[1]=0;
			it->info().m_color[2]=255;
		}
		//else it->info().m_label = 0;// for test manifold
	}
}

void segmentedLine(Delaunay & dt, double b, double c){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		if( (dt.triangle(it)[0].y() - b*dt.triangle(it)[0].x() - c  < 0 ) && 
			(dt.triangle(it)[1].y() - b*dt.triangle(it)[1].x() - c  < 0 ) && 
			(dt.triangle(it)[2].y() - b*dt.triangle(it)[2].x() - c  < 0 ) ){
			it->info().m_label = 3;
			it->info().m_color[0]=0;
			it->info().m_color[1]=153;
			it->info().m_color[2]=0;
		}
	}
}

bool checkSine(double y, double x , double a, double p){
	/*double ysin = 2*a + a*sin(x/p);
	if(y < ysin) return 1;
	else return 0;*/
	double ysin = 1200 - 2*a + a*sin(x/p);
	if(y > ysin) return 1;
	else return 0;
}

void segmentedSine(Delaunay & dt, double a, double p){
	for(Finite_faces_iterator it = dt.finite_faces_begin(); it != dt.finite_faces_end(); it++){
		if( checkSine(dt.triangle(it)[0].y(),dt.triangle(it)[0].x(),a,p) && 
			checkSine(dt.triangle(it)[1].y(),dt.triangle(it)[1].x(),a,p) && 
			checkSine(dt.triangle(it)[2].y(),dt.triangle(it)[2].x(),a,p) ){
			it->info().m_label = 4;
			it->info().m_color[0]=153;
			it->info().m_color[1]=0;
			it->info().m_color[2]=153;
		}
	}
}

void segmentation(Delaunay & dt, int windowWidth, int windowHeight){
	double xc,yc,r; // center (xc,yc) and radio r
	r=(double)( rand()%(int)floor(windowHeight/4)   ) + (windowHeight/6);
	xc=(double)( rand()%windowWidth );
	yc=(double)( rand()%windowHeight );
	r=300;
	xc=500;
	yc=450;
	segmentedCircle(dt,xc,yc,r);
	
	double x,y,side; // left upper point of square (x,y) and area of square = side^2
	side=(double)( rand()%(int)floor(windowHeight/4)   ) + (windowHeight/3);
	x=(double)(rand()%windowWidth - (int)(side/2) );
	y=(double)(rand()%windowHeight + (int)(side/2) );	
	x=900;
	y=950;
	side=580;

	// for non manifold case
	/*x=200;
	y=750;
	side=600;*/
	segmentedSquare(dt,x,y,side);

	double b,c;
	b=((double) rand() / (RAND_MAX));
	//c=(double)( rand()%(windowHeight/2) ) -  windowHeight/4;
	c=windowHeight - (double)( rand()%(windowHeight/4) ) -  windowHeight/4;
	b=0.18;
	c=0;
	segmentedLine(dt,b,c);

	double amplitude,period;
	amplitude=80 + (double)(rand()%80);
	period= 80 + (double)(rand()%100);

	amplitude = 100;
	period = 90;
	segmentedSine(dt,amplitude,period);
}