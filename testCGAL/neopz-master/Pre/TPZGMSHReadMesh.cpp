/**
 * @file
 * @brief Contains the implementation of the TPZGMSHReadMesh methods. 
 */

#include "TPZGMSHReadMesh.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzvec.h"

TPZGMSHReadMesh::TPZGMSHReadMesh(TPZGeoMesh *gmesh){
	fGeoMesh = gmesh;
}

void TPZGMSHReadMesh::ReadMesh2D(const char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist) {
	
	//o GMSH pode n�o retornar os n�s sequencialmente
	//nesse casso deveram ser resequenciados
	TPZStack<long> Indexes;
	Resequence(Indexes,meshfile);
	std::ifstream mesh(meshfile);
	char title[256];
	long nnodes,number;
	mesh >> title;//$NOD
	mesh >> nnodes;
	fGeoMesh->NodeVec().Resize(nnodes);   
	TPZVec<REAL> coord(3);
	long i;
	for(i=0;i<nnodes;i++){
		mesh >> number;
		mesh >> coord[0] >> coord[1] >> coord[2];
		fGeoMesh->NodeVec()[i].Initialize(coord,*fGeoMesh);
	}
	mesh >> title;//$ENDNOD
	mesh >> title;//$ELM
	long numelem;
	mesh >> numelem;
	TPZVec<long> nodes;
	nodes.Resize(4);
	long index;
	long numb,fgt1,nmat,fgt2,nvert,no1,no2,no3,no4;
	for(i=0;i<numelem;i++){
		mesh >> numb >> fgt1 >> nmat >> fgt2 >> nvert;
		mesh >> no1 >> no2;
		nodes[0] = Indexes[no1-1];//no1-1;
		nodes[1] = Indexes[no2-1];//no2-1;
		if(nvert == 2){//s� elementos de contorno
			nodes.Resize(2);
			elembclist.Push(TPZGeoElSide(fGeoMesh->CreateGeoElement(EOned,nodes,1,index),-nmat));
			continue;
		}
		if(nvert == 3){//s� elementos de volume
			nodes.Resize(3);
			mesh >> no3;
			nodes[2] = Indexes[no3-1];//no3-1;
			elemlist.Push(fGeoMesh->CreateGeoElement(ETriangle,nodes,1,index));
		}
		if(nvert == 4){//s� elementos de volume
			nodes.Resize(4);
			mesh >> no3 >> no4;
			nodes[2] = Indexes[no3-1];//no3-1;
			nodes[3] = Indexes[no4-1];//no4-1;
			elemlist.Push(fGeoMesh->CreateGeoElement(EQuadrilateral,nodes,1,index));
		}
	}
	mesh >> title;//$ENDELM
	mesh.close();
}

void TPZGMSHReadMesh::ReadMesh2D2(char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist){
	
	std::ifstream mesh(meshfile);
	char title[256];
	long nnodes,number;
	mesh >> title;//$NOD
	mesh >> nnodes;
	fGeoMesh->NodeVec().Resize(nnodes);   
	TPZVec<REAL> coord(3);
	long i;
	for(i=0;i<nnodes;i++){
		mesh >> number;
		mesh >> coord[0] >> coord[1] >> coord[2];
		fGeoMesh->NodeVec()[i].Initialize(coord,*fGeoMesh);
	}
	mesh >> title;//$ENDNOD
	mesh >> title;//$ELM
	long numelem;
	mesh >> numelem;
	TPZVec<long> nodes;
	nodes.Resize(3);
	long index;
	long numb,fgt1,nmat,fgt2,nvert,no1,no2,no3,no4;
	for(i=0;i<numelem;i++){
		mesh >> numb >> fgt1 >> nmat >> fgt2 >> nvert;
		mesh >> no1;
		mesh >> no2;
		nodes[0] = no1-1;
		nodes[1] = no2-1;
		if(nvert == 2){//s� elementos de contorno
			nodes.Resize(2);
			elembclist.Push(TPZGeoElSide(fGeoMesh->CreateGeoElement(EOned,nodes,1,index),-nmat));
			continue;
		}
		if(nvert == 3){//s� elementos de volume
			nodes.Resize(3);
			mesh >> no3;
			nodes[2] = no3-1;
			elemlist.Push(fGeoMesh->CreateGeoElement(ETriangle,nodes,1,index));
		}
		if(nvert == 4){//s� elementos de volume
			nodes.Resize(4);
			mesh >> no3 >> no4;
			nodes[2] = no3-1;
			nodes[3] = no4-1;
			elemlist.Push(fGeoMesh->CreateGeoElement(EQuadrilateral,nodes,1,index));
		}
	}
	mesh >> title;//$ENDELM
	mesh.close();
}

void TPZGMSHReadMesh::ReadMesh3D(char *meshfile,TPZStack<TPZGeoEl *> &elemlist,TPZStack<TPZGeoElSide> &elembclist){
	
	//o GMSH pode n�o retornar os n�s sequencialmente
	//nesse casso deveram ser resequenciados
	TPZStack<long> Indexes;
	Resequence(Indexes,meshfile);
	std::ifstream mesh(meshfile);
	char title[256];
	long nnodes,number;
	mesh >> title;//$NOD
	mesh >> nnodes;
	fGeoMesh->NodeVec().Resize(nnodes);   
	TPZVec<REAL> coord(3);
	long i;
	for(i=0;i<nnodes;i++){
		mesh >> number;
		mesh >> coord[0] >> coord[1] >> coord[2];
		fGeoMesh->NodeVec()[i].Initialize(coord,*fGeoMesh);
	}
	mesh >> title;//$ENDNOD
	mesh >> title;//$ELM
	long numelem;
	mesh >> numelem;
	TPZVec<long> nodes;
	nodes.Resize(3);
	long index;
	long numb,fgt1,nmat,fgt2,nvert;
	TPZVec<long> nos(8);//m�ximo do cubo
	for(i=0;i<numelem;i++){
		mesh >> numb >> fgt1 >> nmat >> fgt2 >> nvert;
		for(i=0;i<nvert;i++) mesh >> nos[i];
		if(nvert == 3){//tri�ngulos
			for(i=0;i<nvert;i++) nodes[i] = Indexes[nos[i]-1];//nos[i]-1;
			nodes.Resize(3);
			elembclist.Push(TPZGeoElSide(fGeoMesh->CreateGeoElement(ETriangle,nodes,1,index),-nmat));
			continue;
		} else if(nvert == 4 && fgt1 != 4){//quadrilateros
			nodes.Resize(4);
			for(i=0;i<nvert;i++) nodes[i] = Indexes[nos[i]-1];//nos[i]-1;
			elembclist.Push(TPZGeoElSide(fGeoMesh->CreateGeoElement(EQuadrilateral,nodes,1,index),-nmat));
			continue;
		} else if(nvert == 4 && fgt1 == 4){//tetraedros
			nodes.Resize(4);
			for(i=0;i<nvert;i++) nodes[i] = Indexes[nos[i]-1];//nos[i]-1;
			elemlist.Push(fGeoMesh->CreateGeoElement(ETetraedro,nodes,1,index));
			continue;
		} else if(nvert == 8){//hexaedros
			nodes.Resize(8);
			for(i=0;i<nvert;i++) nodes[i] = Indexes[nos[i]-1];//nos[i]-1;
			elemlist.Push(fGeoMesh->CreateGeoElement(ECube,nodes,1,index));
			continue;
		}
	}
	mesh >> title;//$ENDELM
	mesh.close();
}

void TPZGMSHReadMesh::Resequence(TPZStack<long> &Indexes,const char *meshfile){
	
	std::ifstream mesh(meshfile);
	char title[256];
	long nnod,pos;
	mesh >> title;//$NOD
	mesh >> nnod;
	long index = 0,i;
	REAL no1,no2,no3;
	for(i=0;i<nnod;i++){
		index++;
		mesh >> pos >> no1 >> no2 >> no3;
		if(index == pos){
			Indexes.Push(i);//Indexes[index] = i
		} else {
			while( index  < pos ){
				Indexes.Push(-1);//Indexes[index] = i
				index++;
			}
			Indexes.Push(i);//index = pos //Indexes[index] = i
		}
	}
	mesh.close();
}

void TPZGMSHReadMesh::PrintGeoMesh(std::ostream &out){
	
	fGeoMesh->Print(out);
	std::cout << "TPZGMSHReadMesh::PrintGeoMesh the geometric mesh of the NeoPZ was printed\n";
}
