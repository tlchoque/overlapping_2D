
bool pairCompare(const std::pair<int, int>& firstElem, const std::pair<int, int>& secondElem) {
  return firstElem.second < secondElem.second;

}

void update_vertex_info(Delaunay &dt, Vertex &v){
	vector<pair<int,int>> &regions = v->info().m_singular_regions_and_subgraphs;
	vector<int> &labels = v->info().m_regions_around;
	regions.clear();
	labels.clear();

	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;

	bool isInternal = true;
	int faceLabel =  fc->info().m_label;
	do{
		if( faceLabel != fc->info().m_label){
			fc_start = fc;
			isInternal = false;
			break;
		}
		++fc;
	}while (fc!= fc_start);
	if( isInternal ) return;

	int labelcirculator = -2;
	do{
		int label =  fc->info().m_label;
		if( label == labelcirculator ){ fc++;continue;}
		//insert the label
		labelcirculator = label;
		//cout<<"rou coleur: "<<labelcirculator<<endl;
		if( labelcirculator == HOLE) { fc++;continue;}

		if( !dt.is_infinite(fc) && label != -1 ){			
			if( std::find(labels.begin(), labels.end(), label ) != labels.end() ){//found in labels
				vector<pair<int,int>>::iterator ite = 
					std::find_if(regions.begin(), regions.end(),[&label](std::pair<int, char> const& elem) {return elem.first == label;	});
				if( !( ite != regions.end() ) ){//not found in regions
					regions.push_back( make_pair(label , 2 ) );
				}
				else{//found
					int position = std::distance( regions.begin(), ite );		
					regions[position].second++;
				}
			}
			else
				labels.push_back( label );
		}		
		++fc;
	}while (fc!= fc_start);
	std::sort(regions.begin(), regions.end(), pairCompare);
}	

bool is_singular_after_relabeling_2(Delaunay &dt, Vertex &v){
	vector<int> labels ;
	Face_circulator fc_start = dt.incident_faces(v);
	Face_circulator fc = fc_start;

	bool isInternal = true;
	int faceLabel =  fc->info().m_label;
	do{
		if( faceLabel != fc->info().m_label){
			fc_start = fc;
			isInternal = false;
			break;
		}
		++fc;
	}while (fc!= fc_start);
	if( isInternal ) return false;

	int labelcirculator = -2;
	do{
		int label =  fc->info().m_label;
		if( label == labelcirculator ){ fc++;continue;}
		//insert the label
		labelcirculator = label;
		//cout<<"rou coleur: "<<labelcirculator<<endl;
		if( labelcirculator == HOLE) { fc++;continue;}



		if( !dt.is_infinite(fc) && label != -1 ){			
			if( std::find(labels.begin(), labels.end(), label ) != labels.end() ){//found in label, is singulars
				return true;
			}
			else
				labels.push_back( label );
		}		
		++fc;
	}while (fc!= fc_start);
	return false;
}