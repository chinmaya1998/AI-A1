#include <iostream>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <bits/stdc++.h> 
#include <string> 
#include <fstream>
#include <chrono>
#include <algorithm>
using namespace std;

typedef vector<string> vecofvecs;
typedef vector<vector<int>> matrix;

float time_limit, dash_cost; 
int vocab_len, n_genes;

vector<char> vocab;
vecofvecs genes;
matrix mc;


void print_intvector(std::vector<int> v){
	for (std::vector<int>::iterator i = v.begin(); i != v.end(); ++i)
	{
		std::cout << *i << ' ' ;
	}
	std::cout << std::endl;
}

void print_charvector(std::string v){
	for (std::string::iterator i = v.begin(); i != v.end(); ++i){
		std::cout << *i;
	}
	std::cout << std::endl;
}

void print_matrix(matrix mat){
	for (matrix::iterator i = mat.begin(); i != mat.end(); ++i)
	{
		print_intvector(*i);
	}
}

void print_vecofvecs(vecofvecs v){
	for (vecofvecs::iterator i = v.begin(); i != v.end(); ++i){
		print_charvector(*i);
	}
}

vector<vector<int>> shr(vector<vector<int>> v){
	vector<vector<int>> newset;
	for (int i = 0; i < v.size(); ++i){
		v[i].push_back(0);
		newset.push_back(v[i]);
		v[i].pop_back();
		v[i].push_back(1);
		newset.push_back(v[i]);
		v[i].pop_back();
	}
	return newset;
}

vector<vector<int>> gen_steps(int k){
	int n = 0;
	vector<vector<int>> ya = {{0},{1}};
	while( n < (k -1) ){
		ya = shr(ya);
		n = n + 1;
	}
	ya.erase(ya.begin());
	return ya;
}


int random(int n, int m){
    return rand() % (m - n + 1) + n;
}

int random(int m){
    return rand() % (m + 1);
}

int step_check(vector<int> v1, vector<int> v2, vector<int> step){				// such that v2 > v1
	for (int i = 0; i < v1.size(); ++i){
		if( ( v1[i] + step[i] > v2[i] ) ){
			return -1; 	// not possible
		}
	}
	return 1;			// possible
}


std::vector<int> add_vecs(std::vector<int> v1, std::vector<int> v2){
	std::vector<int> v;
	for (int i = 0; i < v1.size(); ++i){
		v.push_back(v1[i] + v2[i]);
	}
	return v;
}

std::vector<int> sub_vecs(std::vector<int> v1, std::vector<int> v2){
	std::vector<int> v;
	for (int i = 0; i < v1.size(); ++i){
		v.push_back(v1[i] - v2[i]);
	}
	return v;
}

matrix state_initialize(vector<int> v, matrix* step_dict){
	matrix path;
	vector<int> ini;
	int temp;
	int k = v.size();
	for (int i = 0; i < k; ++i){
		ini.push_back(0);
	}

	path.push_back(ini);
	vector<int> step;
	int si = (*step_dict).size();
	while(ini != v){

		temp = random(si-1);

		step = (*step_dict)[temp];

		if(step_check(ini,v,step) == 1){
			ini = add_vecs(ini,step);
			path.push_back(ini);
		}
		
	}
	return path;
}

matrix gen_path(vector<int> initial, vector<int> final, matrix* step_dict){
	matrix path;
	int temp;
	int k = final.size();

	path.push_back(initial);
	vector<int> step;
	int si = (*step_dict).size();
	while(initial != final){

		temp = random(si-1);
		step = (*step_dict)[temp];

		if(step_check(initial,final,step) == 1){
			initial = add_vecs(initial,step);
			path.push_back(initial);
		}
	}
	return path;
}

vecofvecs create_genseq(matrix path, vecofvecs genes){
	vector<string> ans;
	vector<int> help;
	for (int i = 0; i < path[0].size(); ++i){
		help.push_back(0);
	}
	vector<int> temp;
	string temp2;
	for (int i = 0; i < path.size() - 1; ++i){
		temp = sub_vecs(path[i+1], path[i]);
		for (int j = 0; j < path[0].size(); ++j){
			if (temp[j] == 0){
				temp2.push_back('-');
			}
			else{
				temp2.push_back(genes[j][help[j]]);
				help[j] = help[j] + 1;
			}
		}
		ans.push_back(temp2);
		temp2.clear();
	}
	return ans;
}


int column_cost(string s){
	int c = 0;
	std::vector<char>::iterator it;
	int i1,i2;
	for (int i = 0; i < s.size(); ++i){
		if(s[i] == '-'){
			c = c + 1;
		}
	}
	c = c*dash_cost;
	// cout << c << endl;
	for (int i = 0; i < s.size(); ++i){
		for (int j = i + 1; j < s.size(); ++j){

			it = std::find (vocab.begin(), vocab.end(), s[i]);
			i1 = it - vocab.begin();
			it = std::find (vocab.begin(), vocab.end(), s[j]);
			i2 = it - vocab.begin();
			c = c + mc[i1][i2];
		}
	}
	// cout << c << endl;
	return c;
}

int total_genes_cost(vecofvecs genes){

	int c = 0;
	// string temp;
	// cout << 'e'<<endl;
	for (int i = 0; i < genes.size(); ++i){
		c = c + column_cost(genes[i]);
	}
	
	return c;
}

vecofvecs remove_dash(vecofvecs genes){
	vecofvecs ret;
	for (int i = 0; i < genes[0].size(); ++i){
		ret.push_back("");
	}
	char temp;
	for (int i = 0; i < genes.size(); ++i){
		for (int j = 0; j < genes[0].size(); ++j){
			temp = genes[i][j];
			if(temp != '-'){
				ret[j].push_back(temp);
			}
		}
	}
	return ret;
}

vecofvecs remove_dash2(vecofvecs genes){
	vecofvecs ret;
	string temp;
	for (int i = 0; i < genes.size(); ++i){
		for (int j = 0; j < genes[0].size(); ++j){
			if(genes[i][j] != '-'){
				temp.push_back(genes[i][j]);
			}
		}
		ret.push_back(temp);
		temp.clear();
	}
	
	return ret;
}

vecofvecs vertical_to_horizontal(vecofvecs genes){
	vecofvecs ret;
	string temp;
	for(int i = 0; i < genes[0].size(); ++i){
		for (int j = 0; j < genes.size(); ++j){
			temp.push_back(genes[j][i]);
		}
		ret.push_back(temp);
		temp.clear();
	}
	return ret;
}

vecofvecs horizontal_to_vertical(vecofvecs genes){
	vecofvecs ret;
	string temp;
	for(int i = 0; i < genes[0].size(); ++i){
		for (int j = 0; j < genes.size(); ++j){
			temp.push_back(genes[j][i]);
		}
		ret.push_back(temp);
		temp.clear();
	}
	return ret;
}

class neighbor{
	matrix path;
	vecofvecs genes;
	int deltacost;
	int postition;
public:
	neighbor(matrix path1, vecofvecs genes1, int cost1, int postition1){
		path = path1;

		// cout << "Final Path " << endl;
		// print_matrix(path);

		genes = create_genseq(path1, genes1);

		// cout << "Final Seq " <<endl;
		// print_vecofvecs(genes);

		deltacost = total_genes_cost(genes) - cost1;

		// cout << "Final Cost: " << total_genes_cost(genes) << endl;

		postition = postition1;

		// cout << "Position: " << postition << endl;		
	}
	int get_deltacost(){return deltacost;}
	matrix get_path(){return path;}
	vecofvecs get_genes(){return genes;}
	int get_position(){return postition;}
};

class solution{
	matrix path;
	vecofvecs genes;
	int cost;
public:
	solution(matrix path1, vecofvecs genes1){
		path = path1;
		genes =  genes1;
		cost = total_genes_cost(genes1);
	}
	int get_cost(){return cost;}
	matrix get_path(){return path;}
	vecofvecs get_genes(){return genes;}
};



int main(int argc, char const *argv[]){
	std::ifstream fin; 
	fin.open(argv[1]);

	std::string temp;
	getline(fin, temp);
	time_limit = stof(temp); 	// getting time limit
	getline(fin, temp);
	vocab_len = stoi(temp);		// getting length of vocabulary
	  
	// getting vocabulary   
	getline(fin, temp);
	std::stringstream iss(temp);
	std::string temp2;
	while(iss >> temp2){
		vocab.push_back(temp2[0]);
	}
	vocab.push_back('-');		// Adding dash to vocab at the end for purpose of cost calculation

	getline(fin, temp);
	n_genes = stoi(temp);		// getting no. of genes

	// getting list of genes
	char temp3;
	string temp_vec;
	for (int i = 0; i < n_genes; ++i){
		getline(fin, temp);
		std::stringstream iss(temp);
		temp_vec.clear();
		while(iss >> temp3){
			temp_vec.push_back(temp3);
		}
		genes.push_back(temp_vec);
	}

	getline(fin, temp);
	dash_cost = stoi(temp);		// getting cost of adding a dash

	// getting vocabulary cost matrix
	int temp4;
	for (int i = 0; i < (vocab_len + 1); ++i){
		std::vector<int> v;
		mc.push_back(v);
	}
	for (int i = 0; i < (vocab_len + 1); ++i){
		getline(fin, temp);
		std::stringstream iss(temp);
		while(iss >> temp4){
			mc[i].push_back(temp4);
		}
	}

	fin.close();




	srand(time(0));
	std::vector<int> f;
	for (int i = 0; i < genes.size(); ++i){
		f.push_back(genes[i].size());
	}

	// int max_len_genes = *max_element(f.begin(), f.end()); 

	matrix step_dict = gen_steps( f.size() );

	matrix pa = state_initialize(f, &step_dict);
	vecofvecs ge = create_genseq(pa,genes);
	// cout << "Initial Cost : " << total_genes_cost(ge) << endl;


	// print_matrix(pa);
	// print_vecofvecs(ge);
	
	// cout << "--------------------------------" <<endl;

	solution sol = solution(pa, ge);

	int step_size_limit = 4;
	int check_neighbors = 100;
	int max_search = 80;
	int restarts = 100;

	for (int q = 0; q < restarts; ++q){

		pa = state_initialize(f, &step_dict);
		ge = create_genseq(pa,genes);

		int climbs = 0;
		int e = 0;

		int search_fresh = 0;

		while(true){
			// cout << search_fresh << endl;

			search_fresh = search_fresh + 1;

			matrix temp_path;
			int temp_cost;
			int r;

			// cout << "Neighbor: " << i << endl;
			int step_size = random(1, step_size_limit);
			r = random(0,pa.size() - step_size - 1);     //step_size

			temp_path = gen_path( pa[r], pa[r+step_size], &step_dict ); //step_size
			vecofvecs temp_genes( &ge[r], &ge[r+step_size]); //step_size
			temp_cost = total_genes_cost(temp_genes);		
			neighbor temp_neighbor = neighbor( temp_path, remove_dash(temp_genes), temp_cost, r+1);

			for (int i = 0; i < check_neighbors; ++i){
				temp_path = gen_path( pa[r], pa[r+step_size], &step_dict ); //step_size
				vecofvecs temp_genes( &ge[r], &ge[r+step_size]); //step_size
				temp_cost = total_genes_cost(temp_genes);
				neighbor trial_neighbor = neighbor( temp_path, remove_dash(temp_genes), temp_cost, r+1);

				if(trial_neighbor.get_deltacost() < temp_neighbor.get_deltacost() ){
					temp_neighbor = trial_neighbor;
				}
			}

			
			e = e + 1;
			if(temp_neighbor.get_deltacost() <= 0){
				if(temp_neighbor.get_deltacost() == 0 and search_fresh > max_search){
					break;
				}
				else if(temp_neighbor.get_deltacost() < 0){
					climbs = climbs + 1;
					search_fresh = 0;
				}
				else{
					search_fresh = search_fresh + 1;
				}

				
				int pos = temp_neighbor.get_position();
				auto it  = pa.begin() + pos;
				auto it1 = pa.begin() + pos + step_size; //step_size
				pa.erase(it, it1);
				auto it2 = ge.begin() + pos - 1;
				auto it3 = ge.begin() + pos + step_size - 1 ; //step_size
				ge.erase(it2, it3);
				auto it4  = pa.begin() + pos;
				auto it5 = ge.begin() + pos - 1; 
				matrix selected_path  = temp_neighbor.get_path();
				vecofvecs selected_genes = temp_neighbor.get_genes();  
				pa.insert(it4, selected_path.begin() + 1, selected_path.end() );
				ge.insert(it5, selected_genes.begin(), selected_genes.end() );
			}


		}
		cout << "Climbs "<< climbs << endl;
		cout << "Total Steps "<< e << endl;
		cout << "Cost "<< total_genes_cost(ge)<< endl;
		cout << "------------------ " << endl;

		if(total_genes_cost(ge) < sol.get_cost()){
			sol = solution(pa,ge);	
		}
		

	}	


	// cout << "--------------------------------" <<endl;
	vecofvecs final_answer = vertical_to_horizontal(sol.get_genes()); 

	// print_vecofvecs(final_answer);
	cout << "Final Cost : " << total_genes_cost(sol.get_genes()) << endl;

	std::ofstream fout; 
	fout.open(argv[2]);
	for (vecofvecs::iterator i = final_answer.begin(); i != final_answer.end(); ++i){
		fout << *i << endl;
	}
	fout.close();


	if(remove_dash2(final_answer) == genes){
		cout << "Gene Sequence Matched" << endl;
	}

	
}
