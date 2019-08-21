#include <iostream>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <bits/stdc++.h> 
#include <string> 
#include <fstream>
#include <chrono>
#include <algorithm>
#include <chrono>

using namespace std;

typedef vector<string> vecofvecs;
typedef vector<vector<int>> matrix;

float time_limit, dash_cost; 
int vocab_len, n_genes;

vector<char> vocab;
vecofvecs genes;
matrix mc;


void print_intvector(std::vector<int> v){
	for (std::vector<int>::iterator i = v.begin(); i != v.end(); ++i){
		std::cout << *i << ' ' ;
	}
	std::cout << std::endl;
}

void print_matrix(matrix mat){
	for (matrix::iterator i = mat.begin(); i != mat.end(); ++i){
		print_intvector(*i);
	}
}

void print_charvector(std::string v){
	for (std::string::iterator i = v.begin(); i != v.end(); ++i){
		std::cout << *i;
	}
	std::cout << std::endl;
}

void print_vecofvecs(vecofvecs v){
	for (vecofvecs::iterator i = v.begin(); i != v.end(); ++i){
		print_charvector(*i);
	}
}

template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b){
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::plus<T>());
    return result;
}

template <typename T>
std::vector<T> operator-(const std::vector<T>& a, const std::vector<T>& b){
    assert(a.size() == b.size());

    std::vector<T> result;
    result.reserve(a.size());

    std::transform(a.begin(), a.end(), b.begin(), 
                   std::back_inserter(result), std::minus<T>());
    return result;
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

vector<int> random_step(vector<int> path_diff){

	int k = path_diff.size();
	vector<int> step;
	step.resize(k,0);
	int ones=0;
	vector<int> t;
	t.reserve(k);
	for (int i = 0; i < path_diff.size(); ++i){
		if(path_diff[i] !=0){
			ones++;
			t.push_back(i);
		}
	}
	int x= random(1,ones);
	for (int i = 0; i < x; ++i){
		int index = rand() % t.size();
		step[t[index]] = 1;
		t.erase(t.begin()+index);
	}
	return step;
}

matrix state_initialize(vector<int> v){
	matrix path;
	vector<int> ini;
	int temp;
	int k = v.size();
	for (int i = 0; i < k; ++i){
		ini.push_back(0);
	}

	path.push_back(ini);
	vector<int> step;
	while(ini != v){
		step = random_step(v - ini);
		ini = ini + step;
		path.push_back(ini);
	}
	return path;
}

matrix gen_path(vector<int> initial, vector<int> final){
	matrix path;
	path.reserve(15);
	int k = final.size();

	path.push_back(initial);
	vector<int> step;
	while(initial != final){
		step = random_step(final - initial);
		initial = initial + step;
		path.push_back(initial);
	}
	return path;
}

vecofvecs create_genseq(matrix path, vecofvecs genes){
	vector<string> ans;
	ans.reserve(path.size());

	vector<int> help;
	help.resize(path[0].size(), 0);

	vector<int> temp;
	string temp2;
	for (int i = 0; i < path.size() - 1; ++i){
		temp = path[i+1] - path[i];
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
	ret.resize(genes[0].size(),"");
	
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
	neighbor(){

	}
	neighbor(matrix path1, vecofvecs genes1, int cost1, int postition1){
		path = path1;
		genes = create_genseq(path1, genes1);
		deltacost = total_genes_cost(genes) - cost1;
		postition = postition1;		
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


// chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
// chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

// cout << "Time difference = " << chrono::duration_cast<chrono::microseconds>(end - begin).count() << "[µs]" << endl;
// cout << "Time difference = " << chrono::duration_cast<chrono::nanoseconds> (end - begin).count() << "[ns]" << endl;



	srand(time(0));
	std::vector<int> f;
	for (int i = 0; i < genes.size(); ++i){
		f.push_back(genes[i].size());
	}

	// matrix step_dict = gen_steps( f.size() );
	matrix pa = state_initialize(f);
	vecofvecs ge = create_genseq(pa,genes);

	solution sol = solution(pa, ge);

	int max_search = 10;
	int restarts = 1;

	for (int q = 0; q < restarts; ++q){
		
		vector<int> info_data1;
		vector<int> info_data2;
		vector<int> info_data3;

		pa = state_initialize(f);
		ge = create_genseq(pa,genes);

		int step_size_limit = 12;
		int check_neighbors = 50000; 
		int positions = 10;

		for (int i = 0; i <= step_size_limit; ++i){
			info_data1.push_back(0);
			info_data2.push_back(0);
			info_data3.push_back(0);
		}

		int climbs = 0;
		int e = 0;
		int search_fresh = 0;


			while(true){

				matrix temp_path;
				int temp_cost;
				int r;

				// if(e>30){
				// 	step_size_limit = 4;
				// 	check_neighbors = 10000;
				// }


				neighbor temp_neighbor = neighbor();
				int best_step_size;

				for (int j = 0; j < positions; ++j){
					
					int step_size = random(2, step_size_limit);
					info_data3[step_size] = info_data3[step_size] + 1;
					r = random(0,pa.size() - step_size - 1);     //step_size

					neighbor temp1_neighbor = neighbor();

					for (int i = 0; i < check_neighbors; ++i){
						temp_path = gen_path( pa[r], pa[r+step_size] ); //step_size //, &step_dict
						vecofvecs temp_genes( &ge[r], &ge[r+step_size]); //step_size
						temp_cost = total_genes_cost(temp_genes);
						neighbor trial_neighbor = neighbor( temp_path, remove_dash(temp_genes), temp_cost, r+1);

						if(i == 0){
							temp1_neighbor = trial_neighbor;
						}
						else if(trial_neighbor.get_deltacost() < temp1_neighbor.get_deltacost() ){
							temp1_neighbor = trial_neighbor;
						}
					}	

					if(j == 0){
						temp_neighbor = temp1_neighbor;
						best_step_size = step_size;
					}
					else if(temp1_neighbor.get_deltacost() < temp_neighbor.get_deltacost() ){
						temp_neighbor = temp1_neighbor;
						best_step_size = step_size;
					}
				}


				e = e + 1;
				if(temp_neighbor.get_deltacost() <= 0){
					if(search_fresh >= max_search){
						break;
					}
					else if(temp_neighbor.get_deltacost() < 0){
						info_data1[best_step_size] = info_data1[best_step_size] + temp_neighbor.get_deltacost();
						info_data2[best_step_size] = info_data2[best_step_size] + 1;
						cout << "Step No.: "<<e<<" Step Size: " << best_step_size << " Delta Cost: "<< temp_neighbor.get_deltacost() << endl;
						climbs = climbs + 1;
						search_fresh = 0;
					}
					else{
						search_fresh = search_fresh + 1;
					}

				
					int pos = temp_neighbor.get_position();
					auto it  = pa.begin() + pos;
					auto it1 = pa.begin() + pos + best_step_size; //step_size
					pa.erase(it, it1);
					auto it2 = ge.begin() + pos - 1;
					auto it3 = ge.begin() + pos + best_step_size - 1 ; //step_size
					ge.erase(it2, it3);
					auto it4  = pa.begin() + pos;
					auto it5 = ge.begin() + pos - 1; 
					matrix selected_path  = temp_neighbor.get_path();
					vecofvecs selected_genes = temp_neighbor.get_genes();  
					pa.insert(it4, selected_path.begin() + 1, selected_path.end() );
					ge.insert(it5, selected_genes.begin(), selected_genes.end() );
				}
				else{
					search_fresh = search_fresh + 1;
				}
				

			}
		



		cout << "Climbs "<< climbs << endl;
		cout << "Total Steps "<< e << endl;
		cout << "Final Cost "<< total_genes_cost(ge)<< endl;
		// for (int i = 0; i < info_data2.size(); ++i){
		// 	// if(info_data2[i]!=0){
		// 	// 	float x = (float)info_data1[i]/(float)info_data2[i];
		// 	// 	cout << "Step Size: " << i << " Average Delta Cost: " << x << endl;
		// 	// }
		// 	// if(info_data3[i]!=0){
		// 	// 	float x = (float)info_data2[i]/(float)info_data3[i];
		// 	// 	cout << "Step Size: " << i << " Proportion: " << x << endl;
		// 	// 	cout << "Step Size: " << i << " Total Neighbors: " << info_data3[i] << endl;
		// 	// }

		// }
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