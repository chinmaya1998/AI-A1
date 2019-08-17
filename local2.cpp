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

typedef std::vector<std::string> vecofvecs;
typedef std::vector<std::vector<int>> matrix;
int max_length = 0;

void print_intvector(std::vector<int> v){
	for (std::vector<int>::iterator i = v.begin(); i != v.end(); ++i)
	{
		std::cout << *i << ' ' ;
	}
	std::cout << std::endl;
}

void print_charvector(std::string v){
	for (std::string::iterator i = v.begin(); i != v.end(); ++i)
	{
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
	for (vecofvecs::iterator i = v.begin(); i != v.end(); ++i)
	{
		print_charvector(*i);
	}
}

int gene_pair_cost(std::string g1, std::string g2, matrix cost, std::string v){

	int c=0;
	std::string::iterator it;
	int i1,i2;
	if(g1.size()!=g2.size()){
		std::cout << "Error - Gene sizes are different!!" <<std::endl;
		return 0;
	}

	else{
		for (int i = 0; i < g1.size(); ++i){
			it = std::find (v.begin(), v.end(), g1[i]);
			i1 = it - v.begin();
			it = std::find (v.begin(), v.end(), g2[i]);
			i2 = it - v.begin();
			c = c + cost[i1][i2];
		}
	}
	return c;
}

int total_genes_cost(vecofvecs genes, matrix cost, std::string v, int dash_cost){

	int c = 0;
	for (int i = 0; i < genes.size(); ++i){
		for (int j = 0; j < genes[i].size(); ++j){
			if(genes[i][j] == '-'){
				c = c + 1;
			}		
		}
	}
	
	c = c*dash_cost;

	for (int i = 0; i < genes.size(); ++i){
		for (int j = i + 1; j < genes.size() ; ++j){
			c = c + gene_pair_cost(genes[i],genes[j],cost,v);	
		}
	}
	return c;
}

class soln_node{

	vecofvecs genes;
	int state_cost;

public:
	soln_node(vecofvecs genes1){
		genes = genes1;
	}

	void assign_cost(matrix cost, std::string v, int dash_cost){
		state_cost = total_genes_cost(genes,cost,v,dash_cost);
	}

	int get_cost(){
		return state_cost;
	}

	vecofvecs get_genes(){
		return genes;
	}

};

int random(int n, int m){
    return rand() % (m - n + 1) + n;
}


// vecofvecs get_neighbors(curr_node){

// }


soln_node simple_hill_climbing(matrix mc, string vocab, int dash_cost, vecofvecs genes){

	// std::cout << 'a' << endl;
	// simple hill-climbing search algorithm (incomplete)
	vecofvecs start_genes(genes.size());
	for(int i = 0; i < genes.size(); i++){
		string s(max_length-genes[i].size(), '-');
		start_genes[i] = genes[i] +s;
		cout << start_genes[i] << endl;
	}
	// std::cout << 'b' << endl;

	soln_node curr_node(start_genes);	// (start_genes to be randomly initialized)
	// std::cout << 'b' << endl;
	soln_node prev_node(start_genes);
	soln_node best_node = prev_node;
	soln_node temp = prev_node;
	// std::cout << 'b' << endl;


	curr_node.assign_cost(mc,vocab,dash_cost);


	int c = curr_node.get_cost();
	cout << c << endl;

	// std::cout << 'b' << endl;
	std::vector<soln_node> neighbors; 	// vector storing neighbor solution nodes
	std::vector<int> costs;				// vector storing costs of neighbor nodes
	std::vector<int> indices;			// vector storing indices of least-cost value neighbor nodes
	string s;
	vecofvecs curr_genes;
	int start = 0;
	int ct = 0;
	while(ct < pow(10, 5)){
		prev_node = curr_node;
		best_node = curr_node;
		curr_genes = curr_node.get_genes();
		for(int i = 0; i < curr_genes.size(); i++){
			s = curr_genes[i];
			start = -1;
			for(int j = 0; j < s.size(); j++){
				if(s[j] == '-'){
					if(j > 0 && start != -1){
						// cout << i << " " << j << " " << start << endl;
						curr_genes[i][j] = curr_genes[i][start];
						curr_genes[i][start] = '-';

						temp = soln_node(curr_genes);

						temp.assign_cost(mc, vocab, dash_cost);
						if(best_node.get_cost() >= temp.get_cost()){
							// cout << c << endl;
							c = temp.get_cost();
							best_node = temp;
						}
						curr_genes[i][start] = curr_genes[i][j];
						curr_genes[i][j] = '-';
					}
					else if(j < s.size()-1){
						curr_genes[i][j] = curr_genes[i][j+1];
						curr_genes[i][j+1] = '-';
						temp = soln_node(curr_genes);
						temp.assign_cost(mc, vocab, dash_cost);
						if(best_node.get_cost() >= temp.get_cost()){
							c = temp.get_cost();
							best_node = temp;
						}
						curr_genes[i][j+1] = curr_genes[i][j];
						curr_genes[i][j] = '-';
					}
				}
				else start = j;
			}
			curr_node = best_node;
		}

		if(curr_node.get_cost() == prev_node.get_cost()) ct++;
		else ct = 0;
		// if(curr_node.get_genes() == prev_node.get_genes()) break;
	}
	cout << curr_node.get_cost() << endl;
	cout << ct << endl;
	return curr_node;
}










int main(int argc, char const *argv[])
{


	float time_limit, dash_cost; 
	int vocab_len, n_genes;

	std::string vocab;
	vecofvecs genes;
	matrix mc;

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
	std::string temp_vec;
	for (int i = 0; i < n_genes; ++i)
	{
		getline(fin, temp);
		std::stringstream iss(temp);
		temp_vec.clear();
		while(iss >> temp3){
			temp_vec.push_back(temp3);
		}
		max_length = max(max_length, int(temp_vec.size()));
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




	// std::cout << "Time Limit is " << time_limit << std::endl;
	// std::cout << "Length of vocabulary is " << vocab_len << std::endl;
	// std::cout << "Vocabulary is ";
	// print_charvector(vocab);

	// std::cout << "Number of genes is " << n_genes << std::endl;

	// std::cout << "Genes are " << std::endl;
	// print_vecofvecs(genes);

	// std::cout << "Cost of adding a dash is " << dash_cost << std::endl;

	// std::cout << "Vocabulary-cost matrix is " << std::endl;
	// print_matrix(mc);


	// std::string g1 = {'-','A','C','T','G','T','G','A'};
	// std::string g2 = {'T','A','C','T','-','-','G','C'};
	// std::string g3 = {'-','A','C','T','G','-','-','A'};
	// vecofvecs mod_genes = {g1,g2,g3};


	// soln_node node1(mod_genes);
	// cout << "a" << endl;
	// node1.assign_cost(mc,vocab,dash_cost);
	// cout << "b" << endl;
	// std::cout << node1.get_cost() << std::endl;
	// cout << "c" << endl;
	// std::cout << total_genes_cost(mod_genes,mc,vocab,dash_cost) << std::endl;



	soln_node sidd = simple_hill_climbing(mc, vocab, dash_cost, genes);

	print_vecofvecs(sidd.get_genes());

	return 0;
}
