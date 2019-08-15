#include <iostream>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <bits/stdc++.h> 
#include <string> 
#include <fstream>
#include <chrono>
#include <algorithm>

typedef std::vector<std::vector<char>> vecofvecs;
typedef std::vector<std::vector<int>> matrix;

void print_intvector(std::vector<int> v){
	for (std::vector<int>::iterator i = v.begin(); i != v.end(); ++i)
	{
		std::cout << *i << ' ' ;
	}
	std::cout << std::endl;
}

void print_charvector(std::vector<char> v){
	for (std::vector<char>::iterator i = v.begin(); i != v.end(); ++i)
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

int gene_pair_cost(std::vector<char> g1, std::vector<char> g2, matrix cost, std::vector<char> v){

	int c=0;
	std::vector<char>::iterator it;
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

int total_genes_cost(vecofvecs genes, matrix cost, std::vector<char> v, int dash_cost){

	int c = 0;
	for (vecofvecs::iterator i = genes.begin(); i != genes.end(); ++i){
		for (int j = 0; j < genes[0].size(); ++j){
			if((*i)[j] == '-'){
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

	void assign_cost(matrix cost, std::vector<char> v, int dash_cost){
		state_cost = total_genes_cost(genes,cost,v,dash_cost);
	}

	int get_cost(){
		return state_cost;
	}

};

int random(int n, int m){
    return rand() % (m - n + 1) + n;
}

int main(int argc, char const *argv[])
{


	float time_limit, dash_cost; 
	int vocab_len, n_genes;

	std::vector<char> vocab;
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
	std::vector<char> temp_vec;
	for (int i = 0; i < n_genes; ++i)
	{
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


	std::vector<char> g1 = {'-','A','C','T','G','T','G','A'};
	std::vector<char> g2 = {'T','A','C','T','-','-','G','C'};
	std::vector<char> g3 = {'-','A','C','T','G','-','-','A'};
	vecofvecs mod_genes = {g1,g2,g3};

	soln_node node1(mod_genes);
	node1.assign_cost(mc,vocab,dash_cost);
	std::cout << node1.get_cost() << std::endl;
	std::cout << total_genes_cost(mod_genes,mc,vocab,dash_cost) << std::endl;


	// simple hill-climbing search algorithm (incomplete)
	soln_node curr_node(start_genes);	// (start_genes to be randomly initialized)
	curr_node.assign_cost(mc,vocab,dash_cost);
	int c = curr_node.get_cost();

	std::vector<soln_node> neighbors; 	// vector storing neighbor solution nodes
	std::vector<int> costs;				// vector storing costs of neighbor nodes
	std::vector<int> indices;			// vector storing indices of least-cost value neighbor nodes
	while(true){

		neighbors = get_neighbors(curr_node); // function which returns the vector of neighbors of current state (to be made)

		for (int i = 0; i < neighbors.size(); ++i){
			neighbors[i].assign_cost(mc,vocab,dash_cost);
			costs.push_back(neighbors[i].get_cost());
		}

		auto i = std::min_element(costs.begin(), costs.end());
		// std::cout << *i << std::endl;
			
		if(c <= *i){					// cost should be strictly decreasing (can be changed)
			break;
		}
		indices.clear();
		for (int i = 0; i < costs.size(); ++i){
			if(costs[i] == *i){
				indices.push_back(i);
			}
		}


		srand(time(0));
		curr_node = neighbors[indices[random(0,indices.size())]];

	}

	




	return 0;
}