#include <vector>  
#include <string>  
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <ctime>
#include <cstdlib>
#include <stddef.h>
#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <atomic>
#include "structures.hpp"
#include "HashTable.hpp"

#define MAX_REPROBE 10 //log2 of max reprobe
#define MAX_COUNT 4294967285 //max value of count

using namespace boost::numeric::ublas;
using namespace Eigen;

/***********************************************************************************************************/

/*function that checks if two nucleotides are equal*/
bool equal(std::vector<Nucleotide> n1, std::vector<Nucleotide> n2)
{
	int lenght = n1.size();
	for (int i=0; i<lenght; i++)
	{
		if(!n1[i].equal(n2[i]))
			return false;
	}
	return true;
}

/***********************************************************************************************************/

/*constructor that initialize the hash table*/
HashTable::HashTable(int k, int L){
	
	this->k=k;
	this->num=0;
	this->stop = false;
	m = pow(2, L);
	matrix.resize(2*k, 2*k);
	inverse.resize(2*k, 2*k);
	
	std::vector<int> new_table_hash (pow(2,L),0);  
	table_hash = new_table_hash;
	std::vector<int> new_table_count (pow(2,L),0);  
	table_count = new_table_count;
	
	srand((unsigned)time(NULL));
	do{
        for(int i=0;i<2*k;i++){
            
			for(int j=0;j<2*k;j++){
                matrix(i,j) = (rand() % 2);
            	}
        }
	}
	while(matrix.determinant() == 0);
	inverse = matrix.inverse();
	
	this->flush();
}

/*method that given a key(vector of Nucleotides) update the count value in the proper hentry of the hash table*/
void HashTable::incrementValue(std::vector<Nucleotide> key){
	
	int key_value = key_vector(key); //key vector 
	int hash = f(key_value); //hash of key vector
	
	int i=0;
	int pos;
	bool done = false;

	while(stop){}

	do{
		pos = (hash + HashTable::reprobe(i)) % m;
		i++;
		done = __sync_bool_compare_and_swap(&table_hash[pos],0,key_value); 
		if(key_value = table_hash[pos])
			done = true;
		if(table_count[pos] >= MAX_COUNT-1)
			done = false;
	}
	while(!done && i < MAX_REPROBE );
	
	bool is_max_reprobe = false;
	if(i>=MAX_REPROBE) 
	#pragma omp critical
	//block in which the content of the hash table is saved on the disk and the hash table is flushed
	{
		stop = true; //bool used to block the other threads

		//the file is created
		std::string file_name = "result";
		file_name.append(std::to_string(num));
		std::ofstream outfile(file_name);

		//the hash table is ordered and then written on the file
		this->order(0, this->m);
		outfile << this->toString(); 
		outfile.close();
		
		//the hash table is flushed
		num++;
		this->flush();
		stop = false;
		is_max_reprobe = true;

		//the key count is now updated
		this->incrementValue(key);
	}
	
	if(is_max_reprobe)
		return;
	
	done = false;
	int oldCount = table_count[pos];
	
	do{
		oldCount = table_count[pos];
		done = __sync_bool_compare_and_swap(&table_count[pos],oldCount,oldCount+1); 
	}while(!done);
	
}

/*method that given a key returns its hash*/
int HashTable::f(int key){
	
	VectorXf key_vector(2*k), result(2*k);
	
	for(int i=2*k-1; i>=0; i--){
		key_vector[i] = (key >> (2*k-1-i)) %2;
	}
	
	result = matrix * key_vector; 
	int hash=0;
	for (int i=0; i<2*k; i++)
		hash += result[2*k-i-1]*pow(2, i);
	return hash;
}
	
/*method that trasforms the key passed as parameter to an integer*/
int HashTable::key_vector(std::vector<Nucleotide> key){
	int key_value = 0;
	for(int i=0; i<k; i++){
		key_value += pow(2, 2*i)*key[i].getBit(0);
		key_value += pow(2, 2*i+1)*key[i].getBit(1);
	}
	return key_value+1;
}

/*method that given an integer as a parameter returns the number of reprobe*/
int HashTable::reprobe(int i){
	return (i*(i+1))/2;
}

/*method that swaps two entry in the hash table*/
void HashTable::swap(int a, int b)
{
    int t_hash = table_hash[a];
	int t_count = table_count[a];
    table_hash[a] = table_hash[b];
	table_count[a] = table_count[b];
    table_hash[b] = t_hash;
	table_count[b] = t_count;
}
 
/*partition method of the quick sort algorithm*/
int HashTable::partition_hash(int low, int high)
{
    int pivot = table_hash[high];    // pivot
    int i = (low - 1);  // Index of smaller element
 
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than or equal to pivot
        if (table_hash[j] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(i, j);
        }
    }
    swap(i + 1, high);
    return (i + 1);
}
 
/*method that order the hashTable according to quicksort algorithm*/
void HashTable::order(int low, int high)
{
    if (low < high)
    {
        //pi is partitioning index, arr[p] is now at right place 
        int pi = partition_hash(low, high);
 
        // Separately sort elements before partition and after partition
        order(low, pi - 1);
		order(pi + 1, high);
		
    }
}

/*method that flush the content of the hash table*/
void HashTable::flush()
{
	for(int i=0; i<m; i++)
	{
		table_count[i] = 0;
		table_hash[i] = 0;
	}
}

/*method that the return the number of file in which there are the contents of the previous hash table*/
int HashTable::getNum()
{
	return this->num;
}

/*method that increments the number of the field num*/
void HashTable::incrementNum()
{
	this->num = this->num+1;
}

/*method that returns a string that represents the hash table*/
std::string HashTable::toString(){
	std::string temp = "";
	for(int i=0; i<m+2; i++)
		if (table_count[i] != 0 && table_hash[i] != 0)
		{
			temp.append(std::to_string(table_hash[i]));
			temp.append(" ");
			temp.append(std::to_string(table_count[i]));
			temp.append("\n");
		}
	return temp;
}

