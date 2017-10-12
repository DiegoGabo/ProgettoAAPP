#include <vector>  
#include <string>  
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <ctime>
#include <cstdlib>
#include <stddef.h>
#include <Eigen/Dense>
#include <omp.h>
#include "structures.hpp"
#include "HashTable.hpp"

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

/*constructor*/
HashEntry::HashEntry(){
	count=0;
}

/*check if key is empty*/
bool HashEntry::isEmpty(){
	if(key.size()==0) return true;
	else return false;
}

/*getter count*/
int HashEntry::getC(){
	return count;
}

/*getter key*/
std::vector<Nucleotide> HashEntry::getK(){
	return key;
}

/*setter count*/
void HashEntry::setC(int count){
	this->count=count;
}

/*setter key*/
void HashEntry::setK(std::vector<Nucleotide> key){
	this->key=key;
}

/*to string*/
std::string HashEntry::toString(){
	std::string temp;
	for(Nucleotide n : key){
		temp.append(n.toString());
	}
	temp.append(" ");
	temp.append(std::to_string(count));
	return temp;
}

/***********************************************************************************************************/

/*constructor that initialize the hash table*/
HashTable::HashTable(int k, int L){
	
	this->k=k;
	this->L=L;
	matrix.resize(2*k, 2*k);
	inverse.resize(2*k, 2*k);

	for(int i=0;i<pow(2,L);i++){
		HashEntry he;
		table.push_back(he);
	}
	srand((unsigned)time(NULL));
	do{
        for(int i=0;i<2*k;i++){
            
			for(int j=0;j<2*k;j++){
                matrix(i,j) = (rand() % 2);
            	}
        }
	}
	while(matrix.determinant() == 0);
}

/*method that given a key(vector of Nucleotides) update the count value in the proper hentry of the hash table*/
void HashTable::incrementValue(std::vector<Nucleotide> key){
	int i=0, pos;
	int hash=HashTable::f(key);
	
	#pragma omp critical
	{
		do{
		int tableLenght = pow(2, L);
		pos = (hash + HashTable::reprobe(i)) % tableLenght;
		i++;
		}while(!(table[pos].isEmpty()) && !(equal(key,table[pos].getK())) && i < 10);
		if(table[pos].isEmpty()) table[pos].setK(key); 
		table[pos].setC(table[pos].getC()+1);
	}
}

/*method that given an integer as a parameter returns the number of reprobe*/
int HashTable::reprobe(int i){
	return (i*(i+1))/2;
}

/*method that given a key returns its hash*/
int HashTable::f(std::vector<Nucleotide> key){
	VectorXf key_vector(2*k), result(2*k);
	for(int i=0; i<k; i++)
	{
		Nucleotide currentNucleotide = key[i];
		key_vector(2*i) = currentNucleotide.getBit(0);
		key_vector(2*i+1) = currentNucleotide.getBit(1);
	}
	
	result = matrix * key_vector; 
	int hash=0;
	for (int i=0; i<2*k; i++)
		hash += result[i]*pow(2, i);
	return hash;
}

/*method that returns a string that represents the hash table*/
std::string HashTable::toString(){
	std::string temp;
	for(HashEntry he : table){
		if(!he.isEmpty())
			temp.append(he.toString() + "\n");
	}
	return temp;
}


