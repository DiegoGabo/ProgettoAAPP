#include <vector>  
#include <string>  
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <ctime>
#include <cstdlib>
#include <stddef.h>
#include <Eigen/Dense>
#include <atomic>
#include "structures.hpp"
#include "HashTable.hpp"

#define L 8 //log of hash table dimension

using namespace boost::numeric::ublas;
using namespace Eigen;

/***********************************************************************************************************/

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
HashEntry::HashEntry() noexcept{
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

/*getter pointer to count*/
int* HashEntry::getPC(){
	return &count;
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

HashTable::HashTable(int k){
	
	this->k=k;
	matrix.resize(2*k, 2*k);
	inverse.resize(2*k, 2*k);

	for(int i=0;i<pow(2,L);i++){
		std::atomic<HashEntry> he;
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

void HashTable::incrementValue(std::vector<Nucleotide> key){
	int i=0, pos;
	int hash=HashTable::f(key);
	int tableLenght = pow(2, L);
	HashEntry empty_he; /*CHIAVE VUOTA E COUNT 0*/
	HashEntry new_he; new_he.setK(key); /*CHIAVE KEY E COUNT 0*/
	bool done = false;
	do{
		pos = (hash + HashTable::reprobe(i)) % tableLenght;
		i++;
		done = std::atomic_compare_exchange_strong(&table[pos],empty_he,new_he); /*CAS*/
	}while(!done || !(equal(key,table[pos].load(std::memory_order_relaxed).getK())) || i>=10); /*esco anche se key non è vuota ma il confronto della chiave da successo*/
	int oldval = table[pos].load(std::memory_order_relaxed).getC();
	done = false;
	do{
		oldval = table[pos].load(std::memory_order_relaxed).getC();
		HashEntry old_he; old_he.setC(oldval); old_he.setK(key); /*CHIAVE KEY E COUNT VECCHIO*/
		new_he.setC(oldval+1); /*CHIAVE KEY E COUNT VECCHIO INCREMENTATO*/
		done = std::atomic_compare_exchange_strong(&table[pos],old_he,new_he); /*CAS*/
	}while(!done);
	/*NB. con i>=10 non funziona. confronto sempre nel secondo cas con un hash entry la cui key è uguaòe a quella corrente => ciclo infinito
	..forse potrei mettere un table[pos].load(std::memory_order_relaxed).getK())) al posto di key*/	
}


int HashTable::reprobe(int i){
	return (i*(i+1))/2;
}

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

std::string HashTable::toString(){
	std::string temp;
	for(HashEntry he : table){
		if(!he.isEmpty())
			temp.append(he.toString() + "\n");
	}
	return temp;
}


