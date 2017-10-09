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

#define L 8 //log2 of hash table dimension
#define MAX_REPROBE 4 //Log2 of max reprobe

// NB funziona tutto se k<=8 e 0<=2k-l<=16-MAX_REPROBE

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

HashTable::HashTable(int k){
	
	this->k=k;
	m = pow(2, L);
	matrix.resize(2*k, 2*k);
	inverse.resize(2*k, 2*k);
	
	std::vector<unsigned short> new_thr (pow(2,L),0);  
	thr = new_thr;
	std::vector<unsigned short> new_tc (pow(2,L),0);  
	thr = new_tc;
	
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
	
	unsigned short kv = key_vector(key); //key vector 
	unsigned short hkv = f(kv); //hash of key vector
	
	unsigned short i=0;
	int pos;
	
	bool done = false;
	unsigned short inv_hkv = 0;
	do{
		pos = (hkv + HashTable::reprobe(i)) % m;
		i++;
		
		//calcolo della nuova chiave da inserire
		unsigned short newKey = hkv;
		newKey >> L;
		newKey << MAX_REPROBE; //reprobe occupa un numero fisso di bit
		newKey += i+1;
		
		done = __sync_bool_compare_and_swap(&thr[pos],0,newKey); 
		
		//calcolo dell'inverso del key vector da confrontare
		inv_hkv = thr[pos] >> MAX_REPROBE;
		inv_hkv << L;
		unsigned short rep_c = thr[pos] % (unsigned short)pow(2,MAX_REPROBE);
		inv_hkv += (pos-reprobe(rep_c-1)) % m;	
		
	}while(!done && kv!=f_rev(inv_hkv) && i<=MAX_REPROBE);
	
	if(i>MAX_REPROBE) return; 
	
	done = false;
	int oldCount = tc[pos];
	do{
		oldCount = tc[pos];
		done = __sync_bool_compare_and_swap(&tc[pos],oldCount,oldCount+1); 
	}while(!done);
	
}


unsigned short HashTable::f(unsigned short kv){
	
	VectorXf key_vector(2*k), result(2*k);
	
	for(int i=2*k-1; i>=0; i++){
		key_vector(i) = (kv >> (2*k-1-i)) %2;
	}
	
	result = matrix * key_vector; 
	unsigned short hash=0;
	for (int i=0; i<2*k; i++)
		hash += result[i]*pow(2, i);
	return hash;
}

unsigned short HashTable::f_rev(unsigned short inv_hkv){
	
	//TODO (uguale ma con matrice inversa?)
	
}

unsigned short HashTable::key_vector(std::vector<Nucleotide> key){
	unsigned short kv = 0;
	for(int i=0; i<k; i++){
		Nucleotide currentNucleotide = key[i];
		kv << 1;
		kv += currentNucleotide.getBit(0);
		kv << 1;
		kv += currentNucleotide.getBit(1);
	}
	return kv;
}

int HashTable::reprobe(unsigned short i){
	return (i*(i+1))/2;
}

//TODO
std::string HashTable::toString(){
	/*std::string temp;
	for(HashEntry he : table){
		if(!he.isEmpty())
			temp.append(he.toString() + "\n");
	}
	return temp;*/
}


