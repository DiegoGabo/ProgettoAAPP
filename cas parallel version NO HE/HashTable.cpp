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
#define MAX_COUNT 65535

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
	
	std::vector<unsigned short> new_table_hash_reprobe (pow(2,L),0);  
	table_hash_reprobe = new_table_hash_reprobe;
	std::vector<unsigned short> new_table_count (pow(2,L),0);  
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
}

void HashTable::incrementValue(std::vector<Nucleotide> key){
	
	unsigned short key_value = key_vector(key); //key vector 
	unsigned short hash = f(key_value); //hash of key vector
	
	unsigned short i=0;
	int pos;
	
	bool done = false;
	unsigned short inv_hash = 0;
	do{
		pos = (hash + HashTable::reprobe(i)) % m;
		i++;
		
		//calcolo della nuova chiave da inserire
		unsigned short newKey = hash;
		newKey >> L;
		newKey << MAX_REPROBE; //reprobe occupa un numero fisso di bit
		newKey += i+1;
		cout << "\nkey = " << std::to_string(newKey);
		
		done = __sync_bool_compare_and_swap(&table_hash_reprobe[pos],0,newKey); 
		
		//calcolo dell'inverso del key vector da confrontare
		inv_hash = table_hash_reprobe[pos] >> MAX_REPROBE;
		inv_hash << L;
		unsigned short rep_c = table_hash_reprobe[pos] % (unsigned short)pow(2,MAX_REPROBE);
		inv_hash += (pos-reprobe(rep_c-1));
		inv_hash = inv_hash % m;
		
		cout << "\treprobe = " << std::to_string(rep_c);
		cout << "\thash = " << std::to_string(inv_hash);
		cout << "\tf_rev = " << std::to_string(f_rev(inv_hash)%m);
		cout << "\tf(frev) = " << std::to_string(f(f_rev(inv_hash))%m);
		
	}while(!done && !(key_value==f_rev(inv_hash) && table_count[pos] < MAX_COUNT-1) && i<=MAX_REPROBE);
	
	if(i>MAX_REPROBE) 
		{
			cout << " > di max reprobe";			
			return; 
		}
	
	done = false;
	
	int oldCount = table_count[pos];
	
	do{
		oldCount = table_count[pos];
		done = __sync_bool_compare_and_swap(&table_count[pos],oldCount,oldCount+1); 
	}while(!done);
	
}


unsigned short HashTable::f(unsigned short key){
	
	VectorXf key_vector(2*k), result(2*k);
	
	for(int i=2*k-1; i>=0; i--){
		key_vector[i] = (key >> (2*k-1-i)) %2;
	}
	
	result = matrix * key_vector; 
	unsigned short hash=0;
	for (int i=0; i<2*k; i++)
		hash += result[i]*pow(2, i);
	return hash;
}

unsigned short HashTable::f_rev(unsigned short inv_hash){
	
	VectorXf key_vector(2*k), result(2*k);
	
	for(int i=2*k-1; i>=0; i--){
		key_vector[i] = (inv_hash >> (2*k-1-i)) %2;
	}
	
	result = inverse * key_vector; 
	unsigned short hash=0;
	for (int i=0; i<2*k; i++)
		hash += result[i]*pow(2, i);
	return hash;
	
}

unsigned short HashTable::key_vector(std::vector<Nucleotide> key){
	unsigned short key_value = 0;
	for(int i=0; i<k; i++){
		Nucleotide currentNucleotide = key[i];
		key_value << 1;
		key_value += currentNucleotide.getBit(0);
		key_value << 1;
		key_value += currentNucleotide.getBit(1);
	}
	return key_value;
}

int HashTable::reprobe(unsigned short i){
	return (i*(i+1))/2;
}

std::string HashTable::toString(){
	std::string temp = "";
	for(int i=0; i<m; i++)
		if (table_count[i] != 0)
		{
			temp.append(std::to_string(table_hash_reprobe[i]));
			temp.append(" ");
			temp.append(std::to_string(table_count[i]));
			temp.append("\n");
		}
	return temp;
}
