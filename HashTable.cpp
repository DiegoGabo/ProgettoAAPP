#include <vector>  
#include <string>  
#include <math.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <ctime>
#include <cstdlib>

#define L 5 //log of hash table dimension

/***********************************************************************************************************/

/*TODO equal*/

/***********************************************************************************************************/

/*a single entry of the hash table*/
class HashEntry
{
	private:
		std::vector<Nucleotide> key;
		int count;
	public:
		HashEntry();
		bool isEmpty();
		int getC();
		void setC(int);
		void setK(std::vector<Nucleotide>);
		std::string toString();
}

/*constructor*/
HashEntry::HashEntry(){
	count=0;
}

/*check if key is empty*/
HashEntry::isEmpty(){
	if(!key) return true;
	else return false;
}

/*getter count*/
int HashEntry::getC(){
	return count;
}

/*setter count*/
void HashEntry::setC(int count){
	this.count=count;
}

/*setter key*/
void HashEntry::setK(std::vector<Nucleotide> key){
	this.key=key;
}

/*to string*/
std::string HashEntry::toString(){
	std::string temp;
	for(Nucleotide n : key){
		temp.append(n.toString());
	}
	temp.append(" "+count);
	return temp;
}

/***********************************************************************************************************/

using namespace boost::numeric::ublas;
class HashTable
{
	private:
		int k;
		std::vector<HashEntry> table;
		matrix<bool> m;
		int f(std::vector<Nucleotide>); /*TODO*/
		int reprobe(int);
	public:
		HashTable(int); 
		IncrementValue(std::vector<Nucleotide>);
		std::string toString();
}

HashTable::HashTable(int k){
	this.k=k;
	for(int i=0;i<pow(2,L);i++){
		HashEntry he();
		table.push_back(he);
	}
	/*
	matrix<bool> temp(2*k,2*k);
	m=temp;
	do{
		for(int i=0;i<2*k;i++){
			for(int j=0;j<2*k;j++){
				srand((unsigned)time(NULL));
				m[i][j] = (bool)(rand() % 2);
			}
		}
	}while(/*INVERTIBILE*/)
	*/
}

HashTable::IncrementValue(std::vector<Nucleotide> key){
	int i=0, pos;
	int hash=HashTable::f(key);
	do{
		pos = hash + HashTable::reprobe(i);
		i++;
	}while(!table[pos].isEmpty() && !equal(key,table[pos]);
	if(table[pos].isEmpty()) table[pos].setK(key); 
	table[pos].setC(table[pos].getC()+1);
}

int HashTable::reprobe(int i){
	return (i*(i+1))/2;
}

int HashTable::f(std::vector<Nucleotide> key){
	return 5;
}

std::string HashTable::toString(){
	std::string temp;
	for(HashEntry he : table){
		temp.append(he.toString + "/n");
	}
	return temp;
}


