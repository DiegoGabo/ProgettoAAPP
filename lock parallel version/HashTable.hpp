#include <Eigen/Dense>
#include <vector>  
#include <string>

using namespace Eigen;

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
		std::vector<Nucleotide> getK();
		void setC(int);
		void setK(std::vector<Nucleotide>);
		std::string toString();
};

/*class used to implement the hash table in which are counted the k-mers*/
class HashTable
{
	private:
		int k, L;
		std::vector<HashEntry> table;
		int num=0;
		MatrixXf matrix, inverse;
		int f(std::vector<Nucleotide>); 
		int reprobe(int);
	public:
		HashTable(int, int); 
		void incrementValue(std::vector<Nucleotide>);
		std::string toString();
		int getNum();
};
