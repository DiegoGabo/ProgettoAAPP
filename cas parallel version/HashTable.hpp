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
		int* getPC();
};

class HashTable
{
	private:
		int k;
		std::vector<HashEntry> table;
		MatrixXf matrix, inverse;
		int f(std::vector<Nucleotide>); /*TODO*/
		int reprobe(int);
	public:
		HashTable(int); 
		void incrementValue(std::vector<Nucleotide>);
		std::string toString();
};
