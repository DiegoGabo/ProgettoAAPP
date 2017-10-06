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
		HashEntry() noexcept;
		bool isEmpty();
		int getC();
		std::vector<Nucleotide> getK();
		void setC(int);
		void setK(std::vector<Nucleotide>);
		std::string toString();
		int* getPC();
};

class Atomicwrap
{
	public:
		std::atomic<HashEntry> ahe;
		Atomicwrap() noexcept;
};

class HashTable
{
	private:
		int k;
		std::vector<Atomicwrap> table;
		MatrixXf matrix, inverse;
		int f(std::vector<Nucleotide>);
		int reprobe(int);
	public:
		HashTable(int); 
		void incrementValue(std::vector<Nucleotide>);
		std::string toString();
};
